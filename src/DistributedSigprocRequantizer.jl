module DistributedSigprocRequantizer

using Distributed, OnlineStats, Blio, Glob, Statistics, Sockets, ProgressBars

export setup_workers, teardown_workers, glob_files, open_files
export make_unit_hist, make_8bit_hist, dist_channel_fits
export dist_fits!, dist_fit!
export get_local_hists, get_global_hist
export write_quantized_files, splice_quantized_file

include("DistributedSumFil.jl")
export DistributedSumFil

function setup_workers(hosts)
    ws = addprocs(hosts .* " " .* string.(getaddrinfo.(hosts));
                  exeflags="--project=$(dirname(dirname(@__FILE__)))")

    @eval Main using Distributed
    @eval Main @everywhere $(ws[1]) using DistributedSigprocRequantizer
    @eval Main @everywhere $(ws[2:end]) using DistributedSigprocRequantizer
    @eval Main @everywhere $(ws[1]) using OnlineStats, Blio, Glob
    @eval Main @everywhere $(ws[2:end]) using OnlineStats, Blio, Glob

    ws
end

function teardown_workers(ws)
    rmprocs(ws)
end

function glob_files(ws, dir, pattern)
    fs = map(w->@spawnat(w, glob(pattern, dir)), ws)
    fnames = fetch.(fs)
    if all(isempty, fnames)
        @warn "no files found on any host"
    elseif any(isempty, fnames)
        @warn "missing files on one or more hosts"
    elseif any(v->length(v)>1, fnames)
        @warn "multiple matching files on one or more hosts"
    end
    fnames
end

function open_files(ws::AbstractVector{Int},
                    wfnames::AbstractVector{<:AbstractVector{String}})
    # Only work with workers on which files were found
    wfns = Iterators.filter(wf->!isempty(wf[2]), zip(ws,wfnames))
    # wwf is workers with found files
    wwf = first.(wfns)
    # fnexist is filenames that were found on workers
    fnexist = first.(last.(wfns))

    fs = map(wfns) do (w,f)
        @spawnat w begin
            Main.fbh, Main.fbd = Filterbank.mmap(f[1])
            Main.fbh
        end
    end
    # fbhs is filterbank headers for files that were found on workers
    fbhs = fetch.(fs)

    # Permute workers with files, filenames, and headers into descending
    # frequency order
    p = sortperm(fbhs, by=h->h[:fch1], rev=true)
    permute!(wwf, p)
    permute!(fnexist, p)
    permute!(fbhs, p)

    wwf, fnexist, fbhs
end

function make_unit_hist(lo, hi; closed=true)
    r = lo-0.5:hi+0.5
    Hist(r, Float32; closed)
end

make_unit_hist(lohi::Tuple; closed=true) = make_unit_hist(first(lohi), last(lohi); closed)
make_unit_hist(e; closed=true) = make_unit_hist(extrema(e); closed)

function make_8bit_hist(lo, hi; closed=true)
    r = range(lo, stop=hi, length=257)
    Hist(r, Float32; closed)
end

make_8bit_hist(lohi::Tuple; closed=true) = make_8bit_hist(first(lohi), last(lohi); closed)
make_8bit_hist(e; closed=true) = make_8bit_hist(extrema(e); closed)

function dist_channel_fits(o, ws; qlen=5_000)
    fs = map(ws) do w
        @spawnat w begin
            ochans = [copy(o) for _ in axes(Main.fbd, 1)]
            nsamps = size(Main.fbd, 3)
            qlen = clamp(qlen, 1, nsamps)
            # Fit first qlen samples for each channel
            fbd_chans = eachslice(@view(Main.fbd[:,:,1:qlen]); dims=1)
            for (oc, fbd_chan) in zip(ochans, fbd_chans)
                fit!(oc, fbd_chan)
            end
            # If more samples exist beyond qlen
            if qlen < nsamps
                # Subsample them
                stride = cld(nsamps-qlen, qlen)
                fbd_chans = eachslice(@view(Main.fbd[:,:,qlen+1:stride:end]); dims=1)
                for (oc, fbd_chan) in zip(ochans, fbd_chans)
                    fit!(oc, fbd_chan)
                end
            end
            ochans
        end
    end
    fetch.(fs)
end

function dist_fits!(o, ws; qlen=5_000)
    fs = map(ws) do w
        @spawnat w begin
            nsamps = size(Main.fbd, 3)
            qlen = clamp(qlen, 1, nsamps)
            # Fit first qlen samples
            fit!(o, @view(Main.fbd[:,:,1:qlen]))
            # If more samples exist beyond qlen
            if qlen < nsamps
                # Subsample them
                stride = cld(nsamps-qlen, qlen)
                fit!(o, @view(Main.fbd[:,:,qlen+1:stride:end]))
            end
            o
        end
    end
    fetch.(fs)
end

function dist_fit!(o, ws; qlen=5_000)
    local_fits = dist_fits!(o, ws; qlen)
    foreach(l->merge!(o, l), local_fits)
    o
end

function get_local_hists(ws; qlen=5_000, closed=true)
    global_extrema = dist_fit!(Extrema(Float32), ws; qlen)
    edges = global_extrema.min-0.5f0:global_extrema.max+0.5f0
    dist_fits!(Hist(edges, Float32; closed), ws; qlen)
end

function get_global_hist(ws; qlen=5_000, closed=true)
    local_hists = get_local_hists(ws; qlen, closed)
    # Merge all local hists into local_hists[1]
    foldl(merge!, local_hists)
    local_hists[1]
end

function quantize(x, lo, hi, d=hi-lo)
    ifelse(x < lo, 0f0,
        ifelse(x > hi, 255f0,
            round((x-lo)/d * 255f0)
        )
    )
end

function output_8bit(fname::AbstractString, fbh, fbd, lo, hi; chunk_size=2^20)
    nspec_per_chunk = fld(chunk_size, fbh[:sample_size])
    vws=(@view(fbd[:,:,t]) for t in Iterators.partition(axes(fbd,3), nspec_per_chunk))
    q8 = similar(first(vws), UInt8)
    # Copy fbh so that we don't modify fbh (so fbh will still reflect fbd)
    outfbh = Filterbank.Header(fbh)
    outfbh[:nbits] = 8
    open(fname, "w") do io
        write(io, outfbh)
        for v in vws
            if size(v) != size(q8)
                q8 = similar(v, UInt8)
            end
            @. q8 = quantize(v, lo, hi)
            write(io, q8)
        end
    end
end

function output_8bit(chan::Union{AbstractChannel,RemoteChannel}, fbh, fbd, lo, hi; chunk_size=2^20)
    nspec_per_chunk = fld(chunk_size, fbh[:sample_size])
    vws=(@view(fbd[:,:,t]) for t in Iterators.partition(axes(fbd,3), nspec_per_chunk))
    q8 = similar(first(vws), UInt8)

    for v in vws
        if size(v) != size(q8)
            q8 = similar(v, UInt8)
        end
        @. q8 = quantize(v, lo, hi)
        try
            put!(chan, q8)
        catch ex
            ex isa InvalidStateException || rethrow()
            return nothing
        end
    end

    # Put zero sized Array into channel to signify EOF
    try
        put!(chan, similar(q8, ntuple(_->0, ndims(q8))))
    catch ex
        ex isa InvalidStateException || rethrow()
    end

    nothing
end

function write_quantized_files(ws, outnames, lo, hi; chunk_size=2^20)
    fs = map(zip(ws, outnames)) do (w,fn)
        @spawnat w begin
            mkpath(dirname(fn))
            output_8bit(fn, Main.fbh, Main.fbd, lo, hi; chunk_size)
        end
    end
    foreach(wait, fs)
end

function splice_quantized_file(ws, fbhs, outname, lo, hi; chunk_size=2^20)

    # Create a remote channel for each worker that will each hold up to 32
    # Array{UInt8,3} instances
    chans = [RemoteChannel(()->Channel{Array{UInt8,3}}(32)) for _ in ws]

    # Spawn quantization on workers
    fs = map(zip(ws, chans)) do (w,ch)
        @spawnat w begin
            output_8bit(ch, Main.fbh, Main.fbd, lo, hi; chunk_size)
        end
    end

    # Copy fbhs[1] so that we don't modify the caller's copy!!!
    outfbh = Filterbank.Header(fbhs[1])

    # Create splicebuf Array into which we will copy chunks from workers
    nspec_per_chunk = fld(chunk_size, outfbh[:sample_size])
    splicebuf = Array{UInt8}(undef, outfbh[:nchans], length(ws), outfbh[:nifs], nspec_per_chunk)

    # Compute number of output spectra
    nspec_total = minimum(get.(fbhs, :nsamps))
    pb = ProgressBar(total=nspec_total, unit="spectra")

    # Update outfbh for spliced 8-bit file.  For now assume all nodes have data
    # that are contiguous in frequency.
    outfbh[:nbits] = 8
    outfbh[:nchans] *= length(ws)

    open(outname, "w") do io
        write(io, outfbh)

        # Loop through chunks from workers
        while true
            # Take data chunks from channels
            chunks = take!.(chans)

            # Get number of spectra in each chunk
            nspecs = size.(chunks, 3)

            # If any chunk has 0 spectra then we are done
            if any(==(0), nspecs)
                @info "got empty chunk, done"
                break
            # Else if all chunks have the nspec_per_chunk spectra, output them
            elseif all(==(nspec_per_chunk), nspecs)
                for (dst, src) in zip(eachslice(splicebuf; dims=2), chunks)
                    copyto!(dst, src)
                end
                write(io, splicebuf)
                update(pb, nspec_per_chunk)
            # Else, find smallest number of spectra, make views for that many
            # spectra of each chunk, output views, and then we are done
            else
                @info "got short chunk, done"
                nmin = minimum(nspecs)
                vw = @view splicebuf[:, :, :, 1:nmin]
                for (dst, src) in zip(eachslice(splicebuf; dims=2), chunks)
                    # Copy first nmin spectra from src chunk
                    copyto!(dst, @view(src[:,:,1:nmin]))
                end
                write(io, vw)
                update(pb, nmin)
                break
            end
        end # while true
    end # open

    # Close channels (will cause workers to finish)
    foreach(close, chans)

    # Wait for workers to finish (should already be done)
    foreach(wait, fs)
end

end # module DistributedSigprocRequantizer
