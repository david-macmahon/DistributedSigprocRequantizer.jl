module DistributedSigprocRequantizer

using Distributed, OnlineStats, Blio, Glob, Statistics, Sockets

export setup_workers, teardown_workers, glob_files, open_files
export dist_fits!, dist_fit!
export get_local_hists, get_global_hist
export write_quantized_files

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

function output_8bit(fname, fbh, fbd, lo, hi; chunk_size=2^20)
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

function write_quantized_files(ws, outnames, lo, hi; chunk_size=2^20)
    fs = map(zip(ws, outnames)) do (w,fn)
        @spawnat w begin
            mkpath(dirname(fn))
            output_8bit(fn, Main.fbh, Main.fbd, lo, hi; chunk_size)
        end
    end
    foreach(wait, fs)
end

end # module DistributedSigprocRequantizer