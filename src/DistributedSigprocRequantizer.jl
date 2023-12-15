module DistributedSigprocRequantizer

using OnlineStats, Blio, Glob, Statistics
using Distributed

export setup_workers, teardown_workers, glob_files, open_files
export get_global_hist, write_quantized_files

include("DistributedSumFil.jl")
export DistributedSumFil

#=
hosts=["blc0$i" for i in 0:7]
dir = "/datax/dibas/AGBT23A_999_37/GUPPI"
pattern = "BLP*/*_guppi_*_0009.rawspec.0001.fil"
=#

function setup_workers(hosts)
    ws = addprocs(hosts; exeflags="--project=$(dirname(dirname(@__FILE__)))")

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
    wfns = Iterators.filter(wf->!isempty(wf[2]), zip(ws,wfnames))
    fs = map(wfns) do (w,f)
        @spawnat w begin
            Main.fbh, Main.fbd = Filterbank.mmap(f[1])
            Main.fbh
        end
    end
    fbhs = fetch.(fs)
    wwf = first.(wfns)
    fnexist = first.(last.(wfns))
    # wwf is workers with files
    # fnexist is filenames that exist on workers with files
    # fbhs is filterbank headers of those files
    wwf, fnexist, fbhs
end

function get_extremas(ws; qlen=10_000)
    fs = map(ws) do w
        @spawnat w begin
            qlen = min(qlen, size(Main.fbd, 3))
            fit!(Extrema(Float32), @view(Main.fbd[:,:,1:qlen]))
        end
    end
    fetch.(fs)
end

function make_hist(extremas; closed=true)
    global_extrema = Extrema(Float32)
    foreach(e->merge!(global_extrema, e), extremas)
    r = floor(global_extrema.min):ceil(global_extrema.max)
    Hist(r, Float32; closed)
end

function fit_hists(ws, h; qlen=10_000)
    fs = map(ws) do w
        @spawnat w begin
            qlen = min(qlen, size(Main.fbd, 3))
            fit!(h, @view(Main.fbd[:,:,1:qlen]))
        end
    end
    fetch.(fs)
end

function get_global_hist(ws; qlen=10_000, closed=true)
    @info "getting local extremas"
    @time local_extremas = get_extremas(ws; qlen)
    global_hist = make_hist(local_extremas; closed)
    @info "fitting local hists"
    # Empty global_hist is copied to each worker
    @time local_hists = fit_hists(ws, global_hist; qlen)
    @info "merging local hists into global hist"
    @time foreach(e->merge!(global_hist, e), local_hists)

    global_hist
end

function quantize(x, lo, hi, d=hi-lo)::UInt8
#try
    UInt8(
        ifelse(x < lo, 0f0,
            ifelse(x > hi, 255f0,
                round((x-lo)/d * 255f0)
            )
        )
    )
#catch
#    @show x lo hi d
#    rethrow()
#end
end

function output_8bit(fname, fbh, fbd, lo, hi; chunk_size=2^20)
    nspec_per_chunk = fld(chunk_size, fbh[:sample_size])
    vws=(@view(fbd[:,:,t]) for t in Iterators.partition(axes(fbd,3), nspec_per_chunk))
    q8 = similar(first(vws), UInt8)
    fbh[:nbits] = 8
    open(fname, "w") do io
        write(io, fbh)
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