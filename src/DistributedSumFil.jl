module DistributedSumFil

using DistributedSigprocRequantizer
using Statistics
using ArgParse

function parse_commandline(args=ARGS)
    aps = ArgParseSettings()

    @add_arg_table! aps begin
        "-d", "--dir"
            help = "base input directory"
            arg_type = String
            required = true
        "-o", "--outdir"
            help = "directory for output files (defaults to DIR)"
            arg_type = String
        "-p", "--pattern"
            help = "glob pattern relative to DIR to match files"
            arg_type = String
            required = true
        "-q", "--qlen"
            help = "number of spectra to use for stats"
            arg_type = Int
            default = 10_000
        "-t", "--tail"
            help = "fraction of samples to ignore at tails"
            action = :append_arg
            arg_type = Float32
            range_tester = x->(eps(Float32) < x < 0.5-eps(Float32))
        "HOSTS"
            help = "one or more hosts containing files to 8 bit"
            arg_type = String
            nargs = '+'
            required = true
    end

    return parse_args(args, aps)
end

function main()::Int
    parsed_args = parse_commandline()
    #@show parsed_args; exit()

    #=
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    =#

    hosts = parsed_args["HOSTS"]
    dir = parsed_args["dir"]
    pattern = parsed_args["pattern"]
    qlen = parsed_args["qlen"]
    tail = parsed_args["tail"]
    outdir = get(parsed_args, "outdir", "dir")

    # If no tail options were specified, use 3f-6
    isempty(tail) && push!(tail, 3f-6)

    @info "setting up workers on $(length(hosts)) hosts"
    ws = setup_workers(hosts)

    @info "globbing for input files"
    fnames = glob_files(ws, dir, pattern)
    if all(isempty, fnames)
        @info "no files found"
        teardown_workers(ws)
        return 1
    end

    @info "opening input files"
    wsfexist, fnexist = open_files(ws, fnames)

    @info "gathering stats from $(length(wsfexist)) hosts"
    global_hist = get_global_hist(ws; qlen)

    # Compute lo/hi thresholds
    lo, hi = quantile.(global_hist, (first(tail), 1-last(tail)))
    @info "using thresholds $lo and $hi"

    outbasenames = replace.(basename.(fnexist), r"((\.\d\d\d\d)?\.fil$)"=>s".8\1")
    outnames = joinpath.(outdir, outbasenames)
    
    @info "quantizing and writing output to $outdir"
    write_quantized_files(wsfexist, outnames, lo, hi)

    @info "tearing down workers"
    teardown_workers(ws)

    @info "done"
    return 0
end

end # module DistributedSumFil