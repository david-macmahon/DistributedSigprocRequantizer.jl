# DistributedSigprocRequantizer.jl

`DistributedSigprocRequantizer.jl` is a package for requantizing related Sigproc
Filterbank files that are spread across multiple computers before splicing them.
Requantizing on multiple computers benefits from parallelization of the
requantization step.  Requantizing to a smaller data type also leads to faster
splicing due to smaller files.

Currently the only requantization supported is from `Float32` data to `UInt8`
data.  Historically, this has been performed by the `sum_fil` utility of the
Sigproc library, but `sum_fil` only operates on a single input file.  If the
data from a given observation have been split across multiple files (e.g. due to
the nature of the data acquisition instrumentation) `sum_fil` requires splicing
the input data files together into a single file prior to requantization.

`DistributedSigprocRequantizer.jl` operates in a distributed manner by launching
worker processes on remote computers where the input files reside.  This allows
the requantization to occur on each input file in parallel with the others.  To
calculate the global offset and scaling factor for the full bandwidth across all
input files, the workers compute statistics for their local input files.  These
local statistics are combined into global statistics for the entire dataset.
Global offset and scale factors are derived from the global statistics in the
same way as the original `sum_fil`.  These values are then used by all the
workers to requantize the data from their local input files and write out new
Filterbank files with `UInt8` data.  The new Filterbank files are suitable for
splicing and will splice approximately four times faster than the original
`Float32` files given their smaller size.

# Installation

The recommended way to install this package is to clone the git repository to
somewhere in/under your home directory.  Then you will need to install the
additional packages used by this package.  One of those packages, `Blio`, is not
yet in the General package registry so you will have to install it from GitHub
manually before using the `instantiate` command to install the remaining
packages automatically.

After cloning this repository, change into the top level directory where this
file and, more importantly, the `Project.toml` file are located.  Then start
Julia with the `--project` flag.  Once you are at the `julia>` prompt, press `]`
to enter the built-in package manager.  The prompt should change to:

```
(DistributedSigprocRequantizer) pkg>
```

At this prompt, enter the following command:

```
add https://github.com/david-macmahon/Blio.jl
```

When that completes, you will be able to instantiate the rest of the packages
automatically by running:

```
instantiate
```

When that finishes you can press the backspace key to exit the package manager
and return to the `julia>` prompt.  The installation is complete and you may
exit Julia.

To test that the package is installed correctly, you can run the following
command from your shell (with the current directory still the same as the
directory containing this file) and you should see the following output:

```
$ ./bin/dist_sum_fil.jl -h
usage: dist_sum_fil.jl -d DIR [-o OUTDIR] -p PATTERN [-q QLEN]
                       [-t TAIL] [-h] HOSTS...

positional arguments:
  HOSTS                 one or more hosts containing files to 8 bit

optional arguments:
  -d, --dir DIR         base input directory
  -o, --outdir OUTDIR   directory for output files (defaults to DIR)
  -p, --pattern PATTERN
                        glob pattern relative to DIR to match files
  -q, --qlen QLEN       number of spectra to use for stats (type:
                        Int64, default: 10000)
  -t, --tail TAIL       fraction of samples to ignore at tails (type:
                        Float32, default: 3.0f-6)
  -h, --help            show this help message and exit
```

# Usage

The main utility of this package is the `dist_sum_fil.jl` script.  It is located
in the `bin` subdirectory.  The script needs to be told the common root
directory that the input files live under.  This directory, which must be
identical across all hosts being used, is specified with the `-d` option.  To
find input files on each host, the workers will use the *glob* pattern specified
with the `-p` option to find input files.  The glob pattern is relative to the
common base directory.  The script also must be given the list of hosts on which
to run worker processes that will actually look for the input files and process
them.  The hosts are simply given as positional parameters.

By default, the 8 bit output files will be created in the same directory as the
input files, but a different output directory may be given with the `-o` option.
Other options include specifying how many input spectra are examined when
determining the offset and scale values (`-q`) as well as what fraction of
samples at the tail ends of input values can be ignored when computing the
offset and scale values (`-t`).  Samples with values in these tail regions will
be clamped to the minimum and maximum values of the output data type.

# Example

```
$ dist_sum_fil.jl -d /datax/dibas/AGBT23B_999_08/GUPPI \
                  -p '*/*_guppi_60277_04433_PSR_B0355+54_0010.rawspec.0001.fil' \
                  -o /datax/users/davidm \
                  blc{0..3}{0..7}
[ Info: setting up workers on 32 hosts
[ Info: globbing for input files
[ Info: opening input files
[ Info: gathering stats from 32 hosts
[ Info: getting local extremas
  2.579612 seconds (85.70 k allocations: 5.777 MiB, 4.76% compilation time)
[ Info: fitting local hists
  3.073296 seconds (94.46 k allocations: 13.153 MiB, 5.52% compilation time)
[ Info: merging local hists into global hist
  0.204158 seconds (131.54 k allocations: 8.749 MiB, 8.07% gc time, 99.86% compilation time)
[ Info: using thresholds 0.5 and 463.5
[ Info: quantizing and writing output to /datax/users/davidm
[ Info: tearing down workers
[ Info: done
```
