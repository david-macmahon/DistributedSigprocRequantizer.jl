#!/bin/bash
#=
PROJDIR=$(dirname $(dirname $(readlink -e "${BASH_SOURCE[0]}")))
exec julia --project=$PROJDIR --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#

using DistributedSigprocRequantizer

exit(@eval Main DistributedSumFil.main())