#!/bin/bash
set -e

DBIN=$(pwd)/bin/

DDATA=results/q5k3/
DORBIT=results/orbit/
mkdir -p ${DDATA}
mkdir -p ${DORBIT}

FORBIT=${DORBIT}/q${q}k${k}_orbit.dat

make -C grp gen_grp_o3
make -C source s2_uniform
make -C source ising_s2_crit

