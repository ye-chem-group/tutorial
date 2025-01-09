#!/bin/bash

out_path="output"
mkdir -p $out_path

mol="h2o"
basis="cc-pvdz"
kmesh="211"
fxyz="geom/$mol.xyz"
falat="geom/$mol.alat"
fchk="$out_path/mf.chk"
flo="$out_path/lo.chk"
fcderi="$out_path/cderi"
fovL="$out_path/ovL"
loname="PM"
mem="4"
thresh_occ="1e-5"
thresh_vir="1e-6"

# HF
python run_hf.py $fxyz $falat $basis $fchk $kmesh $mem $fcderi 2> $out_path/hf.err > $out_path/hf.stdout

# PM localization; sort LOs by unit cell
python run_lo.py $fchk $flo $loname $mem 2> $out_path/lo.err > $out_path/lo.stdout


# Periodic LNO-CC with PM orbitals (orbital-based fragment)
# cell #1: LO index 0 ~ 3
python run_klno.py $fchk $flo $loname 0,4 $thresh_occ $thresh_vir $mem $fcderi $fovL 2> $out_path/lno04.err > $out_path/lno04.stdout
# cell #2: LO index 4 ~ 7
python run_klno.py $fchk $flo $loname 4,8 $thresh_occ $thresh_vir $mem $fcderi $fovL 2> $out_path/lno48.err > $out_path/lno48.stdout

# Periodic LNO-CC with IAOs (atom-based fragments)
# cell #1: atom index 0 ~ 2
python run_klno_iao.py $fchk 0,3 $thresh_occ $thresh_vir $mem $fcderi $fovL 2> $out_path/lno_iao03.err > $out_path/lno_iao03.stdout
# cell #2: atom index 3 ~ 5
python run_klno_iao.py $fchk 3,6 $thresh_occ $thresh_vir $mem $fcderi $fovL 2> $out_path/lno_iao36.err > $out_path/lno_iao36.stdout
