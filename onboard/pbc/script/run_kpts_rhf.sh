#!/bin/bash

#######################
#     LIBRARY PATH    #
#######################
pyscf_path= # fill in yourself

ncore="12"
mem0="4"    # mem/core, GB
mem=`echo $ncore $mem0 | awk '{print $1*$2}'`
walltime="12:00:00" # 12 h

run_path="$PWD"
exe="$run_path/run_kpts_rhf.py"

calc_path="$run_path/run"
mkdir -p $calc_path
cd $calc_path

fchk="mf.chk"
pseudo="none"   # all-electron potential

for alat in 3.567
do
    mkdir -p $alat
    cd $alat

    for zeta in dz
    do
        mkdir -p $zeta
        cd $zeta

        basis="cc-pv${zeta}"

        for kmesh in 111 222 333 444 555 666
        do
            mkdir -p $kmesh
            cd $kmesh

            jobname="${alat}_${zeta}_${kmesh}"

            # python $exe $fchk $basis $alat $pseudo $kmesh $mem 2> err > stdout
            cat > submit.sh << eof
#!/bin/bash

#SBATCH -e sbatch.err
#SBATCH -o sbatch.log
#SBATCH -n 1
#SBATCH -c $ncore
#SBATCH -t $walltime
#SBATCH --mem-per-cpu=${mem0}gb

export OMP_NUM_THREADS=$ncore
export MKL_NUM_THREADS=1

module restore pyscf

export PYTHONPATH=${pyscf_path}:\${PYTHONPATH}

python $exe $fchk $basis $alat $pseudo $kmesh $mem 2> err > stdout
eof
            sbatch -J $jobname submit.sh

            cd .. # $kmesh
        done

        cd .. # $zeta
    done

    cd .. # $alat
done
