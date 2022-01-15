#!/bin/bash

echo " Job begin time: " `date`

export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12
export OMP_STACKSIZE=1000m
#ulimit -s unlimited

mkdir c15 ### results will be saved
touch c15/network.tsv ### the reaction network sampled by reactor

ids=(bisabolyl) ### the list of initial structures
for id in ${ids[@]}
do
        for step in `seq 1 5` ### here 5 MTD simulations will be conducted for every carbocation
        do
            xtb --md --input metadyn.inp ${id}.xyz --etemp 298 >/dev/null ### conduct MTD
            mv xtb.trj c15/xtb_${id}_${step}.xyz
            init=0
            for num in `seq 1 3000` ### the structures will be saved every 5fs so there is a total of 3000 structures in 30ps MTD trajectory
            do
                line_begin=$(($(expr $num \* 42)-41))
                line_end=$(expr $num \* 42)
                sed -n ''"${line_begin}","${line_end}"'p' xtb_${id}_${step}_opt.xyz > temp.xyz
                xtb temp.xyz --opt --charge 1 > temp.out ### structure optimization for every snapshot
                cat xtbopt.xyz >> c15/xtb_${id}_${step}_opt.xyz
                ### get the raw SMILES
                flag=`printf "%02d" \`python charge.py\`` ### locate the most positive carbocation atom
                if [ $flag -eq 0 ];then
                    continue
                fi
                smi=`python xyz2mol.py xtbopt.xyz -o smiles -c 1 --no-graph` ### get SMILES from xyz
                if [[ "$smi" == *"."* ]] || [[ $smi == "" ]]; then
                    continue
                fi
                if [[ "$smi" == *" "* ]];then
                    continue
                fi
                if [ $init -eq 0 ];then
                    init=1
                    site=$flag
                    cp xtbopt.xyz sub.xyz
                    smi_sub=$smi
                    num_sub=$num
                    continue
                fi
                if [ $flag -eq $site ];then
                    cp xtbopt.xyz sub.xyz
                    num_sub=$num
                    continue
                fi
                ### reaction will be detected once the flag has changed
                smi_prod=$smi
                cp xtbopt.xyz prod.xyz
                cp sub.xyz c15/${id}_${step}_${num_sub}.xyz
                cp prod.xyz c15/${id}_${step}_${num}.xyz
                ### estimate reaction with RMSD-PP
                xtb sub.xyz --path prod.xyz --input path.inp > path.out 
                barrier=(`grep 'forward  barrier' path.out | awk '{print $5}'`)
                energy=(`grep 'reaction energy' path.out | awk '{print $5}'`)
                minData=9999.9
                for i in ${barrier[@]}
                do
                    if [[ ${i} < ${minData} ]];then
                        minData=${i}
                    fi
                done
                echo -e ${id}_${step}_${num_sub}"\t"${smi_sub}"\t"${id}_${step}_${num}"\t"${smi_prod}"\t"$minData"\t"${energy[0]} >> c15/network.tsv
                site=$flag
                smi_sub=$smi
                num_sub=$num
            done
        done
done
echo " Job finish time: " `date`
