#!/bin/bash

runs=(4729 4733 4734 4738 4739 4740 4741 4746 4747 4748 4749 4750)

for ((i=0; i<${#runs[@]}; ++i))
do
    export run=${runs[i]}
    mkdir -p /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run${run}/imported_data
    mkdir -p /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run${run}/reco
    mkdir -p /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run${run}/ana

    for FILE in `ls /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run${run}/raw_files/data*.root`
    do
	export FILE=$FILE
	export base_file_name=$(basename $FILE ".root")
	echo $base_file_name
	lar -c import_crt_sharps_data_sbnd.fcl -s $FILE -o /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/imported_data/${base_file_name}_imported_data.root
	lar -c reco_data_crt_sharps_sbnd.fcl -s /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/imported_data/${base_file_name}_imported_data.root -o /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/reco/${base_file_name}_reco.root
	lar -c ana_crt_sharps_sbnd.fcl -s /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/reco/${base_file_name}_reco.root -T /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/ana/${base_file_name}_ana.root
    done
done
