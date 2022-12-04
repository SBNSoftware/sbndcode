for FILE in `ls /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/raw_files/*run4534*.root`
do
    export FILE=$FILE
    export base_file_name=$(basename $FILE ".root")
    echo $base_file_name
    lar -c import_crt_sharps_data_sbnd.fcl -s $FILE -o /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/imported_data/${base_file_name}_imported_data.root
    lar -c reco_data_crt_sharps_sbnd.fcl -s /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/imported_data/${base_file_name}_imported_data.root -o /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/reco/${base_file_name}_reco.root
    lar -c ana_crt_sharps_sbnd.fcl -s /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/reco/${base_file_name}_reco.root -T /pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500/ana/${base_file_name}_ana.root
done
