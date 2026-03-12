export top_dir=/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500

for FILE in `ls ${top_dir}/imported_data/*run45*.root`
do
    export FILE=$FILE
    export base_file_name=$(basename $FILE "_imported_data.root")

    if [[ ! -f ${top_dir}/reco/${base_file_name}_reco.root ]];
    then
	echo $base_file_name
	lar -c reco_data_crt_sharps_sbnd.fcl -s ${FILE} -o ${top_dir}/reco/${base_file_name}_reco.root
    fi
 done
