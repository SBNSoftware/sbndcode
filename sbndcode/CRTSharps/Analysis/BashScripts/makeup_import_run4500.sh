export top_dir=/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500

for FILE in `ls ${top_dir}/raw_files/*run45*.root`
do
    export FILE=$FILE
    export base_file_name=$(basename $FILE ".root")

    if [[ ! -f ${top_dir}/imported_data/${base_file_name}_imported_data.root ]];
    then
	echo $base_file_name
	lar -c import_crt_sharps_data_sbnd.fcl -s ${FILE} -o ${top_dir}/imported_data/${base_file_name}_imported_data.root
    fi
 done
