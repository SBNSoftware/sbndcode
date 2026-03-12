export top_dir=/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/run4500

for FILE in `ls ${top_dir}/reco/*run45*.root`
do
    export FILE=$FILE
    export base_file_name=$(basename $FILE "_reco.root")

    if [[ ! -f ${top_dir}/ana/${base_file_name}_ana.root ]];
    then
	echo $base_file_name
	lar -c ana_crt_sharps_sbnd.fcl -s ${FILE} -o ${top_dir}/ana/${base_file_name}_ana.root
    fi
 done
