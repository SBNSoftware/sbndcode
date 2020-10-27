output_file="file_list.txt"
samweb_definition="official_MCP2.0_prodsingle_electron_pi+_bnblike_forward_reco_sbnd"
nfiles_limit=100
limit_files=true
job_name=""

# initialize file count
file_id=0

for file in $(samweb -e sbnd list-definition-files $samweb_definition); do
	if [ $file_id -le $nfiles_limit ]; then
		file_location=$(samweb -e sbnd get-file-access-url --schema=root $file)
		echo ${file_location/root:\/\/fndca1.fnal.gov:1094\/pnfs\/fnal.gov\/usr/\/pnfs} >> output_file
	fi
	if [ $limit_files == true ]; then
		file_id=$((file_id+1))
	fi
done;