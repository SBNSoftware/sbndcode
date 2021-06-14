mrbsetenv && mrbslp

#lar -c run_Select.fcl -S /pnfs/sbnd/scratch/users/lnguyen/prodcorsika_cosmics_job/v09_09_00/gen/prodcorsika_cosmics/crossing_muons_etime/files.list

#lar -c run_Select.fcl -S /pnfs/sbnd/scratch/users/lnguyen/prodcorsika_cosmics_no_diffusion_job/v09_09_00/gen/prodcorsika_cosmics/no_diffusion/files.list

#lar -c run_Select.fcl $(samweb -e sbnd get-file-access-url --schema root reco2-6c529450-fcb6-4c24-b440-001c46b12a3d.root )

#gdb --args lar -c run_Select.fcl $(samweb -e sbnd get-file-access-url --schema root reco2-6c529450-fcb6-4c24-b440-001c46b12a3d.root )

#lar -c eventdump.fcl -s $(samweb -e sbnd get-file-access-url --schema root reco2-6c529450-fcb6-4c24-b440-001c46b12a3d.root ) -n 1 | grep Hit | grep pandora

#lar -c run_Select.fcl -S /pnfs/sbnd/scratch/users/lnguyen/prodcorsika_cosmics_job/v09_09_00/reco2/prodcorsika_cosmics_3ms_events/etime_3ms/files.list -n -1

#lar -c run_Select.fcl -S /pnfs/sbnd/scratch/users/lnguyen/prodcorsika_cosmics_job/v09_09_00/reco2/prodcorsika_cosmics_3ms_no_diffusion_events/etime_3ms/files.list -n -1

#lar -c run_multi_Select.fcl -S 12k_files.list -n 50
lar -c run_Select_SCE.fcl -S 12k_files.list -n 20
