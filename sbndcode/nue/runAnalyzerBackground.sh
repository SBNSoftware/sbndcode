#!/bin/bash

rm /exp/sbnd/data/users/coackley/BNBEvents/analysed_DL_uboone/*.root
rm /exp/sbnd/data/users/coackley/BNBEvents/analysed_Current/*.root
rm /exp/sbnd/data/users/coackley/BNBEvents/merged.root

count=1
while [ $count -le 300 ]
do
    
	lar -c run_nuE_DLUboone.fcl -s /exp/sbnd/data/users/coackley/BNBEvents/reco2_DL_uboone/CRUMBS"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/BNBEvents/analysed_DL_uboone/"${count}".root
	lar -c run_nuE_Current.fcl -s /exp/sbnd/data/users/coackley/BNBEvents/reco2_Current/CRUMBS/"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/BNBEvents/analysed_Current/"${count}".root
    ((count=count+1))

done

hadd /exp/sbnd/data/users/coackley/BNBEvents/merged.root /exp/sbnd/data/users/coackley/BNBEvents/analysed_DL_uboone/*.root /exp/sbnd/data/users/coackley/BNBEvents/analysed_Current/*.root
