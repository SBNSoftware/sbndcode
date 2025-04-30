#!/bin/bash

rm /exp/sbnd/data/users/coackley/Nu+E/analysed_DL_uboone/CRUMBS/*.root
rm /exp/sbnd/data/users/coackley/Nu+E/analysed_DL_dune/CRUMBS/*.root
rm /exp/sbnd/data/users/coackley/Nu+E/analysed_Current/NoRefinement/CRUMBS/*.root
rm /exp/sbnd/data/users/coackley/Nu+E/analysed_Cheated/CRUMBS/*.root
rm /exp/sbnd/data/users/coackley/Nu+E/merged.root

count=1
while [ $count -le 10 ]
do
    
	lar -c run_nuE_DLUboone.fcl -s /exp/sbnd/data/users/coackley/Nu+E/reco2_DL_uboone/CRUMBS/"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/Nu+E/analysed_DL_uboone/CRUMBS/"${count}".root
	lar -c run_nuE_DLDune.fcl -s /exp/sbnd/data/users/coackley/Nu+E/reco2_DL_dune/CRUMBS/"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/Nu+E/analysed_DL_dune/CRUMBS/"${count}".root
	lar -c run_nuE_Current.fcl -s /exp/sbnd/data/users/coackley/Nu+E/reco2_Current/NoRefinement/CRUMBS/"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/Nu+E/analysed_Current/NoRefinement/CRUMBS/"${count}".root
	lar -c run_nuE_Cheated.fcl -s /exp/sbnd/data/users/coackley/Nu+E/reco2_Cheated/CRUMBS/"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/Nu+E/analysed_Cheated/CRUMBS/"${count}".root
	((count=count+1))

done

hadd /exp/sbnd/data/users/coackley/Nu+E/merged.root /exp/sbnd/data/users/coackley/Nu+E/analysed_DL_uboone/CRUMBS/*.root /exp/sbnd/data/users/coackley/Nu+E/analysed_DL_dune/CRUMBS/*.root /exp/sbnd/data/users/coackley/Nu+E/analysed_Current/NoRefinement/CRUMBS/*.root /exp/sbnd/data/users/coackley/Nu+E/analysed_Cheated/CRUMBS/*.root
