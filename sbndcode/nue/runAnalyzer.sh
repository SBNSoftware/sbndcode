#!/bin/bash

count=1
while [ $count -le 10 ]
do
	lar -c run_nuE.fcl -s /exp/sbnd/data/users/coackley/Nu+E/reco2_Cheated/CRUMBS/"${count}".root -n -1 -T /exp/sbnd/data/users/coackley/Nu+E/analysed_Cheated/CRUMBS/"${count}".root
	((count=count+1))

done
