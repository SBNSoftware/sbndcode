#!/bin/bash

root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/protons_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/pions_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/nu_energy_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/lep_mom_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/cos_lep_theta_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/proton_mom_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/nu_energy_config_fsi.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/lep_mom_config_fsi.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nuSTORM2023/ccinc/cos_lep_theta_config_fsi.txt")'
