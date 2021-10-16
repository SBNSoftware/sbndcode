#!/bin/bash

root -l -b -q PhysicsBookPlots.C'("configs/cc0pi/protons_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/cc0pi/nu_energy_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/cc0pi/lep_mom_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/cc0pi/cos_lep_theta_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/cc0pi/proton_mom_config.txt")'
