#!/bin/bash

root -l -b -q PhysicsBookPlots.C'("configs/protons_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/pions_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/nu_energy_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/lep_mom_config.txt")'
root -l -b -q PhysicsBookPlots.C'("configs/cos_lep_theta_config.txt")'
