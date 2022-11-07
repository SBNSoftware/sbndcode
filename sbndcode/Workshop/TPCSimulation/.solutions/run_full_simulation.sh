# --- run full simulation (gen, g4, detsim)
nevents=10
path="$MRB_SOURCE/sbndconde/sbndcode/Workshop/TPCSimulation"
lar -c $path/sim_tutorial_gen_non0_T0.fcl -n $nevents -o output_non0T0_${nevents}events_gen.root
lar -c $path/standard_g4_sbnd.fcl -s output_non0T0_${nevents}events_gen.root -o output_non0T0_${nevents}events_g4.root
lar -c $path/standard_detsim_sbnd.fcl -s output_non0T0_${nevents}events_g4.root -o output_non0T0_${nevents}events_detsim.root

# --- run event display
lar -c evd_sbnd.fcl -s output_non0T0_${nevents}events_detsim.root
