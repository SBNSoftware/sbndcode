#include "/sbnd/app/users/hlay/plotting_utils/Structs.h"

const std::map<int, int> razzledMap = { { 11, 0 },
                                        { 13, 1 },
                                        { 22, 2 },
                                        { 211, 3 },
                                        { 2212, 4 },
};

const std::map<int, int> razzleMap = { { 11, 0 },
                                       { 22, 1 },
                                       { 1, 2 },
                                       { 0, 3 },
};

const std::map<int, int> dazzleMap = { { 13, 0 },
                                       { 211, 1 },
                                       { 2212, 2 },
                                       { 1, 3 },
                                       { 0, 4 },
};

const std::map<int, TString> pdgStrings = { { 11, "e^{#pm}" },
                                            { 13, "#mu^{#pm}" },
                                            { 22, "#gamma" },
                                            { 211, "#pi^{#pm}" },
                                            { 2212, "p" }
};

const std::vector<PIDTraining> razzle_trainings = { { "nominal_razzle", "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_22_00/PID/Razzle.weights.xml", false, true, false, "Nominal Razzle", kCyan-3 },
                                                    { "updated_razzle", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/razzle/training/dataset/weights/ShowerPIDMVA_BDTG.weights.xml", false, true, false, "Updated Razzle", kMagenta+2 },
                                                    { "razzled", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/razzled/training/Razzled/weights/Razzled_BDTG.weights.xml", true, false, false, "Razzled", kViolet-5 },
};

const std::vector<PIDTraining> dazzle_trainings = { { "nominal_dazzle", "/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/v01_22_00/PID/Dazzle.weights.xml", false, false, true, "Nominal Dazzle", kCyan-3 },
                                                    { "updated_dazzle", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/dazzle/training/dataset/weights/TrackPIDMVA_BDTG.weights.xml", false, false, true, "Updated Dazzle", kMagenta+2 },
                                                    { "razzled", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/razzled/training/Razzled/weights/Razzled_BDTG.weights.xml", true, false, false, "Razzled", kViolet-5 },
};

const std::vector<PIDTraining> razzled_trainings = { { "razzled", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/razzled/training/Razzled/weights/Razzled_BDTG.weights.xml", true, false, false, "Nominal", kViolet-5 },
                                                     { "razzled_lower_track_length_threshold", "/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv2/razzled_lower_track_length_threshold/Razzled/weights/Razzled_BDTG.weights.xml", true, false, false, "Lower Thresholds", kTeal+6 }
};

const std::vector<TString> razzleAxisLabels  = { "", "e^{#pm}", "#gamma", "Track", "Other" };
const std::vector<TString> dazzleAxisLabels  = { "", "#mu^{#pm}", "#pi^{#pm}", "p", "Shower", "Other" };
const std::vector<TString> razzledAxisLabels  = { "", "e^{#pm}", "#mu^{#pm}", "#gamma", "#pi^{#pm}", "p" };
