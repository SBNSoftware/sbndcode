void MacroComputeXSec()
{

    
    // Fiducial volume properties
    int fNTPCs = 2;
    double fFVX = 188.5; // cm
    double fFVY = 380.; // cm
    double fFVZ = 480.; // cm
    double fFVcm3 = fNTPCs * (fFVX * fFVY * fFVZ); // cm^3
    double fFV = fFVcm3 * 1e-6; // m^3

    // LAr properties
    double fLArDensity = 1397.3; // kg/m^3
    int fLArNNucleons = 40; // number of nucleons
    int fLArNProtons = 18; // number of protons
    double fLArAtomicMassU = 39.95; // u
    double fLArAtomicMassKg = fLArAtomicMassU * 1.66054e-27; // kg

    // Number of targets
    double fNNuclei = (fFV * fLArDensity )/ fLArAtomicMassKg;
    double fNTargets = fLArNProtons * fNNuclei;

    // Integrated flux
    double fIntegratedFlux = 4.3725e+10; // cm^-2

    // Decay branching ratio
    double fDecayBR = 0.642; // branching ratio

    // Event selection
    int fNEvents = 22; // number of events
    int fNBackground = 10; // number of background events
    double fEfficiency = 0.075; // efficiency

    // Cross section
    double fXSec = (fNEvents - fNBackground) / (fNTargets * fIntegratedFlux * fEfficiency * fDecayBR); // cm^2

    // Print results
    std::cout << " * The fiducial volume is: " << fFV << " m^3 or " << fFVcm3 << " cm^3" << std::endl;
    std::cout<<" * NNuclei is: " << fNNuclei << " NTargets: " << fNTargets << std::endl;
    std::cout << " * The cross section is: " << fXSec << " cm^2" << std::endl;
    
}