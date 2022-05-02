#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>



int convertToHepevt_argo(std::string fname, std::string lifetime) {
  
  fstream eventfile;
  eventfile.open(fname.c_str());

  // vectors for the different components in the file
  std::vector<string> particle;
  std::vector<double> pdg;
  std::vector<double> px; std::vector<double> py; std::vector<double> pz;
  std::vector<double> x; std::vector<double> y; std::vector<double> z; 
  std::vector<double> mass;
  std::vector<double> E; std::vector<double> t;
  std::vector<double> status;

  std::string line;

  // {E_e^-, px_e^-, py_e^-, pz_e^-, E_e^+, px_e^+, py_e^+, pz_e^+, E_{darkNeutrino}, Vrt_z, weight} in the file.  
  //TRandom3 RandGen;
  TRandom3 RandGen(0);

  double crt_upper_z = - 165;// cm
  while(!eventfile.eof()) {
    
    // random number generation for the particle statrting position for particles for simulations. 
    // assuming that they are even-distributed. 

    // x: [-395.831226cm, 460.95 cm] -- for upperstream panel
    // x: [-145.570605,140.0] - for downstream panel
    // Y: [-390.95, -105.370026]
    double x_lim[2] = {-145., 140.}; //cm
    double y_lim[2] = {-390., -100.};
    double newx = RandGen.Uniform(x_lim[0], x_lim[1]);//(2*RandGen.Rndm()-1)*200; // cm
    double newy = RandGen.Uniform(y_lim[0], y_lim[1]);//RandGen.Rndm()*400; 
    
    int PDG[2]={11,-11}; // read in electron and positron 4-momentum
    double massloc[3]={0.5110000000e-03,0.5110000000e-03};
    for(int ipart=0;ipart<2;ipart++){
      double Eloc,pxloc,pyloc,pzloc;   
      eventfile >>  Eloc >> pxloc >> pyloc >> pzloc; 
      pdg.push_back(PDG[ipart]);  
      mass.push_back(massloc[ipart]); status.push_back(1);   
      px.push_back(pxloc); 
      py.push_back(pyloc);
      pz.push_back(pzloc);
      E.push_back(Eloc);

      t.push_back(0);
      x.push_back(newx); 
      y.push_back(newy); 
    } 
    double Etemp, Vrt_z_temp, weight_temp;
    eventfile >> Etemp >> Vrt_z_temp >> weight_temp;
    for (int i = 0; i<2; i++){
      // crt panel is from -165.0 cm to 781.25 cm. 
      z.push_back(crt_upper_z + Vrt_z_temp*100); //adding vertex position in z-direction to the electron-positron pair.
      // File from the theorist is using [cm] as unit. 
    }

    //dark neutrino.
    pdg.push_back(9999);  
    mass.push_back(0.175000000e-01); status.push_back(0);   
    px.push_back(0); 
    py.push_back(0);
    pz.push_back(0);

    E.push_back(Etemp);
    t.push_back(0);
    x.push_back(weight_temp); //using x_position of dark neutrino to store the weight as in LArSoft doesn't support weight. 
    y.push_back(0); 
    z.push_back(0);  

    
  }

  for (int j=0; j<20; j++) {
    std::cout<<pdg[j]<<" "<<mass[j]<<" "<<px[j]<<" "<<py[j]<<" "<<pz[j]<<" "<<E[j]<<" "<<x[j]<<" "<<y[j]<<" "<<z[j]<<" "<<t[j]<<std::endl;
  }
  
  eventfile.close();

  // open file to write

  //ofstream outfile;
  //outfile.open("BPCevent.hepevt");
  int counter(0); int file(0);

  const int nparticlesinevent=3;
  const int nevents_in_file=300000;

  ofstream outfile;
  std::stringstream ss;

  ofstream outfile_weight; 
  std::stringstream ss1;

  ss  << "DarkNeutrinos_hepevt/DarkNeutrino_"<<lifetime<<"lifetime_file.hepevt";

  outfile.open(ss.str());

  for (int e=0; e<nevents_in_file; e++) {
    outfile<<e<<" " << nparticlesinevent <<std::endl;
    for(int n=0; n<nparticlesinevent; n++) {
      int m = n + e*3;
      outfile<<status[m]<<" "<<pdg[m]<<" 0 0 0 0 "<<px[m]<<" "<<py[m]<<" "<<pz[m]<<" "<<E[m]<<" "<<mass[m]<<" "<<x[m]<<" "<<y[m]<<" "<<z[m]<<" "<<t[m]<<std::endl;
    }
  }
  
  outfile.close();




/*


  for (int k=0; k<10000; k++) {
    ofstream outfile;
    if(counter == 0) {
      
      file++;
      std::stringstream ss;
    
      ss << "BPBeventsApr2/BPBevent_" << file << ".hepevt";

      //ofstream outfile;
      outfile.open(ss.str());  
    }

    outfile<<k<<" 6"<<std::endl;
    for(int n=0; n<6; n++) {
      int m=n+6*k;
      outfile<<status[m]<<" "<<pdg[m]<<" 0 0 0 0 "<<px[m]<<" "<<py[m]<<" "<<pz[m]<<" "<<E[m]<<" "<<mass[m]<<" "<<x[m]<<" "<<y[m]<<" "<<z[m]<<" "<<t[m]<<std::endl;
    }
    counter++;

    if (counter==50) {
      counter =0; 
      outfile.close();  
    }


  }
*/

  return 0;
}
