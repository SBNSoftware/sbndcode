#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <map>
#include "TTree.h"
#include "TFile.h"

void inclAnna_to_hepmc_sbnd_filteredForVisibleCluster(std::string indir , std::string inclxx_file_name , std::string outdir, int eventsperfile=100) {

   //assumes your file is inclxx_file_name.root

   map<int, int> PDGmap = {{1002,1000010020}, {1003,1000010030}, {2003,1000020030}, {2004, 1000020040}};
   map<int, float> EnergyMap = {{1002,64.44}, {1003,76.86}, {2003,169.27}, {2004, 191.71}};

   std::string cmd = "mkdir -p " + outdir;
   gSystem->Exec(cmd.c_str());

   // Input file
   std::string infile = indir + inclxx_file_name + ".root";
   TFile* f = new TFile(infile.c_str());
   
   TTreeReader et_reader("et", f);

   TTreeReader at_reader("at", f);

   //number of events converted
   int event_count=0;

   TTreeReaderArray<float> Ekin(et_reader,"EKin");
   TTreeReaderArray<float> px(et_reader,"px");
   TTreeReaderArray<float> py(et_reader,"py");
   TTreeReaderArray<float> pz(et_reader,"pz");
   TTreeReaderArray<int> pdg(et_reader,"PDGCode");
   TTreeReaderArray<short> A(et_reader,"A");
   TTreeReaderArray<short> Z(et_reader,"Z");
   TTreeReaderValue<short> nParticles(et_reader,"nParticles");

   TTreeReaderValue<float> protonPx(at_reader,"protonPx");
   TTreeReaderValue<float> protonPy(at_reader,"protonPy");
   TTreeReaderValue<float> protonPz(at_reader,"protonPz");
   TTreeReaderValue<float> neutronPx(at_reader,"neutronPx");
   TTreeReaderValue<float> neutronPy(at_reader,"neutronPy");
   TTreeReaderValue<float> neutronPz(at_reader,"neutronPz");
   TTreeReaderValue<float> muonPx(at_reader,"muonPx");
   TTreeReaderValue<float> muonPy(at_reader,"muonPy");
   TTreeReaderValue<float> muonPz(at_reader,"muonPz");
   TTreeReaderValue<float> rx(at_reader,"rx");
   TTreeReaderValue<float> ry(at_reader,"ry");
   TTreeReaderValue<float> rz(at_reader,"rz");
   TTreeReaderValue<float> nuPz(at_reader,"nuPz");


   // Output file
   std::ofstream _hepmc_file;

   //std::cout << outdir + inclxx_file_name + "_" + std::to_string(0) + ".hepmc" << std::endl;

   _hepmc_file.open(outdir + inclxx_file_name + "_" + std::to_string(0) + ".hepmc");




   double GlobalTimeOffset = 3125.;
   double RandomTimeOffset = 1600.;

   srand(time(0));

   while (et_reader.Next() && at_reader.Next()){

      if (event_count % 1000 == 0) std::cout << "At event " << event_count << std::endl;

      if (event_count % eventsperfile == 0) {

         //  std::cout << outdir+inclxx_file_name+"_"+std::to_string(i/50)+".hepmc" << std::endl;

         // Close and open a new file

         _hepmc_file.close();
         _hepmc_file.open(outdir + inclxx_file_name + "_" + std::to_string(event_count / eventsperfile) + ".hepmc");


      }

      /*
      std::cout << *nParticles << std::endl;
      for(int i= 0; i < A.GetSize(); i++){
         std::cout<< A[i] << " " << Z[i] << " " << pdg[i] << std::endl;
      }
      for(int i= 0; i < pdg.GetSize(); i++){
         std::cout<< A[i] << " " << Z[i] << " " << pdg[i] << std::endl;
      }
      std::cout << Ekin.GetSize() << std::endl;
      std::cout << Ekin[0] << std::endl;
      std::cout << Ekin[1] << std::endl;
      std::cout << *protonPx << std::endl;
      std::cout << *rx << std::endl;

      */

      // Filter for nuclear clusters

      bool containsCluster = 0;
      bool containsVisibleCluster = 0;

      for (size_t j = 0; j < *nParticles; j++) {

         if(pdg[j] == 1002 || pdg[j] == 1003 || pdg[j] == 2003 || pdg[j] == 2004){
            containsCluster = 1;
         }

         if(pdg[j] == 1002 || pdg[j] == 1003 || pdg[j] == 2003 || pdg[j] == 2004){
            if(Ekin[j] > EnergyMap[pdg[j]]){
               containsVisibleCluster = 1;
            }
         }

      }

      if(!containsVisibleCluster)
      continue;

      // Neutrino time in the spill
      double nu_time = rand() / double(RAND_MAX) * RandomTimeOffset + GlobalTimeOffset;



      // Save the number of particles for this events (+1 for the neutrino)
      _hepmc_file << event_count << " " << *nParticles + 2 << std::endl;


      // Save the neutrino
      _hepmc_file << 0 << " "
         << 14 << " "
         << -1 << " " << 0 << " " 
         << 0 << " "                         // Should be 1st and 2nd daugther, but putting ccnc and mode
         << 0 << " "                         // Should be 1st and 2nd daugther, but putting ccnc and mode
         << 0. << " "       // NuWro uses MeV, we want GeV
         << 0. << " "       // NuWro uses MeV, we want GeV
         << *nuPz / 1000. << " "       // NuWro uses MeV, we want GeV
         << *nuPz / 1000. << " "         // NuWro uses MeV, we want GeV
         << 0. << " "
         << *rx << " "                       // NuWro uses cm, ok
         << *ry << " "                       // NuWro uses cm, ok
         << *rz << " "                       // NuWro uses cm, ok
         << nu_time << " "
         << std::endl;


      // Save the neutron
      float neutronMass = 939.565;
      float neutronEsq = neutronMass*neutronMass + *neutronPx * *neutronPx + *neutronPy * *neutronPy + *neutronPz * *neutronPz;

      _hepmc_file << 0 << " "
         << 2112 << " "
         << -1 << " " << 0 << " " 
         << 0 << " "                         // Should be 1st and 2nd daugther, but putting ccnc and mode
         << 0 << " "                         // Should be 1st and 2nd daugther, but putting ccnc and mode
         << *neutronPx / 1000. << " "       // NuWro uses MeV, we want GeV
         << *neutronPy / 1000. << " "       // NuWro uses MeV, we want GeV
         << *neutronPz / 1000. << " "       // NuWro uses MeV, we want GeV
         << sqrt(neutronEsq) / 1000. << " "         // NuWro uses MeV, we want GeV
         << neutronMass / 1000. << " "
         << *rx << " "                       // NuWro uses cm, ok
         << *ry << " "                       // NuWro uses cm, ok
         << *rz << " "                       // NuWro uses cm, ok
         << nu_time << " "
         << std::endl;

      // Save the final state particles
      for (size_t j = 0; j < *nParticles; j++) {

         float EkinTemp = Ekin[j];
         float pTempSq = px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j];
         float Etemp = (pTempSq + EkinTemp*EkinTemp)/(2*EkinTemp);
         float mTemp = Etemp - EkinTemp;
         if(pdg[j] == 1002 || pdg[j] == 1003 || pdg[j] == 2003 || pdg[j] == 2004){
            pdg[j] = PDGmap[pdg[j]];
         }

         _hepmc_file << 1 << " "
            << pdg[j] << " "
            << 0 << " " << 0 << " " << 0 << " " << 0 << " "
            << px[j] / 1000. << " " // INCL++ uses MeV, we want GeV
            << py[j] / 1000. << " " // INCL++ uses MeV, we want GeV
            << pz[j] / 1000. << " " // INCL++ uses MeV, we want GeV
            << Etemp / 1000. << " "   // INCL++ uses MeV, we want GeV
            << mTemp / 1000. << " " // INCL++ uses MeV, we want GeV
            << *rx << " "                   // NuWro uses cm, ok
            << *ry << " "                   // NuWro uses cm, ok
            << *rz << " "                   // NuWro uses cm, ok
            << nu_time << " "
            << std::endl;
      }

      event_count++;

   }

   _hepmc_file.close();

   std::cout << event_count  << " events converted to hepmc" << std::endl;

   std::ofstream output("eventcount.meta");
   output << event_count;
   output.close();

}
