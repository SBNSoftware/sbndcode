/**
 * @file    TriggerEmulationService.cc
 * @brief   Service implementation associated with TriggerEmulationService.h and TriggerEmulationService_service.cc
 * @author  Nikki Pallat (palla110@umn.edu)
 * @date    July 7, 2025
 * @see     TriggerEmulationService.h 
 * @ingroup TriggerEmulationService
 *
 */

#include "TriggerEmulationService.h"

//namespace lar {
namespace calib {

  TriggerEmulationService::TriggerEmulationService(fhicl::ParameterSet const& pset, art::ActivityRegistry& amp)
    : fMonWidth(pset.get<int>("MonWidth", 3)),
      fTotalCAENBoards(pset.get<int>("TotalCAENBoards", 10)),
      PMTPerBoard(pset.get<int>("PMTPerBoard", 16)),
      Baseline(pset.get<int>("Baseline", 14250)),
      fMC(pset.get<bool>("MC", false))
  {}

  TriggerEmulationService::~TriggerEmulationService()
  {} // deconstructor

  void TriggerEmulationService::ConstructMonPulse(
      std::vector<raw::OpDetWaveform> fWaveforms,
      int MonThreshold,
      std::vector<int> *MonPulse,
      bool Saving,
      int FlashCounter,
      int *numPairsOverThreshold
  )
  {
      // Loop over the entries in our waveform vector
      // We care about getting the pairing correct

      std::fill(MonPulse->begin(), MonPulse->end(), 0);    

      if (fMC) { // monte carlo
          if (fWaveforms.empty()) {
              std::cout << "Empty waveform vector. Exiting ConstructMonPulse." << std::endl;
              return;
          }

          std::map<int, const raw::OpDetWaveform*> channel_to_waveform;
          for (const auto& wvf : fWaveforms)
              channel_to_waveform[wvf.ChannelNumber()] = &wvf;

          std::vector<int> Pair2 = { 6,   8,  10,  12,  14,  16,  36,  38,  40,  60,  62,  66,  68, 70,  84,  86,  88,  90,  92,  94, 114, 116, 118, 138, 140, 144, 146, 148, 162, 164, 168, 170, 172, 192, 194, 196, 216, 218, 220, 222, 224, 226, 240, 242, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304};
          std::vector<int> Pair1 = { 7,   9,  11,  13,  15,  17,  37,  39,  41,  61,  63,  67,  69, 71,  85,  87,  89,  91,  93,  95, 115, 117, 119, 139, 141, 145, 147, 149, 163, 165, 169, 171, 173, 193, 195, 197, 217, 219, 221, 223, 225, 227, 241, 243, 247, 249, 251, 271, 273, 275, 295, 297, 299, 301, 303, 305};
          std::vector<int> Unpaired = {65,  64, 143, 142, 167, 166, 245, 244};
          std::set<int> used_channels;

          // resize
          int ReadoutSize = fWaveforms[0].size();
          MonPulse->assign(ReadoutSize, 0);

          int countPairs = 0;

          for (size_t i = 0; i < Pair1.size(); ++i) {
              int ch1 = Pair1[i];
              int ch2 = Pair2[i];

              // skip if either waveform is missing
              if (channel_to_waveform.count(ch1) == 0 || channel_to_waveform.count(ch2) == 0) continue;
              // skip if already processed
              if (used_channels.count(ch1) || used_channels.count(ch2)) continue;

              const auto& wvf1 = *channel_to_waveform[ch1];
              const auto& wvf2 = *channel_to_waveform[ch2];

              auto bin1 = ConstructBinaryResponse(wvf1, MonThreshold);
              auto bin2 = ConstructBinaryResponse(wvf2, MonThreshold);

              bool pairOverThreshold = false;

              for (int j = 0; j < ReadoutSize; ++j) {
                  if (bin1[j] || bin2[j]) {
                      (*MonPulse)[j]++;
                      pairOverThreshold = true;
                  }
              }

              if (pairOverThreshold) countPairs++;

              used_channels.insert(ch1);
              used_channels.insert(ch2);
          }

          for (int ch : Unpaired) { // Unpaired channels
              if (used_channels.count(ch)) continue;
              if (channel_to_waveform.count(ch) == 0) continue;

              const auto& wvf = *channel_to_waveform[ch];
              auto bin = ConstructBinaryResponse(wvf, MonThreshold);

              bool pairOverThreshold = false;

              for (int j = 0; j < ReadoutSize; ++j) {
                  if (bin[j]) {
                      (*MonPulse)[j]++;
                      pairOverThreshold = true;
                  }
              }

              if (pairOverThreshold) countPairs++;

          }

          if (numPairsOverThreshold) *numPairsOverThreshold = countPairs;

      } else { // data
          if (fWaveforms.empty()) {
              std::cout << "Empty waveform vector. Exiting ConstructMonPulse." << std::endl;
              return;
          }

          int NumFlash = fWaveforms.size() / (fTotalCAENBoards * PMTPerBoard);
          int FirstReadoutIndex = 0 + FlashCounter*PMTPerBoard + 0*PMTPerBoard*NumFlash;
          int ReadoutSize = fWaveforms[FirstReadoutIndex].size();
 
          for (int CurrentBoard = 0; CurrentBoard < fTotalCAENBoards; ++CurrentBoard) {
              int CAENChannel = 0;
              // Loop over each PMT in a board
              while (CAENChannel < PMTPerBoard) {

                  int ChannelStep = 1;
                  std::vector<bool> BinaryMonContrib(ReadoutSize);

                  if (CAENChannel != 14) {
                      int WaveIndex_Pair1 = CAENChannel + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NumFlash;
                      int WaveIndex_Pair2 = CAENChannel + 1 + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NumFlash;

                      ChannelStep = 2;
  
                      auto const& wvf_Pair1 = fWaveforms[WaveIndex_Pair1];
                      auto const& wvf_Pair2 = fWaveforms[WaveIndex_Pair2];

                      std::vector<bool> BinaryResponse_Pair1 = ConstructBinaryResponse(wvf_Pair1, MonThreshold);
                      std::vector<bool> BinaryResponse_Pair2 = ConstructBinaryResponse(wvf_Pair2, MonThreshold);

                      for (int i=0; i<ReadoutSize; i++) BinaryMonContrib[i] = (BinaryResponse_Pair1[i] || BinaryResponse_Pair2[i]);

                  } else {
                      int WaveIndex = CAENChannel + FlashCounter*PMTPerBoard + CurrentBoard*PMTPerBoard*NumFlash;
                      auto const& wvf_Unpaired = fWaveforms[WaveIndex];
                      std::vector<bool> BinaryResponse_Unpaired = ConstructBinaryResponse(wvf_Unpaired, MonThreshold);
                      for (int i=0; i<ReadoutSize; i++) BinaryMonContrib[i] = (BinaryResponse_Unpaired[i]);
                  }  
                  for (int i=0; i<ReadoutSize; i++)
                  {
                    if (BinaryMonContrib[i]) (*MonPulse)[i] = (*MonPulse)[i]+1;
                  }

                  CAENChannel=CAENChannel+ChannelStep;
              } //loop over channels
          } //loop over boards
      } // data  
  } // ConstructMonPulse

  std::vector<bool> TriggerEmulationService::ConstructBinaryResponse(const raw::OpDetWaveform &wvf, int MonThreshold) 
  {
    std::vector<bool> BinaryResponse(wvf.size());
    int WaveformIndex=1;
    while (WaveformIndex<int(wvf.size()))
    {
      bool CrossedThreshold = (((wvf[WaveformIndex-1]-Baseline)>-MonThreshold) && ((wvf[WaveformIndex]-Baseline)<-MonThreshold));
      if (CrossedThreshold)
      {
        int StartIndex = WaveformIndex;
        int EndIndex = WaveformIndex+4*fMonWidth;
        if (StartIndex>=int(wvf.size() )) StartIndex=wvf.size()-1; //should never happen
        if (EndIndex>=int(wvf.size())) EndIndex=wvf.size()-1; //May happen
        for (int k=StartIndex; k<EndIndex; k++)
        {
          BinaryResponse[k] = true;
        }
        WaveformIndex=EndIndex+1;
      }
      else WaveformIndex=WaveformIndex+1;
    }

    return BinaryResponse; 
  } // ConstructBinaryResponse

} // namespace calib



