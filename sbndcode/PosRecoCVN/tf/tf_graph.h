////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       Graph
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
////              P.Plonski,                      from DUNE, WUT, Sept. 2017
////              T.Cai (tejinc@yorku.ca)         from DUNE, YorkU, March 2022
////              B.N.Nayak (nayakb@uci.edu)      from DUNE, BNL, November 2022
////
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Graph_h
#define Graph_h

#include <memory>
#include <string>
#include <vector>

namespace tensorflow {
  class Session;
  class Tensor;
  struct SavedModelBundle;
}

namespace tf {

  class Graph {
  public:
    static std::unique_ptr<Graph> create(const char* graph_file_name,
                                         const std::vector<std::string>& inputs = {},
                                         const std::vector<std::string>& outputs = {},
                                         bool use_bundle = false,
                                         int ninputs = 1,
                                         int noutputs = 1)
    {
      bool success;
      std::unique_ptr<Graph> ptr(
        new Graph(graph_file_name, inputs, outputs, success, use_bundle, ninputs, noutputs));
      if (success) { return ptr; }
      else {
        return nullptr;
      }
    }

    ~Graph();

    std::vector<std::vector<float>> run(const std::vector<std::vector<float>>& x);

    // process vector of 3D inputs, return vector of 1D outputs; use all inputs
    // if samples = -1, or only the specified number of first samples;
    // can deal with multiple inputs
    std::vector<std::vector<std::vector<float>>> run(
      const std::vector<std::vector<std::vector<std::vector<float>>>>& x,
      long long int samples = -1);
    std::vector<std::vector<std::vector<float>>> run(const std::vector<tensorflow::Tensor>& x);
    std::vector<std::vector<float>> runx(const std::vector<tensorflow::Tensor>& x);
    std::vector<std::vector<float>> runae(const std::vector<tensorflow::Tensor>& x);

  private:
    int n_inputs;
    int n_outputs;
    /// Not-throwing constructor.
    Graph(const char* graph_file_name,
          const std::vector<std::string>& inputs,
          const std::vector<std::string>& outputs,
          bool& success,
          bool use_bundle = false,
          int ninputs = 1,
          int noutputs = 1);

    tensorflow::Session* fSession;
    bool fUseBundle;
    tensorflow::SavedModelBundle* fBundle;
    std::vector<std::string> fInputNames;
    std::vector<std::string> fOutputNames;
  };

} // namespace tf

#endif
