////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       Graph
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//              T.Cai (tejinc@yorku.ca)         from DUNE, YorkU, March 2022
//              B.N.Nayak (nayakb@uci.edu)      from DUNE, BNL, November 2022
//
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tf_graph.h"

#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/platform/env.h"

#include "tensorflow/core/public/session_options.h"

// -------------------------------------------------------------------
tf::Graph::Graph(const char* graph_file_name,
                 const std::vector<std::string>& inputs,
                 const std::vector<std::string>& outputs,
                 bool& success,
                 bool use_bundle,
                 int ninputs,
                 int noutputs)
{
  fUseBundle = use_bundle;
  success = false; // until all is done correctly

  n_inputs = ninputs;
  n_outputs = noutputs;

  // Force tf to only use a single core so it doesn't eat batch farms
  tensorflow::SessionOptions options;
  tensorflow::ConfigProto& config = options.config;
  config.set_inter_op_parallelism_threads(1);
  config.set_intra_op_parallelism_threads(1);
  config.set_use_per_session_threads(false);

  auto status = tensorflow::NewSession(options, &fSession);
  if (!status.ok()) {
    std::cout << status.ToString() << std::endl;
    return;
  }

  if (fUseBundle) {

    fBundle = new tensorflow::SavedModelBundle();
    status = tensorflow::LoadSavedModel(tensorflow::SessionOptions(),
                                        tensorflow::RunOptions(),
                                        graph_file_name,
                                        {tensorflow::kSavedModelTagServe},
                                        fBundle);
    std::cout << "tf_graph loaded SavedModelBundle with status: " << status.ToString() << std::endl;
    if (!status.ok()) return;

    auto sig_map = fBundle->meta_graph_def.signature_def();
    std::string sig_def = "serving_default";
    bool has_default_key = false;
    std::vector<std::string> sig_map_keys;
    for (auto const& p : sig_map) {
      if (p.first == sig_def) has_default_key = true;
      sig_map_keys.push_back(p.first);
    }
    auto model_def = sig_map.at((has_default_key) ? sig_def : sig_map_keys.back());

    // ... Get the input names
    if (inputs.empty()){
      std::cout << "tf_graph using all inputs:" << std::endl;
      for (auto const& p : model_def.inputs()) {
        fInputNames.push_back(p.second.name());
        std::cout << "tf_graph InputName: " << fInputNames.back() << std::endl;
        std::cout << "  key: " << p.first << " value: " << p.second.name() << std::endl;
      }
    }
    else{
      std::cout << "tf_graph using selected inputs:" << std::endl;
      for (const auto& s : inputs) {
        for (auto const& p : model_def.inputs()) {
          if (p.first == s) {
            fInputNames.push_back(p.second.name());
            std::cout << "  key: " << p.first << " value: " << p.second.name() << std::endl;
          }
        }
      }
    }

    // ... Get the output names
    //  .. get all outputs if no specific name provided
    if (outputs.empty()) {
      std::cout << "tf_graph using all outputs:" << std::endl;
      for (auto const& p : model_def.outputs()) {
        fOutputNames.push_back(p.second.name());
        std::cout << "  key: " << p.first << " value: " << p.second.name() << std::endl;
      }
    }
    //  .. or use only the outputs whose keys are specified
    else {
      std::cout << "tf_graph using selected outputs:" << std::endl;
      for (const auto& s : outputs) {
        for (auto const& p : model_def.outputs()) {
          if (p.first == s) {
            fOutputNames.push_back(p.second.name());
            std::cout << "  key: " << p.first << " value: " << p.second.name() << std::endl;
          }
        }
      }
    }
    if (fOutputNames.empty()) {
      std::cout << "tf_graph did not find outputs in SaveModelBundle." << std::endl;
      return;
    }
  }
  else {

    tensorflow::GraphDef graph_def;
    status = tensorflow::ReadBinaryProto(tensorflow::Env::Default(), graph_file_name, &graph_def);
    std::cout << "tf_graph loaded ProtoBuf graph with status: " << status.ToString() << std::endl;
    if (!status.ok()) return;

    size_t ng = graph_def.node().size();
    for (int i = 0; i < n_inputs; ++i) {
      fInputNames.push_back(graph_def.node()[i].name());
    }

    // last node as output if no specific name provided
    if (outputs.empty()) {
      for (int i = n_outputs; i > 0; --i) {
        fOutputNames.push_back(graph_def.node()[ng - i].name());
      }
    }
    else // or last nodes with names containing provided strings
    {
      std::string last, current, basename, name;
      for (size_t n = 0; n < ng; ++n) {
        name = graph_def.node()[n].name();
        auto pos = name.find("/");
        if (pos != std::string::npos) { basename = name.substr(0, pos); }
        else {
          continue;
        }

        bool found = false;
        for (const auto& s : outputs) {
          if (name.find(s) != std::string::npos) {
            found = true;
            break;
          }
        }
        if (found) {
          if (!last.empty() && (basename != current)) { fOutputNames.push_back(last); }
          current = basename;
          last = name;
        }
      }
      if (!last.empty()) { fOutputNames.push_back(last); }
    }
    if (fOutputNames.empty()) {
      std::cout << "Output nodes not found in the graph." << std::endl;
      return;
    }
    status = fSession->Create(graph_def);
    if (!status.ok()) {
      std::cout << status.ToString() << std::endl;
      return;
    }
  }

  success = true; // ok, graph loaded from the file
}

tf::Graph::~Graph()
{
  fSession->Close().IgnoreError();
  delete fSession;
  if (fUseBundle) { delete fBundle; }
}
// -------------------------------------------------------------------

std::vector<std::vector<float>> tf::Graph::run(const std::vector<std::vector<float>>& x)
{
  if (x.empty() || x.front().empty()) { return std::vector<std::vector<float>>(); }

  long long int rows = x.size(), cols = x.front().size();

  std::vector<tensorflow::Tensor> _x;
  _x.push_back(
    tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, rows, cols, 1})));
  auto input_map = _x[0].tensor<float, 4>();

  for (long long int r = 0; r < rows; ++r) {
    const auto& row = x[r];
    for (long long int c = 0; c < cols; ++c) {
      input_map(0, r, c, 0) = row[c];
    }
  }

  auto result = run(_x);
  if (!result.empty()) { return result.front(); }
  else {
    return std::vector<std::vector<float>>();
  }
}
// -------------------------------------------------------------------

std::vector<std::vector<std::vector<float>>> tf::Graph::run(
  const std::vector<std::vector<std::vector<std::vector<float>>>>& x,
  long long int samples)
{
  if ((samples == 0) || x.empty() || x.front().empty() || x.front().front().empty() ||
      x.front().front().front().empty())
    return std::vector<std::vector<std::vector<float>>>();

  if ((samples == -1) || (samples > (long long int)x.size())) { samples = x.size(); }

  long long int rows = x.front().size(), cols = x.front().front().size(),
                depth = x.front().front().front().size();

  std::vector<tensorflow::Tensor> _x;

  // Single-input network
  if (n_inputs == 1) {
    _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT,
                                    tensorflow::TensorShape({samples, rows, cols, depth})));
    auto input_map = _x[0].tensor<float, 4>();
    for (long long int s = 0; s < samples; ++s) {
      const auto& sample = x[s];
      for (long long int r = 0; r < rows; ++r) {
        const auto& row = sample[r];
        for (long long int c = 0; c < cols; ++c) {
          const auto& col = row[c];
          for (long long int d = 0; d < depth; ++d) {
            input_map(s, r, c, d) = col[d];
          }
        }
      }
    }
  }
  // Multi-input network
  else {
    for (int i = 0; i < depth; ++i) {
      _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT,
                                      tensorflow::TensorShape({samples, rows, cols, 1})));
    }

    for (int view = 0; view < depth; ++view) {
      auto input_map = _x[view].tensor<float, 4>();
      for (long long int s = 0; s < samples; ++s) {
        const auto& sample = x[s];
        for (long long int r = 0; r < rows; ++r) {
          const auto& row = sample[r];
          for (long long int c = 0; c < cols; ++c) {
            const auto& col = row[c];
            long long int d = view;
            input_map(s, r, c, 0) = col[d];
          }
        }
      }
    }
  }

  return run(_x);
}

// -------------------------------------------------------------------

std::vector<std::vector<std::vector<float>>> tf::Graph::run(
  const std::vector<tensorflow::Tensor>& x)
{
  std::vector<std::pair<std::string, tensorflow::Tensor>> inputs;
  for (int i = 0; i < n_inputs; ++i) {
    inputs.push_back({fInputNames[i], x[i]});
  }

  std::vector<tensorflow::Tensor> outputs;
  std::vector<std::string> outputNames;
  auto status = (fUseBundle) ?
                  fBundle->GetSession()->Run(inputs, fOutputNames, outputNames, &outputs) :
                  fSession->Run(inputs, fOutputNames, outputNames, &outputs);

  if (status.ok()) {
    size_t samples = 0;

    for (size_t o = 0; o < outputs.size(); ++o) {
      if (o == 0) { samples = outputs[o].dim_size(0); }
      else if ((int)samples != outputs[o].dim_size(0)) {
        throw std::string("TF outputs size inconsistent.");
      }
    }

    std::vector<std::vector<std::vector<float>>> result;
    result.resize(samples, std::vector<std::vector<float>>(outputs.size()));

    for (size_t s = 0; s < samples; ++s) {
      for (size_t o = 0; o < outputs.size(); ++o) {
        size_t n = outputs[o].dim_size(1);
        auto output_map = outputs[o].tensor<float, 2>();

        result[s][o].resize(outputs[o].dim_size(1));

        std::vector<float>& vs = result[s][o];
        for (size_t i = 0; i < n; ++i) {
          vs[i] = output_map(s, i);
        }
      }
    }

    return result;
  }
  else {
    std::cout << status.ToString() << std::endl;
    return std::vector<std::vector<std::vector<float>>>();
  }
}
// -------------------------------------------------------------------

std::vector<std::vector<float>> tf::Graph::runx(const std::vector<tensorflow::Tensor>& x)
{
  std::vector<std::pair<std::string, tensorflow::Tensor>> inputs;
  for (int i = 0; i < n_inputs; ++i) {
    inputs.push_back({fInputNames[i], x[i]});
  }

  std::vector<tensorflow::Tensor> outputs;
  std::vector<std::string> outputNames;
  auto status = (fUseBundle) ?
                  fBundle->GetSession()->Run(inputs, fOutputNames, outputNames, &outputs) :
                  fSession->Run(inputs, fOutputNames, outputNames, &outputs);

  if (status.ok()) {
    size_t samples = 0, nouts = 0;

    for (size_t o = 0; o < outputs.size(); ++o) {
      if (o == 0) { samples = outputs[o].dim_size(0); }
      else if ((int)samples != outputs[o].dim_size(0)) {
        throw std::string("TF outputs size inconsistent.");
      }
      nouts += outputs[o].dim_size(1);
    }

    std::vector<std::vector<float>> result;
    result.resize(samples, std::vector<float>(nouts));

    size_t idx0 = 0;
    for (size_t o = 0; o < outputs.size(); ++o) {
      auto output_map = outputs[o].tensor<float, 2>();

      size_t n = outputs[o].dim_size(1);
      for (size_t s = 0; s < samples; ++s) {
        std::vector<float>& vs = result[s];
        for (size_t i = 0; i < n; ++i) {
          vs[idx0 + i] = output_map(s, i);
        }
      }
      idx0 += n;
    }

    return result;
  }
  else {
    std::cout << status.ToString() << std::endl;
    return std::vector<std::vector<float>>();
  }
}
// -------------------------------------------------------------------

std::vector<std::vector<float>> tf::Graph::runae(const std::vector<tensorflow::Tensor>& x)
{
  std::vector<std::pair<std::string, tensorflow::Tensor>> inputs;
  for (int i = 0; i < n_inputs; ++i) {
    inputs.push_back({fInputNames[i], x[i]});
  }

  std::vector<tensorflow::Tensor> outputs;
  std::vector<std::string> outputNames;
  auto status = (fUseBundle) ?
                  fBundle->GetSession()->Run(inputs, fOutputNames, outputNames, &outputs) :
                  fSession->Run(inputs, fOutputNames, outputNames, &outputs);

  if (status.ok()) {
    size_t samples = 0, npoints = 0;

    if (outputs.size() > 1) { throw std::string("TF runae: detected more than one output."); }

    samples = outputs[0].dim_size(0);
    npoints = outputs[0].dim_size(1);

    std::vector<std::vector<float>> result;
    result.resize(samples, std::vector<float>(npoints));

    auto output_map = outputs[0].tensor<float, 3>();

    for (size_t s = 0; s < samples; ++s) {
      std::vector<float>& vs = result[s];
      for (size_t i = 0; i < npoints; ++i) {
        vs[i] = output_map(s, i, 0);
      }
    }

    return result;
  }
  else {
    std::cout << status.ToString() << std::endl;
    return std::vector<std::vector<float>>();
  }
}
