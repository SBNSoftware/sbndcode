#include "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroBv2/integrated_flux_ray_trace/FluxMap.h"

void IntegratedFluxCalcSD()
{
  const double nom = 1.65974e+13;

  double total_frac_error = 0., total_frac_bias = 0.;
  
  for(auto const& [ name, map ] : univsIntegratedFluxMap)
    {
      std::cout << name << std::endl;

      double total = 0.;

      for(auto const& [ id, value ] : map)
	{
	  total += value;
	}

      const double mean = total / map.size();
      double sdtotal = 0.;

      for(auto const& [ id, value ] : map)
	{
	  sdtotal += (value - mean) * (value - mean);
	}

      const double sd = std::sqrt(sdtotal / (map.size() - 1));
      const double frac_error = sd / mean;

      const double bias = std::abs(mean - nom);
      const double frac_bias = bias / nom;
      
      std::cout << "Mean: " << mean << " SD: " << sd << " ("
		<< frac_error * 100 << "%, " << frac_bias * 100 << "%)" << std::endl;

      if(name != "flux_weights_all")
	{
	  total_frac_error += std::pow(frac_error, 2);
	  total_frac_bias += std::pow(frac_bias, 2);
	}
    }
  std::cout << '\n'
	    << "Total Error: " << std::sqrt(total_frac_error) * 100 << "%, " << std::sqrt(total_frac_bias) * 100 << "%" << std::endl;
}
