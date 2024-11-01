// WARNING: the SP C++ code has a lot of hard coded names for various
// filter components.  Until this is cleaned up, one MUST configure
// the filter competents with matching type names and not change their
// instance names.


local wc = import 'wirecell.jsonnet';

local lf(name, data={}) = {
  type: 'LfFilter',
  name: name,
  data: {
    max_freq: 1 * wc.megahertz,
    tau: 0.0 * wc.megahertz,
  } + data,
};
local hf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1 * wc.megahertz,
    sigma: 0.0 * wc.megahertz,
    power: 2,
    flag: true,
  } + data,
};
// All "wire" filters are Hf with different base values.
local wf(name, data={}) = {
  type: 'HfFilter',
  name: name,
  data: {
    max_freq: 1,  // warning: units
    power: 2,
    flag: false,
    sigma: 0.0,  // caller should provide
  } + data,
};

// Zeus take my eyes! Magic numbers are everywhere!
/**  
 *  Default SP parameters (till May 2019)
 */
// [
//   lf('ROI_tight_lf', { tau: 0.02 * wc.megahertz }),  // 0.02 -> 0.027
//   lf('ROI_tighter_lf', { tau: 0.1 * wc.megahertz }),  // 0.1 -> 0.075
//   lf('ROI_loose_lf', { tau: 0.0025 * wc.megahertz }),  // 0.0025 ->  0.004
// 
//   hf('Gaus_tight'),
//   hf('Gaus_wide', { sigma: 1.11408e-01 * wc.megahertz }),
//   hf('Wiener_tight_U', {
//     sigma: 5.75416e+01 / 800.0 * 2 * wc.megahertz,
//     power: 4.10358e+00,
//   }),
//   hf("Wiener_tight_V", { sigma: 5.99306e+01/800.0*2 * wc.megahertz,
//   	                   power: 4.20820e+00 }),
//   hf('Wiener_tight_W', {
//     sigma: 5.88802e+01 / 800.0 * 2 * wc.megahertz,
//     power: 4.17455e+00,
//   }),
//   hf('Wiener_wide_U', {
//     sigma: 1.78695e+01 / 200.0 * 2 * wc.megahertz,
//     power: 5.33129e+00,
//   }),
//   hf("Wiener_wide_V", { sigma: 1.84666e+01/200.0*2 * wc.megahertz,
//   	                  power: 5.60489e+00 }),
//   hf('Wiener_wide_W', {
//     sigma: 1.83044e+01 / 200.0 * 2 * wc.megahertz,
//     power: 5.44945e+00,
//   }),
// 
//   wf('Wire_ind', { sigma: 1.0 / wc.sqrtpi * 1.4 }),
//   wf('Wire_col', { sigma: 1.0 / wc.sqrtpi * 3.0 }),
// ]

[
  lf('ROI_tight_lf', { tau: 0.0185 * wc.megahertz }),  // 0.02 
  lf('ROI_tighter_lf', { tau: 0.145 * wc.megahertz }),  // 0.1 
  lf('ROI_loose_lf', { tau: 0.00175 * wc.megahertz }),  // 0.0025 

  hf('Gaus_tight'),
  hf('Gaus_wide', { sigma: 0.10 * wc.megahertz }), 


  hf('Wiener_tight_U', {
    sigma: 0.15  * wc.megahertz,
    power: 5.5,
  }),
  hf("Wiener_tight_V", {
    sigma: 0.15 * wc.megahertz,
    power: 5.0 }),
  hf('Wiener_tight_W', {
    sigma: 0.25 * wc.megahertz,
    power: 3.0,
  }),

  hf('Wiener_wide_U', {
    sigma: 0.186765  * wc.megahertz,
    power: 5.05429,
  }),
  hf("Wiener_wide_V", {
    sigma: 0.1936 * wc.megahertz,
    power: 5.77422,
  }),
  hf('Wiener_wide_W', {
    sigma: 0.175722  * wc.megahertz,
    power: 4.37928,
  }),

  wf('Wire_ind', { sigma: 1.0 / wc.sqrtpi * 1.05 }), 
  wf('Wire_col', { sigma: 1.0 / wc.sqrtpi * 3.60 }),
]
