#include "sbndPDMapAlg.h"

#ifndef SBNDPDMAPALG_CXX
#define SBNDPDMAPALG_CXX
 
//------------------------------------------------------------------------------
//--- opdet::sbndPDMapAlg implementation
//------------------------------------------------------------------------------

namespace opdet{

  sbndPDMapAlg::sbndPDMapAlg()
  {
    // Inserting data in std::map
  PDmap[0] = "bar";
  PDmap[1] = "bar";
  PDmap[2] = "bar";
  PDmap[3] = "bar";
  PDmap[4] = "bar";
  PDmap[5] = "bar";
  PDmap[6] = "pmt";
  PDmap[7] = "pmt";
  PDmap[8] = "pmt";
  PDmap[9] = "pmt";
  PDmap[10] = "pmt";
  PDmap[11] = "pmt";
  PDmap[12] = "pmt";
  PDmap[13] = "pmt";
  PDmap[14] = "pmt";
  PDmap[15] = "pmt";
  PDmap[16] = "pmt";
  PDmap[17] = "pmt";
  PDmap[18] = "bar";
  PDmap[19] = "bar";
  PDmap[20] = "bar";
  PDmap[21] = "bar";
  PDmap[22] = "bar";
  PDmap[23] = "bar";
  PDmap[24] = "barepmt";
  PDmap[25] = "barepmt";
  PDmap[26] = "barepmt";
  PDmap[27] = "barepmt";
  PDmap[28] = "barepmt";
  PDmap[29] = "barepmt";
  PDmap[30] = "bar";
  PDmap[31] = "bar";
  PDmap[32] = "bar";
  PDmap[33] = "bar";
  PDmap[34] = "bar";
  PDmap[35] = "bar";
  PDmap[36] = "pmt";
  PDmap[37] = "pmt";
  PDmap[38] = "pmt";
  PDmap[39] = "pmt";
  PDmap[40] = "pmt";
  PDmap[41] = "pmt";
  PDmap[42] = "pmt";
  PDmap[43] = "pmt";
  PDmap[44] = "pmt";
  PDmap[45] = "pmt";
  PDmap[46] = "pmt";
  PDmap[47] = "pmt";
  PDmap[48] = "bar";
  PDmap[49] = "bar";
  PDmap[50] = "bar";
  PDmap[51] = "bar";
  PDmap[52] = "bar";
  PDmap[53] = "bar";
  PDmap[54] = "bar";
  PDmap[55] = "bar";
  PDmap[56] = "bar";
  PDmap[57] = "bar";
  PDmap[58] = "bar";
  PDmap[59] = "bar";
  PDmap[60] = "pmt";
  PDmap[61] = "pmt";
  PDmap[62] = "pmt";
  PDmap[63] = "pmt";
  PDmap[64] = "pmt";
  PDmap[65] = "pmt";
  PDmap[66] = "xarapucaprime";
  PDmap[67] = "xarapucaprime";
  PDmap[68] = "xarapuca";
  PDmap[69] = "xarapuca";
  PDmap[70] = "pmt";
  PDmap[71] = "pmt";
  PDmap[72] = "pmt";
  PDmap[73] = "pmt";
  PDmap[74] = "pmt";
  PDmap[75] = "pmt";
  PDmap[76] = "bar";
  PDmap[77] = "bar";
  PDmap[78] = "bar";
  PDmap[79] = "bar";
  PDmap[80] = "bar";
  PDmap[81] = "bar";
  PDmap[82] = "barepmt";
  PDmap[83] = "barepmt";
  PDmap[84] = "barepmt";
  PDmap[85] = "barepmt";
  PDmap[86] = "barepmt";
  PDmap[87] = "barepmt";
  PDmap[88] = "bar";
  PDmap[89] = "bar";
  PDmap[90] = "bar";
  PDmap[91] = "bar";
  PDmap[92] = "bar";
  PDmap[93] = "bar";
  PDmap[94] = "pmt";
  PDmap[95] = "pmt";
  PDmap[96] = "xarapucaprime";
  PDmap[97] = "xarapucaprime";
  PDmap[98] = "xarapuca";
  PDmap[99] = "xarapuca";
  PDmap[100] = "pmt";
  PDmap[101] = "pmt";
  PDmap[102] = "pmt";
  PDmap[103] = "pmt";
  PDmap[104] = "xarapucaprime";
  PDmap[105] = "xarapucaprime";
  PDmap[106] = "xarapuca";
  PDmap[107] = "xarapuca";
  PDmap[108] = "pmt";
  PDmap[109] = "pmt";
  PDmap[110] = "pmt";
  PDmap[111] = "pmt";
  PDmap[112] = "xarapucaprime";
  PDmap[113] = "xarapucaprime";
  PDmap[114] = "xarapuca";
  PDmap[115] = "xarapuca";
  PDmap[116] = "pmt";
  PDmap[117] = "pmt";
  PDmap[118] = "bar";
  PDmap[119] = "bar";
  PDmap[120] = "bar";
  PDmap[121] = "bar";
  PDmap[122] = "bar";
  PDmap[123] = "bar";
  PDmap[124] = "bar";
  PDmap[125] = "bar";
  PDmap[126] = "bar";
  PDmap[127] = "bar";
  PDmap[128] = "bar";
  PDmap[129] = "bar";
  PDmap[130] = "pmt";
  PDmap[131] = "pmt";
  PDmap[132] = "arapucaT1";
  PDmap[133] = "arapucaT1";
  PDmap[134] = "arapucaT2";
  PDmap[135] = "arapucaT2";
  PDmap[136] = "pmt";
  PDmap[137] = "pmt";
  PDmap[138] = "pmt";
  PDmap[139] = "pmt";
  PDmap[140] = "arapucaT1";
  PDmap[141] = "arapucaT1";
  PDmap[142] = "arapucaT2";
  PDmap[143] = "arapucaT2";
  PDmap[144] = "pmt";
  PDmap[145] = "pmt";
  PDmap[146] = "pmt";
  PDmap[147] = "pmt";
  PDmap[148] = "arapucaT1";
  PDmap[149] = "arapucaT1";
  PDmap[150] = "arapucaT2";
  PDmap[151] = "arapucaT2";
  PDmap[152] = "pmt";
  PDmap[153] = "pmt";
  PDmap[154] = "bar";
  PDmap[155] = "bar";
  PDmap[156] = "bar";
  PDmap[157] = "bar";
  PDmap[158] = "bar";
  PDmap[159] = "bar";
  PDmap[160] = "barepmt";
  PDmap[161] = "barepmt";
  PDmap[162] = "barepmt";
  PDmap[163] = "barepmt";
  PDmap[164] = "barepmt";
  PDmap[165] = "barepmt";
  PDmap[166] = "bar";
  PDmap[167] = "bar";
  PDmap[168] = "bar";
  PDmap[169] = "bar";
  PDmap[170] = "bar";
  PDmap[171] = "bar";
  PDmap[172] = "pmt";
  PDmap[173] = "pmt";
  PDmap[174] = "pmt";
  PDmap[175] = "pmt";
  PDmap[176] = "pmt";
  PDmap[177] = "pmt";
  PDmap[178] = "arapucaT1";
  PDmap[179] = "arapucaT1";
  PDmap[180] = "arapucaT2";
  PDmap[181] = "arapucaT2";
  PDmap[182] = "pmt";
  PDmap[183] = "pmt";
  PDmap[184] = "pmt";
  PDmap[185] = "pmt";
  PDmap[186] = "pmt";
  PDmap[187] = "pmt";
  PDmap[188] = "bar";
  PDmap[189] = "bar";
  PDmap[190] = "bar";
  PDmap[191] = "bar";
  PDmap[192] = "bar";
  PDmap[193] = "bar";
  PDmap[194] = "bar";
  PDmap[195] = "bar";
  PDmap[196] = "bar";
  PDmap[197] = "bar";
  PDmap[198] = "bar";
  PDmap[199] = "bar";
  PDmap[200] = "pmt";
  PDmap[201] = "pmt";
  PDmap[202] = "pmt";
  PDmap[203] = "pmt";
  PDmap[204] = "pmt";
  PDmap[205] = "pmt";
  PDmap[206] = "pmt";
  PDmap[207] = "pmt";
  PDmap[208] = "pmt";
  PDmap[209] = "pmt";
  PDmap[210] = "pmt";
  PDmap[211] = "pmt";
  PDmap[212] = "bar";
  PDmap[213] = "bar";
  PDmap[214] = "bar";
  PDmap[215] = "bar";
  PDmap[216] = "bar";
  PDmap[217] = "bar";
  PDmap[218] = "barepmt";
  PDmap[219] = "barepmt";
  PDmap[220] = "barepmt";
  PDmap[221] = "barepmt";
  PDmap[222] = "barepmt";
  PDmap[223] = "barepmt";
  PDmap[224] = "bar";
  PDmap[225] = "bar";
  PDmap[226] = "bar";
  PDmap[227] = "bar";
  PDmap[228] = "bar";
  PDmap[229] = "bar";
  PDmap[230] = "pmt";
  PDmap[231] = "pmt";
  PDmap[232] = "pmt";
  PDmap[233] = "pmt";
  PDmap[234] = "pmt";
  PDmap[235] = "pmt";
  PDmap[236] = "pmt";
  PDmap[237] = "pmt";
  PDmap[238] = "pmt";
  PDmap[239] = "pmt";
  PDmap[240] = "pmt";
  PDmap[241] = "pmt";
  PDmap[242] = "bar";
  PDmap[243] = "bar";
  PDmap[244] = "bar";
  PDmap[245] = "bar";
  PDmap[246] = "bar";
  PDmap[247] = "bar";
 
  }

  sbndPDMapAlg::~sbndPDMapAlg()
  { }

  bool sbndPDMapAlg::pdType(int ch, std::string pdname)
  {
    if(PDmap[ch]==pdname) return true;
    return false;
  }

  std::string sbndPDMapAlg::pdName(int ch)
  {
    if(ch<(int)PDmap.size()) return PDmap[ch];
    return "There is no such channel";
  }

  int sbndPDMapAlg::size()
  {
    return (int)PDmap.size();
  }
}

#endif
