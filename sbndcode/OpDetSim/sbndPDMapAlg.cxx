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
  PDmap[24] = "dummy";
  PDmap[25] = "dummy";
  PDmap[26] = "barepmt";
  PDmap[27] = "barepmt";
  PDmap[28] = "dummy";
  PDmap[29] = "dummy";
  PDmap[30] = "barepmt";
  PDmap[31] = "barepmt";
  PDmap[32] = "dummy";
  PDmap[33] = "dummy";
  PDmap[34] = "barepmt";
  PDmap[35] = "barepmt";
  PDmap[36] = "bar";
  PDmap[37] = "bar";
  PDmap[38] = "bar";
  PDmap[39] = "bar";
  PDmap[40] = "bar";
  PDmap[41] = "bar";
  PDmap[42] = "pmt";
  PDmap[43] = "pmt";
  PDmap[44] = "pmt";
  PDmap[45] = "pmt";
  PDmap[46] = "pmt";
  PDmap[47] = "pmt";
  PDmap[48] = "pmt";
  PDmap[49] = "pmt";
  PDmap[50] = "pmt";
  PDmap[51] = "pmt";
  PDmap[52] = "pmt";
  PDmap[53] = "pmt";
  PDmap[54] = "bar";
  PDmap[55] = "bar";
  PDmap[56] = "bar";
  PDmap[57] = "bar";
  PDmap[58] = "bar";
  PDmap[59] = "bar";
  PDmap[60] = "bar";
  PDmap[61] = "bar";
  PDmap[62] = "bar";
  PDmap[63] = "bar";
  PDmap[64] = "bar";
  PDmap[65] = "bar";
  PDmap[66] = "pmt";
  PDmap[67] = "pmt";
  PDmap[68] = "pmt";
  PDmap[69] = "pmt";
  PDmap[70] = "pmt";
  PDmap[71] = "pmt";
  PDmap[72] = "arapucaT1";
  PDmap[73] = "arapucaT1";
  PDmap[74] = "arapucaT2";
  PDmap[75] = "arapucaT2";
  PDmap[76] = "pmt";
  PDmap[77] = "pmt";
  PDmap[78] = "pmt";
  PDmap[79] = "pmt";
  PDmap[80] = "pmt";
  PDmap[81] = "pmt";
  PDmap[82] = "bar";
  PDmap[83] = "bar";
  PDmap[84] = "bar";
  PDmap[85] = "bar";
  PDmap[86] = "bar";
  PDmap[87] = "bar";
  PDmap[88] = "dummy";
  PDmap[89] = "dummy";
  PDmap[90] = "barepmt";
  PDmap[91] = "barepmt";
  PDmap[92] = "dummy";
  PDmap[93] = "dummy";
  PDmap[94] = "barepmt";
  PDmap[95] = "barepmt";
  PDmap[96] = "dummy";
  PDmap[97] = "dummy";
  PDmap[98] = "barepmt";
  PDmap[99] = "barepmt";
  PDmap[100] = "bar";
  PDmap[101] = "bar";
  PDmap[102] = "bar";
  PDmap[103] = "bar";
  PDmap[104] = "bar";
  PDmap[105] = "bar";
  PDmap[106] = "pmt";
  PDmap[107] = "pmt";
  PDmap[108] = "arapucaT1";
  PDmap[109] = "arapucaT1";
  PDmap[110] = "arapucaT2";
  PDmap[111] = "arapucaT2";
  PDmap[112] = "pmt";
  PDmap[113] = "pmt";
  PDmap[114] = "pmt";
  PDmap[115] = "pmt";
  PDmap[116] = "arapucaT1";
  PDmap[117] = "arapucaT1";
  PDmap[118] = "arapucaT2";
  PDmap[119] = "arapucaT2";
  PDmap[120] = "pmt";
  PDmap[121] = "pmt";
  PDmap[122] = "pmt";
  PDmap[123] = "pmt";
  PDmap[124] = "arapucaT1";
  PDmap[125] = "arapucaT1";
  PDmap[126] = "arapucaT2";
  PDmap[127] = "arapucaT2";
  PDmap[128] = "pmt";
  PDmap[129] = "pmt";
  PDmap[130] = "bar";
  PDmap[131] = "bar";
  PDmap[132] = "bar";
  PDmap[133] = "bar";
  PDmap[134] = "bar";
  PDmap[135] = "bar";
  PDmap[136] = "bar";
  PDmap[137] = "bar";
  PDmap[138] = "bar";
  PDmap[139] = "bar";
  PDmap[140] = "bar";
  PDmap[141] = "bar";
  PDmap[142] = "pmt";
  PDmap[143] = "pmt";
  PDmap[144] = "arapucaT1";
  PDmap[145] = "arapucaT1";
  PDmap[146] = "arapucaT2";
  PDmap[147] = "arapucaT2";
  PDmap[148] = "pmt";
  PDmap[149] = "pmt";
  PDmap[150] = "pmt";
  PDmap[151] = "pmt";
  PDmap[152] = "arapucaT1";
  PDmap[153] = "arapucaT1";
  PDmap[154] = "arapucaT2";
  PDmap[155] = "arapucaT2";
  PDmap[156] = "pmt";
  PDmap[157] = "pmt";
  PDmap[158] = "pmt";
  PDmap[159] = "pmt";
  PDmap[160] = "arapucaT1";
  PDmap[161] = "arapucaT1";
  PDmap[162] = "arapucaT2";
  PDmap[163] = "arapucaT2";
  PDmap[164] = "pmt";
  PDmap[165] = "pmt";
  PDmap[166] = "bar";
  PDmap[167] = "bar";
  PDmap[168] = "bar";
  PDmap[169] = "bar";
  PDmap[170] = "bar";
  PDmap[171] = "bar";
  PDmap[172] = "dummy";
  PDmap[173] = "dummy";
  PDmap[174] = "barepmt";
  PDmap[175] = "barepmt";
  PDmap[176] = "dummy";
  PDmap[177] = "dummy";
  PDmap[178] = "barepmt";
  PDmap[179] = "barepmt";
  PDmap[180] = "dummy";
  PDmap[181] = "dummy";
  PDmap[182] = "barepmt";
  PDmap[183] = "barepmt";
  PDmap[184] = "bar";
  PDmap[185] = "bar";
  PDmap[186] = "bar";
  PDmap[187] = "bar";
  PDmap[188] = "bar";
  PDmap[189] = "bar";
  PDmap[190] = "pmt";
  PDmap[191] = "pmt";
  PDmap[192] = "pmt";
  PDmap[193] = "pmt";
  PDmap[194] = "pmt";
  PDmap[195] = "pmt";
  PDmap[196] = "arapucaT1";
  PDmap[197] = "arapucaT1";
  PDmap[198] = "arapucaT2";
  PDmap[199] = "arapucaT2";
  PDmap[200] = "pmt";
  PDmap[201] = "pmt";
  PDmap[202] = "pmt";
  PDmap[203] = "pmt";
  PDmap[204] = "pmt";
  PDmap[205] = "pmt";
  PDmap[206] = "bar";
  PDmap[207] = "bar";
  PDmap[208] = "bar";
  PDmap[209] = "bar";
  PDmap[210] = "bar";
  PDmap[211] = "bar";
  PDmap[212] = "bar";
  PDmap[213] = "bar";
  PDmap[214] = "bar";
  PDmap[215] = "bar";
  PDmap[216] = "bar";
  PDmap[217] = "bar";
  PDmap[218] = "pmt";
  PDmap[219] = "pmt";
  PDmap[220] = "pmt";
  PDmap[221] = "pmt";
  PDmap[222] = "pmt";
  PDmap[223] = "pmt";
  PDmap[224] = "pmt";
  PDmap[225] = "pmt";
  PDmap[226] = "pmt";
  PDmap[227] = "pmt";
  PDmap[228] = "pmt";
  PDmap[229] = "pmt";
  PDmap[230] = "bar";
  PDmap[231] = "bar";
  PDmap[232] = "bar";
  PDmap[233] = "bar";
  PDmap[234] = "bar";
  PDmap[235] = "bar";
  PDmap[236] = "dummy";
  PDmap[237] = "dummy";
  PDmap[238] = "barepmt";
  PDmap[239] = "barepmt";
  PDmap[240] = "dummy";
  PDmap[241] = "dummy";
  PDmap[242] = "barepmt";
  PDmap[243] = "barepmt";
  PDmap[244] = "dummy";
  PDmap[245] = "dummy";
  PDmap[246] = "barepmt";
  PDmap[247] = "barepmt";
  PDmap[248] = "bar";
  PDmap[249] = "bar";
  PDmap[250] = "bar";
  PDmap[251] = "bar";
  PDmap[252] = "bar";
  PDmap[253] = "bar";
  PDmap[254] = "pmt";
  PDmap[255] = "pmt";
  PDmap[256] = "pmt";
  PDmap[257] = "pmt";
  PDmap[258] = "pmt";
  PDmap[259] = "pmt";
  PDmap[260] = "pmt";
  PDmap[261] = "pmt";
  PDmap[262] = "pmt";
  PDmap[263] = "pmt";
  PDmap[264] = "pmt";
  PDmap[265] = "pmt";
  PDmap[266] = "bar";
  PDmap[267] = "bar";
  PDmap[268] = "bar";
  PDmap[269] = "bar";
  PDmap[270] = "bar";
  PDmap[271] = "bar";
 
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
    if(ch<272) return PDmap[ch];
    return "There is no such channel";
  }
}

#endif
