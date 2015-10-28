/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "internal.h"

/*Constants for use with OD_DIVU_SMALL().
  See \cite{Rob05} for details on computing these constants.
  @INPROCEEDINGS{Rob05,
    author="Arch D. Robison",
    title="{N}-bit Unsigned Division via {N}-bit Multiply-Add",
    booktitle="Proc. of the 17th IEEE Symposium on Computer Arithmetic
     (ARITH'05)",
    pages="131--139",
    address="Cape Cod, MA",
    month=Jun,
    year=2005
  }*/
uint32_t OD_DIVU_SMALL_CONSTS[OD_DIVU_DMAX][2] = {
  {0xFFFFFFFF,0xFFFFFFFF}, {0xFFFFFFFF,0xFFFFFFFF}, {0xAAAAAAAB,         0},
  {0xFFFFFFFF,0xFFFFFFFF}, {0xCCCCCCCD,         0}, {0xAAAAAAAB,         0},
  {0x92492492,0x92492492}, {0xFFFFFFFF,0xFFFFFFFF}, {0xE38E38E4,         0},
  {0xCCCCCCCD,         0}, {0xBA2E8BA3,         0}, {0xAAAAAAAB,         0},
  {0x9D89D89E,         0}, {0x92492492,0x92492492}, {0x88888889,         0},
  {0xFFFFFFFF,0xFFFFFFFF}, {0xF0F0F0F1,         0}, {0xE38E38E4,         0},
  {0xD79435E5,0xD79435E5}, {0xCCCCCCCD,         0}, {0xC30C30C3,0xC30C30C3},
  {0xBA2E8BA3,         0}, {0xB21642C9,         0}, {0xAAAAAAAB,         0},
  {0xA3D70A3E,         0}, {0x9D89D89E,         0}, {0x97B425ED,0x97B425ED},
  {0x92492492,0x92492492}, {0x8D3DCB09,         0}, {0x88888889,         0},
  {0x84210842,0x84210842}, {0xFFFFFFFF,0xFFFFFFFF}, {0xF83E0F84,         0},
  {0xF0F0F0F1,         0}, {0xEA0EA0EA,0xEA0EA0EA}, {0xE38E38E4,         0},
  {0xDD67C8A6,0xDD67C8A6}, {0xD79435E5,0xD79435E5}, {0xD20D20D2,0xD20D20D2},
  {0xCCCCCCCD,         0}, {0xC7CE0C7D,         0}, {0xC30C30C3,0xC30C30C3},
  {0xBE82FA0C,         0}, {0xBA2E8BA3,         0}, {0xB60B60B6,0xB60B60B6},
  {0xB21642C9,         0}, {0xAE4C415D,         0}, {0xAAAAAAAB,         0},
  {0xA72F053A,         0}, {0xA3D70A3E,         0}, {0xA0A0A0A1,         0},
  {0x9D89D89E,         0}, {0x9A90E7D9,0x9A90E7D9}, {0x97B425ED,0x97B425ED},
  {0x94F2094F,0x94F2094F}, {0x92492492,0x92492492}, {0x8FB823EE,0x8FB823EE},
  {0x8D3DCB09,         0}, {0x8AD8F2FC,         0}, {0x88888889,         0},
  {0x864B8A7E,         0}, {0x84210842,0x84210842}, {0x82082082,0x82082082},
  {0xFFFFFFFF,0xFFFFFFFF}, {0xFC0FC0FD,         0}, {0xF83E0F84,         0},
  {0xF4898D60,         0}, {0xF0F0F0F1,         0}, {0xED7303B6,         0},
  {0xEA0EA0EA,0xEA0EA0EA}, {0xE6C2B449,         0}, {0xE38E38E4,         0},
  {0xE070381C,0xE070381C}, {0xDD67C8A6,0xDD67C8A6}, {0xDA740DA8,         0},
  {0xD79435E5,0xD79435E5}, {0xD4C77B04,         0}, {0xD20D20D2,0xD20D20D2},
  {0xCF6474A9,         0}, {0xCCCCCCCD,         0}, {0xCA4587E7,         0},
  {0xC7CE0C7D,         0}, {0xC565C87C,         0}, {0xC30C30C3,0xC30C30C3},
  {0xC0C0C0C1,         0}, {0xBE82FA0C,         0}, {0xBC52640C,         0},
  {0xBA2E8BA3,         0}, {0xB81702E1,         0}, {0xB60B60B6,0xB60B60B6},
  {0xB40B40B4,0xB40B40B4}, {0xB21642C9,         0}, {0xB02C0B03,         0},
  {0xAE4C415D,         0}, {0xAC769184,0xAC769184}, {0xAAAAAAAB,         0},
  {0xA8E83F57,0xA8E83F57}, {0xA72F053A,         0}, {0xA57EB503,         0},
  {0xA3D70A3E,         0}, {0xA237C32B,0xA237C32B}, {0xA0A0A0A1,         0},
  {0x9F1165E7,0x9F1165E7}, {0x9D89D89E,         0}, {0x9C09C09C,0x9C09C09C},
  {0x9A90E7D9,0x9A90E7D9}, {0x991F1A51,0x991F1A51}, {0x97B425ED,0x97B425ED},
  {0x964FDA6C,0x964FDA6C}, {0x94F2094F,0x94F2094F}, {0x939A85C4,0x939A85C4},
  {0x92492492,0x92492492}, {0x90FDBC09,0x90FDBC09}, {0x8FB823EE,0x8FB823EE},
  {0x8E78356D,0x8E78356D}, {0x8D3DCB09,         0}, {0x8C08C08C,0x8C08C08C},
  {0x8AD8F2FC,         0}, {0x89AE408A,         0}, {0x88888889,         0},
  {0x8767AB5F,0x8767AB5F}, {0x864B8A7E,         0}, {0x85340853,0x85340853},
  {0x84210842,0x84210842}, {0x83126E98,         0}, {0x82082082,0x82082082},
  {0x81020408,0x81020408}, {0xFFFFFFFF,0xFFFFFFFF}, {0xFE03F810,         0},
  {0xFC0FC0FD,         0}, {0xFA232CF3,         0}, {0xF83E0F84,         0},
  {0xF6603D99,         0}, {0xF4898D60,         0}, {0xF2B9D649,         0},
  {0xF0F0F0F1,         0}, {0xEF2EB720,         0}, {0xED7303B6,         0},
  {0xEBBDB2A6,         0}, {0xEA0EA0EA,0xEA0EA0EA}, {0xE865AC7C,         0},
  {0xE6C2B449,         0}, {0xE525982B,         0}, {0xE38E38E4,         0},
  {0xE1FC780F,         0}, {0xE070381C,0xE070381C}, {0xDEE95C4D,         0},
  {0xDD67C8A6,0xDD67C8A6}, {0xDBEB61EF,         0}, {0xDA740DA8,         0},
  {0xD901B204,         0}, {0xD79435E5,0xD79435E5}, {0xD62B80D7,         0},
  {0xD4C77B04,         0}, {0xD3680D37,         0}, {0xD20D20D2,0xD20D20D2},
  {0xD0B69FCC,         0}, {0xCF6474A9,         0}, {0xCE168A77,0xCE168A77},
  {0xCCCCCCCD,         0}, {0xCB8727C1,         0}, {0xCA4587E7,         0},
  {0xC907DA4F,         0}, {0xC7CE0C7D,         0}, {0xC6980C6A,         0},
  {0xC565C87C,         0}, {0xC4372F86,         0}, {0xC30C30C3,0xC30C30C3},
  {0xC1E4BBD6,         0}, {0xC0C0C0C1,         0}, {0xBFA02FE8,0xBFA02FE8},
  {0xBE82FA0C,         0}, {0xBD691047,0xBD691047}, {0xBC52640C,         0},
  {0xBB3EE722,         0}, {0xBA2E8BA3,         0}, {0xB92143FA,0xB92143FA},
  {0xB81702E1,         0}, {0xB70FBB5A,0xB70FBB5A}, {0xB60B60B6,0xB60B60B6},
  {0xB509E68B,         0}, {0xB40B40B4,0xB40B40B4}, {0xB30F6353,         0},
  {0xB21642C9,         0}, {0xB11FD3B8,0xB11FD3B8}, {0xB02C0B03,         0},
  {0xAF3ADDC7,         0}, {0xAE4C415D,         0}, {0xAD602B58,0xAD602B58},
  {0xAC769184,0xAC769184}, {0xAB8F69E3,         0}, {0xAAAAAAAB,         0},
  {0xA9C84A48,         0}, {0xA8E83F57,0xA8E83F57}, {0xA80A80A8,0xA80A80A8},
  {0xA72F053A,         0}, {0xA655C439,0xA655C439}, {0xA57EB503,         0},
  {0xA4A9CF1E,         0}, {0xA3D70A3E,         0}, {0xA3065E40,         0},
  {0xA237C32B,0xA237C32B}, {0xA16B312F,         0}, {0xA0A0A0A1,         0},
  {0x9FD809FE,         0}, {0x9F1165E7,0x9F1165E7}, {0x9E4CAD24,         0},
  {0x9D89D89E,         0}, {0x9CC8E161,         0}, {0x9C09C09C,0x9C09C09C},
  {0x9B4C6F9F,         0}, {0x9A90E7D9,0x9A90E7D9}, {0x99D722DB,         0},
  {0x991F1A51,0x991F1A51}, {0x9868C80A,         0}, {0x97B425ED,0x97B425ED},
  {0x97012E02,0x97012E02}, {0x964FDA6C,0x964FDA6C}, {0x95A02568,0x95A02568},
  {0x94F2094F,0x94F2094F}, {0x94458094,0x94458094}, {0x939A85C4,0x939A85C4},
  {0x92F11384,0x92F11384}, {0x92492492,0x92492492}, {0x91A2B3C5,         0},
  {0x90FDBC09,0x90FDBC09}, {0x905A3863,0x905A3863}, {0x8FB823EE,0x8FB823EE},
  {0x8F1779DA,         0}, {0x8E78356D,0x8E78356D}, {0x8DDA5202,0x8DDA5202},
  {0x8D3DCB09,         0}, {0x8CA29C04,0x8CA29C04}, {0x8C08C08C,0x8C08C08C},
  {0x8B70344A,0x8B70344A}, {0x8AD8F2FC,         0}, {0x8A42F870,0x8A42F870},
  {0x89AE408A,         0}, {0x891AC73B,         0}, {0x88888889,         0},
  {0x87F78088,         0}, {0x8767AB5F,0x8767AB5F}, {0x86D90545,         0},
  {0x864B8A7E,         0}, {0x85BF3761,0x85BF3761}, {0x85340853,0x85340853},
  {0x84A9F9C8,0x84A9F9C8}, {0x84210842,0x84210842}, {0x83993052,0x83993052},
  {0x83126E98,         0}, {0x828CBFBF,         0}, {0x82082082,0x82082082},
  {0x81848DA9,         0}, {0x81020408,0x81020408}, {0x80808081,         0},
  {0xFFFFFFFF,0xFFFFFFFF}, {0xFF00FF01,         0}, {0xFE03F810,         0},
  {0xFD08E551,         0}, {0xFC0FC0FD,         0}, {0xFB188566,         0},
  {0xFA232CF3,         0}, {0xF92FB222,         0}, {0xF83E0F84,         0},
  {0xF74E3FC3,         0}, {0xF6603D99,         0}, {0xF57403D6,         0},
  {0xF4898D60,         0}, {0xF3A0D52D,         0}, {0xF2B9D649,         0},
  {0xF1D48BCF,         0}, {0xF0F0F0F1,         0}, {0xF00F00F0,0xF00F00F0},
  {0xEF2EB720,         0}, {0xEE500EE5,0xEE500EE5}, {0xED7303B6,         0},
  {0xEC979119,         0}, {0xEBBDB2A6,         0}, {0xEAE56404,         0},
  {0xEA0EA0EA,0xEA0EA0EA}, {0xE9396520,         0}, {0xE865AC7C,         0},
  {0xE79372E3,         0}, {0xE6C2B449,         0}, {0xE5F36CB0,0xE5F36CB0},
  {0xE525982B,         0}, {0xE45932D8,         0}, {0xE38E38E4,         0},
  {0xE2C4A689,         0}, {0xE1FC780F,         0}, {0xE135A9CA,         0},
  {0xE070381C,0xE070381C}, {0xDFAC1F75,         0}, {0xDEE95C4D,         0},
  {0xDE27EB2D,         0}, {0xDD67C8A6,0xDD67C8A6}, {0xDCA8F159,         0},
  {0xDBEB61EF,         0}, {0xDB2F171E,         0}, {0xDA740DA8,         0},
  {0xD9BA4257,         0}, {0xD901B204,         0}, {0xD84A598F,         0},
  {0xD79435E5,0xD79435E5}, {0xD6DF43FD,         0}, {0xD62B80D7,         0},
  {0xD578E97D,         0}, {0xD4C77B04,         0}, {0xD417328A,         0},
  {0xD3680D37,         0}, {0xD2BA083C,         0}, {0xD20D20D2,0xD20D20D2},
  {0xD161543E,0xD161543E}, {0xD0B69FCC,         0}, {0xD00D00D0,0xD00D00D0},
  {0xCF6474A9,         0}, {0xCEBCF8BC,         0}, {0xCE168A77,0xCE168A77},
  {0xCD712753,         0}, {0xCCCCCCCD,         0}, {0xCC29786D,         0},
  {0xCB8727C1,         0}, {0xCAE5D85F,0xCAE5D85F}, {0xCA4587E7,         0},
  {0xC9A633FD,         0}, {0xC907DA4F,         0}, {0xC86A7890,0xC86A7890},
  {0xC7CE0C7D,         0}, {0xC73293D8,         0}, {0xC6980C6A,         0},
  {0xC5FE7403,0xC5FE7403}, {0xC565C87C,         0}, {0xC4CE07B0,0xC4CE07B0},
  {0xC4372F86,         0}, {0xC3A13DE6,0xC3A13DE6}, {0xC30C30C3,0xC30C30C3},
  {0xC2780614,         0}, {0xC1E4BBD6,         0}, {0xC152500C,0xC152500C},
  {0xC0C0C0C1,         0}, {0xC0300C03,0xC0300C03}, {0xBFA02FE8,0xBFA02FE8},
  {0xBF112A8B,         0}, {0xBE82FA0C,         0}, {0xBDF59C92,         0},
  {0xBD691047,0xBD691047}, {0xBCDD535E,         0}, {0xBC52640C,         0},
  {0xBBC8408D,         0}, {0xBB3EE722,         0}, {0xBAB65610,0xBAB65610},
  {0xBA2E8BA3,         0}, {0xB9A7862A,0xB9A7862A}, {0xB92143FA,0xB92143FA},
  {0xB89BC36D,         0}, {0xB81702E1,         0}, {0xB79300B8,         0},
  {0xB70FBB5A,0xB70FBB5A}, {0xB68D3134,0xB68D3134}, {0xB60B60B6,0xB60B60B6},
  {0xB58A4855,0xB58A4855}, {0xB509E68B,         0}, {0xB48A39D4,0xB48A39D4},
  {0xB40B40B4,0xB40B40B4}, {0xB38CF9B0,0xB38CF9B0}, {0xB30F6353,         0},
  {0xB2927C2A,         0}, {0xB21642C9,         0}, {0xB19AB5C5,         0},
  {0xB11FD3B8,0xB11FD3B8}, {0xB0A59B42,         0}, {0xB02C0B03,         0},
  {0xAFB321A1,0xAFB321A1}, {0xAF3ADDC7,         0}, {0xAEC33E20,         0},
  {0xAE4C415D,         0}, {0xADD5E632,0xADD5E632}, {0xAD602B58,0xAD602B58},
  {0xACEB0F89,0xACEB0F89}, {0xAC769184,0xAC769184}, {0xAC02B00B,         0},
  {0xAB8F69E3,         0}, {0xAB1CBDD4,         0}, {0xAAAAAAAB,         0},
  {0xAA392F36,         0}, {0xA9C84A48,         0}, {0xA957FAB5,0xA957FAB5},
  {0xA8E83F57,0xA8E83F57}, {0xA8791709,         0}, {0xA80A80A8,0xA80A80A8},
  {0xA79C7B17,         0}, {0xA72F053A,         0}, {0xA6C21DF7,         0},
  {0xA655C439,0xA655C439}, {0xA5E9F6ED,0xA5E9F6ED}, {0xA57EB503,         0},
  {0xA513FD6C,         0}, {0xA4A9CF1E,         0}, {0xA4402910,0xA4402910},
  {0xA3D70A3E,         0}, {0xA36E71A3,         0}, {0xA3065E40,         0},
  {0xA29ECF16,0xA29ECF16}, {0xA237C32B,0xA237C32B}, {0xA1D13986,         0},
  {0xA16B312F,         0}, {0xA105A933,         0}, {0xA0A0A0A1,         0},
  {0xA03C1689,         0}, {0x9FD809FE,         0}, {0x9F747A15,0x9F747A15},
  {0x9F1165E7,0x9F1165E7}, {0x9EAECC8D,0x9EAECC8D}, {0x9E4CAD24,         0},
  {0x9DEB06C9,0x9DEB06C9}, {0x9D89D89E,         0}, {0x9D2921C4,         0},
  {0x9CC8E161,         0}, {0x9C69169B,0x9C69169B}, {0x9C09C09C,0x9C09C09C},
  {0x9BAADE8E,0x9BAADE8E}, {0x9B4C6F9F,         0}, {0x9AEE72FD,         0},
  {0x9A90E7D9,0x9A90E7D9}, {0x9A33CD67,0x9A33CD67}, {0x99D722DB,         0},
  {0x997AE76B,0x997AE76B}, {0x991F1A51,0x991F1A51}, {0x98C3BAC7,0x98C3BAC7},
  {0x9868C80A,         0}, {0x980E4156,0x980E4156}, {0x97B425ED,0x97B425ED},
  {0x975A7510,         0}, {0x97012E02,0x97012E02}, {0x96A8500A,         0},
  {0x964FDA6C,0x964FDA6C}, {0x95F7CC73,         0}, {0x95A02568,0x95A02568},
  {0x9548E498,         0}, {0x94F2094F,0x94F2094F}, {0x949B92DE,         0},
  {0x94458094,0x94458094}, {0x93EFD1C5,0x93EFD1C5}, {0x939A85C4,0x939A85C4},
  {0x93459BE7,         0}, {0x92F11384,0x92F11384}, {0x929CEBF5,         0},
  {0x92492492,0x92492492}, {0x91F5BCB9,         0}, {0x91A2B3C5,         0},
  {0x91500915,0x91500915}, {0x90FDBC09,0x90FDBC09}, {0x90ABCC02,0x90ABCC02},
  {0x905A3863,0x905A3863}, {0x90090090,0x90090090}, {0x8FB823EE,0x8FB823EE},
  {0x8F67A1E4,         0}, {0x8F1779DA,         0}, {0x8EC7AB3A,         0},
  {0x8E78356D,0x8E78356D}, {0x8E2917E1,         0}, {0x8DDA5202,0x8DDA5202},
  {0x8D8BE340,         0}, {0x8D3DCB09,         0}, {0x8CF008CF,0x8CF008CF},
  {0x8CA29C04,0x8CA29C04}, {0x8C55841D,         0}, {0x8C08C08C,0x8C08C08C},
  {0x8BBC50C9,         0}, {0x8B70344A,0x8B70344A}, {0x8B246A88,         0},
  {0x8AD8F2FC,         0}, {0x8A8DCD20,         0}, {0x8A42F870,0x8A42F870},
  {0x89F8746A,         0}, {0x89AE408A,         0}, {0x89645C4F,0x89645C4F},
  {0x891AC73B,         0}, {0x88D180CD,0x88D180CD}, {0x88888889,         0},
  {0x883FDDF0,0x883FDDF0}, {0x87F78088,         0}, {0x87AF6FD6,         0},
  {0x8767AB5F,0x8767AB5F}, {0x872032AC,0x872032AC}, {0x86D90545,         0},
  {0x869222B2,         0}, {0x864B8A7E,         0}, {0x86053C34,0x86053C34},
  {0x85BF3761,0x85BF3761}, {0x85797B91,0x85797B91}, {0x85340853,0x85340853},
  {0x84EEDD36,         0}, {0x84A9F9C8,0x84A9F9C8}, {0x84655D9C,         0},
  {0x84210842,0x84210842}, {0x83DCF94E,         0}, {0x83993052,0x83993052},
  {0x8355ACE4,         0}, {0x83126E98,         0}, {0x82CF7504,         0},
  {0x828CBFBF,         0}, {0x824A4E61,         0}, {0x82082082,0x82082082},
  {0x81C635BC,0x81C635BC}, {0x81848DA9,         0}, {0x814327E4,         0},
  {0x81020408,0x81020408}, {0x80C121B3,         0}, {0x80808081,         0},
  {0x80402010,0x80402010}, {0xFFFFFFFF,0xFFFFFFFF}, {0xFF803FE1,         0},
  {0xFF00FF01,         0}, {0xFE823CA6,         0}, {0xFE03F810,         0},
  {0xFD863087,         0}, {0xFD08E551,         0}, {0xFC8C15B5,         0},
  {0xFC0FC0FD,         0}, {0xFB93E673,         0}, {0xFB188566,         0},
  {0xFA9D9D20,         0}, {0xFA232CF3,         0}, {0xF9A9342D,         0},
  {0xF92FB222,         0}, {0xF8B6A622,0xF8B6A622}, {0xF83E0F84,         0},
  {0xF7C5ED9D,         0}, {0xF74E3FC3,         0}, {0xF6D7054E,         0},
  {0xF6603D99,         0}, {0xF5E9E7FD,         0}, {0xF57403D6,         0},
  {0xF4FE9083,         0}, {0xF4898D60,         0}, {0xF414F9CE,         0},
  {0xF3A0D52D,         0}, {0xF32D1EE0,         0}, {0xF2B9D649,         0},
  {0xF246FACC,         0}, {0xF1D48BCF,         0}, {0xF16288B9,         0},
  {0xF0F0F0F1,         0}, {0xF07FC3E0,0xF07FC3E0}, {0xF00F00F0,0xF00F00F0},
  {0xEF9EA78C,         0}, {0xEF2EB720,         0}, {0xEEBF2F19,         0},
  {0xEE500EE5,0xEE500EE5}, {0xEDE155F4,         0}, {0xED7303B6,         0},
  {0xED05179C,0xED05179C}, {0xEC979119,         0}, {0xEC2A6FA0,0xEC2A6FA0},
  {0xEBBDB2A6,         0}, {0xEB5159A0,         0}, {0xEAE56404,         0},
  {0xEA79D14A,         0}, {0xEA0EA0EA,0xEA0EA0EA}, {0xE9A3D25E,0xE9A3D25E},
  {0xE9396520,         0}, {0xE8CF58AB,         0}, {0xE865AC7C,         0},
  {0xE7FC600F,         0}, {0xE79372E3,         0}, {0xE72AE476,         0},
  {0xE6C2B449,         0}, {0xE65AE1DC,         0}, {0xE5F36CB0,0xE5F36CB0},
  {0xE58C544A,         0}, {0xE525982B,         0}, {0xE4BF37D9,         0},
  {0xE45932D8,         0}, {0xE3F388AF,         0}, {0xE38E38E4,         0},
  {0xE32942FF,         0}, {0xE2C4A689,         0}, {0xE260630B,         0},
  {0xE1FC780F,         0}, {0xE198E520,         0}, {0xE135A9CA,         0},
  {0xE0D2C59A,         0}, {0xE070381C,0xE070381C}, {0xE00E00E0,0xE00E00E0},
  {0xDFAC1F75,         0}, {0xDF4A9369,         0}, {0xDEE95C4D,         0},
  {0xDE8879B3,         0}, {0xDE27EB2D,         0}, {0xDDC7B04D,         0},
  {0xDD67C8A6,0xDD67C8A6}, {0xDD0833CE,         0}, {0xDCA8F159,         0},
  {0xDC4A00DD,         0}, {0xDBEB61EF,         0}, {0xDB8D1428,         0},
  {0xDB2F171E,         0}, {0xDAD16A6B,         0}, {0xDA740DA8,         0},
  {0xDA17006D,0xDA17006D}, {0xD9BA4257,         0}, {0xD95DD300,         0},
  {0xD901B204,         0}, {0xD8A5DEFF,         0}, {0xD84A598F,         0},
  {0xD7EF2152,         0}, {0xD79435E5,0xD79435E5}, {0xD73996E9,         0},
  {0xD6DF43FD,         0}, {0xD6853CC1,         0}, {0xD62B80D7,         0},
  {0xD5D20FDF,         0}, {0xD578E97D,         0}, {0xD5200D52,0xD5200D52},
  {0xD4C77B04,         0}, {0xD46F3235,         0}, {0xD417328A,         0},
  {0xD3BF7BA9,         0}, {0xD3680D37,         0}, {0xD310E6DB,         0},
  {0xD2BA083C,         0}, {0xD2637101,         0}, {0xD20D20D2,0xD20D20D2},
  {0xD1B71759,         0}, {0xD161543E,0xD161543E}, {0xD10BD72C,         0},
  {0xD0B69FCC,         0}, {0xD061ADCA,         0}, {0xD00D00D0,0xD00D00D0},
  {0xCFB8988C,         0}, {0xCF6474A9,         0}, {0xCF1094D4,         0},
  {0xCEBCF8BC,         0}, {0xCE69A00D,         0}, {0xCE168A77,0xCE168A77},
  {0xCDC3B7A9,0xCDC3B7A9}, {0xCD712753,         0}, {0xCD1ED924,         0},
  {0xCCCCCCCD,         0}, {0xCC7B0200,         0}, {0xCC29786D,         0},
  {0xCBD82FC7,         0}, {0xCB8727C1,         0}, {0xCB36600D,         0},
  {0xCAE5D85F,0xCAE5D85F}, {0xCA95906C,         0}, {0xCA4587E7,         0},
  {0xC9F5BE86,         0}, {0xC9A633FD,         0}, {0xC956E803,0xC956E803},
  {0xC907DA4F,         0}, {0xC8B90A96,         0}, {0xC86A7890,0xC86A7890},
  {0xC81C23F5,0xC81C23F5}, {0xC7CE0C7D,         0}, {0xC78031E0,0xC78031E0},
  {0xC73293D8,         0}, {0xC6E5321D,         0}, {0xC6980C6A,         0},
  {0xC64B2278,0xC64B2278}, {0xC5FE7403,0xC5FE7403}, {0xC5B200C6,         0},
  {0xC565C87C,         0}, {0xC519CAE0,0xC519CAE0}, {0xC4CE07B0,0xC4CE07B0},
  {0xC4827EA8,0xC4827EA8}, {0xC4372F86,         0}, {0xC3EC1A06,         0},
  {0xC3A13DE6,0xC3A13DE6}, {0xC3569AE6,         0}, {0xC30C30C3,0xC30C30C3},
  {0xC2C1FF3E,         0}, {0xC2780614,         0}, {0xC22E4507,         0},
  {0xC1E4BBD6,         0}, {0xC19B6A42,         0}, {0xC152500C,0xC152500C},
  {0xC1096CF6,         0}, {0xC0C0C0C1,         0}, {0xC0784B2F,         0},
  {0xC0300C03,0xC0300C03}, {0xBFE80300,         0}, {0xBFA02FE8,0xBFA02FE8},
  {0xBF589280,         0}, {0xBF112A8B,         0}, {0xBEC9F7CE,         0},
  {0xBE82FA0C,         0}, {0xBE3C310C,         0}, {0xBDF59C92,         0},
  {0xBDAF3C64,         0}, {0xBD691047,0xBD691047}, {0xBD231803,         0},
  {0xBCDD535E,         0}, {0xBC97C21E,0xBC97C21E}, {0xBC52640C,         0},
  {0xBC0D38EE,0xBC0D38EE}, {0xBBC8408D,         0}, {0xBB837AB1,         0},
  {0xBB3EE722,         0}, {0xBAFA85A9,0xBAFA85A9}, {0xBAB65610,0xBAB65610},
  {0xBA725820,0xBA725820}, {0xBA2E8BA3,         0}, {0xB9EAF063,         0},
  {0xB9A7862A,0xB9A7862A}, {0xB9644CC4,         0}, {0xB92143FA,0xB92143FA},
  {0xB8DE6B9A,         0}, {0xB89BC36D,         0}, {0xB8594B41,         0},
  {0xB81702E1,         0}, {0xB7D4EA19,0xB7D4EA19}, {0xB79300B8,         0},
  {0xB7514689,         0}, {0xB70FBB5A,0xB70FBB5A}, {0xB6CE5EF9,0xB6CE5EF9},
  {0xB68D3134,0xB68D3134}, {0xB64C31D9,         0}, {0xB60B60B6,0xB60B60B6},
  {0xB5CABD9B,         0}, {0xB58A4855,0xB58A4855}, {0xB54A00B5,0xB54A00B5},
  {0xB509E68B,         0}, {0xB4C9F9A5,         0}, {0xB48A39D4,0xB48A39D4},
  {0xB44AA6E9,0xB44AA6E9}, {0xB40B40B4,0xB40B40B4}, {0xB3CC0706,         0},
  {0xB38CF9B0,0xB38CF9B0}, {0xB34E1884,         0}, {0xB30F6353,         0},
  {0xB2D0D9EF,         0}, {0xB2927C2A,         0}, {0xB25449D7,         0},
  {0xB21642C9,         0}, {0xB1D866D1,0xB1D866D1}, {0xB19AB5C5,         0},
  {0xB15D2F76,         0}, {0xB11FD3B8,0xB11FD3B8}, {0xB0E2A260,0xB0E2A260},
  {0xB0A59B42,         0}, {0xB068BE31,         0}, {0xB02C0B03,         0},
  {0xAFEF818C,         0}, {0xAFB321A1,0xAFB321A1}, {0xAF76EB19,         0},
  {0xAF3ADDC7,         0}, {0xAEFEF982,         0}, {0xAEC33E20,         0},
  {0xAE87AB76,0xAE87AB76}, {0xAE4C415D,         0}, {0xAE10FFA9,         0},
  {0xADD5E632,0xADD5E632}, {0xAD9AF4D0,         0}, {0xAD602B58,0xAD602B58},
  {0xAD2589A4,         0}, {0xACEB0F89,0xACEB0F89}, {0xACB0BCE1,0xACB0BCE1},
  {0xAC769184,0xAC769184}, {0xAC3C8D4A,         0}, {0xAC02B00B,         0},
  {0xABC8F9A0,0xABC8F9A0}, {0xAB8F69E3,         0}, {0xAB5600AC,         0},
  {0xAB1CBDD4,         0}, {0xAAE3A136,         0}, {0xAAAAAAAB,         0},
  {0xAA71DA0D,         0}, {0xAA392F36,         0}, {0xAA00AA01,         0},
  {0xA9C84A48,         0}, {0xA9900FE6,         0}, {0xA957FAB5,0xA957FAB5},
  {0xA9200A92,0xA9200A92}, {0xA8E83F57,0xA8E83F57}, {0xA8B098E0,0xA8B098E0},
  {0xA8791709,         0}, {0xA841B9AD,         0}, {0xA80A80A8,0xA80A80A8},
  {0xA7D36BD8,         0}, {0xA79C7B17,         0}, {0xA765AE44,         0},
  {0xA72F053A,         0}, {0xA6F87FD6,0xA6F87FD6}, {0xA6C21DF7,         0},
  {0xA68BDF79,         0}, {0xA655C439,0xA655C439}, {0xA61FCC16,0xA61FCC16},
  {0xA5E9F6ED,0xA5E9F6ED}, {0xA5B4449D,         0}, {0xA57EB503,         0},
  {0xA54947FE,         0}, {0xA513FD6C,         0}, {0xA4DED52C,0xA4DED52C},
  {0xA4A9CF1E,         0}, {0xA474EB1F,0xA474EB1F}, {0xA4402910,0xA4402910},
  {0xA40B88D0,         0}, {0xA3D70A3E,         0}, {0xA3A2AD39,0xA3A2AD39},
  {0xA36E71A3,         0}, {0xA33A575A,0xA33A575A}, {0xA3065E40,         0},
  {0xA2D28634,         0}, {0xA29ECF16,0xA29ECF16}, {0xA26B38C9,         0},
  {0xA237C32B,0xA237C32B}, {0xA2046E1F,0xA2046E1F}, {0xA1D13986,         0},
  {0xA19E2540,         0}, {0xA16B312F,         0}, {0xA1385D35,         0},
  {0xA105A933,         0}, {0xA0D3150C,         0}, {0xA0A0A0A1,         0},
  {0xA06E4BD4,0xA06E4BD4}, {0xA03C1689,         0}, {0xA00A00A0,0xA00A00A0},
  {0x9FD809FE,         0}, {0x9FA63284,         0}, {0x9F747A15,0x9F747A15},
  {0x9F42E095,0x9F42E095}, {0x9F1165E7,0x9F1165E7}, {0x9EE009EE,0x9EE009EE},
  {0x9EAECC8D,0x9EAECC8D}, {0x9E7DADA9,         0}, {0x9E4CAD24,         0},
  {0x9E1BCAE3,         0}, {0x9DEB06C9,0x9DEB06C9}, {0x9DBA60BB,0x9DBA60BB},
  {0x9D89D89E,         0}, {0x9D596E54,0x9D596E54}, {0x9D2921C4,         0},
  {0x9CF8F2D1,0x9CF8F2D1}, {0x9CC8E161,         0}, {0x9C98ED58,         0},
  {0x9C69169B,0x9C69169B}, {0x9C395D10,0x9C395D10}, {0x9C09C09C,0x9C09C09C},
  {0x9BDA4124,0x9BDA4124}, {0x9BAADE8E,0x9BAADE8E}, {0x9B7B98C0,         0},
  {0x9B4C6F9F,         0}, {0x9B1D6311,0x9B1D6311}, {0x9AEE72FD,         0},
  {0x9ABF9F48,0x9ABF9F48}, {0x9A90E7D9,0x9A90E7D9}, {0x9A624C97,         0},
  {0x9A33CD67,0x9A33CD67}, {0x9A056A31,         0}, {0x99D722DB,         0},
  {0x99A8F74C,         0}, {0x997AE76B,0x997AE76B}, {0x994CF320,0x994CF320},
  {0x991F1A51,0x991F1A51}, {0x98F15CE7,         0}, {0x98C3BAC7,0x98C3BAC7},
  {0x989633DB,0x989633DB}, {0x9868C80A,         0}, {0x983B773B,         0},
  {0x980E4156,0x980E4156}, {0x97E12644,0x97E12644}, {0x97B425ED,0x97B425ED},
  {0x97874039,         0}, {0x975A7510,         0}, {0x972DC45B,         0},
  {0x97012E02,0x97012E02}, {0x96D4B1EF,         0}, {0x96A8500A,         0},
  {0x967C083B,         0}, {0x964FDA6C,0x964FDA6C}, {0x9623C686,0x9623C686},
  {0x95F7CC73,         0}, {0x95CBEC1B,         0}, {0x95A02568,0x95A02568},
  {0x95747844,         0}, {0x9548E498,         0}, {0x951D6A4E,         0},
  {0x94F2094F,0x94F2094F}, {0x94C6C187,         0}, {0x949B92DE,         0},
  {0x94707D3F,         0}, {0x94458094,0x94458094}, {0x941A9CC8,0x941A9CC8},
  {0x93EFD1C5,0x93EFD1C5}, {0x93C51F76,         0}, {0x939A85C4,0x939A85C4},
  {0x9370049C,         0}, {0x93459BE7,         0}, {0x931B4B91,         0},
  {0x92F11384,0x92F11384}, {0x92C6F3AC,0x92C6F3AC}, {0x929CEBF5,         0},
  {0x9272FC48,0x9272FC48}, {0x92492492,0x92492492}, {0x921F64BF,         0},
  {0x91F5BCB9,         0}, {0x91CC2C6C,0x91CC2C6C}, {0x91A2B3C5,         0},
  {0x917952AF,         0}, {0x91500915,0x91500915}, {0x9126D6E5,         0},
  {0x90FDBC09,0x90FDBC09}, {0x90D4B86F,         0}, {0x90ABCC02,0x90ABCC02},
  {0x9082F6B0,         0}, {0x905A3863,0x905A3863}, {0x9031910A,         0},
  {0x90090090,0x90090090}, {0x8FE086E3,         0}, {0x8FB823EE,0x8FB823EE},
  {0x8F8FD7A0,         0}, {0x8F67A1E4,         0}, {0x8F3F82A8,0x8F3F82A8},
  {0x8F1779DA,         0}, {0x8EEF8766,         0}, {0x8EC7AB3A,         0},
  {0x8E9FE542,0x8E9FE542}, {0x8E78356D,0x8E78356D}, {0x8E509BA8,0x8E509BA8},
  {0x8E2917E1,         0}, {0x8E01AA05,         0}, {0x8DDA5202,0x8DDA5202},
  {0x8DB30FC6,0x8DB30FC6}, {0x8D8BE340,         0}, {0x8D64CC5C,         0},
  {0x8D3DCB09,         0}, {0x8D16DF35,0x8D16DF35}, {0x8CF008CF,0x8CF008CF},
  {0x8CC947C5,         0}, {0x8CA29C04,0x8CA29C04}, {0x8C7C057D,         0},
  {0x8C55841D,         0}, {0x8C2F17D2,0x8C2F17D2}, {0x8C08C08C,0x8C08C08C},
  {0x8BE27E39,0x8BE27E39}, {0x8BBC50C9,         0}, {0x8B963829,0x8B963829},
  {0x8B70344A,0x8B70344A}, {0x8B4A451A,         0}, {0x8B246A88,         0},
  {0x8AFEA483,0x8AFEA483}, {0x8AD8F2FC,         0}, {0x8AB355E0,0x8AB355E0},
  {0x8A8DCD20,         0}, {0x8A6858AB,         0}, {0x8A42F870,0x8A42F870},
  {0x8A1DAC60,0x8A1DAC60}, {0x89F8746A,         0}, {0x89D3507D,         0},
  {0x89AE408A,         0}, {0x89894480,         0}, {0x89645C4F,0x89645C4F},
  {0x893F87E8,0x893F87E8}, {0x891AC73B,         0}, {0x88F61A37,0x88F61A37},
  {0x88D180CD,0x88D180CD}, {0x88ACFAEE,         0}, {0x88888889,         0},
  {0x8864298F,         0}, {0x883FDDF0,0x883FDDF0}, {0x881BA59E,         0},
  {0x87F78088,         0}, {0x87D36EA0,         0}, {0x87AF6FD6,         0},
  {0x878B841B,         0}, {0x8767AB5F,0x8767AB5F}, {0x8743E595,         0},
  {0x872032AC,0x872032AC}, {0x86FC9296,0x86FC9296}, {0x86D90545,         0},
  {0x86B58AA8,         0}, {0x869222B2,         0}, {0x866ECD53,0x866ECD53},
  {0x864B8A7E,         0}, {0x86285A23,0x86285A23}, {0x86053C34,0x86053C34},
  {0x85E230A3,0x85E230A3}, {0x85BF3761,0x85BF3761}, {0x859C5060,0x859C5060},
  {0x85797B91,0x85797B91}, {0x8556B8E7,0x8556B8E7}, {0x85340853,0x85340853},
  {0x851169C7,0x851169C7}, {0x84EEDD36,         0}, {0x84CC6290,         0},
  {0x84A9F9C8,0x84A9F9C8}, {0x8487A2D1,         0}, {0x84655D9C,         0},
  {0x84432A1B,0x84432A1B}, {0x84210842,0x84210842}, {0x83FEF802,0x83FEF802},
  {0x83DCF94E,         0}, {0x83BB0C18,         0}, {0x83993052,0x83993052},
  {0x837765F0,0x837765F0}, {0x8355ACE4,         0}, {0x83340520,0x83340520},
  {0x83126E98,         0}, {0x82F0E93D,0x82F0E93D}, {0x82CF7504,         0},
  {0x82AE11DE,         0}, {0x828CBFBF,         0}, {0x826B7E99,0x826B7E99},
  {0x824A4E61,         0}, {0x82292F08,         0}, {0x82082082,0x82082082},
  {0x81E722C2,0x81E722C2}, {0x81C635BC,0x81C635BC}, {0x81A55963,         0},
  {0x81848DA9,         0}, {0x8163D283,         0}, {0x814327E4,         0},
  {0x81228DBF,         0}, {0x81020408,0x81020408}, {0x80E18AB3,         0},
  {0x80C121B3,         0}, {0x80A0C8FB,0x80A0C8FB}, {0x80808081,         0},
  {0x80604836,0x80604836}, {0x80402010,0x80402010}, {0x80200802,0x80200802},
  {0xFFFFFFFF,0xFFFFFFFF}
};

const daala_plane_info od_plane_info_tab[OD_PIX_NB] = {
  {{ 0, 0, 0 }, { 0, 0, 0 }}, /* n/a */
  {{ 0, 0, 0 }, { 0, 0, 0 }}, /* OD_PIX_YUV444 */
  {{ 0, 1, 0 }, { 0, 1, 0 }}, /* OD_PIX_YUV422 */
  {{ 0, 2, 0 }, { 0, 2, 0 }}, /* OD_PIX_YUV411 */
  {{ 0, 1, 1 }, { 0, 1, 1 }}, /* OD_PIX_YUV420 */
  {{ 0, 0, 0 }, { 0, 0, 0 }}, /* OD_PIX_YUV400 */
};

#if defined(OD_ENABLE_ASSERTIONS)
# include <stdio.h>

void od_fatal_impl(const char *_str, const char *_file, int _line) {
  fprintf(stderr, "Fatal (internal) error in %s, line %d: %s\n",
   _file, _line, _str);
  abort();
}
#endif

int od_ilog(uint32_t _v) {
#if defined(OD_CLZ)
  return (OD_CLZ0-OD_CLZ(_v))&-!!_v;
#else
  /*On a Pentium M, this branchless version tested as the fastest on
     1,000,000,000 random 32-bit integers, edging out a similar version with
     branches, and a 256-entry LUT version.*/
  int ret;
  int m;
  ret = !!_v;
  m = !!(_v&0xFFFF0000)<<4;
  _v >>= m;
  ret |= m;
  m = !!(_v&0xFF00)<<3;
  _v >>= m;
  ret |= m;
  m = !!(_v&0xF0)<<2;
  _v >>= m;
  ret |= m;
  m = !!(_v&0xC)<<1;
  _v >>= m;
  ret |= m;
  ret += !!(_v&0x2);
  return ret;
#endif
}

void **od_malloc_2d(size_t _height, size_t _width, size_t _sz) {
  size_t  rowsz;
  size_t  colsz;
  size_t  datsz;
  char *ret;
  colsz = _height*sizeof(void *);
  rowsz = _sz*_width;
  datsz = rowsz*_height;
  /*Alloc array and row pointers.*/
  ret = (char *)malloc(datsz+colsz);
  /*Initialize the array.*/
  if (ret != NULL) {
    size_t   i;
    void **p;
    char *datptr;
    p = (void **)ret;
    i = _height;
    for (datptr = ret + colsz; i-- > 0; p++, datptr += rowsz)
      *p = (void *)datptr;
  }
  return (void **)ret;
}

void **od_calloc_2d(size_t _height, size_t _width, size_t _sz) {
  size_t  colsz;
  size_t  rowsz;
  size_t  datsz;
  char *ret;
  colsz = _height*sizeof(void *);
  rowsz = _sz*_width;
  datsz = rowsz*_height;
  /*Alloc array and row pointers.*/
  ret = (char *)calloc(datsz + colsz, 1);
  /*Initialize the array.*/
  if (ret != NULL) {
    size_t   i;
    void **p;
    char *datptr;
    p = (void **)ret;
    i = _height;
    for (datptr = ret + colsz; i-- > 0; p++, datptr += rowsz)
      *p = (void *)datptr;
  }
  return (void **)ret;
}

void od_free_2d(void *_ptr) {
  free(_ptr);
}

#define BUFFER_INCREMENT (256)

void oggbyte_writeinit(oggbyte_buffer *_b) {
  OD_CLEAR(_b, 1);
  _b->ptr = _b->buf = (unsigned char *)malloc(BUFFER_INCREMENT);
  _b->storage = BUFFER_INCREMENT;
}

void oggbyte_writetrunc(oggbyte_buffer *_b, ptrdiff_t _bytes) {
  OD_ASSERT(_bytes >= 0);
  _b->ptr = _b->buf + _bytes;
}

void oggbyte_write1(oggbyte_buffer *_b, unsigned _value) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte >= _b->storage) {
    _b->buf = (unsigned char *)realloc(_b->buf,
            _b->storage + BUFFER_INCREMENT);
    _b->storage += BUFFER_INCREMENT;
    _b->ptr = _b->buf+endbyte;
  }
  *(_b->ptr++) = (unsigned char)_value;
}

void oggbyte_write4(oggbyte_buffer *_b, uint32_t _value) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr - _b->buf;
  if (endbyte+4 > _b->storage) {
    _b->buf = (unsigned char *)realloc(_b->buf,
            _b->storage + BUFFER_INCREMENT);
    _b->storage += BUFFER_INCREMENT;
    _b->ptr = _b->buf + endbyte;
  }
  *(_b->ptr++) = (unsigned char)_value;
  _value >>= 8;
  *(_b->ptr++) = (unsigned char)_value;
  _value >>= 8;
  *(_b->ptr++) = (unsigned char)_value;
  _value >>= 8;
  *(_b->ptr++) = (unsigned char)_value;
}

void oggbyte_writecopy(oggbyte_buffer *_b, const void *_source,
 ptrdiff_t _bytes) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte+_bytes > _b->storage) {
    _b->storage = endbyte+_bytes+BUFFER_INCREMENT;
    _b->buf = (unsigned char *)realloc((void *)_b->buf, _b->storage);
    _b->ptr = _b->buf+endbyte;
  }
  memmove(_b->ptr, _source, _bytes);
  _b->ptr += _bytes;
}

void oggbyte_reset(oggbyte_buffer *_b) {
  _b->ptr = _b->buf;
}

void oggbyte_writeclear(oggbyte_buffer *_b) {
  free(_b->buf);
  OD_CLEAR(_b, 1);
}

void oggbyte_readinit(oggbyte_buffer *_b, unsigned char *_buf,
 ptrdiff_t _bytes) {
  OD_ASSERT(_bytes >= 0);
  OD_CLEAR(_b, 1);
  _b->buf = _b->ptr = _buf;
  _b->storage = _bytes;
}

int oggbyte_look1(oggbyte_buffer *_b) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr - _b->buf;
  if (endbyte >= _b->storage) return -1;
  else return _b->ptr[0];
}

int oggbyte_look4(oggbyte_buffer *_b, uint32_t *_val) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte > _b->storage-4) {
    if (endbyte < _b->storage) {
      *_val = _b->ptr[0];
      endbyte++;
      if (endbyte < _b->storage) {
        *_val |= (uint32_t)_b->ptr[1]<<8;
        endbyte++;
        if (endbyte < _b->storage) *_val |= (uint32_t)_b->ptr[2]<<16;
      }
    }
    return -1;
  }
  else {
    *_val = _b->ptr[0];
    *_val |= (uint32_t)_b->ptr[1]<<8;
    *_val |= (uint32_t)_b->ptr[2]<<16;
    *_val |= (uint32_t)_b->ptr[3]<<24;
  }
  return 0;
}

void oggbyte_adv1(oggbyte_buffer *_b) {
  _b->ptr++;
}

void oggbyte_adv4(oggbyte_buffer *_b) {
  _b->ptr += 4;
}

int oggbyte_read1(oggbyte_buffer *_b) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr-_b->buf;
  if (endbyte >= _b->storage) return -1;
  else return *(_b->ptr++);
}

int oggbyte_read4(oggbyte_buffer *_b, uint32_t *_val) {
  unsigned char *end;
  end = _b->buf+_b->storage;
  if (_b->ptr+4 > end) {
    if (_b->ptr < end) {
      *_val = *(_b->ptr++);
      if (_b->ptr < end) {
        *_val |= (uint32_t)*(_b->ptr++)<<8;
        if (_b->ptr < end) *_val |= (uint32_t)*(_b->ptr++)<<16;
      }
    }
    return -1;
  }
  else {
    *_val = (*_b->ptr++);
    *_val |= (uint32_t)*(_b->ptr++)<<8;
    *_val |= (uint32_t)*(_b->ptr++)<<16;
    *_val |= (uint32_t)*(_b->ptr++)<<24;
  }
  return 0;
}

int oggbyte_readcopy(oggbyte_buffer *_b, void *_dest, uint32_t _bytes) {
  ptrdiff_t endbyte;
  endbyte = _b->ptr - _b->buf;
  OD_ASSERT(endbyte >= 0);
  OD_ASSERT(endbyte <= _b->storage);
  if ((size_t)(_b->storage - endbyte) < _bytes) return -1;
  memcpy(_dest, _b->ptr, _bytes);
  _b->ptr += _bytes;
  return 0;
}

ptrdiff_t oggbyte_bytes(oggbyte_buffer *_b) {
  return _b->ptr-_b->buf;
}

ptrdiff_t oggbyte_bytes_left(oggbyte_buffer *_b) {
  return _b->storage - oggbyte_bytes(_b);
}

unsigned char *oggbyte_get_buffer(oggbyte_buffer *_b) {
  return _b->buf;
}

const char *daala_version_string(void) {
  return OD_VENDOR_STRING;
}

uint32_t daala_version_number(void) {
  return OD_VERSION_MAJOR<<16|OD_VERSION_MINOR<<8|OD_VERSION_SUB;
}

int daala_packet_isheader(daala_packet *dpkt) {
  return dpkt->bytes > 0 ? dpkt->packet[0] >> 7 : 0;
}

int daala_packet_iskeyframe(daala_packet *dpkt) {
  return dpkt->bytes > 0 ? dpkt->packet[0] & 0x40 : 0;
}
