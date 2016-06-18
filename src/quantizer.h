/*Daala video codec
Copyright (c) 2015 Daala project contributors.  All rights reserved.

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


#if !defined(_quantizer_H)
# define _quantizer_H (1)

/*Fractional_coded_quantizer ~=
   log2(quantizer / (1 << OD_COEFF_SHIFT))*6.307 + 6.235*/
/*Base/scale factor for linear quantizer to fractional coded quantizer
   conversion (6.307 * 2^12) */
#define OD_LOG_QUANTIZER_BASE_Q12 (0x0064EB)
/*Inverse of above scale factor.*/
#define OD_LOG_QUANTIZER_EXP_Q12 (0x000289)
/*Offset for linear quantizer to fractional coded quantizer
   conversion (6.235 * 2^45) */
#define OD_LOG_QUANTIZER_OFFSET_Q45 (0x0000C7851EB851ECLL)

extern const int OD_N_CODED_QUANTIZERS;
int od_quantizer_to_codedquantizer(int q);
int od_codedquantizer_to_quantizer(int cq);

#endif
