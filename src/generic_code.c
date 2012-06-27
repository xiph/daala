/*
    Daala video codec
    Copyright (C) 2012 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "generic_code.h"

/** Takes the base-2 log of E(x)
 *
 * @param [in] ExQ16 expectation of x in Q16
 *
 * @retval 2*log2(ExQ16/2^16)
 */
int logEx(int ExQ16)
{
  int lg;
  int lgQ1;
  int odd;

  lg = od_ilog(ExQ16);
  if (lg<15)
  {
    odd = ExQ16*ExQ16 > 2<<2*lg;
  } else {
    int tmp=ExQ16>>(lg-8);
    odd = tmp*tmp > (1<<15);
  }
  lgQ1 = OD_MAXI(0,2*lg - 33 + odd);
  return lgQ1;
}

/** Updates the probability model based on the encoded/decoded value
 *
 * @param [in,out] model generic prob model
 * @param [in,out] ExQ16 expectation of x
 * @param [in]     x     variable encoded/decoded (used for ExQ16)
 * @param [in]     xs    variable x after shift (used for the model)
 * @param [in]     id    id of the icdf to adapt
 * @param [in]     integration integration period of ExQ16 (leaky average over 1<<integration samples)
 *
 */
void generic_model_update(GenericEncoder *model,int *ExQ16,int x,int xs,int id,int integration)
{
  int i;
  int xenc;
  ogg_uint16_t *cdf;

  cdf=model->cdf[id];
  /* Renormalize if we cannot add increment */
  if (cdf[15]+model->increment>32767){
    for (i=0;i<16;i++){
      /* Second term ensures that the pdf is non-null */
      cdf[i]=(cdf[i]>>1)+i+1;
    }
  }

  /* Update freq count */
  xenc=OD_MINI(15,xs);
  /* This can be easily vectorized */
  for (i=xenc;i<16;i++)
    cdf[i]+=model->increment;

  *ExQ16+=(x<<(16-integration))-(*ExQ16>>integration);
}
