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

#include <stdio.h>

#include "generic_code.h"
#include "entdec.h"
#include "odintrin.h"
#include "pvq_code.h"

/** Encodes a random variable using a "generic" model, assuming that the distribution is
 * one-sided (zero and up), has a single mode, and decays exponentially passed the mode
 *
 * @param [in,out] dec   range decoder
 * @param [in,out] model generic probability model
 * @param [in,out] ExQ16 expectation of x (adapted)
 * @param [in]     integration integration period of ExQ16 (leaky average over 1<<integration samples)
 *
 * @retval decoded variable x
 */
int generic_decode(ec_dec *dec, GenericEncoder *model, int *ExQ16, int integration)
{
  int lgQ1;
  int shift;
  int id;
  unsigned short *icdf;
  int xs;
  int lsb=0;
  int x;

  lgQ1=logEx(*ExQ16);

  shift=OD_MAXI(0,(lgQ1-5)>>1);

  id=OD_MINI(GENERIC_TABLES-1,lgQ1);
  icdf=model->icdf[id];

  xs=ec_dec_icdf16_ft(dec,icdf,model->tot[id]);

  if(xs==15){
    unsigned decay;
    /* Bounds should be OK as long as shift is consistent with Ex */
    decay=decayE[((*ExQ16>>12)+(1<<shift>>1))>>shift];

    xs+=laplace_decode_special(dec,decay,-1);
  }

  if(shift!=0){
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(xs==0);
    lsb=ec_dec_bits(dec,shift-special);
    lsb-=!special<<(shift-1);
  }
  x = (xs<<shift)+lsb;

  generic_model_update(model,ExQ16,x,xs,id,integration);

  /*printf("dec: %d %d %d %d %d %x\n", *ExQ4, x, shift, id, xs, dec->rng);*/
  return x;
}
