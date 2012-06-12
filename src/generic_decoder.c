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

extern const unsigned short expectA[];

int generic_decode(ec_dec *dec, GenericEncoder *model, int *ExQ4)
{
  int lgQ1;
  int shift;
  int id;
  unsigned char *icdf;
  int xs;
  int xdec;
  int i;
  int lsb=0;
  int x;

  lgQ1 = OD_MAXI(0, od_ilog(*ExQ4**ExQ4)-8);

  shift = OD_MAXI(0, (lgQ1-5)>>1);

  id = OD_MINI(GENERIC_TABLES-1, lgQ1);
  icdf = model->icdf[id];

  /* FIXME: need to have an adaptive ft, not 2^15 */
  /*for (i=0;i<16;i++)printf("%d ", icdf[i]);printf("(%d %d %d %d)\n", (int)(model->tot[id]), x, *ExQ4, id);*/
  xs = ec_dec_icdf_ft(dec, icdf, model->tot[id]);

  if (xs == 15)
  {
    int a;
    unsigned decay;
    /* FIXME: Check the bounds here */
    /*fprintf(stderr, "xs15: %d %d\n", *ExQ4, shift);*/
    a=expectA[(*ExQ4+(1<<shift>>1))>>shift];
    decay = 256*exp(-256./a);

    //printf("%x ", dec->rng);
    xs += laplace_decode_special(dec, decay);

    //printf("dec15 %d %d %d %x\n", xs, a, decay, dec->rng);
    /*printf("xs>15: %d %d %f\n", x, xs, *ExQ4/16.);*/
  }

  if (shift!=0)
  {
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(xs==0);
    lsb = ec_dec_bits(dec, shift-special);
    lsb -= !special<<(shift-1);
    /*bits_used += shift-special;*/
  }
  /* Renormalize if we cannot add increment */
  if (model->tot[id]>255-model->increment)
  {
    for (i=0;i<16;i++)
    {
      /* Second term ensures that the pdf is non-null */
      icdf[i]=(icdf[i]>>1) + (15-i);
    }
    model->tot[id]=model->tot[id]/2+16;
  }

  /* Update freq count */
  xdec = OD_MINI(15, xs);
  /* This can be easily vectorized */
  for (i=0;i<xdec;i++)
    icdf[i]+=model->increment;
  model->tot[id] += model->increment;

  x = (xs<<shift)+lsb;
  *ExQ4 += 1 + x - (*ExQ4>>4);
  *ExQ4 = OD_MINI(32767, *ExQ4);

  /*printf("dec: %d %d %d %d %d %x\n", *ExQ4, x, shift, id, xs, dec->rng);*/
  return x;
}
