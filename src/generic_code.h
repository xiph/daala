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

#ifndef GENERIC_ENCODER_H_
#define GENERIC_ENCODER_H_

#include "entenc.h"
#include "entdec.h"

#define GENERIC_TABLES 12

typedef struct {
  unsigned short icdf[GENERIC_TABLES][16];
  unsigned short tot[GENERIC_TABLES];
  int increment;
} GenericEncoder;

void generic_encode(ec_enc *enc, GenericEncoder *model, int x, int *ExQ16, int integration);

int generic_decode(ec_dec *dec, GenericEncoder *model, int *ExQ16, int integration);

static inline int logEx(int ExQ16)
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

static inline void generic_model_update(GenericEncoder *model,int *ExQ16,int x,int xs,int id,int integration)
{
  int i;
  int xenc;
  unsigned short *icdf;

  icdf=model->icdf[id];
  /* Renormalize if we cannot add increment */
  if (model->tot[id]>65535-model->increment){
    for (i=0;i<16;i++){
      /* Second term ensures that the pdf is non-null */
      icdf[i]=(icdf[i]>>1)+(15-i);
    }
    model->tot[id]=model->tot[id]/2+16;
  }

  /* Update freq count */
  xenc=OD_MINI(15,xs);
  /* This can be easily vectorized */
  for (i=0;i<xenc;i++)
    icdf[i]+=model->increment;
  model->tot[id]+=model->increment;

  *ExQ16+=(x<<(16-integration))-(*ExQ16>>integration);
}

#endif
