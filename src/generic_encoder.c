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
#include "entenc.h"
#include "entdec.h"
#include "odintrin.h"
#include "pvq_code.h"

extern const unsigned short expectA[];

double generic_bits=0;

void generic_model_init(GenericEncoder *model)
{
  int i, j;
  model->increment = 2;
  for (i=0;i<GENERIC_TABLES;i++)
  {
    for(j=0;j<16;j++)
    {
      /* FIXME: Come on, we can do better than that! */
      model->icdf[i][j] = (15-j)<<3;
    }
    model->tot[i] = 128;
  }
}

void generic_encode(ec_enc *enc, GenericEncoder *model, int x, int *ExQ4)
{
  int lgQ1;
  int shift;
  int id;
  unsigned char *icdf;
  int xs;
  int xenc;
  int i;

  lgQ1 = OD_MAXI(0, od_ilog(*ExQ4**ExQ4)-8);

  shift = OD_MAXI(0, (lgQ1-5)>>1);

  id = OD_MINI(GENERIC_TABLES-1, lgQ1);
  icdf = model->icdf[id];

  xs = (x+(1<<shift>>1))>>shift;
  xenc = OD_MINI(15, xs);

  /* FIXME: need to have an adaptive ft, not 2^15 */
  /*for (i=0;i<16;i++)printf("%d ", icdf[i]);printf("(%d %d %d %d)\n", (int)(model->tot[id]), x, *ExQ4, id);*/
  ec_enc_icdf_ft(enc, xenc, icdf, model->tot[id]);

  if (xs >= 15)
  {
    int a;
    unsigned decay;
    /* FIXME: Check the bounds here */
    a=expectA[(*ExQ4+(1<<shift>>1))>>shift];
    decay = 256*exp(-256./a);

    //printf ("%d %x ", x, enc->rng);
    laplace_encode_special(enc, xs-15, decay);
    //printf("enc15 %d %d %d %x\n", xs, a, decay, enc->rng);

    /*printf("xs>15: %d %d %f\n", x, xs, *ExQ4/16.);*/
  }

  if (shift!=0)
  {
    /*printf("shift=%d: %d %d %f\n", shift, x, xs, *ExQ4/16.);*/
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(xs==0);
    printf("encbits: %d %d %d %d\n", x, xs, x-(xs<<shift)+!special, shift-special);
    ec_enc_bits(enc, x-(xs<<shift)+!special, shift-special);
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
  xenc = OD_MINI(15, xs);
  /* This can be easily vectorized */
  for (i=0;i<xenc;i++)
    icdf[i]+=model->increment;
  model->tot[id] += model->increment;

  *ExQ4 += x - (*ExQ4>>4);
  //printf("%d %d %d\n", *ExQ4, x, enc->rng);
}


int main(void)
{
  int i;
  ec_enc enc;
  ec_dec dec;
  unsigned char buf[1<<22];
  int nb_vectors=0;
  int ExQ4=16;
  GenericEncoder gen;

  ec_enc_init(&enc, buf, 1<<22);
  generic_model_init(&gen);
  while (1)
  {
    int K;
    scanf("%d ", &K);
    if (feof(stdin))
      break;
    nb_vectors++;
    generic_encode(&enc, &gen, K, &ExQ4);
  }
  printf("DECODE\n");
  /*for (i=0;i<GENERIC_TABLES;i++)
  {
    int j;
    printf("%d ", gen.tot[i]);
    for (j=0;j<16;j++)
      printf("%d ", gen.icdf[i][j]);
    printf("\n");
  }*/
  ExQ4=16;
  ec_dec_init(&dec, buf, 1<<22);
  generic_model_init(&gen);
  for (i=0;i<nb_vectors;i++)
  {
    int K = generic_decode(&dec, &gen, &ExQ4);
    printf ("%d\n", K);
  }
  printf("average bits = %f\n", generic_bits/(float)nb_vectors);
  printf("tell()'s average = %f\n", ec_tell(&enc)/(float)nb_vectors);
  return 0;
}
