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
  unsigned short icdf[GENERIC_TABLES][16]; /**< "inverse" cdf for multiple expectations of x */
  unsigned short tot[GENERIC_TABLES]; /**< total frequency for each icdf */
  int increment; /**< Frequency increment for learning the icdfs */
} GenericEncoder;

void generic_model_init(GenericEncoder *model);

void generic_encode(od_ec_enc *enc, GenericEncoder *model, int x, int *ExQ16, int integration);

int generic_decode(od_ec_dec *dec, GenericEncoder *model, int *ExQ16, int integration);

/** Takes the base-2 log of E(x)
 *
 * @param [in] ExQ16 expectation of x in Q16
 *
 * @retval 2*log2(ExQ16/2^16)
 */
int logEx(int ExQ16);

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
void generic_model_update(GenericEncoder *model,int *ExQ16,int x,int xs,int id,int integration);

#endif
