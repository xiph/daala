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

#ifndef PVQ_CODE_H_
#define PVQ_CODE_H_

#include "entenc.h"
#include "entdec.h"

void pvq_encoder(ec_enc *enc, const int *y,int N,int K,int *num, int *den, int *u);
void pvq_decoder(ec_dec *dec, int *y,int N,int K,int *num, int *den, int *u);


#endif
