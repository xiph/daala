/*
    Daala video codec
    Copyright (C) 2006-2010 Daala project contributors

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


#include "x86int.h"
#include "cpu.h"
#if defined(OD_X86ASM)

void od_state_opt_vtbl_init_x86(od_state *_state){
  od_state_opt_vtbl_init_c(_state);
  _state->cpu_flags=od_cpu_flags_get();
  if(_state->cpu_flags&OD_CPU_X86_SSE2){
    if(_state->info.frame_width+(OD_UMV_PADDING<<1)<<1<0x7FFF){
      /*We can only use this optimization with (signed) 16-bit strides, until
         SSE4 comes out with pmulld for true 32x32 multiplies.*/
      _state->opt_vtbl.mc_predict1imv8=od_mc_predict1imv8_sse2;
    }
    _state->opt_vtbl.mc_predict1fmv8=od_mc_predict1fmv8_sse2;
    _state->opt_vtbl.mc_blend_full8=od_mc_blend_full8_sse2;
    _state->opt_vtbl.mc_blend_full_split8=od_mc_blend_full_split8_sse2;
  }
}

#endif
