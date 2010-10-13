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
