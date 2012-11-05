/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

/*TODO:
 - Develop a real encoding and measure real bits.
 - Compute bits needed for labels during DP (to bias towards using the same
    label).
 - Thresholds for DP.
   + How do we calculate them?
   + How do they propagate between frames (block sizes change)
   + Compute rate change of trailing MVs correctly.*/

/*The frame number to animate.*/
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
#define ANI_FRAME (69)
#endif

/*Flags describing which edge types we allow.
  One of these, at a minimum, must be set.*/
/*Flag indicating we use `B' edges.*/
#define OD_MC_USEB       (1<<0)
/*Flag indicating we use `V' edges.*/
#define OD_MC_USEV       (1<<1)
/*Flag indicating we include the chroma planes in our SAD calculations.*/
#define OD_MC_USE_CHROMA (1<<2)

typedef struct od_mv_node            od_mv_node;
typedef struct od_mv_dp_state        od_mv_dp_state;
typedef struct od_mv_dp_node         od_mv_dp_node;
typedef struct od_mv_err_node        od_mv_err_node;

#include "mc.h"
#include "encint.h"

typedef int          od_offset[2];
typedef int          od_pattern[8];
typedef ogg_uint16_t od_sad4[4];



/*The state information used by the motion estimation process that is not
   required by the decoder.
  Some of this information corresponds to a vertex in the MV mesh.
  Other pieces correspond to a block whose upper-left corner is located at that
   vertex.*/
struct od_mv_node{
  /*The historical motion vectors for EPZS^2, stored at full-pel resolution.
    Indexed by [time][reference_type][component].*/
  int           mvs[3][2][2];
  /*The current estimated rate of this MV.*/
  unsigned      mv_rate:16;
  /*The current estimated rate of the edge labels.*/
  unsigned      lb_rate:4;
  /*The number of blocks influenced by this MV who failed their SAD checks.*/
  unsigned      needs_check:4;
  /*The current size of the block with this MV at its upper-left.*/
  unsigned      log_mvb_sz:2;
  /*The index of the exterior corner of that block.*/
  unsigned      c:2;
  /*The edge splitting index of that block.*/
  unsigned      s:2;
  /*The current distortion of that block.*/
  ogg_int32_t   sad;
  /*The SAD for BMA predictor centered on this node.
    Used for the dynamic thresholds of the initial EPZS^2 pass.*/
  ogg_int32_t   bma_sad;
  /*The location of this node in the grid.*/
  int           vx;
  int           vy;
  /*The change in global distortion for decimating this node.*/
  ogg_int32_t   dd;
  /*The change in global rate for decimating this node.*/
  int           dr;
  /*The position of this node in the heap.*/
  int           heapi;
};

/*The square pattern, the largest we use, has 9 states.*/
#define OD_DP_NSTATES_MAX     (9)
/*Up to 8 blocks can be influenced by this MV and the previous MV.*/
#define OD_DP_NBLOCKS_MAX     (8)
/*Up to 20 MVs can be predicted by this one, but 3 of those are MVs on the
   DP trellis whose value we have yet to determine.*/
#define OD_DP_NPREDICTED_MAX  (17)
/*At most 6 of them can be changed by a subsequent MV on the DP path.*/
#define OD_DP_NCHANGEABLE_MAX (6)

/*One of the trellis states in the dynamic prgram.*/
struct od_mv_dp_state{
  /*The MV to install for this state.*/
  int           mv[2];
  /*The best state in the previous DP node to use with this one, or -1 to
     indicate the start of the path.*/
  int           prevsi;
  /*The total rate change (thus far) produced by choosing this path.*/
  int           dr;
  /*The total distortion change (thus far) produced by choosing this path.*/
  ogg_int32_t   dd;
  /*The new SAD of each block affected by the the DP between this node and the
     previous node.
    These are installed if the path is selected.*/
  ogg_int32_t   block_sads[OD_DP_NBLOCKS_MAX];
  /*The new rate of each MV predicted by this node.
    These are installed if the path is selected.
    These may supersede the rates reported in previous nodes on the path.*/
  int           pred_mv_rates[OD_DP_NPREDICTED_MAX];
  /*The new rate of this MV.*/
  int           mv_rate;
};

/*A node on the dynamic programming path.*/
struct od_mv_dp_node{
  od_mv_grid_pt  *mvg;
  od_mv_node     *mv;
  /*The number of states considered in this node.*/
  int             nstates;
  /*The number of blocks affected by states in this node.*/
  int             nblocks;
  /*The number of MVs predicted by this node.*/
  int             npredicted;
  /*The number of those MVs that are potentially changeable by future DP
     states.*/
  int             npred_changeable;
  /*The original MV used by this node.*/
  int             original_mv[2];
  /*The original edge label used b this node.*/
  unsigned        original_etype:1;
  /*The original rate of this MV.*/
  int             original_mv_rate;
  /*The original MV rates before predictors were changed by this node.
    This only includes the ones that are actually changeable.*/
  int             original_mv_rates[OD_DP_NCHANGEABLE_MAX];
  /*The last node we save/restore in order to perform prediction.*/
  od_mv_dp_node  *min_predictor_node;
  /*The set of trellis states.*/
  od_mv_dp_state  states[OD_DP_NSTATES_MAX];
  /*The blocks influenced by this MV and the previous MV.*/
  od_mv_node     *blocks[OD_DP_NBLOCKS_MAX];
  /*The vertices whose MV we predict.*/
  od_mv_grid_pt  *predicted_mvgs[OD_DP_NPREDICTED_MAX];
  od_mv_node     *predicted_mvs[OD_DP_NPREDICTED_MAX];
};

struct od_mv_est_ctx{
  od_enc_ctx      *enc;
  /*A cache of the SAD values used during decimation.
    Indexed by [log_mvb_sz][vy>>log_mvb_sz][vx>>log_mvb_sz][s], where s is the
     edge split state.
    The SAD of top-level blocks (log_mvb_sz==2) is not stored in this cache,
     since it is only needed once.*/
  od_sad4        **sad_cache[2];
  /*The state of the MV mesh specific to the encoder.*/
  od_mv_node     **mvs;
  /*A temporary copy of the decoder-side MV grid used to save-and-restore the
     MVs when attempting sub-pel refinement.*/
  od_mv_grid_pt  **refine_grid;
  /*Space for storing the Viterbi trellis used for DP refinment.*/
  od_mv_dp_node   *dp_nodes;
  /*The decimation heap.*/
  od_mv_node     **dec_heap;
  /*The number of vertices in the decimation heap.*/
  int              dec_nheap;
  /*The number of undecimated vertices in each row.*/
  unsigned        *row_counts;
  /*The number of undecimated vertices in each column.*/
  unsigned        *col_counts;
  /*The maximum SAD value for accepting set A predictors for each block size.*/
  int              thresh1[3];
  /*The offsets to inflate the second threshold by for each block size.*/
  int              thresh2_offs[3];
  /*The weights used to produce the accelerated MV predictor.*/
  ogg_int32_t      mvapw[2][2];
  /*Flags indicating which MVs have already been tested during the initial
     EPZS^2 pass.*/
  unsigned char    hit_cache[64][64];
  /*The flag used by the current EPZS search iteration.*/
  unsigned         hit_bit;
  /*The Lagrangian multiplier used for R-D optimization.*/
  int              lambda;
  /*Configuration.*/
  /*The flags indicating which feature to use.*/
  int              flags;
  /*The smallest resolution to refine MVs to.*/
  int              mv_res_min;
  /*The deepest level to refine to (inclusive).*/
  int              level_max;
  /*The shallowest level to decimate to (inclusive).*/
  int              level_min;
};



/*The number of bits to reduce chroma SADs by, if used.*/
#define OD_MC_CHROMA_SCALE (2)



/*The subdivision level of a MV in the mesh, given its position (mod 4).*/
static const int OD_MC_LEVEL[4][4]={
  {0,4,2,4},
  {4,3,4,3},
  {2,4,1,4},
  {4,3,4,3}
};

/*Ancestor lists for a vertex.
  These are stored as lists of offsets to the vertices in the domain.
  Level 0 ancestors are not included, as they cannot be decimated.*/
/*Lists for level 2 vertices.*/
static const od_offset OD_ANCESTORS2[2][2]={
  {{ 0,-2},{ 0, 2}},
  {{-2, 0},{ 2, 0}},
};
/*Lists for level 3 vertices.*/
static const od_offset OD_ANCESTORS3[4][5]={
  {{ 1,-1},{-1, 1},{ 1,-3},{-3, 1},{ 1, 1}},
  {{-1,-1},{ 1, 1},{-1,-3},{-1, 1},{ 3, 1}},
  {{-1,-1},{ 1, 1},{-3,-1},{ 1,-1},{ 1, 3}},
  {{ 1,-1},{-1, 1},{-1,-1},{ 3,-1},{-1, 3}},
};
/*Lists for level 4 vertices.*/
static const od_offset OD_ANCESTORS4[8][9]={
  {{ 0,-1},{ 0, 1},{-1,-2},{ 1, 0},{-1, 2},{-3,-2},{ 1,-2},{-3, 2},{ 1, 2}},
  {{ 0,-1},{ 0, 1},{ 1,-2},{-1, 0},{ 1, 2},{-1,-2},{ 3,-2},{-1, 2},{ 3, 2}},
  {{-1, 0},{ 1, 0},{-2,-1},{ 2,-1},{ 0, 1},{-2,-3},{ 2,-3},{-2, 1},{ 2, 1}},
  {{-1, 0},{ 1, 0},{ 0,-1},{-2, 1},{ 2, 1},{ 0,-3},{-4, 1},{ 0, 1},{ 4, 1}},
  {{ 0,-1},{ 0, 1},{ 1,-2},{-1, 0},{ 1, 2},{ 1,-4},{-3, 0},{ 1, 0},{ 1, 4}},
  {{ 0,-1},{ 0, 1},{-1,-2},{ 1, 0},{-1, 2},{-1,-4},{-1, 0},{ 3, 0},{-1, 4}},
  {{-1, 0},{ 1, 0},{ 0,-1},{-2, 1},{ 2, 1},{-2,-1},{ 2,-1},{-2, 3},{ 2, 3}},
  {{-1, 0},{ 1, 0},{-2,-1},{ 2,-1},{ 0, 1},{-4,-1},{ 0,-1},{ 4,-1},{ 0, 3}},
};
/*The number of ancestors in each list in the grid pattern.*/
static const int OD_NANCESTORS[4][4]={
  {0,9,2,9},
  {9,5,9,5},
  {2,9,0,9},
  {9,5,9,5}
};
/*The lists for each vertex in the grid pattern.*/
static const od_offset *OD_ANCESTORS[4][4]={
  {NULL,            OD_ANCESTORS4[0],OD_ANCESTORS2[0],OD_ANCESTORS4[1]},
  {OD_ANCESTORS4[2],OD_ANCESTORS3[0],OD_ANCESTORS4[3],OD_ANCESTORS3[1]},
  {OD_ANCESTORS2[1],OD_ANCESTORS4[4],NULL            ,OD_ANCESTORS4[5]},
  {OD_ANCESTORS4[6],OD_ANCESTORS3[2],OD_ANCESTORS4[7],OD_ANCESTORS3[3]}
};



/*Computes the SAD of the input image against the given predictor.*/
static ogg_int32_t od_state_sad8(od_state *_state,const unsigned char *_p,
 int _pystride,int _pxstride,int _pli,int _x,int _y,int _log_blk_sz){
  od_img_plane        *iplane;
  const unsigned char *p;
  unsigned char       *src;
  int                  clipx;
  int                  clipy;
  int                  clipw;
  int                  cliph;
  int                  w;
  int                  h;
  int                  i;
  int                  j;
  ogg_int32_t          ret;
  iplane=_state->input.planes+_pli;
  /*Compute the block dimensions in the target image plane.*/
  _x>>=iplane->xdec;
  _y>>=iplane->ydec;
  w=1<<_log_blk_sz-iplane->xdec;
  h=1<<_log_blk_sz-iplane->ydec;
  /*Clip the block against the active picture region.*/
  clipx=(_state->info.pic_x>>iplane->xdec)-_x;
  if(clipx>0){
    w-=clipx;
    _p+=clipx*_pxstride;
    _x+=clipx;
  }
  clipy=(_state->info.pic_y>>iplane->ydec)-_y;
  if(clipy>0){
    h-=clipy;
    _p+=clipy*_pystride;
    _y+=clipy;
  }
  clipw=(_state->info.pic_x+_state->info.pic_width+
   (1<<iplane->xdec)-1>>iplane->xdec)-_x;
  w=OD_MINI(w,clipw);
  cliph=(_state->info.pic_y+_state->info.pic_height+
   (1<<iplane->ydec)-1>>iplane->ydec)-_y;
  h=OD_MINI(h,cliph);
  /*fprintf(stderr,"[%i,%i]x[%i,%i]\n",_x,_y,width,height);*/
  /*Compute the SAD.*/
  src=iplane->data+_y*iplane->ystride+_x*iplane->xstride;
  ret=0;
  for(j=0;j<h;j++){
    p=_p;
    for(i=0;i<w;i++){
      ret+=abs(p[0]-src[i]);
      p+=_pxstride;
    }
    src+=iplane->ystride;
    _p+=_pystride;
  }
  return ret;
}



static void od_mv_est_init(od_mv_est_ctx *_est,od_enc_ctx *_enc){
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  _est->enc=_enc;
  nhmvbs=_enc->state.nhmbs+1<<2;
  nvmvbs=_enc->state.nvmbs+1<<2;
  _est->sad_cache[1]=(od_sad4 **)od_malloc_2d(nvmvbs>>1,nhmvbs>>1,
   sizeof(_est->sad_cache[1][0][0]));
  _est->sad_cache[0]=(od_sad4 **)od_malloc_2d(nvmvbs,nhmvbs,
   sizeof(_est->sad_cache[1][0][0]));
  _est->mvs=(od_mv_node **)od_calloc_2d(nvmvbs+1,nhmvbs+1,
   sizeof(_est->mvs[0][0]));
  _est->refine_grid=(od_mv_grid_pt **)od_malloc_2d(nvmvbs+1,nhmvbs+1,
   sizeof(_est->refine_grid[0][0]));
  _est->dp_nodes=(od_mv_dp_node *)_ogg_malloc(
   sizeof(od_mv_dp_node)*(OD_MAXI(nhmvbs,nvmvbs)+1));
  _est->row_counts=(unsigned *)_ogg_malloc(sizeof(unsigned)*(nvmvbs+1));
  _est->col_counts=(unsigned *)_ogg_malloc(sizeof(unsigned)*(nhmvbs+1));
  for(vy=0;vy<=nvmvbs;vy++)for(vx=0;vx<=nhmvbs;vx++){
    _est->mvs[vy][vx].vx=vx;
    _est->mvs[vy][vx].vy=vy;
    _est->mvs[vy][vx].heapi=-1;
    _enc->state.mv_grid[vy][vx].valid=1;
  }
  _est->dec_heap=(od_mv_node **)_ogg_malloc(
   (nvmvbs+1)*(nhmvbs+1)*sizeof(_est->dec_heap[0]));
  _est->hit_bit=0;
  /*TODO: Allow configuration.*/
  _est->mv_res_min=0;
  _est->flags=OD_MC_USEB|OD_MC_USEV|OD_MC_USE_CHROMA;
  _est->level_max=4;
  _est->level_min=0;
}

static void od_mv_est_clear(od_mv_est_ctx *_est){
  _ogg_free(_est->dec_heap);
  _ogg_free(_est->col_counts);
  _ogg_free(_est->row_counts);
  _ogg_free(_est->dp_nodes);
  od_free_2d(_est->refine_grid);
  od_free_2d(_est->mvs);
  od_free_2d(_est->sad_cache[0]);
  od_free_2d(_est->sad_cache[1]);
}



/*STAGE 1: INITIAL MV ESTIMATES (via EPZS^2).*/



/*The amount to right shift the minimum error by when inflating it for
   computing the second maximum SAD threshold.*/
#define OD_MC_THRESH2_SCALE_BITS (3)

/*The vector offsets in the X direction for each search site in the various
   patterns.*/
static const int OD_SITE_DX[13]={-1,0,1,-1,0,1,-1,0,1,-2,0,2,0};
/*The vector offsets in the Y direction for each search site in the various
   patterns.*/
static const int OD_SITE_DY[13]={-1,-1,-1,0,0,0,1,1,1,0,-2,0,2};

/*The number of sites to search of each boundary condition in the square
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32==dx
  2:      dx==31
  4: -32==dy
  8:      dy==31*/
static const int OD_SQUARE_NSITES[11]={8,5,5,0,5,3,3,0,5,3,3};
/*The list of sites to search for each boundary condition in the square
   pattern.*/
static const od_pattern OD_SQUARE_SITES[11]={
  /* -32<dx<31,   -32<dy<31*/
  {0,1,2,3,5,6,7,8},
  /*-32==dx,      -32<dy<31*/
  {1,2,5,7,8},
  /*     dx==31,  -32<dy<31*/
  {0,1,3,6,7},
  /*-32==dx==31,  -32<dy<31*/
  {-1},
  /* -32<dx<31,  -32==dy*/
  {3,5,6,7,8},
  /*-32==dx,     -32==dy*/
  {5,7,8},
  /*     dx==31, -32==dy*/
  {3,6,7},
  /*-32==dx==31, -32==dy*/
  {-1},
  /* -32<dx<31,       dy==31*/
  {0,1,2,3,5},
  /*-32==dx,          dy==31*/
  {1,2,5},
  /*     dx==31,      dy==31*/
  {0,1,3}
};

/*The number of sites to search of each boundary condition in the diamond
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32==dx
  2:      dx==31
  4: -32==dy
  8:      dy==31*/
static const int OD_DIAMOND_NSITES[11]={4,3,3,0,3,2,2,0,3,2,2};
/*The list of sites to search for each boundary condition in the square
   pattern.*/
static const od_pattern OD_DIAMOND_SITES[11]={
  /* -32<dx<31,   -32<dy<31*/
  {1,3,5,7},
  /*-32==dx,      -32<dy<31*/
  {1,5,7},
  /*     dx==31,  -32<dy<31*/
  {1,3,7},
  /*-32==dx==31,  -32<dy<31*/
  {-1},
  /* -32<dx<31,  -32==dy*/
  {3,5,7},
  /*-32==dx,     -32==dy*/
  {5,7},
  /*     dx==31, -32==dy*/
  {3,7},
  /*-32==dx==31, -32==dy*/
  {-1},
  /* -32<dx<31,       dy==31*/
  {1,3,5},
  /*-32==dx,          dy==31*/
  {1,5},
  /*     dx==31,      dy==31*/
  {1,3}
};

/*The number of sites to search of each boundary condition in the horizontal
   hex pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32==dx
  2:      dx==31
  4: -32==dy
  8:      dy==31*/
static const int OD_HHEX_NSITES[11]={6,3,3,0,4,2,2,0,4,2,2};
/*The list of sites to search for each boundary condition in the horizontal
   hex pattern.*/
static const od_pattern OD_HHEX_SITES[11]={
  /* -32<dx<31,   -32<dy<31*/
  {0,2,6,8,9,11},
  /*-32==dx,      -32<dy<31*/
  {2,8,11},
  /*     dx==31,  -32<dy<31*/
  {0,6,9},
  /*-32==dx==31,  -32<dy<31*/
  {-1},
  /* -32<dx<31,  -32==dy*/
  {6,8,9,11},
  /*-32==dx,     -32==dy*/
  {8,11},
  /*     dx==31, -32==dy*/
  {6,9},
  /*-32==dx==31, -32==dy*/
  {-1},
  /* -32<dx<31,       dy==31*/
  {0,2,9,11},
  /*-32==dx,          dy==31*/
  {2,11},
  /*     dx==31,      dy==31*/
  {0,9}
};

/*The number of sites to search of each boundary condition in the vertical hex
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32==dx
  2:      dx==31
  4: -32==dy
  8:      dy==31*/
static const int OD_VHEX_NSITES[11]={6,4,4,0,3,2,2,0,3,2,2};
/*The list of sites to search for each boundary condition in the vertical hex
   pattern.*/
static const od_pattern OD_VHEX_SITES[11]={
  /* -32<dx<31,   -32<dy<31*/
  {0,2,6,8,10,12},
  /*-32==dx,      -32<dy<31*/
  {2,8,10,12},
  /*     dx==31,  -32<dy<31*/
  {0,6,10,12},
  /*-32==dx==31,  -32<dy<31*/
  {-1},
  /* -32<dx<31,  -32==dy*/
  {6,8,12},
  /*-32==dx,     -32==dy*/
  {8,12},
  /*     dx==31, -32==dy*/
  {6,12},
  /*-32==dx==31, -32==dy*/
  {-1},
  /* -32<dx<31,       dy==31*/
  {0,2,10},
  /*-32==dx,          dy==31*/
  {2,10},
  /*     dx==31,      dy==31*/
  {0,10}
};

#if 0
/*The search state indicating we found a local minimum.*/
#define OD_SEARCH_STATE_DONE (1)

/*The number of sites in the pattern to use for each state.*/
static const int *const OD_SEARCH_NSITES[1]={
  OD_SQUARE_NSITES
};

/*The sites in the pattern to use for each state.*/
static const od_pattern *const OD_SEARCH_SITES[1]={
  OD_SQUARE_SITES
};

/*The successor state given the current state and the terminating site.*/
static const int OD_SEARCH_STATES[1][13]={
  /*Just use a square pattern for the whole search.*/
  { 0, 0, 0, 0, 1, 0, 0, 0, 0,-1,-1,-1,-1}
};
#else
/*The state machine used to select which patterns to use in the gradient
   descent BMA search.
  These are based on the procedure described in
  @ARTICLE{TP06,
    author="Tsung-Han Tsai and Yu-Nan Pan",
    title="A Novel 3-D Predict Hexagon Search Algorithm for Fast Block Motion
     Estimation on H.264 Video Coding",
    journal="{IEEE} Transactions on Circuits and Systems for Video Technology",
    volume=16,
    number=12,
    pages="1542--1549",
    month=Dec,
    year=2006
  }*/

/*The search state indicating we found a local minimum.*/
#define OD_SEARCH_STATE_DONE (6)

/*The number of sites in the pattern to use for each state.*/
static const int *const OD_SEARCH_NSITES[6]={
  OD_DIAMOND_NSITES,
  OD_DIAMOND_NSITES,
  OD_DIAMOND_NSITES,
  OD_HHEX_NSITES,
  OD_VHEX_NSITES,
  OD_DIAMOND_NSITES
};

/*The sites in the pattern to use for each state.*/
static const od_pattern *const OD_SEARCH_SITES[6]={
  OD_DIAMOND_SITES,
  OD_DIAMOND_SITES,
  OD_DIAMOND_SITES,
  OD_HHEX_SITES,
  OD_VHEX_SITES,
  OD_DIAMOND_SITES
};

/*The successor state given the current state and the terminating site.*/
static const int OD_SEARCH_STATES[6][13]={
  /*Start with a small diamond in the first step.*/
  {-1, 2,-1, 1, 6, 1,-1, 2,-1,-1,-1,-1,-1},
  /*Use a small diamond for the second step, too, but remember if we took a
     horizontal step in the first step...*/
  {-1, 3,-1, 3, 6, 3,-1, 3,-1,-1,-1,-1,-1},
  /*...or a vertical one.*/
  {-1, 4,-1, 4, 6, 4,-1, 4,-1,-1,-1,-1,-1},
  /*Then switch to a horizontal hex pattern for all remaining steps...*/
  { 3,-1, 3,-1, 5,-1, 3,-1, 3, 3,-1, 3,-1},
  /*...or a vertical hex pattern, depending.*/
  { 4,-1, 4,-1, 5,-1, 4,-1, 4,-1, 4,-1, 4},
  /*And revert back to a small diamond for the last step.*/
  {-1, 6,-1, 6, 6, 6,-1, 6,-1,-1,-1,-1,-1}
};
#endif



/*Clear the cache of motion vectors we've examined.*/
static void od_mv_est_clear_hit_cache(od_mv_est_ctx *_est){
  if(_est->hit_bit++==0)memset(_est->hit_cache,0,sizeof(_est->hit_cache));
  else _est->hit_bit&=UCHAR_MAX;
}

/*Test if a motion vector has been examined.*/
static int od_mv_est_is_hit(od_mv_est_ctx *_est,int _mvx,int _mvy){
  return _est->hit_cache[_mvy+32][_mvx+32]==_est->hit_bit;
}

/*Mark a motion vector examined.*/
static void od_mv_est_set_hit(od_mv_est_ctx *_est,int _mvx,int _mvy){
  _est->hit_cache[_mvy+32][_mvx+32]=(unsigned char)_est->hit_bit;
}

/*Gets the predictor for a given MV node at the given MV resolution.*/
static void od_state_get_predictor(od_state *_state,int _pred[2],
 int _vx,int _vy,int _level,int _mv_res){
  int nhmvbs;
  int nvmvbs;
  nhmvbs=_state->nhmbs+1<<2;
  nvmvbs=_state->nvmbs+1<<2;
  if(_vx<2||_vy<2||_vx>nhmvbs-2||_vy>nvmvbs-2)_pred[0]=_pred[1]=0;
  else{
    od_mv_grid_pt *cneighbors[4];
    int            a[4][2];
    int            mvb_sz;
    int            ncns;
    int            ci;
    mvb_sz=1<<(4-_level>>1);
    ncns=4;
    if(_level==0){
      cneighbors[0]=_state->mv_grid[_vy-4]+_vx-4;
      cneighbors[1]=_state->mv_grid[_vy-4]+_vx;
      cneighbors[2]=_state->mv_grid[_vy-4]+_vx+4;
      cneighbors[3]=_state->mv_grid[_vy]+_vx-4;
    }
    else{
      if(_level&1){
        cneighbors[0]=_state->mv_grid[_vy-mvb_sz]+_vx-mvb_sz;
        cneighbors[1]=_state->mv_grid[_vy-mvb_sz]+_vx+mvb_sz;
        cneighbors[2]=_state->mv_grid[_vy+mvb_sz]+_vx-mvb_sz;
        cneighbors[3]=_state->mv_grid[_vy+mvb_sz]+_vx+mvb_sz;
      }
      else{
        cneighbors[0]=_state->mv_grid[_vy-mvb_sz]+_vx;
        cneighbors[1]=_state->mv_grid[_vy]+_vx-mvb_sz;
        /*NOTE: Only one of these candidates can be excluded at a time, so
           there will always be at least 3.*/
        if(_vx+mvb_sz>(_vx+3&~3))ncns--;
        else cneighbors[2]=_state->mv_grid[_vy]+_vx+mvb_sz;
        if(_vy+mvb_sz>(_vy+3&~3))ncns--;
        else cneighbors[ncns-1]=_state->mv_grid[_vy+mvb_sz]+_vx;
      }
    }
    for(ci=0;ci<ncns;ci++){
      a[ci][0]=cneighbors[ci]->mv[0];
      a[ci][1]=cneighbors[ci]->mv[1];
    }
    /*Median-of-4.*/
    if(ncns>3){
      /*fprintf(stderr,"Median of 4:\n");
      fprintf(stderr,"(%i,%i) (%i,%i) (%i,%i) (%i,%i)\n",
       a[0][0],a[0][1],a[1][0],a[1][1],a[2][0],a[2][1],a[3][0],a[3][1]);*/
/*
Sorting network for 4 elements:
0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
0001 0010 0011 0100 0101 0110 0111 1001 1010 1011 1101
0:1
0010 0010 0011 0100 0110 0110 0111 1010 1010 1011 1110
0010 0011 0100 0110 0111 1010 1011
2:3
0010 0011 1000 1010 1011 1100 1110
0010 0011 1010 1011
0:2
0010 0110 1010 1110
0010 0110 1010
1:3
1000 1100 1010
1010
This last compare is unneeded for a median:
1:2
1100
*/
      OD_SORT2I(a[0][0],a[1][0]);
      OD_SORT2I(a[0][1],a[1][1]);
      OD_SORT2I(a[2][0],a[3][0]);
      OD_SORT2I(a[2][1],a[3][1]);
      OD_SORT2I(a[0][0],a[2][0]);
      OD_SORT2I(a[0][1],a[2][1]);
      OD_SORT2I(a[1][0],a[3][0]);
      OD_SORT2I(a[1][1],a[3][1]);
      /*fprintf(stderr,"(%i,%i) (%i,%i) (%i,%i) (%i,%i)\n",
       a[0][0],a[0][1],a[1][0],a[1][1],a[2][0],a[2][1],a[3][0],a[3][1]);*/
      _pred[0]=OD_DIV_POW2_RE(a[1][0]+a[2][0],_mv_res+1);
      _pred[1]=OD_DIV_POW2_RE(a[1][1]+a[2][1],_mv_res+1);
    }
    /*Median-of-3.*/
    else{
      /*fprintf(stderr,"Median of 3:\n");
      fprintf(stderr,"(%i,%i) (%i,%i) (%i,%i)\n",
       a[0][0],a[0][1],a[1][0],a[1][1],a[2][0],a[2][1]);*/
      OD_SORT2I(a[0][0],a[1][0]);
      OD_SORT2I(a[0][1],a[1][1]);
      OD_SORT2I(a[1][0],a[2][0]);
      OD_SORT2I(a[1][1],a[2][1]);
      OD_SORT2I(a[0][0],a[1][0]);
      OD_SORT2I(a[0][1],a[1][1]);
      /*fprintf(stderr,"(%i,%i) (%i,%i) (%i,%i)\n",
       a[0][0],a[0][1],a[1][0],a[1][1],a[2][0],a[2][1]);*/
      _pred[0]=OD_DIV_POW2_RE(a[1][0],_mv_res);
      _pred[1]=OD_DIV_POW2_RE(a[1][1],_mv_res);
    }
  }
}

static const int OD_MV_EST_RATE[256]={
  1, 4, 4, 6, 6, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8,10,
 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,12,
 12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,
 12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,
 14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
 14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
 14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
 14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,18
};

/*Estimate the number of bits that will be used to encode the given MV and its
   predictor.*/
static int od_mv_est_bits(int _dx,int _dy,int _predx,int _predy){
  int dx;
  int dy;
  dx=OD_MINI(abs(_dx-_predx),255);
  dy=OD_MINI(abs(_dy-_predy),255);
  _dx=OD_MINI(abs(_dx),255);
  _dy=OD_MINI(abs(_dy),255);
  return 1+OD_MINI(OD_MV_EST_RATE[_dx]+OD_MV_EST_RATE[_dy],
   OD_MV_EST_RATE[dx]+OD_MV_EST_RATE[dy]);
}

/*Computes the SAD of a whole-pel BMA block with the given parameters.*/
static ogg_int32_t od_mv_est_bma_sad8(od_mv_est_ctx *_est,int _ref,int _bx,
 int _by,int _mvx,int _mvy,int _log_mvb_sz){
  od_state     *state;
  od_img_plane *iplane;
  ogg_int32_t   ret;
  int           refi;
  int           mvx;
  int           mvy;
  int           bx;
  int           by;
  int           dx;
  int           dy;
  state=&_est->enc->state;
  refi=state->ref_imgi[_ref];
  iplane=state->ref_imgs[refi].planes+0;
  mvx=OD_DIV_POW2_RE(_mvx<<1,iplane->xdec);
  mvy=OD_DIV_POW2_RE(_mvy<<1,iplane->ydec);
  bx=_bx+(1<<iplane->xdec)-1&~((1<<iplane->xdec)-1);
  by=_by+(1<<iplane->ydec)-1&~((1<<iplane->ydec)-1);
  dx=(bx<<1>>iplane->xdec)+mvx;
  dy=(by<<1>>iplane->ydec)+mvy;
  ret=od_state_sad8(state,iplane->data+dy*iplane->ystride+dx,
   iplane->ystride<<1,2,0,bx,by,_log_mvb_sz+2);
  if(_est->flags&OD_MC_USE_CHROMA){
    int pli;
    for(pli=1;pli<state->input.nplanes;pli++){
      iplane=state->ref_imgs[refi].planes+pli;
      mvx=OD_DIV_POW2_RE(_mvx<<1,iplane->xdec);
      mvy=OD_DIV_POW2_RE(_mvy<<1,iplane->ydec);
      bx=_bx+(1<<iplane->xdec)-1&~((1<<iplane->xdec)-1);
      by=_by+(1<<iplane->ydec)-1&~((1<<iplane->ydec)-1);
      dx=(bx<<1>>iplane->xdec)+mvx;
      dy=(by<<1>>iplane->ydec)+mvy;
      ret+=od_state_sad8(state,iplane->data+dy*iplane->ystride+dx,
       iplane->ystride<<1,2,pli,bx,by,_log_mvb_sz+2)>>OD_MC_CHROMA_SCALE;
    }
  }
  return ret;
}

/*Computes the SAD of a block with the given parameters.*/
static ogg_int32_t od_mv_est_sad8(od_mv_est_ctx *_est,int _ref,
 int _vx,int _vy,int _c,int _s,int _log_mvb_sz){
  unsigned char __attribute((aligned(16))) pred[16][16];
  od_state      *state;
  ogg_int32_t    ret;
  state=&_est->enc->state;
  od_state_pred_block_from_setup(state,pred[0],sizeof(pred[0]),_ref,0,
   _vx,_vy,_c,_s,_log_mvb_sz);
  ret=od_state_sad8(state,pred[0],sizeof(pred[0]),1,0,_vx-2<<2,_vy-2<<2,
   _log_mvb_sz+2);
  if(_est->flags&OD_MC_USE_CHROMA){
    int pli;
    for(pli=1;pli<state->input.nplanes;pli++){
      od_state_pred_block_from_setup(state,pred[0],sizeof(pred[0]),_ref,pli,
       _vx,_vy,_c,_s,_log_mvb_sz);
      ret+=od_state_sad8(state,pred[0],sizeof(pred[0]),1,pli,_vx-2<<2,_vy-2<<2,
       _log_mvb_sz+2)>>OD_MC_CHROMA_SCALE;
    }
  }
  return ret;
}

/*Checks to make sure our current mv_rate and sad values are correct.
  This is used for debugging only.*/
void od_mv_est_check_rd_block_state(od_mv_est_ctx *_est,int _ref,
 int _vx,int _vy,int _log_mvb_sz){
  od_state      *state;
  int half_mvb_sz;
  state=&_est->enc->state;
  half_mvb_sz=1<<_log_mvb_sz-1;
  if(_log_mvb_sz>0&&state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
    od_mv_est_check_rd_block_state(_est,_ref,_vx,_vy,_log_mvb_sz-1);
    od_mv_est_check_rd_block_state(_est,_ref,
     _vx+half_mvb_sz,_vy,_log_mvb_sz-1);
    od_mv_est_check_rd_block_state(_est,_ref,
     _vx,_vy+half_mvb_sz,_log_mvb_sz-1);
    od_mv_est_check_rd_block_state(_est,_ref,
     _vx+half_mvb_sz,_vy+half_mvb_sz,_log_mvb_sz-1);
  }
  else{
    od_mv_node  *block;
    ogg_int32_t  sad;
    int          c;
    int          s;
    block=_est->mvs[_vy]+_vx;
    if(block->log_mvb_sz!=_log_mvb_sz){
      fprintf(stderr,
       "Failure at node (%i,%i): log_mvb_sz should be %i (is %i)\n",
       _vx,_vy,_log_mvb_sz,block->log_mvb_sz);
    }
    if(_log_mvb_sz<2){
      int mask;
      mask=(1<<_log_mvb_sz+1)-1;
      c=!!(_vx&mask);
      if(_vy&mask)c=3-c;
      if(block->c!=c){
        fprintf(stderr,"Failure at node (%i,%i): c should be %i (is %i)\n",
         _vx,_vy,c,block->c);
      }
      s=state->mv_grid[_vy+(OD_VERT_DY[c+1&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+1&3]<<_log_mvb_sz)].valid|
       state->mv_grid[_vy+(OD_VERT_DY[c+3&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+3&3]<<_log_mvb_sz)].valid<<1;
    }
    else{
      c=0;
      s=3;
    }
    if(block->s!=s){
      fprintf(stderr,"Failure at node (%i,%i): s should be %i (is %i)\n",
       _vx,_vy,s,block->s);
    }
    sad=od_mv_est_sad8(_est,_ref,_vx,_vy,c,s,_log_mvb_sz);
    if(block->sad!=sad){
      fprintf(stderr,"Failure at node (%i,%i): sad should be %i (is %i)\n",
       _vx,_vy,sad,block->sad);
    }
  }
}

/*Checks to make sure our current mv_rate and sad values are correct.
  This is used for debugging only.*/
void od_mv_est_check_rd_state(od_mv_est_ctx *_est,int _ref,int _mv_res){
  od_state *state;
  int       nhmvbs;
  int       nvmvbs;
  int       vx;
  int       vy;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  for(vy=0;vy<nvmvbs;vy+=4){
    for(vx=0;vx<nhmvbs;vx+=4){
      od_mv_est_check_rd_block_state(_est,_ref,vx,vy,2);
    }
  }
  for(vy=0;vy<nvmvbs;vy++)for(vx=0;vx<nhmvbs;vx++){
    od_mv_grid_pt *mvg;
    od_mv_node    *mv;
    int            pred[2];
    int            mv_rate;
    mvg=state->mv_grid[vy]+vx;
    if(!mvg->valid)continue;
    mv=_est->mvs[vy]+vx;
    if(vx>=2&&vx<=nhmvbs-2&&vy>=2&&vy<=nvmvbs-2){
      od_state_get_predictor(state,pred,vx,vy,OD_MC_LEVEL[vy&3][vx&3],_mv_res);
      mv_rate=od_mv_est_bits(mvg->mv[0]>>_mv_res,mvg->mv[1]>>_mv_res,
       pred[0],pred[1]);
    }
    else pred[0]=pred[1]=mv_rate=0;
    if(mv_rate!=mv->mv_rate){
      fprintf(stderr,"Failure at node (%i,%i): mv_rate should be %i (is %i)\n",
       vx,vy,mv_rate,mv->mv_rate);
      fprintf(stderr,"Predictor was: (%i,%i)   MV was: (%i,%i)\n",
       pred[0],pred[1],mvg->mv[0]>>_mv_res,mvg->mv[1]>>_mv_res);
    }
  }
}

#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
static const unsigned char OD_YCbCr_MVCAND[3]={210, 16,214};
#endif

static void od_mv_est_init_mv(od_mv_est_ctx *_est,int _ref,int _vx,int _vy){
  od_state      *state;
  od_mv_grid_pt *mvg;
  od_mv_node    *mv;
  od_mv_node    *cneighbors[4];
  od_mv_node    *pneighbors[4];
  ogg_int32_t    t2;
  ogg_int32_t    best_sad;
  ogg_int32_t    best_cost;
  int            best_rate;
  int            cands[6][2];
  int            best_vec[2];
  int            a[4][2];
  int            level;
  int            log_mvb_sz;
  int            mvb_sz;
  int            bx;
  int            by;
  int            ncns;
  int            mvxmin;
  int            mvxmax;
  int            mvymin;
  int            mvymax;
  int            candx;
  int            candy;
  int            predx;
  int            predy;
  int            ci;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
  int            animating;
  int            x0;
  int            y0;
#endif
  /*fprintf(stderr,"Initial search for MV (%i,%i):\n",_vx,_vy);*/
  state=&_est->enc->state;
  mv=_est->mvs[_vy]+_vx;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  mvg=state->mv_grid[_vy]+_vx;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
  animating=daala_granule_basetime(state,state->cur_time)==ANI_FRAME;
  if(animating){
    od_state_mc_predict(state,_ref);
    od_state_fill_vis(state);
    x0=(_vx-2<<3)+(OD_UMV_PADDING<<1);
    y0=(_vy-2<<3)+(OD_UMV_PADDING<<1);
  }
#endif
  /*fprintf(stderr,"Level %i (%ix%i block)\n",level,mvb_sz<<2,mvb_sz<<2);*/
  bx=_vx-2<<2;
  by=_vy-2<<2;
  mvxmin=OD_MAXI(bx-(mvb_sz<<2)-32,-16)-(bx-(mvb_sz<<2));
  mvxmax=OD_MINI(bx+(mvb_sz<<2)+32,state->info.frame_width+16)-
   (bx+(mvb_sz<<2))-1;
  mvymin=OD_MAXI(by-(mvb_sz<<2)-32,-16)-(by-(mvb_sz<<2));
  mvymax=OD_MINI(by+(mvb_sz<<2)+32,state->info.frame_height+16)-
   (by+(mvb_sz<<2))-1;
  /*fprintf(stderr,"(%i,%i): Search range: [%i,%i]x[%i,%i]\n",
   bx,by,mvxmin,mvymin,mvxmax,mvymax);*/
  bx-=mvb_sz<<1;
  by-=mvb_sz<<1;
  ncns=4;
  if(level==0){
    cneighbors[0]=_est->mvs[_vy-4]+_vx-4;
    cneighbors[1]=_est->mvs[_vy-4]+_vx;
    cneighbors[2]=_est->mvs[_vy-4]+_vx+4;
    cneighbors[3]=_est->mvs[_vy]+_vx-4;
    pneighbors[0]=_est->mvs[_vy-4]+_vx;
    pneighbors[1]=_est->mvs[_vy]+_vx-4;
    pneighbors[2]=_est->mvs[_vy]+_vx+4;
    pneighbors[3]=_est->mvs[_vy+4]+_vx;
  }
  else{
    if(level&1){
      pneighbors[0]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      pneighbors[1]=_est->mvs[_vy-mvb_sz]+_vx+mvb_sz;
      pneighbors[2]=_est->mvs[_vy+mvb_sz]+_vx-mvb_sz;
      pneighbors[3]=_est->mvs[_vy+mvb_sz]+_vx+mvb_sz;
      memcpy(cneighbors,pneighbors,sizeof(cneighbors));
    }
    else{
      pneighbors[0]=_est->mvs[_vy-mvb_sz]+_vx;
      pneighbors[1]=_est->mvs[_vy]+_vx-mvb_sz;
      pneighbors[2]=_est->mvs[_vy]+_vx+mvb_sz;
      pneighbors[3]=_est->mvs[_vy+mvb_sz]+_vx;
      cneighbors[0]=pneighbors[0];
      cneighbors[1]=pneighbors[1];
      /*NOTE: Only one of these candidatss can be excluded at a time, so
         there will always be at least 3.*/
      if(_vx+mvb_sz>(_vx+3&~3))ncns--;
      else cneighbors[2]=pneighbors[2];
      if(_vy+mvb_sz>(_vy+3&~3))ncns--;
      else cneighbors[ncns-1]=pneighbors[3];
    }
  }
  /*Spatially correlated predictors (from the current frame):*/
  for(ci=0;ci<ncns;ci++){
    a[ci][0]=cneighbors[ci]->mvs[0][_ref][0];
    a[ci][1]=cneighbors[ci]->mvs[0][_ref][1];
    cands[ci][0]=OD_CLAMPI(mvxmin,a[ci][0],mvxmax);
    cands[ci][1]=OD_CLAMPI(mvymin,a[ci][1],mvymax);
  }
  /*Compute the median predictor:*/
  if(ncns>3){
    /*Median-of-4.*/
    OD_SORT2I(a[0][0],a[1][0]);
    OD_SORT2I(a[0][1],a[1][1]);
    OD_SORT2I(a[2][0],a[3][0]);
    OD_SORT2I(a[2][1],a[3][1]);
    OD_SORT2I(a[0][0],a[2][0]);
    OD_SORT2I(a[0][1],a[2][1]);
    OD_SORT2I(a[1][0],a[3][0]);
    OD_SORT2I(a[1][1],a[3][1]);
    predx=a[1][0]+a[2][0];
    predy=a[1][1]+a[2][1];
    candx=OD_CLAMPI(mvxmin,OD_DIV2(predx),mvxmax);
    candy=OD_CLAMPI(mvymin,OD_DIV2(predy),mvymax);
  }
  else{
    /*Median-of-3.*/
    OD_SORT2I(a[0][0],a[1][0]);
    OD_SORT2I(a[0][1],a[1][1]);
    OD_SORT2I(a[1][0],a[2][0]);
    OD_SORT2I(a[1][1],a[2][1]);
    OD_SORT2I(a[0][0],a[1][0]);
    OD_SORT2I(a[0][1],a[1][1]);
    predx=a[1][0]<<1;
    predy=a[1][1]<<1;
    candx=OD_CLAMPI(mvxmin,a[1][0],mvxmax);
    candy=OD_CLAMPI(mvymin,a[1][1],mvymax);
  }
  od_mv_est_clear_hit_cache(_est);
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
  if(animating){
    od_img_draw_line(&state->vis_img,x0,y0,x0+(candx<<1),y0+(candy<<1),
     OD_YCbCr_MVCAND);
  }
#endif
  best_sad=od_mv_est_bma_sad8(_est,_ref,bx,by,candx,candy,log_mvb_sz);
  best_rate=od_mv_est_bits(candx<<1,candy<<1,predx,predy);
  best_cost=(best_sad<<OD_LAMBDA_SCALE)+best_rate*_est->lambda;
  /*fprintf(stderr,"Median predictor: (%i,%i)   Error: %i\n",candx,candy,best_err);*/
  od_mv_est_set_hit(_est,candx,candy);
  best_vec[0]=candx;
  best_vec[1]=candy;
  /*fprintf(stderr,"Threshold: %i\n",_est->thresh1[log_mvb_sz]);*/
  if(best_sad>_est->thresh1[log_mvb_sz]){
    ogg_int32_t sad;
    ogg_int32_t cost;
    int         rate;
    /*Compute the early termination threshold for set B.*/
    t2=mv->bma_sad;
    for(ci=0;ci<ncns;ci++){
      int log_cnb_sz;
      log_cnb_sz=4-OD_MC_LEVEL[cneighbors[ci]->vy&3][cneighbors[ci]->vx&3]>>1;
      t2=OD_MINI(t2,cneighbors[ci]->bma_sad>>(log_cnb_sz-log_mvb_sz<<1));
    }
    t2=t2+(t2>>OD_MC_THRESH2_SCALE_BITS)+_est->thresh2_offs[log_mvb_sz];
    /*Constant velocity predictor:*/
    cands[ncns][0]=OD_CLAMPI(mvxmin,OD_DIV8(mv->mvs[1][_ref][0]),mvxmax);
    cands[ncns][1]=OD_CLAMPI(mvymin,OD_DIV8(mv->mvs[1][_ref][1]),mvymax);
    ncns++;
    /*Zero predictor.*/
    cands[ncns][0]=0;
    cands[ncns][1]=0;
    ncns++;
    /*Examine the candidates in Set B.*/
    for(ci=0;ci<ncns;ci++){
      candx=cands[ci][0];
      candy=cands[ci][1];
      /*fprintf(stderr,"Set B predictor %i: (%i,%i) ",ci,candx,candy);*/
      if(od_mv_est_is_hit(_est,candx,candy)){/*fprintf(stderr,"...Skipping.\n");*/continue;}
      od_mv_est_set_hit(_est,candx,candy);
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
      if(animating){
        od_img_draw_line(&state->vis_img,x0,y0,x0+(candx<<1),y0+(candy<<1),
         OD_YCbCr_MVCAND);
      }
#endif
      sad=od_mv_est_bma_sad8(_est,_ref,bx,by,candx,candy,log_mvb_sz);
      rate=od_mv_est_bits(candx<<1,candy<<1,predx,predy);
      cost=(sad<<OD_LAMBDA_SCALE)+rate*_est->lambda;
      /*fprintf(stderr,"   Error: %i\n",err);*/
      if(cost<best_cost){
        best_sad=sad;
        best_rate=rate;
        best_cost=cost;
        best_vec[0]=candx;
        best_vec[1]=candy;
      }
    }
    /*fprintf(stderr,"Threshold: %i\n",t2);*/
    if(best_sad>t2){
      /*Constant velocity predictors from the previous frame:*/
      for(ci=0;ci<4;ci++){
        cands[ci][0]=
         OD_CLAMPI(mvxmin,OD_DIV8(pneighbors[ci]->mvs[1][_ref][0]),mvxmax);
        cands[ci][1]=
         OD_CLAMPI(mvymin,OD_DIV8(pneighbors[ci]->mvs[1][_ref][1]),mvymax);
      }
      /*The constant acceleration predictor:*/
      cands[4][0]=OD_CLAMPI(mvxmin,OD_DIV_ROUND_POW2(
       mv->mvs[1][_ref][0]*_est->mvapw[_ref][0]-
       mv->mvs[2][_ref][0]*_est->mvapw[_ref][1],16,0x8000),mvxmax);
      cands[4][1]=OD_CLAMPI(mvymin,OD_DIV_ROUND_POW2(
       mv->mvs[1][_ref][1]*_est->mvapw[_ref][0]-
       mv->mvs[2][_ref][1]*_est->mvapw[_ref][1],16,0x8000),mvymax);
      /*Examine the candidates in Set C.*/
      for(ci=0;ci<5;ci++){
        candx=cands[ci][0];
        candy=cands[ci][1];
        /*fprintf(stderr,"Set C predictor %i: (%i,%i) ",ci,candx,candy);*/
        if(od_mv_est_is_hit(_est,candx,candy)){/*fprintf(stderr,"...Skipping.\n");*/continue;}
        /*if(od_mv_est_is_hit(_est,candx,candy))continue;*/
        od_mv_est_set_hit(_est,candx,candy);
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
        if(animating){
          od_img_draw_line(&state->vis_img,x0,y0,x0+(candx<<1),y0+(candy<<1),
           OD_YCbCr_MVCAND);
        }
#endif
        sad=od_mv_est_bma_sad8(_est,_ref,bx,by,candx,candy,log_mvb_sz);
        rate=od_mv_est_bits(candx<<1,candy<<1,predx,predy);
        cost=(sad<<OD_LAMBDA_SCALE)+rate*_est->lambda;
        /*fprintf(stderr,"   Error: %i\n",err);*/
        if(cost<best_cost){
          best_sad=sad;
          best_rate=rate;
          best_cost=cost;
          best_vec[0]=candx;
          best_vec[1]=candy;
        }
      }
      /*Use the same threshold for Set C as in Set B.*/
      /*fprintf(stderr,"Threshold: %i\n",t2);*/
      if(best_sad>t2){
        const int *pattern;
        int        mvstate;
        int        best_site;
        int        nsites;
        int        sitei;
        int        site;
        int        b;
        /*Gradient descent pattern search.*/
        mvstate=0;
        for(;;){
          best_site=4;
          b=(best_vec[0]<=mvxmin)|(best_vec[0]>=mvxmax)<<1|
           (best_vec[1]<=mvymin)<<2|(best_vec[1]>=mvymax)<<3;
          pattern=OD_SEARCH_SITES[mvstate][b];
          nsites=OD_SEARCH_NSITES[mvstate][b];
          for(sitei=0;sitei<nsites;sitei++){
            site=pattern[sitei];
            candx=best_vec[0]+OD_SITE_DX[site];
            candy=best_vec[1]+OD_SITE_DY[site];
            /*For the large search patterns, our simple mechanism to move
               bounds checking out of the inner loop doesn't work (it would
               need 2 more bits, or 4 times as much table storage, and require
               4 extra compares, when there are often fewer than 4 sites).
              If the displacement is larger than +/-1 in any direction (which
               happens when site>8), check the bounds explicitly.*/
            if(site>8&&(candx<mvxmin||candx>mvxmax||
             candy<mvymin||candy>mvymax)){
              continue;
            }
            /*fprintf(stderr,"Pattern search %i: (%i,%i) ",site,candx,candy);*/
            if(od_mv_est_is_hit(_est,candx,candy)){/*fprintf(stderr,"...Skipping.\n");*/continue;}
            od_mv_est_set_hit(_est,candx,candy);
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
            if(animating){
              od_img_draw_line(&state->vis_img,x0,y0,
               x0+(candx<<1),y0+(candy<<1),OD_YCbCr_MVCAND);
            }
#endif
            sad=od_mv_est_bma_sad8(_est,_ref,bx,by,candx,candy,log_mvb_sz);
            rate=od_mv_est_bits(candx<<1,candy<<1,predx,predy);
            cost=(sad<<OD_LAMBDA_SCALE)+rate*_est->lambda;
            /*fprintf(stderr,"   Error: %i\n",err);*/
            if(cost<best_cost){
              best_sad=sad;
              best_rate=rate;
              best_cost=cost;
              best_site=site;
            }
          }
          mvstate=OD_SEARCH_STATES[mvstate][best_site];
          best_vec[0]+=OD_SITE_DX[best_site];
          best_vec[1]+=OD_SITE_DY[best_site];
          if(mvstate==OD_SEARCH_STATE_DONE)break;
        }
      }
    }
  }
  /*fprintf(stderr,"Finished. Best vector: (%i,%i)  Best error %i\n",
   best_vec[0],best_vec[1],best_err);*/
  mv->mvs[0][_ref][0]=best_vec[0];
  mv->mvs[0][_ref][1]=best_vec[1];
  mvg->mv[0]=best_vec[0]<<3;
  mvg->mv[1]=best_vec[1]<<3;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
  if(animating){
    char             iter_label[16];
    const od_offset *anc;
    od_mv_grid_pt   *amvg;
    int              nanc;
    int              ai;
    int              ax;
    int              ay;
    mvg->valid=1;
    nanc=OD_NANCESTORS[_vy&3][_vx&3];
    anc=OD_ANCESTORS[_vy&3][_vx&3];
    for(ai=0;ai<nanc;ai++){
      ax=_vx+anc[ai][0];
      if(ax<0||ax>(state->nhmbs+1<<2))continue;
      ay=_vy+anc[ai][1];
      if(ay<0||ay>(state->nvmbs+1<<2))continue;
      amvg=state->mv_grid[ay]+ax;
      amvg->valid=1;
    }
    sprintf(iter_label,"ani%08i",state->ani_iter++);
    od_state_dump_img(state,&state->vis_img,iter_label);
  }
#endif
  mv->bma_sad=best_sad;
  mv->mv_rate=best_rate;
  /*fprintf(stderr,"Initialized MV (%2i,%2i): (%3i,%3i), SAD: %i\n",
   _vx,_vy,best_vec[0],best_vec[1],best_sad);*/
  /*od_state_get_predictor(state,a[0],_vx,_vy,level,2);
  if(a[0][0]!=predx||a[0][1]!=predy){
    fprintf(stderr,"Failure in MV predictor init: (%i,%i)!=(%i,%i)\n",
     a[0][0],a[0][1],predx,predy);
  }
  mv->mv_rate=od_mv_est_bits(mvg->mv[0]>>2,mvg->mv[1]>>2,a[0][0],a[0][1]);
  if(mv->mv_rate!=best_rate){
    fprintf(stderr,"Failure in MV rate init: %i!=%i\n",mv->mv_rate,best_rate);
  }*/
}

static void od_mv_est_init_mvs(od_mv_est_ctx *_est,int _ref){
  od_state *state;
  int       nhmvbs;
  int       nvmvbs;
  int       vx;
  int       vy;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  /*Move the motion vector predictors back a frame.*/
  for(vy=2;vy<=nvmvbs-2;vy++)for(vx=2;vx<=nhmvbs-2;vx++){
    od_mv_node *mv;
    mv=_est->mvs[vy]+vx;
    memmove(mv->mvs+1,mv->mvs+0,sizeof(mv->mvs[0])<<1);
  }
  /*We initialize MVs a MVB at a time for cache coherency.
    Proceeding level-by-level would involve less branching and less complex
     code, but the SADs dominate.*/
  for(vy=4;vy<=nvmvbs;vy+=4){
    for(vx=4;vx<=nhmvbs;vx+=4){
      int b;
      /*Keep track of what edges we're on.
        We only need to know about the bottom/right to start with, since the
         top/left are already initialized (or do not need to be).*/
      b=(vx<nhmvbs)<<1|(vy<nvmvbs)<<3;
      /*Initialization MUST proceed in order by level to ensure the necessary
         predictors are available.
        Order within a level does not matter, except for level 0.
        Level 0 is the only level to use predictors outside the current MVB,
         and must proceed in raster order.*/
      /*Level 0 vertex.*/
      if((b&0xA)==0xA)od_mv_est_init_mv(_est,_ref,vx,vy);
      if(_est->level_max<1)continue;
      /*Level 1 vertex.*/
      od_mv_est_init_mv(_est,_ref,vx-2,vy-2);
      if(_est->level_max<2)continue;
      /*Level 2 vertices.*/
      if(b&2)od_mv_est_init_mv(_est,_ref,vx,vy-2);
      if(b&8)od_mv_est_init_mv(_est,_ref,vx-2,vy);
      if(_est->level_max<3)continue;
      /*Level 3 vertices.*/
      /*Add in flags for the top/left edges.*/
      b|=(vx>4)|(vy>4)<<2;
      if(b&4){
        if(b&1)od_mv_est_init_mv(_est,_ref,vx-3,vy-3);
        if(b&2)od_mv_est_init_mv(_est,_ref,vx-1,vy-3);
      }
      if(b&8){
        if(b&1)od_mv_est_init_mv(_est,_ref,vx-3,vy-1);
        if(b&2)od_mv_est_init_mv(_est,_ref,vx-1,vy-1);
      }
      if(_est->level_max<4)continue;
      /*Level 4 vertices.*/
      if(b&1)od_mv_est_init_mv(_est,_ref,vx-3,vy-2);
      if(b&2)od_mv_est_init_mv(_est,_ref,vx-1,vy-2);
      if(b&4){
        od_mv_est_init_mv(_est,_ref,vx-2,vy-3);
        if(b&2)od_mv_est_init_mv(_est,_ref,vx,vy-3);
      }
      if(b&8){
        od_mv_est_init_mv(_est,_ref,vx-2,vy-1);
        if(b&1)od_mv_est_init_mv(_est,_ref,vx-3,vy);
        if(b&2){
          od_mv_est_init_mv(_est,_ref,vx,vy-1);
          od_mv_est_init_mv(_est,_ref,vx-1,vy);
        }
      }
    }
  }
}


/*STAGE 2: DECIMATION.*/



/*Merging domains.
  These are stored as lists of offsets to the vertices in the domain.
  Note that vertices in the merging domain must appear in order from finest
   scale (largest level) to coarsest (smallest level).
  Each list ends with the vertex (0,0), the actual vertex be decimated.*/
/*Level 4 vertex:
            4
*/
static const od_offset OD_MERGEDOM4[1]={
  {0,0},
};

/*Level 3 vertex:
            4
          4-3-4
            4
*/
static const od_offset OD_MERGEDOM3[5]={
  { 0,-1},{-1, 0},{ 1, 0},{ 0, 1},{ 0, 0}
};

/*Level 2 vertex:
          4   4
          |   |
        4-3-4-3-4
          | | |
          4-2-4
          | | |
        4-3-4-3-4
          |   |
          4   4
*/
static const od_offset OD_MERGEDOM2[17]={
  {-1,-2},{ 1,-2},{-2,-1},{ 0,-1},{ 2,-1},{-1, 0},{ 1, 0},{-2, 1},
  { 0, 1},{ 2, 1},{-1, 2},{ 1, 2},{-1,-1},{ 1,-1},{-1, 1},{ 1, 1},
  { 0, 0}
};

/*Level 1 vertex:
          4   4
          |   |
        4-3-4-3-4
          | | |
      4   4-2-4   4
      |   | | |   |
    4-3-4-3-4-3-4-3-4
      | | | | | | |
      4-2-4-1-4-2-4
      | | | | | | |
    4-3-4-3-4-3-4-3-4
      |   | | |   |
      4   4-2-4   4
          | | |
        4-3-4-3-4
          |   |
          4   4
*/
static const od_offset OD_MERGEDOM1[49]={
  {-1,-4},{ 1,-4},{-2,-3},{ 0,-3},{ 2,-3},{-3,-2},{-1,-2},{ 1,-2},
  { 3,-2},{-4,-1},{-2,-1},{ 0,-1},{ 2,-1},{ 4,-1},{-3, 0},{-1, 0},
  { 1, 0},{ 3, 0},{-4, 1},{-2, 1},{ 0, 1},{ 2, 1},{ 4, 1},{-3, 2},
  {-1, 2},{ 1, 2},{ 3, 2},{-2, 3},{ 0, 3},{ 2, 3},{-1, 4},{ 1, 4},
  {-1,-3},{ 1,-3},{-3,-1},{-1,-1},{ 1,-1},{ 3,-1},{-3, 1},{-1, 1},
  { 1, 1},{ 3, 1},{-1, 3},{ 1, 3},{ 0,-2},{-2, 0},{ 2, 0},{ 0, 2},
  { 0, 0}
};

/*The merging domain for a vertex, indexed by level-1.*/
static const od_offset *OD_MERGEDOM[4]={
  OD_MERGEDOM1,
  OD_MERGEDOM2,
  OD_MERGEDOM3,
  OD_MERGEDOM4
};

/*Error support regions.
  These are the blocks whose SAD will change after decimating a vertex at a
   given level, assuming no other vertices in the mesh have been decimated.
  Vertices in the figures at a higher level than the one removed illustrate one
   possible configuration; there may be others.*/
struct od_mv_err_node{
  int dx;
  int dy;
  int log_mvb_sz;
};

/*Level 4 support:
          4-3-4
          |/|\|
          2-.-1
          |\|/|
          4-3-4
*/
static const od_mv_err_node OD_ERRDOM4[4]={
  {-1,-1,0},{ 0,-1,0},{-1, 0,0},{ 0, 0,0}
};

/*Level 3 support:
          4-3-4
          |/|\|
        4-0-.-2-4
        |/|   |\|
        3-.   .-3
        |\|   |/|
        4-2-.-1-4
          |\|/|
          4-3-4
*/
static const od_mv_err_node OD_ERRDOM3[9]={
  {-1,-2,0},{ 0,-2,0},{-2,-1,0},{ 1,-1,0},
  {-2, 0,0},{ 1, 0,0},{-1, 1,0},{ 0, 1,0},
  {-1,-1,1}
};

/*Level 2 support:
        4-3-4-3-4
        |/|\|/|\|
      4-2-.-1-.-2-4
      |/|  /|\  |\|
      3-. / | \ .-3
      |\|/  |  \|/|
      4-0---.---0-4
      |/|\  |  /|\|
      3-. \ | / .-3
      |\|  \|/  |/|
      4-2-.-1-.-2-4
        |\|/|\|/|
        4-3-4-3-4
*/
static const od_mv_err_node OD_ERRDOM2[20]={
  {-2,-3,0},{-1,-3,0},{ 0,-3,0},{ 1,-3,0},
  {-3,-2,0},{ 2,-2,0},{-3,-1,0},{ 2,-1,0},
  {-3, 0,0},{ 2, 0,0},{-3, 1,0},{ 2, 1,0},
  {-2, 2,0},{-1, 2,0},{ 0, 2,0},{ 1, 2,0},
  {-2,-2,1},{ 0,-2,1},{-2, 0,1},{ 0, 0,1}
};

/*Level 1 support:
        4-3-4-3-4
        |/|\|/|\|
      4-2-.-1-.-2-4
      |/|  /|\  |\|
    4-3-. / | \ .-3-4
    |/| |/  |  \| |\|
  4-2-.-0---.---0-.-2-4
  |/|  /|       |\  |\|
  3-. / |       | \ .-3
  |\|/  |       |  \|/|
  4-1---.       .---1-4
  |/|\  |       |  /|\|
  3-. \ |       | / .-3
  |\|  \|       |/  |/|
  4-2-.-0---.---0-.-2-4
    |\| |\  |  /| |/|
    4-3-. \ | / .-3-4
      |\|  \|/  |/|
      4-2-.-1-.-2-4
        |\|/|\|/|
        4-3-4-3-4
*/
static const od_mv_err_node OD_ERRDOM1[37]={
  {-2,-5,0},{-1,-5,0},{ 0,-5,0},{ 1,-5,0},
  {-3,-4,0},{ 2,-4,0},{-4,-3,0},{-3,-3,0},
  { 2,-3,0},{ 3,-3,0},{-5,-2,0},{ 4,-2,0},
  {-5,-1,0},{ 4,-1,0},{-5, 0,0},{ 4, 0,0},
  {-5, 1,0},{ 4, 1,0},{-4, 2,0},{-3, 2,0},
  { 2, 2,0},{ 3, 2,0},{-3, 3,0},{ 2, 3,0},
  {-2, 4,0},{-1, 4,0},{ 0, 4,0},{ 1, 4,0},
  {-2,-4,1},{ 0,-4,1},{-4,-2,1},{ 2,-2,1},
  {-4, 0,1},{ 2, 0,1},{-2, 2,1},{ 0, 2,1},
  {-2,-2,2}
};

/*The number of blocks in each decimated error domain.*/
static const int OD_NERRDOM[4]={37,20,9,4};
/*The error domain for a vertex, indexed by level-1.*/
static const od_mv_err_node *OD_ERRDOM[4]={
  OD_ERRDOM1,
  OD_ERRDOM2,
  OD_ERRDOM3,
  OD_ERRDOM4
};

/*Returns a negative value, 0, or a positive value, depending on whether
  -_dd1/_dr1 is less, equal or greater than -_dd2/_dr2.*/
static int od_mv_dddr_cmp(ogg_int32_t _dd1,int _dr1,
 ogg_int32_t  _dd2,int _dr2){
  ogg_int64_t diff;
  /*dr==0 and dd!=0 should not be possible, but we check for it anyway just in
     case, to prevent a bug from trashing the whole optimization process.*/
  if(_dr1==0)return _dr2==0?OD_SIGNI(_dd2-_dd1):(OD_SIGNI(_dd1)<<1)-1;
  else if(_dr2==0)return (OD_SIGNI(-_dd2)<<1)+1;
  diff=_dd2*(ogg_int64_t)_dr1-_dd1*(ogg_int64_t)_dr2;
  return OD_SIGNI(diff);
}

/*Compare two nodes on the decimation heap.*/
static int od_mv_dec_cmp(od_mv_node *_n1,od_mv_node *_n2){
  return od_mv_dddr_cmp(_n1->dd,_n1->dr,_n2->dd,_n2->dr);
}

/*Swap the two nodes on the decimation heap at indices _p and _q.*/
static void od_mv_dec_heap_swap(od_mv_node **_heap,int _p,int _q){
  od_mv_node *t;
  _heap[_p]->heapi=_q;
  _heap[_q]->heapi=_p;
  t=_heap[_p];
  _heap[_p]=_heap[_q];
  _heap[_q]=t;
}

/*Convert the list of nodes to be decimated to a heap.*/
static void od_mv_dec_heapify(od_mv_est_ctx *_est){
  od_mv_node **heap;
  int          l;
  int          r;
  int          i;
  heap=_est->dec_heap;
  l=_est->dec_nheap>>1;
  r=_est->dec_nheap-1;
  for(i=l;i-->0;){
    int p;
    p=i;
    do{
      int q;
      q=(p<<1)+1;
      if(q<r&&od_mv_dec_cmp(heap[q],heap[q+1])>=0)q++;
      if(od_mv_dec_cmp(heap[p],heap[q])<=0)break;
      od_mv_dec_heap_swap(heap,p,q);
      p=q;
    }
    while(p<l);
  }
}

/*Restore the heap structure at the given index by moving it down the heap.*/
static void od_mv_dec_heap_down(od_mv_est_ctx *_est,int _heapi){
  od_mv_node **heap;
  int          l;
  int          r;
  int          p;
  heap=_est->dec_heap;
  l=_est->dec_nheap>>1;
  r=_est->dec_nheap-1;
  p=_heapi;
  while(p<l){
    int q;
    q=(p<<1)+1;
    if(q<r&&od_mv_dec_cmp(heap[q],heap[q+1])>=0)q++;
    if(od_mv_dec_cmp(heap[p],heap[q])<=0)break;
    od_mv_dec_heap_swap(heap,p,q);
    p=q;
  }
}

/*Restore the heap structure at the given index by moving it up the heap.*/
static void od_mv_dec_heap_up(od_mv_est_ctx *_est,int _heapi){
  od_mv_node **heap;
  int          p;
  heap=_est->dec_heap;
  p=_heapi;
  while(p>0){
    int q;
    q=p;
    p=(q+1>>1)-1;
    if(od_mv_dec_cmp(heap[p],heap[q])<=0)break;
    od_mv_dec_heap_swap(heap,p,q);
  }
}

/*Retrieve the item at the top of the heap.
  Returns NULL if there are no more nodes to decimate.*/
static od_mv_node *od_mv_dec_heap_delhead(od_mv_est_ctx *_est){
  od_mv_node *ret;
  if(_est->dec_nheap<=0)return NULL;
  ret=_est->dec_heap[0];
  ret->heapi=-1;
  if(--_est->dec_nheap>0){
    _est->dec_heap[0]=_est->dec_heap[_est->dec_nheap];
    _est->dec_heap[0]->heapi=0;
    od_mv_dec_heap_down(_est,0);
  }
  return ret;
}

static void od_mv_dec_heap_del(od_mv_est_ctx *_est,od_mv_node *_node){
  int heapi;
  heapi=_node->heapi;
  if(heapi>=0){
    _node->heapi=-1;
    _est->dec_nheap--;
    if(_est->dec_nheap>heapi){
      _est->dec_heap[heapi]=_est->dec_heap[_est->dec_nheap];
      _est->dec_heap[heapi]->heapi=heapi;
      if(od_mv_dec_cmp(_node,_est->dec_heap[heapi])>=0){
        od_mv_dec_heap_up(_est,heapi);
      }
      else od_mv_dec_heap_down(_est,heapi);
    }
    else _est->dec_heap[_est->dec_nheap]=NULL;
  }
}

/*Sets the dd and dr values of the given node, restoring the heap structure
   afterwards.*/
static void od_mv_dec_update(od_mv_est_ctx *_est,od_mv_node *_node,
 int _dd,int _dr){
  int diff;
  diff=od_mv_dddr_cmp(_dd,_dr,_node->dd,_node->dr);
  _node->dd=_dd;
  _node->dr=_dr;
  if(_node->heapi>=0){
    if(diff<=0)od_mv_dec_heap_up(_est,_node->heapi);
    else od_mv_dec_heap_down(_est,_node->heapi);
  }
}

static void od_mv_est_init_nodes(od_mv_est_ctx *_est){
  od_state      *state;
  od_mv_node    *mv_row;
  od_mv_grid_pt *grid;
  int            nhmvbs;
  int            nvmvbs;
  int            etype;
  int            ebits;
  int            vx;
  int            vy;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  if(_est->flags&OD_MC_USEV){
    if(_est->flags&OD_MC_USEB){
      etype=0;
      ebits=3;
    }
    else{
      etype=1;
      ebits=0;
    }
  }
  else etype=ebits=0;
  for(vy=0;vy<=nvmvbs;vy++){
    mv_row=_est->mvs[vy];
    grid=state->mv_grid[vy];
    for(vx=0;vx<=nhmvbs;vx++){
      int level;
      level=OD_MC_LEVEL[vy&3][vx&3];
      if(level<=_est->level_max){
        /*While we're here, reset the MV state.*/
        grid[vx].valid=1;
        grid[vx].right=etype;
        grid[vx].down=etype;
        _est->row_counts[vy]++;
        _est->col_counts[vx]++;
        mv_row[vx].dr=-(mv_row[vx].mv_rate+
         /*Inbetween the level limits, vertices require on average 2 bits to
            indicate the presence of children.
           TODO: Fix a more exact representation.*/
         ((_est->level_min<=level&&level<_est->level_max)<<1));
        /*Vertices on even levels require new edge labels.
          Vertices on the frame border require one fewer, while those outside
           the border require none.*/
        if(!(level&1)&&vx>=2&&vx<=nhmvbs-2&&vy>=2&&vy<=nvmvbs-2){
          if(vx>2&&vx<nhmvbs-2&&vy>2&&vy<nvmvbs-2)mv_row[vx].dr-=ebits;
          else mv_row[vx].dr-=ebits+1>>1;
        }
        /*TODO: Similarly fix-up child flags.*/
      }
      else grid[vx].valid=0;
    }
  }
}

/*Computes the SAD of all blocks at all scales with all possible edge
   splittings, using OBMC.
  These are what will drive the error of the adaptive subdivision process.*/
static void od_mv_est_calc_sads(od_mv_est_ctx *_est,int _ref){
  od_state *state;
  int       nhmvbs;
  int       nvmvbs;
  int       vx;
  int       vy;
  int       c;
  int       s;
  state=&_est->enc->state;
  /*TODO: Interleaved evaluation would probably provide better cache
     coherency.*/
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  if(_est->level_max>=3){
    for(vy=0;vy<nvmvbs;vy++){
      od_mv_node *mv_row;
      mv_row=_est->mvs[vy];
      for(vx=0;vx<nhmvbs;vx++){
        c=(vx&1)^((vy&1)<<1|vy&1);
        /*While we're here, fill in the block's setup state.*/
        mv_row[vx].c=c;
        mv_row[vx].log_mvb_sz=0;
        if(_est->level_max>=4){
          for(s=0;s<4;s++){
            _est->sad_cache[0][vy][vx][s]=
            (ogg_uint16_t)od_mv_est_sad8(_est,_ref,vx,vy,c,s,0);
          }
          mv_row[vx].s=3;
          mv_row[vx].sad=_est->sad_cache[0][vy][vx][3];
        }
        else{
          mv_row[vx].s=0;
          mv_row[vx].sad=od_mv_est_sad8(_est,_ref,vx,vy,c,0,0);
        }
      }
    }
  }
  nhmvbs>>=1;
  nvmvbs>>=1;
  if(_est->level_max>=1){
    if(_est->level_min<3)for(vy=0;vy<nvmvbs;vy++){
      od_mv_node *mv_row;
      mv_row=_est->mvs[vy<<1];
      for(vx=0;vx<nhmvbs;vx++){
        c=(vx&1)^((vy&1)<<1|vy&1);
        if(_est->level_max>=2){
          for(s=0;s<4;s++){
            _est->sad_cache[1][vy][vx][s]=
             (ogg_uint16_t)od_mv_est_sad8(_est,_ref,vx<<1,vy<<1,c,s,1);
          }
          if(_est->level_max<=2){
            mv_row[vx<<1].c=c;
            mv_row[vx<<1].s=3;
            mv_row[vx<<1].log_mvb_sz=1;
            mv_row[vx<<1].sad=_est->sad_cache[1][vy][vx][3];
          }
        }
        else{
          mv_row[vx<<1].c=c;
          mv_row[vx<<1].s=0;
          mv_row[vx<<1].log_mvb_sz=1;
          mv_row[vx<<1].sad=od_mv_est_sad8(_est,_ref,vx<<1,vy<<1,c,0,1);
        }
      }
    }
  }
  else{
    nhmvbs>>=1;
    nvmvbs>>=1;
    for(vy=0;vy<nvmvbs;vy++){
      od_mv_node *mv_row;
      mv_row=_est->mvs[vy<<2];
      for(vx=0;vx<nhmvbs;vx++){
        mv_row[vx<<2].c=0;
        mv_row[vx<<2].s=3;
        mv_row[vx<<2].log_mvb_sz=2;
        mv_row[vx<<2].sad=od_mv_est_sad8(_est,_ref,vx<<2,vy<<2,0,3,2);
      }
    }
  }
}

static void od_mv_est_init_du(od_mv_est_ctx *_est,int _ref,int _vx,int _vy){
  od_state             *state;
  od_mv_node           *dec;
  od_mv_node           *merge;
  const od_mv_err_node *errdom;
  int                   nerrdom;
  const od_offset      *mergedom;
  int                   nhmvbs;
  int                   nvmvbs;
  int                   level;
  int                   dlev;
  int                   log_mvb_sz_min;
  int                   log_mvb_sz;
  int                   di;
  int                   vx;
  int                   vy;
  int                   dx;
  int                   dy;
  /*fprintf(stderr,"Computing du's for (%i,%i)\n",_vx,_vy);*/
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  dec=_est->mvs[_vy]+_vx;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  dlev=_est->level_max<=2;
  log_mvb_sz_min=5-_est->level_max>>1;
  errdom=OD_ERRDOM[level-1+(dlev<<1)];
  nerrdom=OD_NERRDOM[level-1+(dlev<<1)];
  mergedom=OD_MERGEDOM[level-1+(dlev<<1)];
  dec->dd=0;
  /*Subtract off the error before decimation.*/
  for(di=0;di<nerrdom;di++){
    vx=_vx+(errdom[di].dx<<dlev);
    vy=_vy+(errdom[di].dy<<dlev);
    if(vx>=0&&vy>=0&&vx<nhmvbs&&vy<nvmvbs){
      int mvb_sz;
      log_mvb_sz=errdom[di].log_mvb_sz+dlev;
      if(log_mvb_sz<log_mvb_sz_min)continue;
      mvb_sz=1<<log_mvb_sz-dlev;
      for(dy=0;dy<mvb_sz;dy++){
        for(dx=0;dx<mvb_sz;dx++){
          dec->dd-=_est->mvs[vy+(dy<<dlev)][vx+(dx<<dlev)].sad;
          /*dec->dd-=_est->sad_cache[dlev][(vy>>dlev)+dy][(vx>>dlev)+dx][undecs];*/
          /*fprintf(stderr,"Added error (%i,%i) [%ix%i]: %i\n",
           vx+(dx<<dlev),vy+(dy<<dlev),4<<dlev,4<<dlev,dec->dd);*/
        }
      }
    }
    /*else fprintf(stderr,"(%i,%i) outside [%i,%i]x[%i,%i]\n",vx,vy,0,0,nhmvbs,nvmvbs);*/
  }
  /*fprintf(stderr,"Subtracted initial error: %i\n",dec->dd);*/
  /*Decimate the vertices in the merging domain.
    Also sum up the rate changes while we do it.*/
  for(di=0;;di++){
    vx=_vx+(mergedom[di][0]<<dlev);
    if(vx<0||vx>nhmvbs)continue;
    vy=_vy+(mergedom[di][1]<<dlev);
    if(vy<0||vy>nvmvbs)continue;
    if(OD_MC_LEVEL[vy&3][vx&3]>_est->level_max)continue;
    state->mv_grid[vy][vx].valid=0;
    merge=_est->mvs[vy]+vx;
    if(merge==dec)break;
    dec->dr+=merge->dr;
    /*fprintf(stderr,"Merged vertex (%2i,%2i), dr: %i\n",vx,vy,dec->dr);*/
  }
  /*fprintf(stderr,"Merged vertex (%2i,%2i)\n",vx,vy);*/
  /*fprintf(stderr,"Decimated vertices in merging domain.\n");*/
  /*Add in the error after decimation.*/
  for(di=0;di<nerrdom;di++){
    vx=_vx+(errdom[di].dx<<dlev);
    vy=_vy+(errdom[di].dy<<dlev);
    if(vx>=0&&vy>=0&&vx<nhmvbs&&vy<nvmvbs){
      log_mvb_sz=errdom[di].log_mvb_sz+dlev;
      if(log_mvb_sz<log_mvb_sz_min)continue;
      else if(log_mvb_sz<2){
        int mask;
        int c;
        int s;
        mask=(1<<log_mvb_sz+1)-1;
        c=!!(vx&mask);
        if(vy&mask)c=3-c;
        s=state->mv_grid[vy+(OD_VERT_DY[c+1&3]<<log_mvb_sz)][
         vx+(OD_VERT_DX[c+1&3]<<log_mvb_sz)].valid|
         state->mv_grid[vy+(OD_VERT_DY[c+3&3]<<log_mvb_sz)][
         vx+(OD_VERT_DX[c+3&3]<<log_mvb_sz)].valid<<1;
        dec->dd+=
         _est->sad_cache[log_mvb_sz][vy>>log_mvb_sz][vx>>log_mvb_sz][s];
        /*fprintf(stderr,"Added error (%i,%i) [%ix%i] {%i,%i}: %i\n",vx,vy,1<<log_mvb_sz+2,1<<log_mvb_sz+2,c,s,dec->dd);*/
      }
      else{
        /*Cache the SAD for top-level blocks in the dd field, which is
           otherwise unused (since they cannot be decimated).*/
        _est->mvs[vy][vx].dd=od_mv_est_sad8(_est,_ref,vx,vy,0,3,2);
        dec->dd+=_est->mvs[vy][vx].dd;
        /*fprintf(stderr,"Added error (%i,%i) [%ix%i]: %i\n",
         vx,vy,1<<log_mvb_sz+2,1<<log_mvb_sz+2,dec->dd);*/
      }
    }
  }
  /*fprintf(stderr,"Total merging error: %i\n",dec->dd);*/
  /*Restore the vertices in the merging domain.*/
  for(di=0;;di++){
    vx=_vx+(mergedom[di][0]<<dlev);
    if(vx<0||vx>nhmvbs)continue;
    vy=_vy+(mergedom[di][1]<<dlev);
    if(vy<0||vy>nvmvbs)continue;
    if(OD_MC_LEVEL[vy&3][vx&3]>_est->level_max)continue;
    state->mv_grid[vy][vx].valid=1;
    if(vx==_vx&&vy==_vy)break;
  }
  /*fprintf(stderr,"Restored vertices in merging domain.\n");*/
  /*Add this node to the heap.*/
  dec->heapi=_est->dec_nheap;
  _est->dec_heap[_est->dec_nheap++]=dec;
}

static void od_mv_est_init_dus(od_mv_est_ctx *_est,int _ref){
  od_state *state;
  int       nhmvbs;
  int       nvmvbs;
  int       vx;
  int       vy;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  memset(_est->row_counts,0,sizeof(_est->col_counts[0])*(nvmvbs+1));
  memset(_est->col_counts,0,sizeof(_est->col_counts[0])*(nhmvbs+1));
  od_mv_est_init_nodes(_est);
  /*fprintf(stderr,"Finished MV bits.\n");*/
  od_mv_est_calc_sads(_est,_ref);
  /*fprintf(stderr,"Finished SADs.\n");*/
  /*Clear the merge heap.*/
  _est->dec_nheap=0;
  _est->dec_heap[0]=NULL;
  /*The initialization is destructive to dr, and so must proceed by level from
     top to bottom.*/
  if(_est->level_max>=1){
    /*Level 1 vertices.*/
    if(_est->level_min<1){
      for(vy=2;vy<=nvmvbs;vy+=4){
        for(vx=2;vx<=nhmvbs;vx+=4)od_mv_est_init_du(_est,_ref,vx,vy);
      }
    }
    if(_est->level_max>=2){
      /*Level 2 vertices.*/
      if(_est->level_min<2){
        for(vy=0;;vy+=2){
          for(vx=2;vx<=nhmvbs;vx+=4)od_mv_est_init_du(_est,_ref,vx,vy);
          vy+=2;
          if(vy>nvmvbs)break;
          for(vx=0;vx<=nhmvbs;vx+=4)od_mv_est_init_du(_est,_ref,vx,vy);
        }
      }
      if(_est->level_max>=3){
        if(_est->level_min<3){
          /*Level 3 vertices.*/
          for(vy=1;vy<=nvmvbs;vy+=2){
            for(vx=1;vx<=nhmvbs;vx+=2)od_mv_est_init_du(_est,_ref,vx,vy);
          }
        }
        if(_est->level_max>=4){
          /*Level 4 vertices.*/
          if(_est->level_min<4){
            for(vy=0;;vy++){
              for(vx=1;vx<=nhmvbs;vx+=2)od_mv_est_init_du(_est,_ref,vx,vy);
              vy++;
              if(vy>nvmvbs)break;
              for(vx=0;vx<=nhmvbs;vx+=2)od_mv_est_init_du(_est,_ref,vx,vy);
            }
          }
        }
      }
    }
  }
  /*Make the node list into a proper heap.*/
  od_mv_dec_heapify(_est);
}

static void od_mv_est_decimate(od_mv_est_ctx *_est,int _ref){
  od_mv_node *dec;
  od_state   *state;
  int         nhmvbs;
  int         nvmvbs;
  int         dlev;
  int         vx;
  int         vy;
  od_mv_est_init_dus(_est,_ref);
  od_mv_est_check_rd_state(_est,_ref,2);
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  dlev=_est->level_max<=2;
  for(;;){
    const od_offset *mergedom;
    int              level;
    int              di;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
    if(daala_granule_basetime(state,state->cur_time)==ANI_FRAME){
      char iter_label[16];
      od_state_mc_predict(state,_ref);
      od_state_fill_vis(state);
      sprintf(iter_label,"ani%08i",state->ani_iter++);
      od_state_dump_img(state,&state->vis_img,iter_label);
    }
#endif
    dec=od_mv_dec_heap_delhead(_est);
    /*Stop if we've fully decimated the mesh, or if this decimation would not
       improve R-D performance at the current lambda.*/
    if(dec==NULL||dec->dr*_est->lambda+(dec->dd<<OD_LAMBDA_SCALE)>0)break;
    level=OD_MC_LEVEL[dec->vy&3][dec->vx&3];
    /*fprintf(stderr,"Merging node (%2i,%2i), level %i, dd %5i, dr %5i, dopt %5i:\n",
     dec->vx,dec->vy,level,dec->dd,dec->dr,
     dec->dr*_est->lambda+(dec->dd<<OD_LAMBDA_SCALE));*/
    mergedom=OD_MERGEDOM[level-1+(dlev<<1)];
    for(di=0;;di++){
      od_mv_node      *merge;
      od_mv_node      *ancestor;
      od_mv_node      *block;
      const od_offset *anc;
      int              nanc;
      int              ai;
      int              ax;
      int              ay;
      int              bx;
      int              by;
      int              log_mvb_sz;
      int              mask;
      /*Don't decimate vertices outside of the mesh.*/
      vx=dec->vx+(mergedom[di][0]<<dlev);
      if(vx<0||vx>nhmvbs)continue;
      vy=dec->vy+(mergedom[di][1]<<dlev);
      if(vy<0||vy>nvmvbs)continue;
      merge=_est->mvs[vy]+vx;
      /*Don't decimate vertices that have already been decimated.*/
      if(!state->mv_grid[vy][vx].valid){/*fprintf(stderr,"Skipping node (%i,%i) (already merged).\n",vx,vy);*/continue;}
      /*fprintf(stderr,"Merging node (%2i,%2i), dd %5i, dr %5i:\n",vx,vy,
       merge->dd,merge->dr);*/
      /*Update the deltas for this vertex in the merging domain.
        The simple rule applied below handles overlapped domains with an
         inclusion-exclusion approach.
        See Balmelli 2001 for details.*/
      nanc=OD_NANCESTORS[vy&3][vx&3];
      anc=OD_ANCESTORS[vy&3][vx&3];
      for(ai=0;ai<nanc;ai++){
        ax=vx+anc[ai][0];
        if(ax<0||ax>nhmvbs)continue;
        ay=vy+anc[ai][1];
        if(ay<0||ay>nvmvbs)continue;
        ancestor=_est->mvs[ay]+ax;
        od_mv_dec_update(_est,ancestor,
         ancestor->dd-merge->dd,ancestor->dr-merge->dr);
        /*fprintf(stderr,"Updated ancestor (%2i,%2i) of (%2i,%2i): dd %5i, dr %5i\n",
         ax,ay,vx,vy,ancestor->dd,ancestor->dr);*/
      }
      state->mv_grid[vy][vx].valid=0;
      od_mv_dec_heap_del(_est,merge);
      _est->row_counts[vy]--;
      _est->col_counts[vx]--;
      level=OD_MC_LEVEL[vy&3][vx&3];
      log_mvb_sz=4-level>>1;
      /*Account for quadrilaterals which may have only partially belonged to
         the merging domain (e.g., that would not have belonged were we using
         triangles).*/
      if(!(level&1)){
        static const int OD_CDX[4]={-1,1,-1,1};
        static const int OD_CDY[4]={-1,-1,1,1};
        int k;
        mask=(1<<log_mvb_sz+1)-1;
        for(k=0;k<4;k++){
          int cx;
          int cy;
          int ddd;
          int s;
          cx=vx+(OD_CDX[k]<<log_mvb_sz);
          if(cx<0||cx>nhmvbs)continue;
          cy=vy+(OD_CDY[k]<<log_mvb_sz);
          if(cy<0||cy>nvmvbs)continue;
          bx=vx+(OD_ERRDOM4[k].dx<<log_mvb_sz);
          by=vy+(OD_ERRDOM4[k].dy<<log_mvb_sz);
          block=_est->mvs[by]+bx;
          by>>=log_mvb_sz;
          bx>>=log_mvb_sz;
          if(!state->mv_grid[cy][cx].valid){
            block->s=0;
            block->sad=_est->sad_cache[log_mvb_sz][by][bx][0];
            /*If the opposing corner has already been decimated, the remaining
               adjustments have already been made.*/
            continue;
          }
          /*s is the split state of the error block with (vx,vy) decimated, and
             (cx,cy) undecimated.*/
          s=1<<(((k+3&3)>>1)^!!(vx&mask));
          block->s=s;
          block->sad=_est->sad_cache[log_mvb_sz][by][bx][s];
          /*Replace the old decimation error change with the new one.*/
          ddd=_est->sad_cache[log_mvb_sz][by][bx][0]-
           _est->sad_cache[log_mvb_sz][by][bx][s^3]+
           _est->sad_cache[log_mvb_sz][by][bx][3]-
           _est->sad_cache[log_mvb_sz][by][bx][s];
          /*fprintf(stderr,"Checking opposing corner (%2i,%2i): ddd %i\n",
           cx,cy,ddd);*/
          /*This happens in regions of constant motion.*/
          if(ddd==0)continue;
          ancestor=_est->mvs[cy]+cx;
          od_mv_dec_update(_est,ancestor,ancestor->dd+ddd,ancestor->dr);
          /*fprintf(stderr,"Updated corner (%2i,%2i): dd %5i, dr %5i\n",
           cx,cy,ancestor->dd,ancestor->dr);*/
          /*Update the opposing corner's ancestors, which also, of
             necessity, must contain the affected quadrilateral, and must
             not have been decimated yet.*/
          nanc=OD_NANCESTORS[cy&3][cx&3];
          anc=OD_ANCESTORS[cy&3][cx&3];
          for(ai=0;ai<nanc;ai++){
            ax=cx+anc[ai][0];
            if(ax<0||ax>nhmvbs)continue;
            ay=cy+anc[ai][1];
            if(ay<0||ay>nvmvbs)continue;
            ancestor=_est->mvs[ay]+ax;
            od_mv_dec_update(_est,ancestor,ancestor->dd+ddd,ancestor->dr);
            /*fprintf(stderr,"Updated ancestor (%2i,%2i): dd %5i, dr %5i\n",
             ax,ay,ancestor->dd,ancestor->dr);*/
          }
          /*Add back in the components that do not apply to the interior
             corner.*/
          ddd=-ddd;
          if(vx&mask)cx=vx;
          else cy=vy;
          /*fprintf(stderr,"Checking interior corner (%2i,%2i): ddd %i\n",
           cx,cy,ddd);*/
          ancestor=_est->mvs[cy]+cx;
          od_mv_dec_update(_est,ancestor,ancestor->dd+ddd,ancestor->dr);
          /*fprintf(stderr,"Updated corner (%2i,%2i): dd %5i, dr %5i\n",
           cx,cy,ancestor->dd,ancestor->dr);*/
          /*And update all the interior corner's ancestors, which also, of
             necessity, must contain the affected quadrilateral, and must not
             have been decimated yet.*/
          nanc=OD_NANCESTORS[cy&3][cx&3];
          anc=OD_ANCESTORS[cy&3][cx&3];
          for(ai=0;ai<nanc;ai++){
            ax=cx+anc[ai][0];
            if(ax<0||ax>nhmvbs)continue;
            ay=cy+anc[ai][1];
            if(ay<0||ay>nvmvbs)continue;
            ancestor=_est->mvs[ay]+ax;
            od_mv_dec_update(_est,ancestor,ancestor->dd+ddd,ancestor->dr);
            /*fprintf(stderr,"Updated ancestor (%2i,%2i): dd %5i, dr %5i\n",
             ax,ay,ancestor->dd,ancestor->dr);*/
          }
        }
      }
      /*Otherwise, we eliminated several smaller blocks.
        Update the SAD and block setup for the larger block that took their
         place.*/
      else{
        int c;
        bx=vx-(1<<log_mvb_sz);
        by=vy-(1<<log_mvb_sz);
        log_mvb_sz++;
        mask=(1<<log_mvb_sz+1)-1;
        c=!!(bx&mask);
        if(by&mask)c=3-c;
        block=_est->mvs[by]+bx;
        block->log_mvb_sz=log_mvb_sz;
        block->c=c;
        block->s=3;
        if(log_mvb_sz<2){
          block->sad=
           _est->sad_cache[log_mvb_sz][by>>log_mvb_sz][bx>>log_mvb_sz][3];
        }
        /*At the top level, we cached the SAD in the dd field.*/
        else block->sad=block->dd;
      }
      /*If we just decimated our target vertex, stop.*/
      if(merge==dec)break;
    }
  }
  od_mv_est_check_rd_state(_est,_ref,2);
  /*fprintf(stderr,"Finished merging.\n");*/
  /*if(dec!=NULL){
    fprintf(stderr,"Node (%i,%i) dd %i, dr %i, dopt %i: not enough.\n",
     dec->vx,dec->vy,dec->dd,dec->dr,
     dec->dr*_est->lambda+(dec->dd<<OD_LAMBDA_SCALE));
  }*/
  /*if(state->mv_grid[31][1].valid){
    dec=_est->mvs[31]+1;
    fprintf(stderr,"(%i,%i) remains. dd: %5i, dr: %2i, dopt: %6i.\n",
     dec->vx,dec->vy,dec->dd,dec->dr,
     dec->dr*_est->lambda+(dec->dd<<OD_LAMBDA_SCALE));
  }*/
}



/*STAGE 3: Iterated Dynamic Programming.*/



/*The list of MVs that can be predicted by a level 0 MV, excluding those not
   yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED0[17]={
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { 2,-2},{ 1,-1},{ 2, 2},{ 1, 1},{ 0, 4},{ 4, 4},
  /*The remaining ones are not.*/
  {-2,-2},{ 0,-2},{-1,-1},{ 0,-1},{-1, 0},{-2, 0},
  {-1, 1},{ 0, 1},{-2, 2},{ 0, 2},{-4, 4}
};
/*The list of MVs that can be predicted by a level 1 MV, excluding those
   not yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED1[10]={
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { 1,-1},{ 1, 1},
  /*The remaining ones are not.*/
  { 0,-2},{-1,-1},{ 0,-1},{-2, 0},{-1, 0},{-1, 1},{ 0, 1},{ 0, 2}
};
/*The list of MVs that can be predicted by a level 2 MV, excluding those
   not yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED2[7]={
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { 1,-1},{ 1, 1},
  /*The remaining ones are not.*/
  {-1,-1},{ 0,-1},{-1, 0},{-1, 1},{ 0, 1}
};
/*The list of MVs that can be predicted by a level 3 MV, excluding those
   not yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED3[3]={
  /*These predicted MVs are NOT changeable by future MVs in the DP path.*/
  { 0,-1},{-1, 0},{ 0, 1}
};

/*The list of MVs that can be predicted by a level 0 MV, excluding those not
   yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED0[17]={
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { 2, 2},{-2, 2},{-1, 1},{ 1, 1},{ 4, 4},
  /*The remaining ones are not.*/
  {-2,-2},{ 0,-2},{ 2,-2},{-1,-1},{ 0,-1},{ 1,-1},
  {-2, 0},{-1, 0},{ 1, 0},{ 2, 0},{ 4, 0},{-4, 4}
};
/*The list of MVs that can be predicted by a level 1 MV, excluding those
   not yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED1[10]={
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {-1, 1},{ 1, 1},
  /*The remaining ones are not.*/
  { 0,-2},{-1,-1},{ 0,-1},{ 1,-1},{-2, 0},{-1, 0},{ 1, 0},{ 2, 0}
};
/*The list of MVs that can be predicted by a level 2 MV, excluding those
   not yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED2[7]={
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {-1, 1},{ 1, 1},
  /*The remaining ones are not.*/
  {-1,-1},{ 0,-1},{ 1,-1},{-1, 0},{ 1, 0}
};
/*The list of MVs that can be predicted by a level 3 MV, excluding those
   not yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED3[3]={
  /*These predicted MVs are NOT changeable by future MVs in the DP path.*/
  { 0,-1},{-1, 0},{ 1, 0}
};

/*The number of predicted MVs in each list.*/
static const int OD_NPREDICTED[5]={17,10,7,3,0};
/*The number of changeable predicted MVs in each list.*/
static const int OD_NROW_PRED_CHANGEABLE[4]={6,2,2,0};
/*The number of changeable predicted MVs in each list.*/
static const int OD_NCOL_PRED_CHANGEABLE[4]={5,2,2,0};
/*The lists of offsets to predicted MVs for each level.*/
static const od_offset *const OD_ROW_PREDICTED[4]={
  OD_ROW_PREDICTED0,
  OD_ROW_PREDICTED1,
  OD_ROW_PREDICTED2,
  OD_ROW_PREDICTED3
};
/*The lists of offsets to predicted MVs for each level.*/
static const od_offset *const OD_COL_PREDICTED[4]={
  OD_COL_PREDICTED0,
  OD_COL_PREDICTED1,
  OD_COL_PREDICTED2,
  OD_COL_PREDICTED3
};

/*The amount of history to restore in the trellis state to ensure predicted MVs
   are evaluated correctly in row refinement.*/
static const int OD_ROW_PRED_HIST_SIZE[5]={8,4,2,2,1};
/*The amount of history to restore in the trellis state to ensure predicted MVs
   are evaluated correctly in column refinement.*/
static const int OD_COL_PRED_HIST_SIZE[5]={8,4,2,2,1};



/*Returns the boundary case indicating which motion vector range edges the
   current motion vector is abutting.
  _vx:         The horizontal position of the node.
  _vy:         The vertical position of the node.
  _dx:         The horizontal component of the motion vector.
  _dy:         The vertical component of the motion vector.
  _dsz:        The amount the vector is being adjusted by.
  _log_blk_sz: The log base 2 of the maximum size of a block the vector can
                belong to.
  Return: A set of flags indicating the boundary conditions, after the
   documentation at OD_SQUARE_SITES.*/
static int od_mv_est_get_boundary_case(od_state *_state,int _vx,int _vy,
 int _dx,int _dy,int _dsz,int _log_blk_sz){
  int mvxmin;
  int mvxmax;
  int mvymin;
  int mvymax;
  int blk_sz;
  int bx;
  int by;
  blk_sz=1<<_log_blk_sz;
  bx=_vx-2<<2;
  by=_vy-2<<2;
  mvxmin=OD_MAXI(bx-blk_sz-32,-16)-(bx-blk_sz)<<3;
  mvxmax=(OD_MINI(bx+blk_sz+32,_state->info.frame_width+16)-(bx+blk_sz)<<3)-
   _dsz;
  mvymin=OD_MAXI(by-blk_sz-32,-16)-(by-blk_sz)<<3;
  mvymax=(OD_MINI(by+blk_sz+32,_state->info.frame_height+16)-(by+blk_sz)<<3)-
   _dsz;
  return (_dx<=mvxmin)|(_dx>=mvxmax)<<1|(_dy<=mvymin)<<2|(_dy>=mvymax)<<3;
}

/*Computes the SAD of the specified block.*/
static ogg_int32_t od_mv_est_block_sad8(od_mv_est_ctx *_est,int _ref,
 od_mv_node *_block){
  /*ogg_int32_t    ret;*/
  /*fprintf(stderr,"Adding SAD (%3i,%3i) [%2ix%2i]: ",
   _block->vx-2<<2,_block->vy-2<<2,
   4<<_block->log_mvb_sz,4<<_block->log_mvb_sz);*/
  return /*ret=*/od_mv_est_sad8(_est,_ref,_block->vx,_block->vy,
   _block->c,_block->s,_block->log_mvb_sz);
  /*fprintf(stderr,"%6i\n",ret);
  return ret;*/
}

/*Gets the change in SAD for the blocks affected by the given DP node, using
   the current state of the grid.*/
static ogg_int32_t od_mv_dp_get_sad_change8(od_mv_est_ctx *_est,int _ref,
 od_mv_dp_node *_dp,ogg_int32_t _block_sads[8]){
  int         bi;
  ogg_int32_t dd;
  dd=0;
  for(bi=0;bi<_dp->nblocks;bi++){
    od_mv_node *block;
    block=_dp->blocks[bi];
    _block_sads[bi]=od_mv_est_block_sad8(_est,_ref,block);
    /*fprintf(stderr,"SAD change for block (%i,%i) [%ix%i]: %i-%i=%i\n",
     block->vx,block->vy,1<<block->log_mvb_sz+2,1<<block->log_mvb_sz+2,
     _block_sads[bi],block->sad,_block_sads[bi]-block->sad);*/
    dd+=_block_sads[bi]-block->sad;
  }
  return dd;
}

/*Computes a rate adjustment for the predictors changed by following the given
   trellis path.
  As a side effect, enough of the trellis needed to evaluate that change is
   loaded into the MV grid.
  _dp:            The current DP node.
  _cur_mv_rate:   Returns the MV rate for the motion vector set by the current
                   DP node.
  _pred_mv_rates: Returns the MV rate for the motion vectors predicted by the
                   MV set by the current DP node.
  _prevsi:        The state index to follow in the previous DP node.
  _mv_res:        The motion vector resolution (0=1/8th pel to 2=1/2 pel).
  Return: The change in rate for the preceding MVs.*/
static int od_mv_dp_get_rate_change(od_state *_state,od_mv_dp_node *_dp,
 int *_cur_mv_rate,int _pred_mv_rates[17],int _prevsi,int _mv_res){
  od_mv_node    *mv;
  od_mv_grid_pt *mvg;
  int            nhmvbs;
  int            nvmvbs;
  int            pred[2];
  int            pi;
  int            dr;
  /*Move the state from the current trellis path into the grid.*/
  if(_dp->min_predictor_node!=NULL){
    int            pred_sis[8];
    int            pred_si;
    int            npreds;
    od_mv_dp_node *pred_dp;
    npreds=_dp-_dp->min_predictor_node;
    /*if(npreds>8)fprintf(stderr,"Too far back!\n");
    fprintf(stderr,"Restoring ");*/
    /*First, follow the trellis path backwards to find the state used in each
       node.*/
    pred_si=pred_sis[npreds-1]=_prevsi;
    for(pi=2;pi<=npreds;pi++){
      pred_dp=_dp-pi;
      pred_si=pred_dp[1].states[pred_si].prevsi;
      if(pred_si>=pred_dp[0].nstates)pred_si-=pred_dp[0].nstates;
      pred_sis[npreds-pi]=pred_si;
    }
    /*Then restore that state going FORWARDS.*/
    for(pred_dp=_dp->min_predictor_node;pred_dp<_dp;pred_dp++){
      pred_si=pred_sis[pred_dp-_dp->min_predictor_node];
      /*Restore the state for this MV itself.*/
      pred_dp->mv->mv_rate=pred_dp->states[pred_si].mv_rate;
      mvg=pred_dp->mvg;
      mvg->mv[0]=pred_dp->states[pred_si].mv[0];
      mvg->mv[1]=pred_dp->states[pred_si].mv[1];
      /*fprintf(stderr,"(%i,%i:%i)->(%i,%i) ",
       pred_dp->mv->vx,pred_dp->mv->vy,pred_si,mvg->mv[0],mvg->mv[1]);*/
      /*Restore the state for the MVs this one predicted.*/
      for(pi=0;pi<pred_dp->npred_changeable;pi++){
        pred_dp->predicted_mvs[pi]->mv_rate=
         pred_dp->states[pred_si].pred_mv_rates[pi];
      }
    }
    /*fprintf(stderr,"\n");*/
  }
  nhmvbs=_state->nhmbs+1<<2;
  nvmvbs=_state->nvmbs+1<<2;
  /*Compute the new rate for the current MV.*/
  mv=_dp->mv;
  if(mv->vx<2||mv->vx>nhmvbs-2||mv->vy<2||mv->vy>nvmvbs-2)*_cur_mv_rate=dr=0;
  else{
    od_state_get_predictor(_state,pred,mv->vx,mv->vy,
     OD_MC_LEVEL[mv->vy&3][mv->vx&3],_mv_res);
    mvg=_dp->mvg;
    *_cur_mv_rate=od_mv_est_bits(mvg->mv[0]>>_mv_res,mvg->mv[1]>>_mv_res,
     pred[0],pred[1]);
    /*fprintf(stderr,"Current MV rate: %i-%i=%i\n",
     *_cur_mv_rate,mv->mv_rate,*_cur_mv_rate-mv->mv_rate);*/
    dr=*_cur_mv_rate-mv->mv_rate;
    /*Compute the new rates for the MVs this one predicts.*/
    /*fprintf(stderr,
     "Calculating predicted pred_mv_rates for node (%i,%i):\n",
     _dp->mv->vx,_dp->mv->vy);*/
    for(pi=0;pi<_dp->npredicted;pi++){
      mv=_dp->predicted_mvs[pi];
      mvg=_dp->predicted_mvgs[pi];
      od_state_get_predictor(_state,pred,mv->vx,mv->vy,
       OD_MC_LEVEL[mv->vy&3][mv->vx&3],_mv_res);
      _pred_mv_rates[pi]=od_mv_est_bits(
       mvg->mv[0]>>_mv_res,mvg->mv[1]>>_mv_res,pred[0],pred[1]);
      /*fprintf(stderr,"Calculated predicted mv_rate of %i for (%i,%i)\n",
       _pred_mv_rates[pi],mv->vx,mv->vy);
      fprintf(stderr,"Predictor was: (%i,%i)   MV was: (%i,%i)\n",
       pred[0],pred[1],mvg->mv[0]>>_mv_res,mvg->mv[1]>>_mv_res);
      fprintf(stderr,"Predicted MV (%i,%i) rate: %i-%i=%i\n",
       mv->vx,mv->vy,_pred_mv_rates[pi],mv->mv_rate,
       _pred_mv_rates[pi]-mv->mv_rate);*/
      dr+=_pred_mv_rates[pi]-mv->mv_rate;
    }
  }
  return dr;
}

#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
static const unsigned char OD_YCbCr_BEDGE[3]= { 41,240,110};
static const unsigned char OD_YCbCr_VEDGE[3]= {145, 54, 34};
static const unsigned char OD_YCbCr_VBEDGE[3]={170,166, 16};

static void od_mv_dp_animate_state(od_state *_state,int _ref,
 od_mv_dp_node *_dp,int _has_gap){
  od_mv_dp_node *dp;
  char           iter_label[16];
  int            active_states[OD_DP_NSTATES_MAX<<1];
  int            prev_active_states[OD_DP_NSTATES_MAX<<1];
  int            nactive_states;
  int            nprev_active_states;
  int            state;
  int            si;
  int            x0;
  int            y0;
  od_state_mc_predict(_state,_ref);
  od_state_fill_vis(_state);
  /*Now, draw the current state of the DP.*/
  /*First draw the candidate edge labels for the active trellis paths.*/
  for(si=0;si<_dp->nstates;si++){
    prev_active_states[si<<1]=si;
    prev_active_states[si<<1|1]=si+_dp->nstates;
  }
  nprev_active_states=_dp->nstates<<1;
  nactive_states=0;
  dp=_dp;
  do{
    int has_vedge;
    int has_bedge;
    if(nactive_states>0){
      x0=(dp[0].mv->vx-2<<3)+(OD_UMV_PADDING<<1);
      y0=(dp[0].mv->vy-2<<3)+(OD_UMV_PADDING<<1);
      has_vedge=has_bedge=0;
      for(si=0;si<nprev_active_states;si++){
        if(prev_active_states[si]<dp[0].nstates)has_bedge=1;
        else has_vedge=1;
      }
      if(has_vedge||has_bedge){
        int mvb_sz;
        int x1;
        int y1;
        x1=(dp[1].mv->vx-2<<3)+(OD_UMV_PADDING<<1);
        y1=(dp[1].mv->vy-2<<3)+(OD_UMV_PADDING<<1);
        od_img_draw_line(&_state->vis_img,x0,y0,x1,y1,
         has_vedge?has_bedge?OD_YCbCr_VBEDGE:OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
        if(dp[1].mv->vx-dp[0].mv->vx>1){
          mvb_sz=dp[1].mv->vx-dp[0].mv->vx;
          if(!_has_gap||dp+1!=_dp)mvb_sz>>=1;
          if(!_state->mv_grid[dp[0].mv->vy][dp[0].mv->vx+mvb_sz].valid){
            if(dp[0].mv->vy>=mvb_sz&&
             _state->mv_grid[dp[0].mv->vy-mvb_sz][dp[0].mv->vx+mvb_sz].valid){
              od_img_draw_line(&_state->vis_img,
               x0+(mvb_sz<<3),y0-(mvb_sz<<3),x0+(mvb_sz<<3),y1,
               has_vedge?has_bedge?OD_YCbCr_VBEDGE:
               OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
            }
            if(dp[0].mv->vy<=(_state->nvmbs+1<<2)-mvb_sz&&
             _state->mv_grid[dp[0].mv->vy+mvb_sz][dp[0].mv->vx+mvb_sz].valid){
              od_img_draw_line(&_state->vis_img,
               x0+(mvb_sz<<3),y0+(mvb_sz<<3),x0+(mvb_sz<<3),y1,
               has_vedge?has_bedge?OD_YCbCr_VBEDGE:
               OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
            }
          }
        }
        else if(dp[1].mv->vy-dp[0].mv->vy>1){
          mvb_sz=dp[1].mv->vy-dp[0].mv->vy;
          if(!_has_gap||dp+1!=_dp)mvb_sz>>=1;
          if(!_state->mv_grid[dp[0].mv->vy+mvb_sz][dp[0].mv->vx].valid){
            if(dp[0].mv->vx>=mvb_sz&&
             _state->mv_grid[dp[0].mv->vy+mvb_sz][dp[0].mv->vx-mvb_sz].valid){
              od_img_draw_line(&_state->vis_img,
               x0-(mvb_sz<<3),y0+(mvb_sz<<3),x1,y0+(mvb_sz<<3),
               has_vedge?has_bedge?OD_YCbCr_VBEDGE:
               OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
            }
            if(dp[0].mv->vx<=(_state->nhmbs+1<<2)-mvb_sz&&
             _state->mv_grid[dp[0].mv->vy+mvb_sz][dp[0].mv->vx+mvb_sz].valid){
              od_img_draw_line(&_state->vis_img,
               x0+(mvb_sz<<3),y0+(mvb_sz<<3),x1,y0+(mvb_sz<<3),
               has_vedge?has_bedge?OD_YCbCr_VBEDGE:
               OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
            }
          }
        }
      }
    }
    memcpy(active_states,prev_active_states,
     sizeof(active_states[0])*nprev_active_states);
    nactive_states=nprev_active_states;
    /*Follow the chain backwards to find the new active states.*/
    nprev_active_states=0;
    for(si=0;si<nactive_states;si++){
      int sj;
      state=active_states[si];
      if(state>=dp[0].nstates)state-=dp[0].nstates;
      state=dp[0].states[state].prevsi;
      for(sj=0;sj<nprev_active_states&&prev_active_states[sj]!=state;sj++);
      if(sj>=nprev_active_states){
        prev_active_states[nprev_active_states++]=state;
      }
    }
  }
  while((dp--)->states[0].prevsi>=0);
  /*Now, draw all the candidate MVs in the active trellis paths.
    These two steps used to be together; now they're apart.
    Sorry for the mess that caused.*/
  /*Redraw the MVs, so they appear over the edge labels above.*/
  od_state_draw_mvs(_state);
  for(si=0;si<_dp->nstates;si++){
    prev_active_states[si<<1]=si;
    prev_active_states[si<<1|1]=si+_dp->nstates;
  }
  nprev_active_states=_dp->nstates<<1;
  nactive_states=0;
  dp=_dp;
  do{
    if(nactive_states>0){
      x0=(dp[0].mv->vx-2<<3)+(OD_UMV_PADDING<<1);
      y0=(dp[0].mv->vy-2<<3)+(OD_UMV_PADDING<<1);
      if(!_has_gap||dp+1!=_dp){
        x0=(dp[1].mv->vx-2<<3)+(OD_UMV_PADDING<<1);
        y0=(dp[1].mv->vy-2<<3)+(OD_UMV_PADDING<<1);
        for(si=0;si<nactive_states;si++){
          state=active_states[si];
          if(state>=dp[1].nstates)state-=dp[1].nstates;
          od_img_draw_line(&_state->vis_img,x0,y0,
           x0+OD_DIV_ROUND_POW2(dp[1].states[state].mv[0],2,2),
           y0+OD_DIV_ROUND_POW2(dp[1].states[state].mv[1],2,2),
           OD_YCbCr_MVCAND);
        }
      }
    }
    memcpy(active_states,prev_active_states,
     sizeof(active_states[0])*nprev_active_states);
    nactive_states=nprev_active_states;
    /*Follow the chain backwards to find the new active states.*/
    nprev_active_states=0;
    for(si=0;si<nactive_states;si++){
      int sj;
      state=active_states[si];
      if(state>=dp[0].nstates)state-=dp[0].nstates;
      state=dp[0].states[state].prevsi;
      for(sj=0;sj<nprev_active_states&&prev_active_states[sj]!=state;sj++);
      if(sj>=nprev_active_states){
        prev_active_states[nprev_active_states++]=state;
      }
    }
  }
  while((dp--)->states[0].prevsi>=0);
  /*Draw the first state's MV's.*/
  x0=(dp[1].mv->vx-2<<3)+(OD_UMV_PADDING<<1);
  y0=(dp[1].mv->vy-2<<3)+(OD_UMV_PADDING<<1);
  for(si=0;si<nactive_states;si++){
    state=active_states[si];
    if(state>=dp[1].nstates)state-=dp[1].nstates;
    od_img_draw_line(&_state->vis_img,x0,y0,
     x0+OD_DIV_ROUND_POW2(dp[1].states[state].mv[0],2,2),
     y0+OD_DIV_ROUND_POW2(dp[1].states[state].mv[1],2,2),
     OD_YCbCr_MVCAND);
  }
  sprintf(iter_label,"ani%08i",_state->ani_iter++);
  od_state_dump_img(_state,&_state->vis_img,iter_label);
}
#endif

/*Row refinement.*/

static void od_mv_dp_row_init(od_mv_est_ctx *_est,od_mv_dp_node *_dp,
 int _vx,int _vy,od_mv_dp_node *_prev_dp){
  od_state      *state;
  int            nhmvbs;
  int            nvmvbs;
  state=&_est->enc->state;
  _dp->mv=_est->mvs[_vy]+_vx;
  _dp->mvg=state->mv_grid[_vy]+_vx;
  _dp->original_mv[0]=_dp->mvg->mv[0];
  _dp->original_mv[1]=_dp->mvg->mv[1];
  _dp->original_etype=_dp->mvg->right;
  _dp->original_mv_rate=_dp->mv->mv_rate;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  if(_vx<2||_vx>nhmvbs-2||_vy<2||_vy>nvmvbs-2){
    /*Strictly speaking, we may be used to predict others, but since our MV
       can't possibly change, neither can their rate.*/
    _dp->npredicted=_dp->npred_changeable=0;
    /*No one else is used to predict us, or any other MV we predict.
      However, we may still need to load the previous MV into the grid to
       estimate our SADs properly.*/
    _dp->min_predictor_node=_prev_dp;
  }
  else{
    int level;
    int pred_hist;
    int npred;
    int nchangeable;
    int pi;
    /*Get the list of MVs we help predict.*/
    level=OD_MC_LEVEL[_vy&3][_vx&3];
    /*fprintf(stderr,"Initializing node (%i,%i) [%i,%i] at level %i:\n",
     _vx,_vy,_vx-2<<2,_vy-2<<2,level);*/
    npred=nchangeable=0;
    for(pi=0;pi<OD_NPREDICTED[level];pi++){
      int px;
      int py;
      px=_vx+OD_ROW_PREDICTED[level][pi][0];
      if(px<2||px>nhmvbs-2)continue;
      py=_vy+OD_ROW_PREDICTED[level][pi][1];
      if(py<2||py>nvmvbs-2)continue;
      if(state->mv_grid[py][px].valid){
        /*fprintf(stderr,"Adding (%i,%i) as a PREDICTED MV.\n",px,py);*/
        _dp->predicted_mvgs[npred]=state->mv_grid[py]+px;
        _dp->predicted_mvs[npred]=_est->mvs[py]+px;
        if(pi<OD_NROW_PRED_CHANGEABLE[level]){
          /*fprintf(stderr,"It is CHANGEABLE.\n");*/
          _dp->original_mv_rates[npred]=_est->mvs[py][px].mv_rate;
          nchangeable++;
        }
        npred++;
      }
    }
    _dp->npredicted=npred;
    _dp->npred_changeable=nchangeable;
    /*Now, figure out the earliest DP node that influences our own prediction,
       or that of one of the other MVs we predict.*/
    pred_hist=OD_ROW_PRED_HIST_SIZE[level];
    /*fprintf(stderr,"Marking history up to %i back: %i>=%i\n",
     pred_hist,_prev_dp!=NULL?_prev_dp->mv->vx:-1,_vx-pred_hist);*/
    if(_prev_dp!=NULL&&_prev_dp->mv->vx>=_vx-pred_hist){
      od_mv_dp_node *dp_pred;
      for(dp_pred=_prev_dp;dp_pred->mv->vx>_vx-pred_hist&&
       dp_pred->states[0].prevsi>=0;dp_pred--);
      /*fprintf(stderr,"Stopped at (%i,%i) (%i<=%i? %i) (%i<0? %i)\n",
       dp_pred->mv->vx,dp_pred->mv->vy,dp_pred->mv->vx,_vx-pred_hist,
       dp_pred->mv->vx<=_vx-pred_hist,
       dp_pred->states[0].prevsi,dp_pred->states[0].prevsi<0);*/
      if(dp_pred->mv->vx<_vx-pred_hist){dp_pred++;/*fprintf(stderr,"Too far, incrementing to (%i,%i).\n",dp_pred->mv->vx,dp_pred->mv->vy);*/}
      _dp->min_predictor_node=dp_pred;
      /*fprintf(stderr,"State will be restored back to (%i,%i).\n",
       dp_pred->mv->vx,dp_pred->mv->vy);*/
    }
    else _dp->min_predictor_node=NULL;
  }
}

static void od_mv_dp_first_row_block_setup(od_mv_est_ctx *_est,
 od_mv_dp_node *_dp,int _vx,int _vy){
  od_state *state;
  int       nvmvbs;
  int       level;
  int       log_mvb_sz;
  int       mvb_sz;
  int       nblocks;
  state=&_est->enc->state;
  nvmvbs=state->nvmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(_vx>2){
    if(level>=3){
      if(_vy>=mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      if(_vy<=nvmvbs-mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_sz;
    }
    else{
      int half_mvb_sz;
      int mvb_off;
      half_mvb_sz=mvb_sz>>1;
      if(_vy>=mvb_sz){
        if(state->mv_grid[_vy-half_mvb_sz][_vx-half_mvb_sz].valid){
          if(level>0||
           !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
            mvb_off=half_mvb_sz;
          }
          else mvb_off=half_mvb_sz>>1;
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-mvb_off;
          if(!state->mv_grid[_vy-mvb_off][_vx].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx-mvb_off;
          }
          if(!state->mv_grid[_vy][_vx-mvb_off].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-(mvb_off<<1);
          }
        }
        else _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      }
      if(_vy<=nvmvbs-mvb_sz){
        if(state->mv_grid[_vy+half_mvb_sz][_vx-half_mvb_sz].valid){
          if(level>0||
           !state->mv_grid[_vy+(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
            mvb_off=half_mvb_sz;
          }
          else mvb_off=half_mvb_sz>>1;
          _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_off;
          if(!state->mv_grid[_vy+mvb_off][_vx].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_off]+_vx-mvb_off;
          }
          if(!state->mv_grid[_vy][_vx-mvb_off].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-(mvb_off<<1);
          }
        }
        else _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_sz;
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_prev_row_block_setup(od_mv_est_ctx *_est,
 od_mv_dp_node *_dp,int _vx,int _vy){
  od_state *state;
  int       nvmvbs;
  int       level;
  int       prev_level;
  int       log_mvb_sz;
  int       prev_log_mvb_sz;
  int       mvb_sz;
  int       nblocks;
  state=&_est->enc->state;
  nvmvbs=state->nvmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  prev_level=OD_MC_LEVEL[_vy&3][_vx-mvb_sz&3];
  prev_log_mvb_sz=4-prev_level>>1;
  nblocks=0;
  if(level>=3){
    if(_vy>=mvb_sz){
      _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      if(prev_log_mvb_sz>log_mvb_sz&&
       !state->mv_grid[_vy-mvb_sz][_vx-mvb_sz].valid){
        _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_sz<<1)]+_vx-mvb_sz;
      }
    }
    if(_vy<=nvmvbs-mvb_sz){
      _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_sz;
      if(prev_log_mvb_sz>log_mvb_sz&&
       !state->mv_grid[_vy+mvb_sz][_vx-mvb_sz].valid){
        _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_sz]+_vx-mvb_sz;
      }
    }
  }
  else{
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vy>=mvb_sz){
      if(state->mv_grid[_vy-half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-mvb_off;
        if(!state->mv_grid[_vy-mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx-mvb_off;
        }
        if(!state->mv_grid[_vy][_vx-mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-(mvb_off<<1);
          if(!state->mv_grid[_vy-mvb_off][_vx-(mvb_off<<1)].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx-(mvb_off<<1);
          }
        }
      }
      else{
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
        if(prev_log_mvb_sz>log_mvb_sz&&
         !state->mv_grid[_vy-mvb_sz][_vx-mvb_sz].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_sz<<1)]+_vx-mvb_sz;
        }
      }
    }
    if(_vy<=nvmvbs-mvb_sz){
      if(state->mv_grid[_vy+half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy+(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_off;
        if(!state->mv_grid[_vy+mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_off]+_vx-mvb_off;
        }
        if(!state->mv_grid[_vy][_vx-mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-(mvb_off<<1);
          if(!state->mv_grid[_vy+mvb_off][_vx-(mvb_off<<1)].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_off]+_vx-(mvb_off<<1);
          }
        }
      }
      else{
        _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_sz;
        if(prev_log_mvb_sz>log_mvb_sz&&
         !state->mv_grid[_vy+mvb_sz][_vx-mvb_sz].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_sz]+_vx-mvb_sz;
        }
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_last_row_block_setup(od_mv_est_ctx *_est,
 od_mv_dp_node *_dp,int _vx,int _vy){
  od_state *state;
  int       nvmvbs;
  int       level;
  int       log_mvb_sz;
  int       mvb_sz;
  int       nblocks;
  state=&_est->enc->state;
  nvmvbs=state->nvmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(level>=3){
    if(_vy>=mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx;
    if(_vy<=nvmvbs-mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy]+_vx;
  }
  else{
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vy>=mvb_sz){
      if(state->mv_grid[_vy-half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx;
        if(!state->mv_grid[_vy][_vx+mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx+mvb_off;
        }
        if(!state->mv_grid[_vy-mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx;
        }
      }
      else _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx;
    }
    if(_vy<=nvmvbs-mvb_sz){
      if(state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy+(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx;
        if(!state->mv_grid[_vy][_vx+mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx+mvb_off;
        }
        if(!state->mv_grid[_vy+mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_off]+_vx;
        }
      }
      else _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx;
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_restore_row_state(od_mv_dp_node *_dp){
  od_mv_grid_pt *mvg;
  int            pi;
  do{
    /*Restore the state for this MV itself.*/
    _dp->mv->mv_rate=_dp->original_mv_rate;
    mvg=_dp->mvg;
    mvg->mv[0]=_dp->original_mv[0];
    mvg->mv[1]=_dp->original_mv[1];
    mvg->right=_dp->original_etype;
    for(pi=0;pi<_dp->npred_changeable;pi++){
      /*Restore the state for the MVs this one predicted.*/
      _dp->predicted_mvs[pi]->mv_rate=_dp->original_mv_rates[pi];
    }
  }
  while((_dp--)->states[0].prevsi>=0);
}

static void od_mv_dp_install_row_state(od_mv_dp_node *_dp,int _prevsi){
  od_mv_dp_node *dp;
  od_mv_grid_pt *mvg;
  int            nextsi;
  int            si;
  int            pi;
  int            bi;
  /*We must install the state going FORWARDS, since the pred_mv_rates may have
     changed several times over the course of the trellis.
    Therefore, first we reverse all of the prevsi pointers to make them act
     like nextsi pointers.
    We _can_ update the edge type flags here, however, and it is much more
     convenient to do so while going backwards than forwards.*/
  nextsi=-1;
  for(dp=_dp,si=_prevsi;si>=0;si=_prevsi){
    dp--;
    /*fprintf(stderr,"Node %i, prevsi: %i nextsi: %i\n",_dp-dp,_prevsi,nextsi);*/
    if(si>=dp->nstates){
      dp->mvg->right=1;
      si-=dp->nstates;
    }
    else dp->mvg->right=0;
    _prevsi=dp->states[si].prevsi;
    dp->states[si].prevsi=nextsi;
    nextsi=si;
  }
  /*Now we traverse forward installing the rest of the state.*/
  for(si=nextsi;dp<_dp;dp++){
    /*fprintf(stderr,"Installing state %i for (%i,%i):\n",
     si,dp->mv->vx,dp->mv->vy);*/
    /*Install the state for this MV itself.*/
    /*fprintf(stderr,"Installing current mv_rate for (%i,%i): %i\n",
     dp->mv->vx,dp->mv->vy,dp->states[si].mv_rate);*/
    dp->mv->mv_rate=dp->states[si].mv_rate;
    mvg=dp->mvg;
    mvg->mv[0]=dp->states[si].mv[0];
    mvg->mv[1]=dp->states[si].mv[1];
    /*Install the new block SADs.*/
    for(bi=0;bi<dp->nblocks;bi++){
      dp->blocks[bi]->sad=dp->states[si].block_sads[bi];
    }
    /*Install the state for the MVs this one predicted.*/
    for(pi=0;pi<dp->npredicted;pi++){
      /*fprintf(stderr,"Installing predicted mv_rate for (%i,%i): %i\n",
       dp->predicted_mvs[pi]->vx,dp->predicted_mvs[pi]->vy,
       dp->states[si].pred_mv_rates[pi]);*/
      dp->predicted_mvs[pi]->mv_rate=dp->states[si].pred_mv_rates[pi];
    }
    si=dp->states[si].prevsi;
  }
}

static ogg_int32_t od_mv_est_refine_row(od_mv_est_ctx *_est,int _ref,int _vy,
 int _log_dsz,int _mv_res,const int *_pattern_nsites,
 const od_pattern *_pattern){
  od_state       *state;
  od_mv_grid_pt  *grid;
  od_mv_grid_pt  *pmvg;
  od_mv_grid_pt  *mvg;
  od_mv_dp_node  *dp_node;
  od_mv_dp_state *cstate;
  od_mv_dp_state *pstate;
  ogg_int32_t     dcost;
  int             nhmvbs;
  int             nvmvbs;
  int             level;
  int             log_mvb_sz;
  int             mvb_sz;
  int             labels_only;
  int             nsites;
  int             sitei;
  int             site;
  int             curx;
  int             cury;
  int             vx;
  int             b;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  labels_only=_vy<2||_vy>nvmvbs-2;
  grid=state->mv_grid[_vy];
  dcost=0;
  /*fprintf(stderr,"Refining row %i (%i)...\n",_vy,_vy-2<<2);*/
  for(vx=0;;vx++){
    ogg_int32_t block_sads[18][8];
    ogg_int32_t best_cost;
    ogg_int32_t cost;
    ogg_int32_t best_dd;
    ogg_int32_t dd;
    int         cur_mv_rates[9];
    int         pred_mv_rates[9][17];
    int         best_dr;
    int         dr;
    int         best_si;
    int         si;
    int         has_gap;
    for(;vx<=nhmvbs&&!grid[vx].valid;vx++);
    if(vx>nhmvbs)break;
    level=OD_MC_LEVEL[_vy&3][vx&3];
    /*fprintf(stderr,"Starting DP at vertex %i (%i), level %i\n",
     vx,vx-2<<2,level);*/
    log_mvb_sz=4-level>>1;
    mvb_sz=1<<log_mvb_sz;
    mvg=grid+vx;
    curx=mvg->mv[0];
    cury=mvg->mv[1];
    dp_node=_est->dp_nodes;
    od_mv_dp_row_init(_est,dp_node,vx,_vy,NULL);
    od_mv_dp_first_row_block_setup(_est,dp_node,vx,_vy);
    /*fprintf(stderr,"TESTING block SADs:\n");*/
    od_mv_dp_get_sad_change8(_est,_ref,dp_node,block_sads[0]);
    /*Compute the set of states for the first node.*/
    if(vx>=2&&!labels_only){
      b=od_mv_est_get_boundary_case(state,vx,_vy,curx,cury,
       1<<_log_dsz,log_mvb_sz+2);
      nsites=_pattern_nsites[b];
    }
    else b=nsites=0;
    for(sitei=0,site=4;;sitei++){
      cstate=dp_node[0].states+sitei;
      cstate->mv[0]=curx+(OD_SITE_DX[site]<<_log_dsz);
      cstate->mv[1]=cury+(OD_SITE_DY[site]<<_log_dsz);
      cstate->prevsi=-1;
      mvg->mv[0]=cstate->mv[0];
      mvg->mv[1]=cstate->mv[1];
      cstate->dr=od_mv_dp_get_rate_change(state,dp_node,
       &cstate->mv_rate,cstate->pred_mv_rates,-1,_mv_res);
      cstate->dd=od_mv_dp_get_sad_change8(_est,_ref,dp_node,
       cstate->block_sads);
      /*fprintf(stderr,"State: %i (%g,%g)  dr: %i  dd: %i  dopt: %i\n",
       sitei,0.125*cstate->mv[0],0.125*cstate->mv[1],cstate->dr,cstate->dd,
       cstate->dr*_est->lambda+(cstate->dd<<OD_LAMBDA_SCALE));*/
      if(sitei>=nsites)break;
      site=_pattern[b][sitei];
    }
    dp_node[0].nstates=nsites+1;
    has_gap=0;
    pmvg=mvg;
    for(;vx<nhmvbs;){
      /*Find the next available MV to advance to.*/
      if(level&1){
        if(!grid[vx+mvb_sz].valid){
          /*fprintf(stderr,"Gap found at %i (%i), stopping\n",vx,vx-2<<2);*/
          has_gap=1;
          break;
        }
        else if(level>=3)vx++;
        else if(!grid[vx+1].valid)vx+=mvb_sz;
        else vx++;
      }
      else if(level>=4)vx++;
      else if(!grid[vx+(mvb_sz>>1)].valid)vx+=mvb_sz;
      else if(level>=2||!grid[vx+1].valid)vx+=mvb_sz>>1;
      else vx++;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
      if(daala_granule_basetime(state,state->cur_time)==ANI_FRAME){
        od_mv_dp_restore_row_state(dp_node);
        od_mv_dp_animate_state(state,_ref,dp_node,0);
      }
#endif
      level=OD_MC_LEVEL[_vy&3][vx&3];
      /*fprintf(stderr,"Continuing DP at vertex %i (%i), level %i\n",
       vx,vx-2<<2,level);*/
      log_mvb_sz=4-level>>1;
      mvb_sz=1<<log_mvb_sz;
      mvg=grid+vx;
      curx=mvg->mv[0];
      cury=mvg->mv[1];
      od_mv_dp_row_init(_est,dp_node+1,vx,_vy,dp_node);
      od_mv_dp_prev_row_block_setup(_est,dp_node+1,vx,_vy);
      /*fprintf(stderr,"TESTING block SADs:\n");
      pmvg->mv[0]=dp_node[0].original_mv[0];
      pmvg->mv[0]=dp_node[0].original_mv[0];
      od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,block_sads[0]);*/
      /*Compute the set of states for this node.*/
      if(vx<=nhmvbs-2&&!labels_only){
        b=od_mv_est_get_boundary_case(state,vx,_vy,curx,cury,
         1<<_log_dsz,log_mvb_sz+2);
        nsites=_pattern_nsites[b];
      }
      else nsites=0;
      for(sitei=0,site=4;;sitei++){
        cstate=dp_node[1].states+sitei;
        cstate->mv[0]=curx+(OD_SITE_DX[site]<<_log_dsz);
        cstate->mv[1]=cury+(OD_SITE_DY[site]<<_log_dsz);
        best_si=pmvg->right?dp_node[0].nstates:0;
        best_dr=dp_node[0].states[0].dr;
        best_dd=dp_node[0].states[0].dd;
        best_cost=INT_MAX;
        mvg->mv[0]=cstate->mv[0];
        mvg->mv[1]=cstate->mv[1];
        for(si=0;si<dp_node[0].nstates;si++){
          pstate=dp_node[0].states+si;
          /*Get the rate change for this state using previous state si.
            This automatically loads the required bits of the trellis path into
             the grid, like the previous MV.*/
          cstate->dr=od_mv_dp_get_rate_change(state,dp_node+1,
           cur_mv_rates+si,pred_mv_rates[si],si,_mv_res);
          /*Test against the previous state with a B edge.*/
          if(_est->flags&OD_MC_USEB){
            pmvg->right=0;
            dr=pstate->dr+cstate->dr;
            /*Account for label mismatch.*/
            /*if(pstate->prevsi>=0&&pstate->prevsi>=(dp_node-1)->nstates)dr+=2;*/
            dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
             block_sads[si]);
            cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
            /*fprintf(stderr,
             "State: %2i (%7g,%7g) P.State: %iB  dr: %3i  dd: %6i  dopt: %7i\n",
             sitei,0.125*cstate->mv[0],0.125*cstate->mv[1],si,dr,dd,cost);*/
            if(cost<best_cost){
              best_si=si;
              best_cost=cost;
              best_dd=dd;
              best_dr=dr;
            }
          }
          if(_est->flags&OD_MC_USEV){
            /*Test against the previous state with a V edge.*/
            pmvg->right=1;
            dr=pstate->dr+cstate->dr;
            /*Account for label mismatch.*/
            /*if(pstate->prevsi>=0&&pstate->prevsi<(dp_node-1)->nstates)rate+=2;*/
            dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
             block_sads[si+dp_node[0].nstates]);
            cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
            /*fprintf(stderr,
             "State: %2i (%7g,%7g) P.State: %iV  dr: %3i  dd: %6i  dopt: %7i\n",
             sitei,0.125*cstate->mv[0],0.125*cstate->mv[1],si,dr,dd,cost);*/
            if(cost<best_cost){
              best_si=si+dp_node[0].nstates;
              best_cost=cost;
              best_dd=dd;
              best_dr=dr;
            }
          }
        }
        /*fprintf(stderr,"State: %2i  Best P.State: %i%c\n",
         sitei,best_si%dp_node[0].nstates,best_si>dp_node[0].nstates?'V':'B');*/
        cstate->prevsi=best_si;
        cstate->dr=best_dr;
        cstate->dd=best_dd;
        memcpy(cstate->block_sads,block_sads[best_si],
         sizeof(block_sads[0][0])*dp_node[1].nblocks);
        if(best_si<dp_node[0].nstates){
          cstate->mv_rate=cur_mv_rates[best_si];
          memcpy(cstate->pred_mv_rates,pred_mv_rates[best_si],
           sizeof(pred_mv_rates[0][0])*dp_node[1].npredicted);
        }
        else{
          cstate->mv_rate=cur_mv_rates[best_si-dp_node[0].nstates];
          memcpy(cstate->pred_mv_rates,
           pred_mv_rates[best_si-dp_node[0].nstates],
           sizeof(pred_mv_rates[0][0])*dp_node[1].npredicted);
        }
        if(sitei>=nsites)break;
        site=_pattern[b][sitei];
      }
      dp_node[1].nstates=nsites+1;
      dp_node++;
      pmvg=mvg;
    }
    /*fprintf(stderr,"Finished DP at vertex %i (%i)\n",
     dp_node[0].mv->vx,dp_node[0].mv->vx-2<<2);*/
    best_si=pmvg->right?dp_node[0].nstates:0;
    best_cost=INT_MAX;
    /*TODO: Once we stop optimizing at arbitrary places, we'll need to
       compute the rate change of MVs we didn't get to.*/
    dp_node[1].npredicted=dp_node[1].npred_changeable=0;
    if(dp_node[0].mv->vx<nhmvbs-2){
      od_mv_dp_last_row_block_setup(_est,dp_node+1,dp_node[0].mv->vx,_vy);
      /*fprintf(stderr,"TESTING block SADs:\n");
      pmvg->mv[0]=dp_node[0].original_mv[0];
      pmvg->mv[1]=dp_node[0].original_mv[1];
      od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,block_sads[0]);*/
      for(si=0;si<dp_node[0].nstates;si++){
        pstate=dp_node[0].states+si;
        pmvg->mv[0]=pstate->mv[0];
        pmvg->mv[1]=pstate->mv[1];
        /*Test against the state with a following B edge.*/
        if(_est->flags&OD_MC_USEB){
          pmvg->right=0;
          dr=pstate->dr;
          /*Account for label mismatch.*/
          /*if(!has_gap&&pstate->prevsi>=0&&pstate->prevsi>=(dp_node-1)->nstates){
            rate+=2;
          }*/
          dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
           block_sads[si]);
          cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: --  P.State: %iB  dr: %3i  dd: %6i  dopt: %7i\n",
           si,dr,dd,cost);*/
          if(cost<best_cost){
            best_si=si;
            best_cost=cost;
          }
        }
        /*Test against the state with a following V edge.
          If we hit a gap, then the edge label does not actually affect
           anything, so we can skip these tests if we did the ones above.*/
        if((_est->flags&OD_MC_USEV)&&(!has_gap||!(_est->flags&OD_MC_USEB))){
          pmvg->right=1;
          dr=pstate->dr;
          /*Account for label mismatch.*/
          /*if(pstate->prevsi>=0&&pstate->prevsi<(dp_node-1)->nstates)rate+=2;*/
          dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
           block_sads[si+dp_node[0].nstates]);
          cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: --  P.State: %iV  dr: %3i  dd: %6i  dopt: %7i\n",
           si,dr,dd,cost);*/
          if(cost<best_cost){
            best_si=si+dp_node[0].nstates;
            best_cost=cost;
          }
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else{
      dp_node[1].nblocks=0;
      for(si=0;si<dp_node[0].nstates;si++){
        pstate=dp_node[0].states+si;
        cost=pstate->dr*_est->lambda+(pstate->dd<<OD_LAMBDA_SCALE);
        if(cost<best_cost){
          best_si=si;
          best_cost=cost;
        }
      }
      if(pmvg->right)best_si+=dp_node[0].nstates;
    }
    /*fprintf(stderr,"Best P.State: %i%c dopt: %7i\n",
     best_si%dp_node[0].nstates,best_si>dp_node[0].nstates?'V':'B',best_cost);*/
    if(best_cost>0){
      /*Our optimal path is worse than what we started with!
        Restore the original state and give up.*/
      /*fprintf(stderr,"Best cost (%7i) > 0! Optimization failed.\n",best_cost);*/
      od_mv_dp_restore_row_state(dp_node);
    }
    else{
      int bi;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
      if(daala_granule_basetime(state,state->cur_time)==ANI_FRAME){
        char iter_label[16];
        od_mv_dp_restore_row_state(dp_node);
        od_mv_dp_animate_state(state,_ref,dp_node,0);
        od_mv_dp_install_row_state(dp_node+1,best_si);
        od_state_mc_predict(state,_ref);
        od_state_fill_vis(state);
        sprintf(iter_label,"ani%08i",state->ani_iter++);
        od_state_dump_img(state,&state->vis_img,iter_label);
      }
#endif
      /*Update the state along the optimal path.*/
      od_mv_dp_install_row_state(dp_node+1,best_si);
      /*Store the SADs from this last node, too.*/
      for(bi=0;bi<dp_node[1].nblocks;bi++){
        dp_node[1].blocks[bi]->sad=block_sads[best_si][bi];
      }
      dcost+=best_cost;
    }
  }
  od_mv_est_check_rd_state(_est,_ref,_mv_res);
  return dcost;
}

/*Column refinement.*/

static void od_mv_dp_col_init(od_mv_est_ctx *_est,od_mv_dp_node *_dp,
 int _vx,int _vy,od_mv_dp_node *_prev_dp){
  od_state      *state;
  int            nhmvbs;
  int            nvmvbs;
  state=&_est->enc->state;
  _dp->mv=_est->mvs[_vy]+_vx;
  _dp->mvg=state->mv_grid[_vy]+_vx;
  _dp->original_mv[0]=_dp->mvg->mv[0];
  _dp->original_mv[1]=_dp->mvg->mv[1];
  _dp->original_etype=_dp->mvg->down;
  _dp->original_mv_rate=_dp->mv->mv_rate;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  if(_vx<2||_vx>nhmvbs-2||_vy<2||_vy>nvmvbs-2){
    /*Strictly speaking, we may be used to predict others, but since our MV
       can't possibly change, neither can their rate.*/
    _dp->npredicted=_dp->npred_changeable=0;
    /*No one else is used to predict us, or any other MV we predict.
      However, we may still need to load the previous MV into the grid to
       estimate our SADs properly.*/
    _dp->min_predictor_node=_prev_dp;
  }
  else{
    int level;
    int pred_hist;
    int npred;
    int nchangeable;
    int pi;
    /*Get the list of MVs we help predict.*/
    level=OD_MC_LEVEL[_vy&3][_vx&3];
    /*fprintf(stderr,"Initializing node (%i,%i) [%i,%i] at level %i:\n",
     _vx,_vy,_vx-2<<2,_vy-2<<2,level);*/
    npred=nchangeable=0;
    for(pi=0;pi<OD_NPREDICTED[level];pi++){
      int px;
      int py;
      px=_vx+OD_COL_PREDICTED[level][pi][0];
      if(px<2||px>nhmvbs-2)continue;
      py=_vy+OD_COL_PREDICTED[level][pi][1];
      if(py<2||py>nvmvbs-2)continue;
      if(state->mv_grid[py][px].valid){
        /*fprintf(stderr,"Adding (%i,%i) as a PREDICTED MV.\n",px,py);*/
        _dp->predicted_mvgs[npred]=state->mv_grid[py]+px;
        _dp->predicted_mvs[npred]=_est->mvs[py]+px;
        if(pi<OD_NCOL_PRED_CHANGEABLE[level]){
          /*fprintf(stderr,"It is CHANGEABLE.\n");*/
          _dp->original_mv_rates[npred]=_est->mvs[py][px].mv_rate;
          nchangeable++;
        }
        npred++;
      }
    }
    _dp->npredicted=npred;
    _dp->npred_changeable=nchangeable;
    /*Now, figure out the earliest DP node that influences our own prediction,
       or that of one of the other MVs we predict.*/
    pred_hist=OD_ROW_PRED_HIST_SIZE[level];
    if(_prev_dp!=NULL&&_prev_dp->mv->vy>=_vy-pred_hist){
      od_mv_dp_node *dp_pred;
      for(dp_pred=_prev_dp;dp_pred->mv->vy>_vy-pred_hist&&
       dp_pred->states[0].prevsi>=0;dp_pred--);
      if(dp_pred->mv->vy<_vy-pred_hist)dp_pred++;
      _dp->min_predictor_node=dp_pred;
      /*fprintf(stderr,"State will be restored back to (%i,%i).\n",
       dp_pred->mv->vx,dp_pred->mv->vy);*/
    }
    else _dp->min_predictor_node=NULL;
  }
}

static void od_mv_dp_first_col_block_setup(od_mv_est_ctx *_est,
 od_mv_dp_node *_dp,int _vx,int _vy){
  od_state *state;
  int       nhmvbs;
  int       level;
  int       log_mvb_sz;
  int       mvb_sz;
  int       nblocks;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(_vy>2){
    if(level>=3){
      if(_vx>=mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      if(_vx<=nhmvbs-mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx;
    }
    else{
      int half_mvb_sz;
      int mvb_off;
      half_mvb_sz=mvb_sz>>1;
      if(_vx>=mvb_sz){
        if(state->mv_grid[_vy-half_mvb_sz][_vx-half_mvb_sz].valid){
          if(level>0||
           !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
            mvb_off=half_mvb_sz;
          }
          else mvb_off=half_mvb_sz>>1;
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-mvb_off;
          if(!state->mv_grid[_vy][_vx-mvb_off].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-(mvb_off<<1);
          }
          if(!state->mv_grid[_vy-mvb_off][_vx].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx-mvb_off;
          }
        }
        else _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      }
      if(_vx<=nhmvbs-mvb_sz){
        if(state->mv_grid[_vy-half_mvb_sz][_vx+half_mvb_sz].valid){
          if(level>0||
           !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
            mvb_off=half_mvb_sz;
          }
          else mvb_off=half_mvb_sz>>1;
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx;
          if(!state->mv_grid[_vy][_vx+mvb_off].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx+mvb_off;
          }
          if(!state->mv_grid[_vy-mvb_off][_vx].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx;
          }
        }
        else _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx;
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_prev_col_block_setup(od_mv_est_ctx *_est,
 od_mv_dp_node *_dp,int _vx,int _vy){
  od_state *state;
  int       nhmvbs;
  int       level;
  int       prev_level;
  int       log_mvb_sz;
  int       prev_log_mvb_sz;
  int       mvb_sz;
  int       nblocks;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  prev_level=OD_MC_LEVEL[_vy-mvb_sz&3][_vx&3];
  prev_log_mvb_sz=4-prev_level>>1;
  nblocks=0;
  if(level>=3){
    if(_vx>=mvb_sz){
      _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
      if(prev_log_mvb_sz>log_mvb_sz&&
       !state->mv_grid[_vy-mvb_sz][_vx-mvb_sz].valid){
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-(mvb_sz<<1);
      }
    }
    if(_vx<=nhmvbs-mvb_sz){
      _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx;
      if(prev_log_mvb_sz>log_mvb_sz&&
       !state->mv_grid[_vy-mvb_sz][_vx+mvb_sz].valid){
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx+mvb_sz;
      }
    }
  }
  else{
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vx>=mvb_sz){
      if(state->mv_grid[_vy-half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-mvb_off;
        if(!state->mv_grid[_vy][_vx-mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx-(mvb_off<<1);
        }
        if(!state->mv_grid[_vy-mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx-mvb_off;
          if(!state->mv_grid[_vy-(mvb_off<<1)][_vx-mvb_off].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx-(mvb_off<<1);
          }
        }
      }
      else{
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-mvb_sz;
        if(prev_log_mvb_sz>log_mvb_sz&&
         !state->mv_grid[_vy-mvb_sz][_vx-mvb_sz].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx-(mvb_sz<<1);
        }
      }
    }
    if(_vx<=nhmvbs-mvb_sz){
      if(state->mv_grid[_vy-half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy-(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx;
        if(!state->mv_grid[_vy][_vx+mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_off]+_vx+mvb_off;
        }
        if(!state->mv_grid[_vy-mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx;
          if(!state->mv_grid[_vy-(mvb_off<<1)][_vx+mvb_off].valid){
            _dp->blocks[nblocks++]=_est->mvs[_vy-(mvb_off<<1)]+_vx+mvb_off;
          }
        }
      }
      else{
        _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx;
        if(prev_log_mvb_sz>log_mvb_sz&&
         !state->mv_grid[_vy-mvb_sz][_vx+mvb_sz].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy-mvb_sz]+_vx+mvb_sz;
        }
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_last_col_block_setup(od_mv_est_ctx *_est,
 od_mv_dp_node *_dp,int _vx,int _vy){
  od_state *state;
  int       nhmvbs;
  int       level;
  int       log_mvb_sz;
  int       mvb_sz;
  int       nblocks;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(level>=3){
    if(_vx>=mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_sz;
    if(_vx<=nhmvbs-mvb_sz)_dp->blocks[nblocks++]=_est->mvs[_vy]+_vx;
  }
  else{
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vx>=mvb_sz){
      if(state->mv_grid[_vy+half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy+(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_off;
        if(!state->mv_grid[_vy][_vx-mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-(mvb_off<<1);
        }
        if(!state->mv_grid[_vy+mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_off]+_vx-mvb_off;
        }
      }
      else _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx-mvb_sz;
    }
    if(_vx<=nhmvbs-mvb_sz){
      if(state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !state->mv_grid[_vy+(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          mvb_off=half_mvb_sz;
        }
        else mvb_off=half_mvb_sz>>1;
        _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx;
        if(!state->mv_grid[_vy][_vx+mvb_off].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx+mvb_off;
        }
        if(!state->mv_grid[_vy+mvb_off][_vx].valid){
          _dp->blocks[nblocks++]=_est->mvs[_vy+mvb_off]+_vx;
        }
      }
      else _dp->blocks[nblocks++]=_est->mvs[_vy]+_vx;
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_restore_col_state(od_mv_dp_node *_dp){
  od_mv_grid_pt *mvg;
  int            pi;
  do{
    /*Restore the state for this MV itself.*/
    _dp->mv->mv_rate=_dp->original_mv_rate;
    mvg=_dp->mvg;
    mvg->mv[0]=_dp->original_mv[0];
    mvg->mv[1]=_dp->original_mv[1];
    mvg->down=_dp->original_etype;
    for(pi=0;pi<_dp->npred_changeable;pi++){
      /*Restore the state for the MVs this one predicted.*/
      _dp->predicted_mvs[pi]->mv_rate=_dp->original_mv_rates[pi];
    }
  }
  while((_dp--)->states[0].prevsi>=0);
}

static void od_mv_dp_install_col_state(od_mv_dp_node *_dp,int _prevsi){
  od_mv_dp_node *dp;
  od_mv_grid_pt *mvg;
  int            nextsi;
  int            si;
  int            pi;
  int            bi;
  /*We must install the state going FORWARDS, since the pred_mv_rates may have
     changed several times over the course of the trellis.
    Therefore, first we reverse all of the prevsi pointers to make them act
     like nextsi pointers.
    We _can_ update the edge type flags here, however, and it is much more
     convenient to do so while going backwards than forwards.*/
  nextsi=-1;
  for(dp=_dp,si=_prevsi;si>=0;si=_prevsi){
    dp--;
    /*fprintf(stderr,"Node %i, prevsi: %i nextsi: %i\n",_dp-dp,_prevsi,nextsi);*/
    if(si>=dp->nstates){
      dp->mvg->down=1;
      si-=dp->nstates;
    }
    else dp->mvg->down=0;
    _prevsi=dp->states[si].prevsi;
    dp->states[si].prevsi=nextsi;
    nextsi=si;
  }
  /*Now we traverse forward installing the rest of the state.*/
  for(si=nextsi;dp<_dp;dp++){
    /*Install the state for this MV itself.*/
    dp->mv->mv_rate=dp->states[si].mv_rate;
    mvg=dp->mvg;
    mvg->mv[0]=dp->states[si].mv[0];
    mvg->mv[1]=dp->states[si].mv[1];
    /*Install the new block SADs.*/
    for(bi=0;bi<dp->nblocks;bi++){
      dp->blocks[bi]->sad=dp->states[si].block_sads[bi];
    }
    /*Install the state for the MVs this one predicted.*/
    for(pi=0;pi<dp->npredicted;pi++){
      dp->predicted_mvs[pi]->mv_rate=dp->states[si].pred_mv_rates[pi];
    }
    si=dp->states[si].prevsi;
  }
}

static ogg_int32_t od_mv_est_refine_col(od_mv_est_ctx *_est,int _ref,int _vx,
 int _log_dsz,int _mv_res,const int *_pattern_nsites,
 const od_pattern *_pattern){
  od_state        *state;
  od_mv_grid_pt  **grid;
  od_mv_grid_pt   *pmvg;
  od_mv_grid_pt   *mvg;
  od_mv_dp_node   *dp_node;
  od_mv_dp_state  *cstate;
  od_mv_dp_state  *pstate;
  ogg_int32_t      dcost;
  int              nhmvbs;
  int              nvmvbs;
  int              level;
  int              log_mvb_sz;
  int              mvb_sz;
  int              labels_only;
  int              nsites;
  int              sitei;
  int              site;
  int              curx;
  int              cury;
  int              vy;
  int              b;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  labels_only=_vx<2||_vx>nhmvbs-2;
  grid=state->mv_grid;
  dcost=0;
  /*fprintf(stderr,"Refining column %i (%i)...\n",_vx,_vx-2<<2);*/
  for(vy=0;;vy++){
    ogg_int32_t block_sads[18][8];
    ogg_int32_t best_cost;
    ogg_int32_t cost;
    ogg_int32_t best_dd;
    ogg_int32_t dd;
    int         cur_mv_rates[9];
    int         pred_mv_rates[9][17];
    int         best_dr;
    int         dr;
    int         best_si;
    int         si;
    int         has_gap;
    for(;vy<=nvmvbs&&!grid[vy][_vx].valid;vy++);
    if(vy>nvmvbs)break;
    level=OD_MC_LEVEL[vy&3][_vx&3];
    /*fprintf(stderr,"Starting DP at vertex %i (%i), level %i\n",
     vy,vy-2<<2,level);*/
    log_mvb_sz=4-level>>1;
    mvb_sz=1<<log_mvb_sz;
    mvg=grid[vy]+_vx;
    curx=mvg->mv[0];
    cury=mvg->mv[1];
    dp_node=_est->dp_nodes;
    od_mv_dp_col_init(_est,dp_node,_vx,vy,NULL);
    od_mv_dp_first_col_block_setup(_est,dp_node,_vx,vy);
    /*fprintf(stderr,"TESTING block SADs:\n");
    od_mv_dp_get_sad_change8(_est,_ref,dp_node,block_sads[0]);*/
    /*Compute the set of states for the first node.*/
    if(vy>=2&&!labels_only){
      b=od_mv_est_get_boundary_case(state,_vx,vy,curx,cury,
       1<<_log_dsz,log_mvb_sz+2);
      nsites=_pattern_nsites[b];
    }
    else b=nsites=0;
    for(sitei=0,site=4;;sitei++){
      cstate=dp_node[0].states+sitei;
      cstate->mv[0]=curx+(OD_SITE_DX[site]<<_log_dsz);
      cstate->mv[1]=cury+(OD_SITE_DY[site]<<_log_dsz);
      cstate->prevsi=-1;
      mvg->mv[0]=cstate->mv[0];
      mvg->mv[1]=cstate->mv[1];
      cstate->dr=od_mv_dp_get_rate_change(state,dp_node,
       &cstate->mv_rate,cstate->pred_mv_rates,-1,_mv_res);
      cstate->dd=od_mv_dp_get_sad_change8(_est,_ref,dp_node,
       cstate->block_sads);
      /*fprintf(stderr,"State: %i  dr: %i  dd: %i  dopt: %i\n",
       sitei,cstate->dr,cstate->dd,
       cstate->dr*_est->lambda+(cstate->dd<<OD_LAMBDA_SCALE));*/
      if(sitei>=nsites)break;
      site=_pattern[b][sitei];
    }
    dp_node[0].nstates=nsites+1;
    has_gap=0;
    pmvg=mvg;
    for(;vy<nvmvbs;){
      /*Find the next available MV to advance to.*/
      if(level&1){
        if(!grid[vy+mvb_sz][_vx].valid){
          /*fprintf(stderr,"Gap found at %i (%i), stopping\n",vy,vy-2<<2);*/
          has_gap=1;
          break;
        }
        else if(level>=3)vy++;
        else if(!grid[vy+1][_vx].valid)vy+=mvb_sz;
        else vy++;
      }
      else if(level>=4)vy++;
      else if(!grid[vy+(mvb_sz>>1)][_vx].valid)vy+=mvb_sz;
      else if(level>=2||!grid[vy+1][_vx].valid)vy+=mvb_sz>>1;
      else vy++;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
      if(daala_granule_basetime(state,state->cur_time)==ANI_FRAME){
        od_mv_dp_restore_col_state(dp_node);
        od_mv_dp_animate_state(state,_ref,dp_node,0);
      }
#endif
      level=OD_MC_LEVEL[vy&3][_vx&3];
      /*fprintf(stderr,"Continuing DP at vertex %i (%i), level %i\n",
       vy,vy-2<<2,level);*/
      log_mvb_sz=4-level>>1;
      mvb_sz=1<<log_mvb_sz;
      mvg=grid[vy]+_vx;
      curx=mvg->mv[0];
      cury=mvg->mv[1];
      od_mv_dp_col_init(_est,dp_node+1,_vx,vy,dp_node);
      od_mv_dp_prev_col_block_setup(_est,dp_node+1,_vx,vy);
      /*fprintf(stderr,"TESTING block SADs:\n");
      pmvg->mv[0]=dp_node[0].original_mv[0];
      pmvg->mv[0]=dp_node[0].original_mv[0];
      od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,block_sads[0]);*/
      /*Compute the set of states for this node.*/
      if(vy<=nvmvbs-2&&!labels_only){
        b=od_mv_est_get_boundary_case(state,_vx,vy,curx,cury,
         1<<_log_dsz,log_mvb_sz+2);
        nsites=_pattern_nsites[b];
      }
      else nsites=0;
      for(sitei=0,site=4;;sitei++){
        cstate=dp_node[1].states+sitei;
        cstate->mv[0]=curx+(OD_SITE_DX[site]<<_log_dsz);
        cstate->mv[1]=cury+(OD_SITE_DY[site]<<_log_dsz);
        best_si=pmvg->down?dp_node[0].nstates:0;
        best_dr=dp_node[0].states[0].dr;
        best_dd=dp_node[0].states[0].dd;
        best_cost=INT_MAX;
        mvg->mv[0]=cstate->mv[0];
        mvg->mv[1]=cstate->mv[1];
        for(si=0;si<dp_node[0].nstates;si++){
          pstate=dp_node[0].states+si;
          /*Get the rate change for this state using previous state si.
            This automatically loads the required bits of the trellis path into
             the grid, like the previous MV.*/
          cstate->dr=od_mv_dp_get_rate_change(state,dp_node+1,
            cur_mv_rates+si,pred_mv_rates[si],si,_mv_res);
          /*Test against the previous state with a B edge.*/
          if(_est->flags&OD_MC_USEB){
            pmvg->down=0;
            dr=pstate->dr+cstate->dr;
            /*Account for label mismatch.*/
            /*if(pstate->prevsi>=0&&pstate->prevsi>=(dp_node-1)->nstates)dr+=2;*/
            dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
             block_sads[si]);
            cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
            /*fprintf(stderr,"State: %2i  P.State: %iB  dr: %3i  dd: %6i  dopt: %7i\n",
             sitei,si,dr,dd,cost);*/
            if(cost<best_cost){
              best_si=si;
              best_cost=cost;
              best_dd=dd;
              best_dr=dr;
            }
          }
          if(_est->flags&OD_MC_USEV){
            /*Test against the previous state with a V edge.*/
            pmvg->down=1;
            dr=pstate->dr+cstate->dr;
            /*Account for label mismatch.*/
            /*if(pstate->prevsi>=0&&pstate->prevsi<(dp_node-1)->nstates)rate+=2;*/
            dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
             block_sads[si+dp_node[0].nstates]);
            cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
            /*fprintf(stderr,"State: %2i  P.State: %iV  dr: %3i  dd: %6i  dopt: %7i\n",
             sitei,si,dr,dd,cost);*/
            if(cost<best_cost){
              best_si=si+dp_node[0].nstates;
              best_cost=cost;
              best_dd=dd;
              best_dr=dr;
            }
          }
        }
        /*fprintf(stderr,"State: %2i  Best P.State: %i%c\n",
         sitei,best_si%dp_node[0].nstates,best_si>dp_node[0].nstates?'V':'B');*/
        cstate->prevsi=best_si;
        cstate->dr=best_dr;
        cstate->dd=best_dd;
        memcpy(cstate->block_sads,block_sads[cstate->prevsi],
         sizeof(block_sads[0][0])*dp_node[1].nblocks);
        if(best_si<dp_node[0].nstates){
          cstate->mv_rate=cur_mv_rates[best_si];
          memcpy(cstate->pred_mv_rates,pred_mv_rates[best_si],
           sizeof(pred_mv_rates[0][0])*dp_node[1].npredicted);
        }
        else{
          cstate->mv_rate=cur_mv_rates[best_si-dp_node[0].nstates];
          memcpy(cstate->pred_mv_rates,
           pred_mv_rates[best_si-dp_node[0].nstates],
           sizeof(pred_mv_rates[0][0])*dp_node[1].npredicted);
        }
        if(sitei>=nsites)break;
        site=_pattern[b][sitei];
      }
      dp_node[1].nstates=nsites+1;
      dp_node++;
      pmvg=mvg;
    }
    /*fprintf(stderr,"Finished DP at vertex %i (%i)\n",
     dp_node[0].mv->vy,dp_node[0].mv->vy-2<<2);*/
    best_si=pmvg->down?dp_node[0].nstates:0;
    best_cost=INT_MAX;
    /*TODO: Once we stop optimizing at arbitrary places, we'll need to
       compute the rate change of MVs we didn't get to.*/
    dp_node[1].npredicted=dp_node[1].npred_changeable=0;
    if(dp_node[0].mv->vy<nvmvbs-2){
      od_mv_dp_last_col_block_setup(_est,dp_node+1,_vx,dp_node[0].mv->vy);
      /*fprintf(stderr,"TESTING block SADs:\n");
      pmvg->mv[0]=dp_node[0].original_mv[0];
      pmvg->mv[1]=dp_node[0].original_mv[1];
      od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,block_sads[0]);*/
      for(si=0;si<dp_node[0].nstates;si++){
        pstate=dp_node[0].states+si;
        pmvg->mv[0]=pstate->mv[0];
        pmvg->mv[1]=pstate->mv[1];
        /*Test against the state with a following B edge.*/
        if(_est->flags&OD_MC_USEB){
          pmvg->down=0;
          dr=pstate->dr;
          /*Account for label mismatch.*/
          /*if(!has_gap&&pstate->prevsi>=0&&pstate->prevsi>=(dp_node-1)->nstates){
            rate+=2;
          }*/
          dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
           block_sads[si]);
          cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: --  P.State: %iB  dr: %3i  dd: %6i  dopt: %7i\n",
           si,dr,dd,cost);*/
          if(best_si<0||cost<best_cost){
            best_si=si;
            best_cost=cost;
          }
        }
        /*Test against the state with a following V edge.
          If we hit a gap, then the edge label does not actually affect
           anything, so we can skip these tests if we did the ones above.*/
        if((_est->flags&OD_MC_USEV)&&(!has_gap||!(_est->flags&OD_MC_USEB))){
          pmvg->down=1;
          dr=pstate->dr;
          /*Account for label mismatch.*/
          /*if(pstate->prevsi>=0&&pstate->prevsi<(dp_node-1)->nstates)rate+=2;*/
          dd=pstate->dd+od_mv_dp_get_sad_change8(_est,_ref,dp_node+1,
           block_sads[si+dp_node[0].nstates]);
          cost=dr*_est->lambda+(dd<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: --  P.State: %iV  dr: %3i  dd: %6i  dopt: %7i\n",
           si,dr,dd,cost);*/
          if(cost<best_cost){
            best_si=si+dp_node[0].nstates;
            best_cost=cost;
          }
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else{
      dp_node[1].nblocks=0;
      for(si=0;si<dp_node[0].nstates;si++){
        dp_node[1].nblocks=0;
        pstate=dp_node[0].states+si;
        cost=pstate->dr*_est->lambda+(pstate->dd<<OD_LAMBDA_SCALE);
        if(best_si<0||cost<best_cost){
          best_si=si;
          best_cost=cost;
        }
      }
      if(pmvg->down)best_si+=dp_node[0].nstates;
    }
    /*fprintf(stderr,"Best P.State: %i%c dopt: %7i\n",
     best_si%dp_node[0].nstates,best_si>dp_node[0].nstates?'V':'B',best_cost);*/
    if(best_cost>0){
      /*Our optimal path is worse than what we started with!
        Restore the original state and give up.*/
      /*fprintf(stderr,"Best cost (%7i) > 0! Optimization failed.\n",best_cost);*/
      od_mv_dp_restore_col_state(dp_node);
    }
    else{
      int bi;
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
      if(daala_granule_basetime(state,state->cur_time)==ANI_FRAME){
        char iter_label[16];
        od_mv_dp_restore_col_state(dp_node);
        od_mv_dp_animate_state(state,_ref,dp_node,0);
        od_mv_dp_install_col_state(dp_node+1,best_si);
        od_state_mc_predict(state,_ref);
        od_state_fill_vis(state);
        sprintf(iter_label,"ani%08i",state->ani_iter++);
        od_state_dump_img(state,&state->vis_img,iter_label);
      }
#endif
      /*Update the state along the optimal path.*/
      od_mv_dp_install_col_state(dp_node+1,best_si);
      /*Store the SADs from this last node, too.*/
      for(bi=0;bi<dp_node[1].nblocks;bi++){
        dp_node[1].blocks[bi]->sad=block_sads[best_si][bi];
      }
      dcost+=best_cost;
    }
  }
  od_mv_est_check_rd_state(_est,_ref,_mv_res);
  return dcost;
}

#if 0
static void od_mv_dp_first_col_block_setup(od_state *_state,od_mv_dp_node *_dp,
 int _vx,int _vy){
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  nhmvbs=_state->nhmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(level>=3){
    if(_vx>=mvb_sz){
      od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
       _vx-mvb_sz,_vy-mvb_sz,log_mvb_sz);
    }
    if(_vx<=nhmvbs-mvb_sz){
      od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
       _vx,_vy-mvb_sz,log_mvb_sz);
    }
  }
  else{
    int half_mvb_sz;
    int log_mvb_off;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vx>=mvb_sz){
      if(_state->mv_grid[_vy-half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !_state->mv_grid[_vy-(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          log_mvb_off=log_mvb_sz-1;
          mvb_off=half_mvb_sz;
        }
        else{
          mvb_off=half_mvb_sz>>1;
          log_mvb_off=log_mvb_sz-2;
        }
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_off,_vy-mvb_off,log_mvb_off);
        /*No need to test _state->mv_grid[_vy-mvb_off][_vx].
          If it was valid, we wouldn't be here.*/
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_off,_vy-(mvb_off<<1),log_mvb_off);
        if(!_state->mv_grid[_vy][_vx-mvb_off].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx-(mvb_off<<1),_vy-mvb_off,log_mvb_off);
        }
      }
      else{
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_sz,_vy-mvb_sz,log_mvb_sz);
      }
    }
    if(_vx<=nhmvbs-mvb_sz){
      if(_state->mv_grid[_vy-half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !_state->mv_grid[_vy-(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          log_mvb_off=log_mvb_sz-1;
          mvb_off=half_mvb_sz;
        }
        else{
          mvb_off=half_mvb_sz>>1;
          log_mvb_off=log_mvb_sz-2;
        }
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy-mvb_off,log_mvb_off);
        /*No need to test _state->mv_grid[_vy-mvb_off][_vx].
          If it was valid, we wouldn't be here.*/
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy-(mvb_off<<1),log_mvb_off);
        if(!_state->mv_grid[_vy][_vx+mvb_off].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx+mvb_off,_vy-mvb_off,log_mvb_off);
        }
      }
      else{
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy-mvb_sz,log_mvb_sz);
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_prev_col_block_setup(od_state *_state,od_mv_dp_node *_dp,
 int _vx,int _vy){
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  nhmvbs=_state->nhmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(level>=3){
    if(_vx>=mvb_sz){
      od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
       _vx-mvb_sz,_vy-mvb_sz,log_mvb_sz);
    }
    if(_vx<=nhmvbs-mvb_sz){
      od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
       _vx,_vy-mvb_sz,log_mvb_sz);
    }
  }
  else{
    int half_mvb_sz;
    int log_mvb_off;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vx>=mvb_sz){
      if(_state->mv_grid[_vy-half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !_state->mv_grid[_vy-(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          log_mvb_off=log_mvb_sz-1;
          mvb_off=half_mvb_sz;
        }
        else{
          mvb_off=half_mvb_sz>>1;
          log_mvb_off=log_mvb_sz-2;
        }
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_off,_vy-mvb_off,log_mvb_off);
        if(!_state->mv_grid[_vy][_vx-mvb_off].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx-(mvb_off<<1),_vy-mvb_off,log_mvb_off);
        }
        if(!_state->mv_grid[_vy-mvb_off][_vx].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx-mvb_off,_vy-(mvb_off<<1),log_mvb_off);
          if(!_state->mv_grid[_vy-(mvb_off<<1)][_vx-mvb_off].valid){
            od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
             _vx-(mvb_off<<1),_vy-(mvb_off<<1),log_mvb_off);
          }
        }
      }
      else{
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_sz,_vy-mvb_sz,log_mvb_sz);
      }
    }
    if(_vx<=nhmvbs-mvb_sz){
      if(_state->mv_grid[_vy-half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !_state->mv_grid[_vy-(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          log_mvb_off=log_mvb_sz-1;
          mvb_off=half_mvb_sz;
        }
        else{
          mvb_off=half_mvb_sz>>1;
          log_mvb_off=log_mvb_sz-2;
        }
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy-mvb_off,log_mvb_off);
        if(!_state->mv_grid[_vy][_vx+mvb_off].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx+mvb_off,_vy-mvb_off,log_mvb_off);
        }
        if(!_state->mv_grid[_vy-mvb_off][_vx].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx,_vy-(mvb_off<<1),log_mvb_off);
          if(!_state->mv_grid[_vy-(mvb_off<<1)][_vx+mvb_off].valid){
            od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
             _vx+mvb_off,_vy-(mvb_off<<1),log_mvb_off);
          }
        }
      }
      else{
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy-mvb_sz,log_mvb_sz);
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_dp_last_col_block_setup(od_state *_state,od_mv_dp_node *_dp,
 int _vx,int _vy){
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  nhmvbs=_state->nhmbs+1<<2;
  level=OD_MC_LEVEL[_vy&3][_vx&3];
  log_mvb_sz=4-level>>1;
  mvb_sz=1<<log_mvb_sz;
  nblocks=0;
  if(level>=3){
    if(_vx>=mvb_sz){
      od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
       _vx-mvb_sz,_vy,log_mvb_sz);
    }
    if(_vx<=nhmvbs-mvb_sz){
      od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
       _vx,_vy,log_mvb_sz);
    }
  }
  else{
    int half_mvb_sz;
    int log_mvb_off;
    int mvb_off;
    half_mvb_sz=mvb_sz>>1;
    if(_vx>=mvb_sz){
      if(_state->mv_grid[_vy+half_mvb_sz][_vx-half_mvb_sz].valid){
        if(level>0||
         !_state->mv_grid[_vy+(half_mvb_sz>>1)][_vx-(half_mvb_sz>>1)].valid){
          log_mvb_off=log_mvb_sz-1;
          mvb_off=half_mvb_sz;
        }
        else{
          mvb_off=half_mvb_sz>>1;
          log_mvb_off=log_mvb_sz-2;
        }
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_off,_vy,log_mvb_off);
        if(!_state->mv_grid[_vy][_vx-mvb_off].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx-(mvb_off<<1),_vy,log_mvb_off);
        }
        if(!_state->mv_grid[_vy+mvb_off][_vx].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx-mvb_off,_vy+mvb_off,log_mvb_off);
        }
      }
      else{
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx-mvb_sz,_vy,log_mvb_sz);
      }
    }
    if(_vx<=nhmvbs-mvb_sz){
      if(_state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
        if(level>0||
         !_state->mv_grid[_vy+(half_mvb_sz>>1)][_vx+(half_mvb_sz>>1)].valid){
          log_mvb_off=log_mvb_sz-1;
          mvb_off=half_mvb_sz;
        }
        else{
          mvb_off=half_mvb_sz>>1;
          log_mvb_off=log_mvb_sz-2;
        }
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy,log_mvb_off);
        if(!_state->mv_grid[_vy][_vx+mvb_off].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx+mvb_off,_vy,log_mvb_off);
        }
        if(!_state->mv_grid[_vy+mvb_off][_vx].valid){
          od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
           _vx,_vy+mvb_off,log_mvb_off);
        }
      }
      else{
        od_mv_dp_setup_block(_state,_dp->blocks+nblocks++,
         _vx,_vy,log_mvb_sz);
      }
    }
  }
  _dp->nblocks=nblocks;
}

static void od_mv_est_refine_col(od_mv_est_ctx *_est,int _ref,int _vx,
 int _log_dsz,const int *_pattern_nsites,const od_pattern *_pattern){
  od_state        *state;
  od_mv_grid_pt  **grid;
  od_mv_grid_pt   *pmvg;
  od_mv_grid_pt   *mvg;
  od_mv_dp_node   *dp_node;
  od_mv_dp_state  *cstate;
  od_mv_dp_state  *pstate;
  int              pred[2];
  int              nhmvbs;
  int              nvmvbs;
  int              mv_res;
  int              level;
  int              log_mvb_sz;
  int              mvb_sz;
  int              labels_only;
  int              nsites;
  int              sitei;
  int              site;
  int              prevx;
  int              prevy;
  int              petype;
  int              curx;
  int              cury;
  int              vy;
  int              b;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  labels_only=_vx<2||_vx>nhmvbs-2;
  mv_res=OD_MINI(_log_dsz,2);
  /*fprintf(stderr,"Refining col %i (%i)...\n",_vx,_vx-2<<2);*/
  grid=state->mv_grid;
  for(vy=0;;vy++){
    ogg_int32_t best_cost;
    ogg_int32_t cost;
    ogg_int32_t best_sad;
    ogg_int32_t sad;
    int         best_rate;
    int         rate;
    int         best_si;
    int         si;
    int         has_gap;
    for(;vy<=nvmvbs&&!grid[vy][_vx].valid;vy++);
    if(vy>nvmvbs)break;
    level=OD_MC_LEVEL[vy&3][_vx&3];
    /*fprintf(stderr,"Starting DP at vertex %i (%i), level %i\n",
     vy,vy-2<<2,level);*/
    log_mvb_sz=4-level>>1;
    mvb_sz=1<<log_mvb_sz;
    mvg=grid[vy]+_vx;
    curx=mvg->mv[0];
    cury=mvg->mv[1];
    od_state_get_predictor(state,pred,_vx,vy,level,mv_res);
    dp_node=_est->dp_nodes;
    dp_node[0].vi=vy;
    if(vy>=2&&!labels_only){
      b=od_mv_est_get_boundary_case(state,_vx,vy,curx,cury,
       1<<_log_dsz,log_mvb_sz+2);
      nsites=_pattern_nsites[b];
    }
    else nsites=0;
    /*Compute the set of states for the first node.*/
    if(vy>2&&!labels_only)od_mv_dp_first_col_block_setup(state,dp_node,_vx,vy);
    else dp_node->nblocks=0;
    for(sitei=0,site=4;;sitei++){
      cstate=dp_node[0].states+sitei;
      cstate->mv[0]=curx+(OD_SITE_DX[site]<<_log_dsz);
      cstate->mv[1]=cury+(OD_SITE_DY[site]<<_log_dsz);
      cstate->prevsi=-1;
      if(vy<=2||labels_only){
        cstate->sad=0;
        cstate->rate=0;
      }
      else{
        cstate->rate=od_mv_est_bits(
         cstate->mv[0]>>mv_res,cstate->mv[1]>>mv_res,pred[0],pred[1]);
        mvg->mv[0]=cstate->mv[0];
        mvg->mv[1]=cstate->mv[1];
        cstate->sad=od_mv_dp_sad8(_est,_ref,dp_node);
      }
      /*fprintf(stderr,"State: %i  Rate: %i  SAD: %i\n",
       sitei,cstate->rate,cstate->sad);*/
      if(sitei>=nsites)break;
      site=_pattern[b][sitei];
    }
    dp_node[0].nstates=nsites+1;
    mvg->mv[0]=curx;
    mvg->mv[1]=cury;
    has_gap=0;
    pmvg=mvg;
    for(;vy<nvmvbs;){
      od_mv_dp_node *pred_node;
      od_mv_grid_pt *pred_mvg;
      int            predx;
      int            predy;
      prevx=curx;
      prevy=cury;
      petype=pmvg->down;
      if(level&1){
        vy+=mvb_sz;
        if(vy>nvmvbs||!grid[vy][_vx].valid){
          /*fprintf(stderr,"Gap found at %i (%i), stopping\n",vy,vy-2<<2);*/
          has_gap=1;
          break;
        }
      }
      else if(level>=4)vy++;
      else if(!grid[vy+(mvb_sz>>1)][_vx].valid)vy+=mvb_sz;
      else if(level>=2||!grid[vy+1][_vx].valid)vy+=mvb_sz>>1;
      else vy++;
      level=OD_MC_LEVEL[vy&3][_vx&3];
      /*fprintf(stderr,"Continuing DP at vertex %i (%i), level %i\n",
       vy,vy-2<<2,level);*/
      log_mvb_sz=4-level>>1;
      mvb_sz=1<<log_mvb_sz;
      mvg=grid[vy]+_vx;
      curx=mvg->mv[0];
      cury=mvg->mv[1];
      if((level&1)||vy<2||vy>nvmvbs-2||labels_only)pred_node=NULL;
      else pred_node=od_mv_dp_get_pred_node(dp_node,vy-mvb_sz);
      if(pred_node==NULL){
        od_state_get_predictor(state,pred,_vx,vy,level,mv_res);
      }
      else{
        pred_mvg=grid[pred_node->vi]+_vx;
        predx=pred_mvg->mv[0];
        predy=pred_mvg->mv[1];
      }
      dp_node[1].vi=vy;
      if(vy<=nvmvbs-2&&!labels_only){
        b=od_mv_est_get_boundary_case(state,_vx,vy,curx,cury,
         1<<_log_dsz,log_mvb_sz+2);
        nsites=_pattern_nsites[b];
      }
      /*The first node of every 4th column has its motion vector fixed.
        Do not move it.*/
      else nsites=0;
      /*Compute the set of states for the first node.*/
      od_mv_dp_prev_col_block_setup(state,dp_node+1,_vx,vy);
      for(sitei=0,site=4;;sitei++){
        cstate=dp_node[1].states+sitei;
        cstate->mv[0]=curx+(OD_SITE_DX[site]<<_log_dsz);
        cstate->mv[1]=cury+(OD_SITE_DY[site]<<_log_dsz);
        cstate->prevsi=-1;
        /*If no previous node is used as a predictor for this vector, compute
           this MV's rate, once.*/
        if(pred_node==NULL){
          cstate->rate=od_mv_est_bits(
           cstate->mv[0]>>mv_res,cstate->mv[1]>>mv_res,pred[0],pred[1]);
        }
        mvg->mv[0]=cstate->mv[0];
        mvg->mv[1]=cstate->mv[1];
        for(si=0;si<dp_node[0].nstates;si++){
          pstate=dp_node[0].states+si;
          pmvg->mv[0]=pstate->mv[0];
          pmvg->mv[1]=pstate->mv[1];
          if(pred_node!=NULL){
            od_mv_dp_node *pnode;
            int            pred_si;
            /*Find the state in the predictor that would be active if we chose
               this state in the current node.*/
            pred_si=si;
            for(pnode=dp_node;pnode!=pred_node;pnode--){
              pred_si=pnode->states[pred_si].prevsi;
            }
            /*Compute the new predictor and MV rate.*/
            pred_mvg->mv[0]=pred_node->states[pred_si].mv[0];
            pred_mvg->mv[1]=pred_node->states[pred_si].mv[1];
            od_state_get_predictor(state,pred,_vx,vy,level,mv_res);
            cstate->rate=od_mv_est_bits(
             cstate->mv[0]>>mv_res,cstate->mv[1]>>mv_res,pred[0],pred[1]);
          }
          /*Test against the previous state with a B edge.*/
          pmvg->down=0;
          rate=pstate->rate+cstate->rate;
          /*Account for label mismatch.*/
          if(pstate->prevsi>=0&&pstate->prevsi>=(dp_node-1)->nstates)rate+=2;
          sad=pstate->sad+od_mv_dp_sad8(_est,_ref,dp_node+1);
          cost=rate*_est->lambda+(sad<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: %2i  P.State: %iB  Rate: %3i  SAD: %6i  Cost: %7i\n",
           sitei,si,rate,sad,cost);*/
          if(cstate->prevsi<0||cost<best_cost){
            cstate->prevsi=si;
            best_cost=cost;
            best_sad=sad;
            best_rate=rate;
          }
          /*Test against the previous state with a V edge.*/
          pmvg->down=1;
          rate=pstate->rate+cstate->rate;
          /*Account for label mismatch.*/
          if(pstate->prevsi>=0&&pstate->prevsi<(dp_node-1)->nstates)rate+=2;
          sad=pstate->sad+od_mv_dp_sad8(_est,_ref,dp_node+1);
          cost=rate*_est->lambda+(sad<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: %2i  P.State: %iV  Rate: %3i  SAD: %6i  Cost: %7i\n",
           sitei,si,rate,sad,cost);*/
          if(cost<best_cost){
            cstate->prevsi=si+dp_node[0].nstates;
            best_cost=cost;
            best_sad=sad;
            best_rate=rate;
          }
        }
        /*fprintf(stderr,"State: %2i  Best P.State: %i%c\n",
         sitei,cstate->prevsi%dp_node[0].nstates,
         cstate->prevsi>dp_node[0].nstates?'V':'B');*/
        cstate->rate=best_rate;
        cstate->sad=best_sad;
        if(sitei>=nsites)break;
        site=_pattern[b][sitei];
      }
      dp_node[1].nstates=nsites+1;
      dp_node++;
      /*Restore the state we were optimizing.*/
      if(pred_node!=NULL){
        pred_mvg->mv[0]=predx;
        pred_mvg->mv[1]=predy;
      }
      pmvg->mv[0]=prevx;
      pmvg->mv[1]=prevy;
      pmvg->down=petype;
      mvg->mv[0]=curx;
      mvg->mv[1]=cury;
      pmvg=mvg;
    }
    /*fprintf(stderr,"Finished DP at vertex %i (%i)\n",
     dp_node[0].vi,dp_node[0].vi-2<<2);*/
    best_si=-1;
    if(dp_node[0].vi<nvmvbs-2){
      od_mv_dp_last_col_block_setup(state,dp_node+1,_vx,dp_node[0].vi);
      for(si=0;si<dp_node[0].nstates;si++){
        pstate=dp_node[0].states+si;
        pmvg->mv[0]=pstate->mv[0];
        pmvg->mv[1]=pstate->mv[1];
        /*Test against the state with a following B edge.*/
        pmvg->down=0;
        rate=pstate->rate;
        /*Account for label mismatch.*/
        if(!has_gap&&pstate->prevsi>=0&&pstate->prevsi>=(dp_node-1)->nstates){
          rate+=2;
        }
        sad=pstate->sad+od_mv_dp_sad8(_est,_ref,dp_node+1);
        cost=rate*_est->lambda+(sad<<OD_LAMBDA_SCALE);
        /*fprintf(stderr,"State: --  P.State: %iB  Rate: %3i  SAD: %6i  Cost: %7i\n",
         si,pstate->rate,sad,cost);*/
        if(best_si<0||cost<best_cost){
          best_si=si;
          best_cost=cost;
        }
        /*Test against the state with a following V edge.
          If we hit a gap, then the edge label does not actually affect anything,
           so we can skip these tests.*/
        if(!has_gap){
          pmvg->down=1;
          rate=pstate->rate;
          /*Account for label mismatch.*/
          if(pstate->prevsi>=0&&pstate->prevsi<(dp_node-1)->nstates)rate+=2;
          sad=pstate->sad+od_mv_dp_sad8(_est,_ref,dp_node+1);
          cost=rate*_est->lambda+(sad<<OD_LAMBDA_SCALE);
          /*fprintf(stderr,"State: --  P.State: %iV  Rate: %3i  SAD: %6i  Cost: %7i\n",
           si,pstate->rate,sad,cost);*/
          if(cost<best_cost){
            best_si=si+dp_node[0].nstates;
            best_cost=cost;
          }
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else for(si=0;si<dp_node[0].nstates;si++){
      pstate=dp_node[0].states+si;
      cost=pstate->rate*_est->lambda+(pstate->sad<<OD_LAMBDA_SCALE);
      if(best_si<0||cost<best_cost){
        best_si=si;
        best_cost=cost;
      }
    }
    /*fprintf(stderr,"Best P.State: %i%c\n",
     best_si%dp_node[0].nstates,best_si>dp_node[0].nstates?'V':'B');*/
    /*Update the MV state along the optimal path.*/
    for(;;){
      if(best_si>=dp_node[0].nstates){
        pmvg->down=1;
        best_si-=dp_node[0].nstates;
      }
      else pmvg->down=0;
      pstate=dp_node[0].states+best_si;
      pmvg->mv[0]=pstate->mv[0];
      pmvg->mv[1]=pstate->mv[1];
      best_si=pstate->prevsi;
      if(best_si<0)break;
      dp_node--;
      pmvg=grid[dp_node[0].vi]+_vx;
    }
  }
}
#endif

static ogg_int32_t od_mv_est_refine(od_mv_est_ctx *_est,int _ref,int _log_dsz,
 int _mv_res,const int *_pattern_nsites,const od_pattern *_pattern){
  od_state    *state;
  ogg_int32_t  dcost;
  int          nhmvbs;
  int          nvmvbs;
  int          vx;
  int          vy;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  fprintf(stderr,
   "Refining with displacements of %0g and 1/%i pel MV resolution.\n",
   (1<<_log_dsz)*0.125,1<<3-_mv_res);
  dcost=0;
  for(vy=0;vy<=nvmvbs;vy++)if(_est->row_counts[vy]){
    dcost+=od_mv_est_refine_row(_est,_ref,vy,_log_dsz,_mv_res,
     _pattern_nsites,_pattern);
  }
  for(vx=0;vx<=nhmvbs;vx++)if(_est->col_counts[vx]){
    dcost+=od_mv_est_refine_col(_est,_ref,vx,_log_dsz,_mv_res,
     _pattern_nsites,_pattern);
  }
  return dcost;
}



/*STAGE 4: Sub-pel Refinement.*/



/*Stores the full-pel MVs for use by EPZS^2 in the next frame before sub-pel
   refinement.*/
void od_mv_est_update_fullpel_mvs(od_mv_est_ctx *_est,int _ref){
  od_state *state;
  int       nhmvbs;
  int       nvmvbs;
  int       vx;
  int       vy;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  for(vy=2;vy<=nvmvbs-2;vy++)for(vx=2;vx<=nhmvbs-2;vx++){
    od_mv_grid_pt *mvg;
    od_mv_node    *mv;
    mvg=state->mv_grid[vy]+vx;
    if(!mvg->valid)continue;
    mv=_est->mvs[vy]+vx;
    mv->mvs[0][_ref][0]=mvg->mv[0]>>3;
    mv->mvs[0][_ref][1]=mvg->mv[1]>>3;
  }
}

/*Sets the mv_rate of each node in the mesh, using the given MV resolution.
  Returns the change in rate.*/
int od_mv_est_update_mv_rates(od_mv_est_ctx *_est,int _mv_res){
  od_state *state;
  int       nhmvbs;
  int       nvmvbs;
  int       vx;
  int       vy;
  int       dr;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  dr=0;
  for(vy=2;vy<=nvmvbs-2;vy++)for(vx=2;vx<=nhmvbs-2;vx++){
    od_mv_grid_pt *mvg;
    od_mv_node    *mv;
    int            pred[2];
    mvg=state->mv_grid[vy]+vx;
    if(!mvg->valid)continue;
    mv=_est->mvs[vy]+vx;
    od_state_get_predictor(state,pred,vx,vy,OD_MC_LEVEL[vy&3][vx&3],_mv_res);
    dr-=mv->mv_rate;
    mv->mv_rate=od_mv_est_bits(mvg->mv[0]>>_mv_res,mvg->mv[1]>>_mv_res,
     pred[0],pred[1]);
    dr+=mv->mv_rate;
  }
  return dr;
}


od_mv_est_ctx *od_mv_est_alloc(od_enc_ctx *_enc){
  od_mv_est_ctx *ret;
  ret=(od_mv_est_ctx *)_ogg_malloc(sizeof(od_mv_est_ctx));
  od_mv_est_init(ret,_enc);
  return ret;
}

void od_mv_est_free(od_mv_est_ctx *_est){
  if(_est!=NULL){
    od_mv_est_clear(_est);
    _ogg_free(_est);
  }
}

void od_mv_subpel_refine(od_mv_est_ctx *_est,int _ref,int _cost_thresh){
  od_state       *state;
  od_mv_grid_pt **grid;
  ogg_int32_t     dcost;
  ogg_int32_t     subpel_cost;
  int             cost_thresh;
  int             nhmvbs;
  int             nvmvbs;
  int             mv_res;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  cost_thresh=_cost_thresh;
  /*Save the fullpell MVs now for use by EPZS^2 on the next frame.
    We could also try rounding the results after refinement, I guess.
    I'm not sure it makes much difference*/
  od_mv_est_update_fullpel_mvs(_est,_ref);
  do dcost=od_mv_est_refine(_est,_ref,2,2,OD_DIAMOND_NSITES,OD_DIAMOND_SITES);
  while(dcost<cost_thresh);
  for(mv_res=2;mv_res-->_est->mv_res_min;){
    subpel_cost=od_mv_est_update_mv_rates(_est,mv_res)*_est->lambda;
    /*If the rate penalty for refining is small, bump the termination threshold
       down to make sure we actually get a decent improvement.
      We make sure not to let it get too small, however, so we're not here all
       day (a motion field of all (0,0)'s would have a rate penalty of 0!).*/
    cost_thresh=OD_MAXI(cost_thresh,-OD_MAXI(subpel_cost,16<<OD_LAMBDA_SCALE));
    memcpy(_est->refine_grid[0],state->mv_grid[0],
     sizeof(state->mv_grid[0][0])*(nhmvbs+1)*(nvmvbs+1));
    do{
      dcost=od_mv_est_refine(_est,_ref,mv_res,mv_res,
       OD_DIAMOND_NSITES,OD_DIAMOND_SITES);
      subpel_cost+=dcost;
    }
    while(dcost<cost_thresh);
    if(subpel_cost>0){
      fprintf(stderr,"1/%i refinement FAILED:    dopt %7i\n",
       1<<3-mv_res,subpel_cost);
      grid=_est->refine_grid;
      _est->refine_grid=state->mv_grid;
      state->mv_grid=grid;
      break;
    }
    else fprintf(stderr,"1/%i refinement SUCCEEDED: dopt %7i\n",
     1<<3-mv_res,subpel_cost);
  }
}

void od_mv_est(od_mv_est_ctx *_est,int _ref,int _lambda){
  od_state     *state;
  od_img_plane *iplane;
  ogg_int32_t   dcost;
  int           cost_thresh;
  int           nhmvbs;
  int           nvmvbs;
  int           pli;
  state=&_est->enc->state;
  nhmvbs=state->nhmbs+1<<2;
  nvmvbs=state->nvmbs+1<<2;
  iplane=state->input.planes+0;
  /*If the luma plane is decimated for some reason, then our distortions will
     be smaller, so scale lambda appropriately.*/
  _est->lambda=_lambda>>iplane->xdec+iplane->ydec;
  /*Compute termination thresholds for EPZS^2.*/
  _est->thresh1[0]=16>>iplane->xdec+iplane->ydec;
  _est->thresh1[1]=64>>iplane->xdec+iplane->ydec;
  _est->thresh1[2]=256>>iplane->xdec+iplane->ydec;
  /*If we're using the chroma planes, then our distortions will be larger.
    Compensate by increasing lambda and the termination thresholds.*/
  if(_est->flags&OD_MC_USE_CHROMA){
    for(pli=1;pli<state->input.nplanes;pli++){
      iplane=state->input.planes+pli;
      _est->lambda+=_lambda>>iplane->xdec+iplane->ydec+OD_MC_CHROMA_SCALE;
      _est->thresh1[0]+=16>>iplane->xdec+iplane->ydec+OD_MC_CHROMA_SCALE;
      _est->thresh1[1]+=64>>iplane->xdec+iplane->ydec+OD_MC_CHROMA_SCALE;
      _est->thresh1[2]+=256>>iplane->xdec+iplane->ydec+OD_MC_CHROMA_SCALE;
    }
  }
  _est->thresh2_offs[0]=_est->thresh1[0]>>1;
  _est->thresh2_offs[1]=_est->thresh1[1]>>1;
  _est->thresh2_offs[2]=_est->thresh1[2]>>1;
  /*Accelerated predictor weights.*/
  _est->mvapw[_ref][0]=0x20000;
  _est->mvapw[_ref][1]=0x10000;
  /*TODO: Constant velocity predictor weight.*/
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
  /*Set some initial state.
    This would get reset eventually by the algorithm in a more convenient
     place, but is needed earlier by the visualization.*/
  if(daala_granule_basetime(&_est->enc->state,_est->enc->state.cur_time)==
   ANI_FRAME){
    int vx;
    int vy;
    for(vy=0;vy<nvmvbs;vy++){
      od_mv_grid_pt *grid;
      grid=_est->enc->state.mv_grid[vy];
      for(vx=0;vx<nhmvbs;vx++){
        grid[vx].valid=0;
        grid[vx].right=0;
        grid[vx].down=0;
        grid[vx].mv[0]=0;
        grid[vx].mv[1]=1;
      }
    }
  }
#endif
  od_mv_est_init_mvs(_est,_ref);
  od_mv_est_decimate(_est,_ref);
  /*This threshold is somewhat arbitrary.
    Chen and Willson use 6000 (with SSD as an error metric).
    We would like something more dependent on the frame size.
    For CIF, there are a maximum of 6992 vertices in the mesh, which is pretty
     close to 6000.
    With a SAD error metric like we use, the square root of 6000 would be a
     more appropriate value, however that gives a PSNR improvement of less than
     0.01 dB, and requires almost twice as many iterations to achieve.*/
  cost_thresh=-nhmvbs*nvmvbs<<OD_LAMBDA_SCALE;
#if 0
  /*Logarithmic search.*/
  do{
    dcost=0;
    dcost+=od_mv_est_refine(_est,_ref,5,2,OD_SQUARE_NSITES,OD_SQUARE_SITES);
    dcost+=od_mv_est_refine(_est,_ref,4,2,OD_SQUARE_NSITES,OD_SQUARE_SITES);
    dcost+=od_mv_est_refine(_est,_ref,3,2,OD_SQUARE_NSITES,OD_SQUARE_SITES);
  }
  while(dcost<cost_thresh);
#else
  /*Diamond search.
    This appears to give the same quality as the logarithmic search, but at
     nearly 10 times the speed.*/
  do dcost=od_mv_est_refine(_est,_ref,3,2,OD_DIAMOND_NSITES,OD_DIAMOND_SITES);
  while(dcost<cost_thresh);
#endif
  od_mv_subpel_refine(_est,_ref,cost_thresh);
  {
    int vx;
    int vy;
    static int l0flags[4][4];
    static int lhflags[2][8];
    static int lvflags[2][8];
    int pred;
    int flags;
    /*memset(l0flags,0,sizeof(l0flags));
    memset(lhflags,0,sizeof(lhflags));
    memset(lvflags,0,sizeof(lvflags));*/
    for(vy=0;vy<=nvmvbs;vy+=4)for(vx=0;vx<=nhmvbs;vx+=4){
      if(vy>0){
        if(vx>0){
          pred=state->mv_grid[vy][vx-4].right|state->mv_grid[vy-4][vx].down<<1;
        }
        else{
          pred=state->mv_grid[vy-4][vx].down|state->mv_grid[vy-4][vx].down<<1;
        }
      }
      else if(vx>0){
        pred=state->mv_grid[vy][vx-4].right|state->mv_grid[vy][vx-4].right<<1;
      }
      else pred=0;
      flags=state->mv_grid[vy][vx].right|state->mv_grid[vy][vx].down<<1;
      l0flags[pred][flags]++;
    }
    for(vy=0;vy<=nvmvbs;vy+=4)for(vx=2;vx<=nhmvbs;vx+=4){
      if(!state->mv_grid[vy][vx].valid)continue;
      pred=state->mv_grid[vy][vx-2].right;
      flags=state->mv_grid[vy][vx].right<<2;
      if(vy>0)flags|=state->mv_grid[vy-2][vx].down;
      if(vy<nvmvbs)flags|=state->mv_grid[vy][vx].down<<1;
      lhflags[pred][flags]++;
    }
    for(vy=2;vy<=nvmvbs;vy+=4)for(vx=0;vx<=nhmvbs;vx+=4){
      if(!state->mv_grid[vy][vx].valid)continue;
      pred=state->mv_grid[vy-2][vx].down;
      flags=state->mv_grid[vy][vx].down<<2;
      if(vx>0)flags|=state->mv_grid[vy][vx-2].right;
      if(vx<nhmvbs)flags|=state->mv_grid[vy][vx].right<<1;
      lvflags[pred][flags]++;
    }
    for(vy=0;vy<=nvmvbs;vy+=2)for(vx=1;vx<=nhmvbs;vx+=2){
      if(!state->mv_grid[vy][vx].valid)continue;
      pred=state->mv_grid[vy][vx-1].right;
      flags=state->mv_grid[vy][vx].right<<2;
      if(vy>0)flags|=state->mv_grid[vy-1][vx].down;
      if(vy<nvmvbs)flags|=state->mv_grid[vy][vx].down<<1;
      lhflags[pred][flags]++;
    }
    for(vy=1;vy<=nvmvbs;vy+=2)for(vx=0;vx<=nhmvbs;vx+=2){
      if(!state->mv_grid[vy][vx].valid)continue;
      pred=state->mv_grid[vy-1][vx].down;
      flags=state->mv_grid[vy][vx].down<<2;
      if(vx>0)flags|=state->mv_grid[vy][vx-1].right;
      if(vx<nhmvbs)flags|=state->mv_grid[vy][vx].right<<1;
      lvflags[pred][flags]++;
    }
    for(pred=0;pred<4;pred++){
      fprintf(stderr,"l0flags[%i|%i] = %6i:%6i:%6i:%6i\n",pred&1,pred&2,
       l0flags[pred][0],l0flags[pred][1],l0flags[pred][2],l0flags[pred][3]);
    }
    for(pred=0;pred<2;pred++){
      fprintf(stderr,"lhflags[%i] = %6i:%6i:%6i:%6i:%6i:%6i:%6i:%6i\n",pred,
       lhflags[pred][0],lhflags[pred][1],lhflags[pred][2],lhflags[pred][3],
       lhflags[pred][4],lhflags[pred][5],lhflags[pred][6],lhflags[pred][7]);
    }
    for(pred=0;pred<2;pred++){
      fprintf(stderr,"lvflags[%i] = %6i:%6i:%6i:%6i:%6i:%6i:%6i:%6i\n",pred,
       lvflags[pred][0],lvflags[pred][1],lvflags[pred][2],lvflags[pred][3],
       lvflags[pred][4],lvflags[pred][5],lvflags[pred][6],lvflags[pred][7]);
    }
  }
}
