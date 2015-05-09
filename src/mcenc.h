/*Daala video codec
Copyright (c) 2014 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#if !defined(_mcenc_H)
# define _mcenc_H (1)

/*Flag indicating we include the chroma planes in our SAD calculations.*/
# define OD_MC_USE_CHROMA (1 << 0)

/* The maximum search range for BMA. Also controls hit cache size. */
#define OD_MC_SEARCH_RANGE (64)

typedef struct od_mv_node od_mv_node;
typedef struct od_mv_dp_state od_mv_dp_state;
typedef struct od_mv_dp_node od_mv_dp_node;

# include "mc.h"
# include "encint.h"

typedef ogg_uint16_t od_sad4[4];

/*The state information used by the motion estimation process that is not
   required by the decoder.
  Some of this information corresponds to a vertex in the MV mesh.
  Other pieces correspond to a block whose upper-left corner is located at that
   vertex.*/
struct od_mv_node {
  /*The historical motion vectors for EPZS^2, stored at full-pel resolution.
    Indexed by [time][reference_type][component].*/
  int bma_mvs[3][2][2];
  /*The current estimated rate of this MV.*/
  unsigned mv_rate:16;
  /*The current size of the block with this MV at its upper-left.*/
  unsigned log_mvb_sz:2;
  /*The index of the exterior corner of that block.*/
  unsigned oc:2;
  /*The edge splitting index of that block.*/
  unsigned s:2;
  /*The current distortion of that block.*/
  ogg_int32_t sad;
  /*The SAD for BMA predictor centered on this node.
    Used for the dynamic thresholds of the initial EPZS^2 pass.*/
  ogg_int32_t bma_sad;
  /*The location of this node in the grid.*/
  int vx;
  int vy;
  /*The change in global distortion for decimating this node.*/
  ogg_int32_t dd;
  /*The change in global rate for decimating this node.*/
  int dr;
  /*The position of this node in the heap.*/
  int heapi;
};

/*The square pattern, the largest we use, has 9 states.*/
# define OD_DP_NSTATES_MAX (9)
/*Up to 8 blocks can be influenced by this MV and the previous MV.*/
# define OD_DP_NBLOCKS_MAX (8)
/*Up to 28 MVs can be predicted by this one, but 4 of those are MVs on the
   DP trellis whose value we have yet to determine.*/
# define OD_DP_NPREDICTED_MAX (24)
/*At most 8 of them can be changed by a subsequent MV on the DP path.*/
# define OD_DP_NCHANGEABLE_MAX (8)

/*One of the trellis states in the dynamic program.*/
struct od_mv_dp_state {
  /*The MV to install for this state.*/
  int mv[2];
  /*The best state in the previous DP node to use with this one, or -1 to
     indicate the start of the path.*/
  int prevsi;
  /*The total rate change (thus far) produced by choosing this path.*/
  int dr;
  /*The total distortion change (thus far) produced by choosing this path.*/
  ogg_int32_t dd;
  /*The new SAD of each block affected by the the DP between this node and the
     previous node.
    These are installed if the path is selected.*/
  ogg_int32_t block_sads[OD_DP_NBLOCKS_MAX];
  /*The new rate of each MV predicted by this node.
    These are installed if the path is selected.
    These may supersede the rates reported in previous nodes on the path.*/
  int pred_mv_rates[OD_DP_NPREDICTED_MAX];
  /*The new rate of this MV.*/
  int mv_rate;
};

/*A node on the dynamic programming path.*/
struct od_mv_dp_node {
  od_mv_grid_pt *mvg;
  od_mv_node *mv;
  /*The number of states considered in this node.*/
  int nstates;
  /*The number of blocks affected by states in this node.*/
  int nblocks;
  /*The number of MVs predicted by this node.*/
  int npredicted;
  /*The number of those MVs that are potentially changeable by future DP
     states.*/
  int npred_changeable;
  /*The original MV used by this node.*/
  int original_mv[2];
  /*The original rate of this MV.*/
  int original_mv_rate;
  /*The original MV rates before predictors were changed by this node.
    This only includes the ones that are actually changeable.*/
  int original_mv_rates[OD_DP_NCHANGEABLE_MAX];
  /*The last node we save/restore in order to perform prediction.*/
  od_mv_dp_node *min_predictor_node;
  /*The set of trellis states.*/
  od_mv_dp_state states[OD_DP_NSTATES_MAX];
  /*The blocks influenced by this MV and the previous MV.*/
  od_mv_node *blocks[OD_DP_NBLOCKS_MAX];
  /*The vertices whose MV we predict.*/
  od_mv_grid_pt *predicted_mvgs[OD_DP_NPREDICTED_MAX];
  od_mv_node *predicted_mvs[OD_DP_NPREDICTED_MAX];
};

struct od_mv_est_ctx {
  od_enc_ctx *enc;
  /*A cache of the SAD values used during decimation.
    Indexed by [log_mvb_sz][vy >> log_mvb_sz][vx >> log_mvb_sz][s], where s is
     the edge split state.
    The SAD of top-level blocks (log_mvb_sz == OD_LOG_MVB_DELTA0) is not stored
     in this cache, since it is only needed once.*/
  od_sad4 **sad_cache[OD_LOG_MVB_DELTA0];
  /*The state of the MV mesh specific to the encoder.*/
  od_mv_node **mvs;
  /*A temporary copy of the decoder-side MV grid used to save-and-restore the
     MVs when attempting sub-pel refinement.*/
  od_mv_grid_pt **refine_grid;
  /*Space for storing the Viterbi trellis used for DP refinement.*/
  od_mv_dp_node *dp_nodes;
  /*The decimation heap.*/
  od_mv_node **dec_heap;
  /*The number of vertices in the decimation heap.*/
  int dec_nheap;
  /*The number of undecimated vertices in each row.*/
  unsigned *row_counts;
  /*The number of undecimated vertices in each column.*/
  unsigned *col_counts;
  /*The maximum SAD value for accepting set A predictors for each block size.*/
  int thresh1[OD_NMVBSIZES];
  /*The offsets to inflate the second threshold by for each block size.*/
  int thresh2_offs[OD_NMVBSIZES];
  /*The weights used to produce the accelerated MV predictor.*/
  ogg_int32_t mvapw[2][2];
  /*Flags indicating which MVs have already been tested during the initial
     EPZS^2 pass.*/
  unsigned char hit_cache[OD_MC_SEARCH_RANGE*2][OD_MC_SEARCH_RANGE*2];
  /*The flag used by the current EPZS search iteration.*/
  unsigned hit_bit;
  /*The Lagrangian multiplier used for R-D optimization.*/
  int lambda;
  /*Rate estimations (in units of OD_BITRES).*/
  int mv_small_rate_est[5][16];
  /*Configuration.*/
  /*The flags indicating which feature to use.*/
  int flags;
  /*The smallest resolution to refine MVs to.*/
  int mv_res_min;
  /*The deepest level to refine to (inclusive).*/
  int level_max;
  /*The shallowest level to decimate to (inclusive).*/
  int level_min;
  /*This function is set switchable between SAD and SATD based on the current
     stage of ME/MC. At present, SAD function is called for stage 1, 2, and 3,
     and SATD functions are called for stage 4 (i.e. sub-pel refine).*/
  int (*compute_distortion)(od_enc_ctx *enc, const unsigned char *p,
   int pystride, int pxstride, int pli, int x, int y, int log_blk_sz);
};

#endif
