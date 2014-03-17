/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_od_defs_H)
# define _od_defs_H (0)

#include "intra_fit_tools.h"
#include "../src/filter.h"

#if defined(_OPENMP)
# include <omp.h>
# define NUM_PROCS (4)
# define OD_OMP_GET_THREAD (omp_get_thread_num())
# define OD_OMP_SET_THREADS(_x) (omp_set_num_threads(_x))
#else
# pragma GCC warning "Compiler enviroment does not support OpenMP. The training tools will be very slow."
# define NUM_PROCS (1)
# define OD_OMP_GET_THREAD (0)
# define OD_OMP_SET_THREADS(_x)
#endif

#define FAST_MATH (1)

#define FILTER_BITS (4)

#define USE_TYPE3 (1)
#define NN_SEARCH (1)

#define APPLY_FILTER (1)
#define APPLY_DCT    (1)
#define APPLY_PRED   (1)
#define APPLY_PCA    (1)

#define INPUT_SCALE_BITS (4)
#define INPUT_SCALE (1<<INPUT_SCALE_BITS)

#define ZERO_MEAN   (1)
#define MASK_BLOCKS (0)
#define TF_BLOCKS   (1)

#endif
