/*Daala video codec
Copyright (c) 2016 Daala project contributors.  All rights reserved.

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

#include <math.h>
#if !defined(M_PI)
# define M_PI (3.141592653589793238462643)
#endif
#include <stdint.h>
#include <stdio.h>

#define OD_MINI(a,b) ((a) < (b) ? (a) : (b))
#define OD_MULT16_16_Q15(a, b) (((a)*(int)(b)) >> 15)
#define OD_MULT16_16_P15(a, b) (((a)*(int)(b) + 16384) >> 15)
#define OD_MULT16_16_Q16(a, b) (((a)*(int)(b)) >> 16)

static const int32_t C0[4] = {1073758164, -7654, 16573, -2529};

static int16_t od_pvq_cos_pi_2(int16_t x) {
  int16_t x2;
  x2 = OD_MULT16_16_Q15(x, x);
  return OD_MINI(32767, (C0[0] - x*x + x2*(C0[1] + OD_MULT16_16_Q16(x2,
   C0[2] + OD_MULT16_16_Q16(C0[3], x2)))) >> 15);
}

static const int32_t C[4] = {((1 << 30) + (1 << 14)), -7651, 16554, -2504};

static int16_t od_pvq_cos2(int16_t x) {
  int16_t x2;

  x2 = OD_MULT16_16_Q15(x, x);
  return OD_MINI(32767, (C[0] - x*x + x2*
   (C[1] + OD_MULT16_16_Q16(x2, C[2] + OD_MULT16_16_Q16(C[3], x2)))) >> 15);
}

#define NS (10)

int main()
{
   int i;
   int32_t best[4];
   double best_dist;
   int truth[32768];
   double dtruth[32768];
   int32_t c[4];
   best_dist = 1e10;
   for (i = 0; i < 32768; i++) truth[i] = floor(.5 + 32768*cos(i*M_PI/65536.));
   for (i = 0; i < 32768; i++) dtruth[i] = 32768*cos(i*M_PI/65536.);
   for (c[0] = C0[0] - NS; c[0] <= C0[0] + NS; c[0]++) {
     for (c[1] = C0[1] - NS; c[1] <= C0[1] + NS; c[1]++) {
       for (c[2] = C0[2] - NS; c[2] <= C0[2] + NS; c[2]++) {
         for (c[3] = C0[3] - NS; c[3] <= C0[3] + NS; c[3]++) {
           double dist;
           dist = 0;
           for (i = 0; i < 32768; i++) {
             double err;
             err = od_pvq_cos2(i) - dtruth[i];
             dist += err*err;
           }
           if (dist < best_dist) {
             for (i = 0; i < 4; i++) best[i] = c[i];
             best_dist = dist;
           }
         }
       }
     }
   }
   for (i = 0; i < 4; i++) c[i] = best[i];
   fprintf(stderr, "%d %d %d %d\n", c[0], c[1], c[2], c[3]);
   for (i = 0; i < 32768; i++) {
     int cos1;
     int cos2;
     cos1 = od_pvq_cos_pi_2(i);
     cos2 = od_pvq_cos2(i);
     printf("%d %f %d %d\n", truth[i], dtruth[i], cos2, cos1);
   }
   return 0;
}
