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

#include <stdlib.h>
#include <limits.h>
#include "int_search.h"
#include "../src/internal.h"

int int_simplex_max(double *_y,int _dims,isearch_obj _obj,const double *_aux,
 const int *_lb,const int *_ub,int *_x){
  double        best_obj;
  int           dim;
  int           last_dim;
  int           steps;
  OD_ASSERT(_dims>0&&!!_x&&!!_lb&&!!_ub&&!!_obj);
  last_dim=_dims;
  best_obj=_obj(_aux,_x);
  dim=steps=0;
  do{
    double y;
    int i;
    int dir;
    dir=i=0;
    do{
      if(_x[i]>_lb[i]&&(i!=last_dim)){
        _x[i]--;
        y=_obj(_aux,_x);
        if(y>best_obj){
          best_obj=y;
          dim=i;
          dir=-1;
        }
        _x[i]++;
      }
      if(_x[i]<_ub[i]&&(-i!=last_dim)){
        _x[i]++;
        y=_obj(_aux,_x);
        if(y>best_obj){
          best_obj=y;
          dim=i;
          dir=1;
        }
        _x[i]--;
      }
    }while(++i<_dims);
    if(dir==0)break;
    steps++;
    _x[dim]+=dir;
    last_dim=dir<0?-dim:dim;
    dir=0;
  }while(1);
  if(_y)*_y=best_obj;
  return steps;
}
