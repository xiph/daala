#if !defined(_cpack_linalg_matidx_H)
# define _cpack_linalg_matidx_H (1)

/*Returns the total number of elements in an (_n)x(_n) lower-triangular
   matrix stored in row-major order.*/
#define LT_SZ(_n)         (((_n)*((_n)+1)>>1))
/*Returns the index of element (_i,_j) in a lower-triangular matrix stored in
   row-major order ((_i)>=(_j)).*/
#define LT_IDX(_i,_j)     (LT_SZ(_i)+(_j))

/*Returns the total number of elements in an (_n)x(_n) strictly
   lower-triangular matrix stored in row-major order.*/
#define SLT_SZ(_n)        ((((_n)-1)*(_n)>>1))
/*Returns the index of element (_i,_j) in a strictly lower-triangular matrix
   stored in row-major order ((_i)>(_j)).*/
#define SLT_IDX(_i,_j)    (SLT_SZ(_i)+(_j))

/*Returns the total number of elements in an (_m)x(_n) upper-triangular
   matrix, (_m)<=(_n), stored in row-major order.*/
#define UT_SZ(_m,_n)      (LT_SZ(_n)-LT_SZ((_n)-(_m)))
/*Returns the index of element (_i,_j) in an upper-triangular matrix stored in
   row-major order ((_i)<=(_j)).*/
#define UT_IDX(_i,_j,_n)  (UT_SZ(_i,_n)+(_j)-(_i))
/*This version is fewer operations, but the compiler should be able to optimize
   the above version better when _n is fixed.
#define UT_IDX(_i,_j,_n)  (((_i)*(((_n)<<1)-(_i)-1)>>1)+(_j))*/

/*Returns the total number of elements in an (_m)x(_n) strictly
   upper-triangular matrix, (_m)<=(_n), stored in row-major order.*/
#define SUT_SZ(_m,_n)     (SLT_SZ(_n)-SLT_SZ((_n)-(_m)))
/*Returns the index of element (_i,_j) in a strictly upper-triangular matrix
   stored in row-major order ((_i)<(_j)).*/
#define SUT_IDX(_i,_j,_n) (SUT_SZ(_i,_n)+(_j)-(_i)-1)

/*Swaps elements _a and _b, storing _a in _t as a temporary.*/
#define CP_SWAP(_a,_b,_t) \
  do{ \
    (_t)=(_a); \
    (_a)=(_b); \
    (_b)=(_t); \
  } \
  while(0)

#endif
