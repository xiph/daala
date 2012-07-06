#if !defined(_intra_H)
# define _intra_H (1)
# include "filter.h"



# define OD_INTRA_NMODES (10)

# define OD_INTRA_NCONTEXTS (8)

extern double OD_INTRA_PRED_WEIGHTS_4x4[OD_INTRA_NMODES][4][4][2*4][2*4];



/*Applies intra prediction to a 4x4 block of coefficients at _c, using
   UL, U, and L blocks of reconstructed 4x4 coefficients.
  On input:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    original coefficients (before quantization).
   The other blocks contain reconstructed (post quantization and unprediction)
    coefficients.
  On output:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    the input coefficients with the prediction subtracted.
  Return: The intra prediction mode used (0...OD_INTRA_NMODES-1).*/
int od_intra_pred4x4_apply(od_coeff *_c,int _stride);

/*Unapplies intra prediction to a 4x4 block of coefficients at _c, using
   UL, U, and L blocks of reconstructed 4x4 coefficients.
  On input:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    unquantized coefficients.
   The other blocks contain reconstructed (post quantization and unprediction)
    coefficients.
  On output:
   {{_c[0],_c[1],_c[2],_c[3]},{_c[_stride],_c[stride+1],...}} contains
    the input coefficients with the prediction added.
  Return: The intra prediction mode used (0...OD_INTRA_NMODES-1).*/
void od_intra_pred4x4_unapply(od_coeff *_c,int _stride,int _mode);

void od_intra_pred_cdf(ogg_uint16_t cdf[OD_INTRA_NMODES],
    const unsigned char probs[OD_INTRA_NMODES][OD_INTRA_NCONTEXTS],
    const ogg_uint16_t p0[OD_INTRA_NMODES],int left,int upleft,int up);

int od_intra_pred_search(ogg_uint16_t p0[OD_INTRA_NMODES],
    const ogg_uint16_t cdf[OD_INTRA_NMODES],
    const ogg_uint16_t dist[OD_INTRA_NMODES], ogg_uint16_t lambda, int left,
    int upleft, int up);

#endif
