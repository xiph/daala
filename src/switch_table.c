#include <ogg/os_types.h>

/* Initial CDFs for block size encoding. These are the result of
   ad hoc "training" based on akiyo and parkjoy. They can probably
   be much better, although the influence on actual performance is
   likely to be small. */
const ogg_uint16_t range_cdf_init[7] = {
  575, 2594, 7146, 7803, 13036, 22902, 32768
};

const ogg_uint16_t split16_cdf_init[2][16] = {
  {282, 2224, 3949, 6344, 8968, 13421, 14345, 16271,
  17931, 18631, 26412, 27753, 29187, 30851, 32486, 32768},
  {282, 964, 1903, 3495, 4235, 6688, 7709, 11817,
  12513, 13887, 16583, 21650, 23152, 27718, 32486, 32768}
};
