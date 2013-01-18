
#include <ogg/ogg.h>

typedef struct {
  int vstride;
  ogg_int32_t Sxx[17][17];
  ogg_int32_t Sx[17][17];
  ogg_int16_t *var;
  ogg_int16_t *var_1;

} od_switch_decision_state;

/* Allocate memory to process an image of width w */
void od_switch_init(od_switch_decision_state *st, int w);

/* Free state memory */
void od_switch_free(od_switch_decision_state *st);

/* Needs to be called at the beginning of a frame, with img pointing to the top-left of the first superblock of the current frame */
void od_switch_start_frame(od_switch_decision_state *st, const unsigned char *img, int s);

/* Needs to be called at the beginning of a frame, with img pointing to the top-left of the first superblock of the current line */
void od_switch_start_line(od_switch_decision_state *st);

/* Computes switching decision for macro block pointed to by img, returns a decision for each 8x8 block in the superblock */
void od_switch_process_superblock(od_switch_decision_state *st, const unsigned char *img, int id, int s, int blocksize[4][4]);

int switch_decision(const unsigned char *img, int w, int h, int stride);
