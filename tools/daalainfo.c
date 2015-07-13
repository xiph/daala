/* daalainfo - based on Ogginfo
 *
 * A tool to describe ogg file contents and metadata.
 *
 * Copyright 2002-2005 Michael Smith <msmith@xiph.org>
 * Licensed under the GNU GPL, distributed with this program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <math.h>

#include <ogg/ogg.h>

#include "../include/daala/codec.h"
#include "../include/daala/daaladec.h"

#define CHUNK 4500

#ifdef _WIN32
#define I64FORMAT "I64d"
#else
#define I64FORMAT PRId64
#endif

/* TODO:
 *
 * - detect violations of muxing constraints
 * - detect granulepos 'gaps' (possibly vorbis-specific). (seperate from
 *   serial-number gaps)
 */

typedef struct _stream_processor {
  void (*process_page)(struct _stream_processor *, ogg_page *);
  void (*process_end)(struct _stream_processor *);
  int isillegal;
  int constraint_violated;
  int shownillegal;
  int isnew;
  long seqno;
  int lostseq;
  int start;
  int end;
  int num;
  char *type;
  ogg_uint32_t serial; /* must be 32 bit unsigned */
  ogg_stream_state os;
  void *data;
} stream_processor;

typedef struct {
  stream_processor *streams;
  int allocated;
  int used;
  int in_headers;
} stream_set;

typedef struct {
  daala_info di;
  daala_comment dc;
  daala_setup_info *ds;
  ogg_int64_t bytes;
  ogg_int64_t lastgranulepos;
  ogg_int64_t firstgranulepos;
  int doneheaders;
  ogg_int64_t framenum_expected;
} misc_daala_info;

#define CONSTRAINT_PAGE_AFTER_EOS   1
#define CONSTRAINT_MUXING_VIOLATED  2

static stream_set *create_stream_set(void) {
  stream_set *set = calloc(1, sizeof(stream_set));
  set->streams = calloc(5, sizeof(stream_processor));
  set->allocated = 5;
  set->used = 0;
  return set;
}

static void check_xiph_comment(stream_processor *stream, int i,
 const char *comment, int comment_length) {
  char *sep;
  int j;
  int broken;
  unsigned char *val;
  int bytes;
  int remaining;
  broken = 0;
  sep = strchr(comment, '=');
  if (sep == NULL) {
    fprintf(stderr, "WARNING: Comment %d in stream %d has invalid format, "
     "does not contain '=': \"%s\"\n", i, stream->num, comment);
    return;
  }
  for (j = 0; j < sep - comment; j++) {
    if (comment[j] < 0x20 || comment[j] > 0x7D) {
      fprintf(stderr, "WARNING: Invalid comment fieldname in "
       "comment %d (stream %d): \"%s\"\n", i, stream->num, comment);
      broken = 1;
      break;
    }
  }
  if (broken) {
    return;
  }
  val = (unsigned char *)comment;
  j = sep - comment + 1;
  while (j < comment_length) {
    remaining = comment_length - j;
    if ((val[j] & 0x80) == 0) {
      bytes = 1;
    }
    else if ((val[j] & 0x40) == 0x40) {
      if ((val[j] & 0x20) == 0) {
        bytes = 2;
      }
      else if ((val[j] & 0x10) == 0) {
        bytes = 3;
      }
      else if((val[j] & 0x08) == 0) {
        bytes = 4;
      }
      else if((val[j] & 0x04) == 0) {
        bytes = 5;
      }
      else if((val[j] & 0x02) == 0) {
        bytes = 6;
      }
      else {
        fprintf(stderr, "WARNING: Illegal UTF-8 sequence in "
         "comment %d (stream %d): length marker wrong\n",
         i, stream->num);
        broken = 1;
        break;
      }
    }
    else {
      fprintf(stderr, "WARNING: Illegal UTF-8 sequence in comment "
       "%d (stream %d): length marker wrong\n", i, stream->num);
      broken = 1;
      break;
    }
    if(bytes > remaining) {
      fprintf(stderr, "WARNING: Illegal UTF-8 sequence in comment "
       "%d (stream %d): too few bytes\n", i, stream->num);
      broken = 1;
      break;
    }
    switch (bytes) {
      case 1:
        /* No more checks needed */
        break;
      case 2:
        if ((val[j+1] & 0xC0) != 0x80) {
          broken = 1;
        }
        if ((val[j] & 0xFE) == 0xC0) {
          broken = 1;
        }
        break;
      case 3:
        if(!((val[j] == 0xE0 && val[j+1] >= 0xA0 && val[j+1] <= 0xBF &&
                (val[j+2] & 0xC0) == 0x80) ||
              (val[j] >= 0xE1 && val[j] <= 0xEC &&
               (val[j+1] & 0xC0) == 0x80 &&
               (val[j+2] & 0xC0) == 0x80) ||
              (val[j] == 0xED && val[j+1] >= 0x80 &&
               val[j+1] <= 0x9F &&
               (val[j+2] & 0xC0) == 0x80) ||
              (val[j] >= 0xEE && val[j] <= 0xEF &&
               (val[j+1] & 0xC0) == 0x80 &&
               (val[j+2] & 0xC0) == 0x80))) {
          broken = 1;
        }
        if (val[j] == 0xE0 && (val[j+1] & 0xE0) == 0x80) {
          broken = 1;
        }
        break;
      case 4:
        if (!((val[j] == 0xF0 && val[j+1] >= 0x90 &&
                val[j+1] <= 0xBF &&
                (val[j+2] & 0xC0) == 0x80 &&
                (val[j+3] & 0xC0) == 0x80) ||
              (val[j] >= 0xF1 && val[j] <= 0xF3 &&
               (val[j+1] & 0xC0) == 0x80 &&
               (val[j+2] & 0xC0) == 0x80 &&
               (val[j+3] & 0xC0) == 0x80) ||
              (val[j] == 0xF4 && val[j+1] >= 0x80 &&
               val[j+1] <= 0x8F &&
               (val[j+2] & 0xC0) == 0x80 &&
               (val[j+3] & 0xC0) == 0x80))) {
          broken = 1;
        }
        if (val[j] == 0xF0 && (val[j+1] & 0xF0) == 0x80) {
          broken = 1;
        }
        break;
        /* 5 and 6 aren't actually allowed at this point */
      case 5:
        broken = 1;
        break;
      case 6:
        broken = 1;
        break;
    }
    if (broken) {
      char *simple = malloc(comment_length + 1);
      char *seq = malloc(comment_length * 3 + 1);
      static char hex[] = {'0', '1', '2', '3', '4', '5', '6', '7',
        '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
      int c, c1 = 0, c2 = 0;
      for (c = 0; c < comment_length; c++) {
        seq[c1++] = hex[((unsigned char)comment[i]) >> 4];
        seq[c1++] = hex[((unsigned char)comment[i]) & 0xf];
        seq[c1++] = ' ';
        if (comment[i] < 0x20 || comment[i] > 0x7D) {
          simple[c2++] = '?';
        }
        else {
          simple[c2++] = comment[i];
        }
      }
      seq[c1] = 0;
      simple[c2] = 0;
      fprintf(stderr, "WARNING: Illegal UTF-8 sequence in comment "
       "%d (stream %d): invalid sequence \"%s\": %s\n", i,
       stream->num, simple, seq);
      broken = 1;
      free(simple);
      free(seq);
      break;
    }
    j += bytes;
  }
}

static void daala_process(stream_processor *stream, ogg_page *page) {
  ogg_packet packet;
  misc_daala_info *inf;
  int i, header=0;
  int res;
  inf = stream->data;
  ogg_stream_pagein(&stream->os, page);
  if (inf->doneheaders < 3) {
    header = 1;
  }
  while(1) {
    res = ogg_stream_packetout(&stream->os, &packet);
    if (res < 0) {
      fprintf(stderr, "WARNING: discontinuity in stream (%d)\n", stream->num);
      continue;
    }
    else if (res == 0) {
      break;
    }
    if (inf->doneheaders < 3) {
      if (daala_decode_header_in(&inf->di, &inf->dc, &inf->ds, &packet) < 0) {
        fprintf(stderr, "WARNING: Could not decode Daala header packet - "
         "invalid Daala stream (%d)\n", stream->num);
        continue;
      }
      inf->doneheaders++;
      if (inf->doneheaders == 3) {
        if (ogg_page_granulepos(page) != 0
         || ogg_stream_packetpeek(&stream->os, NULL) == 1) {
          fprintf(stderr, "WARNING: Daala stream %d does not have headers "
           "correctly framed. Terminal header page contains "
           "additional packets or has non-zero granulepos\n",
           stream->num);
        }
        printf("Daala headers parsed for stream %d, "
         "information follows...\n", stream->num);
        printf("Version: %d.%d.%d\n", inf->di.version_major,
         inf->di.version_minor, inf->di.version_sub);
        printf("Vendor: %s\n", inf->dc.vendor);
        printf("Width: %d\n", inf->di.pic_width);
        printf("Height: %d\n", inf->di.pic_height);
        if (inf->di.timebase_numerator == 0
         || inf->di.timebase_denominator == 0) {
          fprintf(stderr, "Invalid zero framerate\n");
        }
        else {
          printf("Framerate %d/%d (%.02f fps)\n",
           inf->di.timebase_numerator, inf->di.timebase_denominator,
           (float)inf->di.timebase_numerator/
           (float)inf->di.timebase_denominator);
        }
        if(inf->di.pixel_aspect_numerator == 0
         || inf->di.pixel_aspect_denominator == 0) {
          printf("Aspect ratio undefined\n");
        }
        else {
          float frameaspect;
          frameaspect = (float)inf->di.pic_width/(float)inf->di.pic_height*
           (float)inf->di.pixel_aspect_numerator/
           (float)inf->di.pixel_aspect_denominator;
          printf("Pixel aspect ratio %d:%d (%f:1)\n",
           inf->di.pixel_aspect_numerator, inf->di.pixel_aspect_denominator,
           (float)inf->di.pixel_aspect_numerator/
           (float)inf->di.pixel_aspect_denominator);
          if (fabs(frameaspect - 4.0/3.0) < 0.02) {
            printf("Frame aspect 4:3\n");
          }
          else if(fabs(frameaspect - 16.0/9.0) < 0.02) {
            printf("Frame aspect 16:9\n");
          }
          else {
            printf("Frame aspect %f:1\n", frameaspect);
          }
        }
        if (inf->dc.comments > 0) {
          printf("User comments section follows...\n");
        }
        for (i = 0; i < inf->dc.comments; i++) {
          char *comment;
          comment = inf->dc.user_comments[i];
          check_xiph_comment(stream, i, comment, inf->dc.comment_lengths[i]);
        }
      }
    }
    else {
      ogg_int64_t framenum;
      ogg_int64_t iframe, pframe;
      ogg_int64_t gp;
      gp = packet.granulepos;
      if (gp > 0) {
        iframe = gp >> inf->di.keyframe_granule_shift;
        pframe = gp - (iframe << inf->di.keyframe_granule_shift);
        framenum = iframe + pframe;
        if (inf->framenum_expected >= 0
         && inf->framenum_expected != framenum) {
          fprintf(stderr, "WARNING: Expected frame %" I64FORMAT
           ", got %" I64FORMAT "\n", inf->framenum_expected,
           framenum);
        }
        inf->framenum_expected = framenum + 1;
      }
      else if (inf->framenum_expected >= 0) {
        inf->framenum_expected++;
      }
    }
  }
  if (!header) {
    ogg_int64_t gp;
    gp = ogg_page_granulepos(page);
    if (gp > 0) {
      if (gp < inf->lastgranulepos) {
        fprintf(stderr, "WARNING: granulepos in stream %d decreases from %"
         I64FORMAT " to %" I64FORMAT "\n", stream->num, inf->lastgranulepos,
         gp);
      }
      inf->lastgranulepos = gp;
    }
    if (inf->firstgranulepos < 0) { /* Not set yet */
    }
    inf->bytes += page->header_len + page->body_len;
  }
}

static void daala_end(stream_processor *stream) {
  misc_daala_info *inf;
  long minutes, seconds, milliseconds;
  double bitrate, time;
  ogg_int64_t iframe, pframe;
  int new_gp;
  inf = stream->data;
  new_gp = inf->di.version_major > 3
   || (inf->di.version_major == 3 && (inf->di.version_minor > 2
   || (inf->di.version_minor == 2 && inf->di.version_sub > 0)));
  /* This should be lastgranulepos - startgranulepos, or something like that*/
  iframe = inf->lastgranulepos >> inf->di.keyframe_granule_shift;
  pframe = inf->lastgranulepos - (iframe << inf->di.keyframe_granule_shift);
  /* The granule position starts at 0 for stream version 3.2.0, but starts at
     1 for version 3.2.1 and above. In the former case, we need to add one
     to the final granule position to get the frame count. */
  time = (double)(iframe + pframe + !new_gp)/
   ((float)inf->di.timebase_numerator/(float)inf->di.timebase_denominator);
  minutes = (long)time/60;
  seconds = (long)time - minutes*60;
  milliseconds = (long)((time - minutes*60 - seconds)*1000);
  bitrate = inf->bytes*8/time/1000.0;
  printf("Daala stream %d:\n"
   "\tTotal data length: %" I64FORMAT " bytes\n"
   "\tPlayback length: %ldm:%02ld.%03lds\n"
   "\tAverage bitrate: %f kb/s\n",
   stream->num, inf->bytes, minutes, seconds, milliseconds, bitrate);
  daala_comment_clear(&inf->dc);
  daala_info_clear(&inf->di);
  daala_setup_free(inf->ds);
  free(stream->data);
}

static void process_null(stream_processor *stream, ogg_page *page) {
  /* This is for invalid streams. */
  (void)stream;
  (void)page;
}

static void process_other(stream_processor *stream, ogg_page *page ) {
  ogg_packet packet;
  ogg_stream_pagein(&stream->os, page);
  while (ogg_stream_packetout(&stream->os, &packet) > 0) {
    /* Should we do anything here? Currently, we don't */
  }
}

static void free_stream_set(stream_set *set) {
  int i;
  for (i = 0; i < set->used; i++) {
    if (!set->streams[i].end) {
      fprintf(stderr, "WARNING: EOS not set on stream %d\n",
       set->streams[i].num);
      if (set->streams[i].process_end) {
        set->streams[i].process_end(&set->streams[i]);
      }
    }
    ogg_stream_clear(&set->streams[i].os);
  }
  free(set->streams);
  free(set);
}

static int streams_open(stream_set *set) {
  int i;
  int res;
  res = 0;
  for (i = 0; i < set->used; i++) {
    if (!set->streams[i].end) {
      res++;
    }
  }
  return res;
}

static void null_start(stream_processor *stream) {
  stream->process_end = NULL;
  stream->type = "invalid";
  stream->process_page = process_null;
}

static void other_start(stream_processor *stream, char *type) {
  if (type) {
    stream->type = type;
  }
  else {
    stream->type = "unknown";
  }
  stream->process_page = process_other;
  stream->process_end = NULL;
}

static void daala_start(stream_processor *stream) {
  misc_daala_info *info;
  stream->type = "daala";
  stream->process_page = daala_process;
  stream->process_end = daala_end;
  stream->data = calloc(1, sizeof(misc_daala_info));
  info = stream->data;
  info->framenum_expected = -1;
}

static stream_processor *find_stream_processor(stream_set *set,
 ogg_page *page) {
  ogg_uint32_t serial;
  int i, invalid, constraint;
  stream_processor *stream;
  serial = ogg_page_serialno(page);
  invalid = constraint = 0;
  for (i = 0; i < set->used; i++) {
    if (serial == set->streams[i].serial) {
      /* We have a match! */
      stream = &(set->streams[i]);
      set->in_headers = 0;
      /* if we have detected EOS, then this can't occur here. */
      if (stream->end) {
        stream->isillegal = 1;
        stream->constraint_violated = CONSTRAINT_PAGE_AFTER_EOS;
        return stream;
      }
      stream->isnew = 0;
      stream->start = ogg_page_bos(page);
      stream->end = ogg_page_eos(page);
      stream->serial = serial;
      return stream;
    }
  }
  /* If there are streams open, and we've reached the end of the
   * headers, then we can't be starting a new stream.
   * XXX: might this sometimes catch ok streams if EOS flag is missing,
   * but the stream is otherwise ok?
   */
  if (streams_open(set) && !set->in_headers) {
    constraint = CONSTRAINT_MUXING_VIOLATED;
    invalid = 1;
  }
  set->in_headers = 1;
  if (set->allocated < set->used) {
    stream = &set->streams[set->used];
  }
  else {
    set->allocated += 5;
    set->streams = realloc(set->streams, sizeof(stream_processor)*
     set->allocated);
    stream = &set->streams[set->used];
  }
  set->used++;
  stream->num = set->used; /* We count from 1 */
  stream->isnew = 1;
  stream->isillegal = invalid;
  stream->constraint_violated = constraint;
  {
    int res;
    ogg_packet packet;
    /* We end up processing the header page twice, but that's ok. */
    ogg_stream_init(&stream->os, serial);
    ogg_stream_pagein(&stream->os, page);
    res = ogg_stream_packetout(&stream->os, &packet);
    if (res <= 0) {
      fprintf(stderr, "WARNING: Invalid header page, no packet found\n");
      null_start(stream);
    }
    else if (packet.bytes >= 6
     && memcmp(packet.packet, "\x80""daala", 6) == 0) {
      daala_start(stream);
    }
    else {
      other_start(stream, NULL);
    }
    res = ogg_stream_packetout(&stream->os, &packet);
    if (res > 0) {
      fprintf(stderr, "WARNING: Invalid header page in stream %d, "
       "contains multiple packets\n", stream->num);
    }
    /* re-init, ready for processing */
    ogg_stream_clear(&stream->os);
    ogg_stream_init(&stream->os, serial);
  }
  stream->start = ogg_page_bos(page);
  stream->end = ogg_page_eos(page);
  stream->serial = serial;
  if (stream->serial == 0 || stream->serial == UINT32_MAX) {
    printf("Note: Stream %d has serial number %d, which is legal but may "
        "cause problems with some tools.\n", stream->num, stream->serial);
  }
  return stream;
}

static int get_next_page(FILE *f, ogg_sync_state *sync, ogg_page *page,
 ogg_int64_t *written) {
  int ret;
  char *buffer;
  int bytes;
  while ((ret = ogg_sync_pageseek(sync, page)) <= 0) {
    if (ret < 0) {
      /* unsynced, we jump over bytes to a possible capture - we don't need to
         read more just yet */
      fprintf(stderr, "WARNING: Hole in data (%d bytes) found at approximate "
       "offset %" I64FORMAT " bytes. Corrupted Ogg.\n", -ret, *written);
      continue;
    }
    /* zero return, we didn't have enough data to find a whole page, read */
    buffer = ogg_sync_buffer(sync, CHUNK);
    bytes = fread(buffer, 1, CHUNK, f);
    if (bytes <= 0) {
      ogg_sync_wrote(sync, 0);
      return 0;
    }
    ogg_sync_wrote(sync, bytes);
    *written += bytes;
  }
  return 1;
}

static void process_file(char *filename) {
  FILE *file;
  ogg_sync_state sync;
  ogg_page page;
  stream_set *processors;
  int gotpage;
  ogg_int64_t written;
  file = fopen(filename, "rb");
  processors = create_stream_set();
  written = 0;
  gotpage = 0;
  if (!file) {
    fprintf(stderr, "Error opening input file \"%s\": %s\n", filename,
     strerror(errno));
    return;
  }
  printf("Processing file \"%s\"...\n\n", filename);
  ogg_sync_init(&sync);
  while (get_next_page(file, &sync, &page, &written)) {
    stream_processor *p;
    p = find_stream_processor(processors, &page);
    gotpage = 1;
    if (!p) {
      fprintf(stderr, "Could not find a processor for stream, bailing\n");
      return;
    }
    if (p->isillegal && !p->shownillegal) {
      char *constraint;
      switch (p->constraint_violated) {
        case CONSTRAINT_PAGE_AFTER_EOS:
          constraint = "Page found for stream after EOS flag";
          break;
        case CONSTRAINT_MUXING_VIOLATED:
          constraint = "Ogg muxing constraints violated, new "
            "stream before EOS of all previous streams";
          break;
        default:
          constraint = "Error unknown.";
      }
      fprintf(stderr, "WARNING: illegally placed page(s) for logical stream "
          "%d\nThis indicates a corrupt Ogg file: %s.\n", p->num, constraint);
      p->shownillegal = 1;
      /* If it's a new stream, we want to continue processing this page
       * anyway to suppress additional spurious errors
       */
      if (!p->isnew)
        continue;
    }
    if (p->isnew) {
      printf("New logical stream (#%d, serial: %08x): type %s\n",
       p->num, p->serial, p->type);
      if (!p->start) {
        fprintf(stderr, "WARNING: stream start flag not set on stream %d\n",
         p->num);
      }
    }
    else if(p->start) {
      fprintf(stderr, "WARNING: stream start flag found in mid-stream "
       "on stream %d\n", p->num);
    }
    if (p->seqno++ != ogg_page_pageno(&page)) {
      if (!p->lostseq) {
        fprintf(stderr, "WARNING: sequence number gap in stream %d. Got page "
         "%ld when expecting page %ld. Indicates missing data.\n",
         p->num, ogg_page_pageno(&page), p->seqno - 1);
      }
      p->seqno = ogg_page_pageno(&page);
      p->lostseq = 1;
    }
    else {
      p->lostseq = 0;
    }
    if (!p->isillegal) {
      p->process_page(p, &page);
      if (p->end) {
        if (p->process_end) {
          p->process_end(p);
        }
        printf("Logical stream %d ended\n", p->num);
        p->isillegal = 1;
        p->constraint_violated = CONSTRAINT_PAGE_AFTER_EOS;
      }
    }
  }
  if (!gotpage) {
    fprintf(stderr, "ERROR: No Ogg data found in file \"%s\".\n"
     "Input probably not Ogg.\n", filename);
  }
  free_stream_set(processors);
  ogg_sync_clear(&sync);
  fclose(file);
}

static void version(void) {
  printf("daalainfo from %s\n", daala_version_string());
}

static void usage(void) {
  version();
  printf(" by the Xiph.Org Foundation (http://www.xiph.org/)\n\n");
  printf("Usage: daalainfo [flags] file1.ogg [file2.ogx ... fileN.ogv]\n"
   "Flags supported:\n"
   "\t-h Show this help message\n");
  printf("\t-V Output version information and exit\n");
}

int main(int argc, char **argv) {
  int f, ret;
  if (argc < 2) {
    printf("Usage: daalainfo [flags] file1.ogg [file2.ogx ... fileN.ogv]\n"
        "\n"
        "daalainfo is a tool for printing information about Ogg Daala files\n"
        "and for diagnosing problems with them.\n"
        "Full help shown with \"daalainfo -h\".\n");
    return 1;
  }
  while ((ret = getopt(argc, argv, "hqvV")) >= 0) {
    switch (ret) {
      case 'h':
        usage();
        return 0;
      case 'V':
        version();
        return 0;
    }
  }
  if (optind >= argc) {
    fprintf(stderr, "No input files specified. \"daalainfo -h\" for help\n");
    return 1;
  }
  ret = 0;
  for (f = optind; f < argc; f++) {
    process_file(argv[f]);
  }
  return ret;
}
