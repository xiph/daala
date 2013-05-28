
noinst_PROGRAMS = dump_video encoder_example player_example

AM_CFLAGS = $(OGG_CFLAGS)

dump_video_SOURCES = dump_video.c
dump_video_LDADD = $(DAALA_DEC_LIBS) $(OGG_LIBS)

encoder_example_SOURCES = encoder_example.c
encoder_example_LDADD = $(DAALA_ENC_LIBS) $(OGG_LIBS)

player_example_SOURCES = player_example.c
player_example_LDADD = $(DAALA_ENC_LIBS) $(OGG_LIBS)

