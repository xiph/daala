#!/bin/bash -e
# continuous integration test script
# run this from the top-level source directory

./autogen.sh
CC=clang CFLAGS='-Wall -O1 -g -fsanitize=address -fno-sanitize-recover' LDFLAGS='-fsanitize=address -fno-sanitize-recover' ./configure
make clean all tools
./examples/encoder_example -k 4 ${VIDEOS}/claire_qcif-2frames.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o /dev/null
rm -f out.$$.ogv
./examples/encoder_example -k 4 ${VIDEOS}/tos444.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o /dev/null
rm -f out.$$.ogv
./examples/encoder_example -k 4 ${VIDEOS}/nes420.y4m -o out.$$.ogv
./examples/dump_video out.$$.ogv -o out.$$.y4m
./tools/dump_psnr ${VIDEOS}/nes420.y4m out.$$.y4m
./tools/dump_psnrhvs ${VIDEOS}/nes420.y4m out.$$.y4m
./tools/dump_ssim ${VIDEOS}/nes420.y4m out.$$.y4m
./tools/dump_fastssim -c ${VIDEOS}/nes420.y4m out.$$.y4m
rm -f out.$$.ogv out.$$.y4m
