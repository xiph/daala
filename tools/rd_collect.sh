#!/bin/bash
set -e

CODECS="<daala|libaom|libaom-rt|vp8|vp9|x264|x265|libjpeg|mozjpeg|theora|webp|bpg|rav1e|svt-av1>"

if [ $# == 0 ]; then
  echo "usage: DAALA_ROOT=<build_dir> $0 $CODECS *.y4m"
  exit 1
fi

export CODEC=$1
shift

if [ $# == 0 ] || [[ ! $CODECS =~ [\<\|]$CODEC[\|\>] ]]; then
  echo "usage: DAALA_ROOT=<build_dir> $0 $CODECS *.y4m"
  exit 1
fi

if [ -z $DAALA_ROOT ]; then
  DAALA_ROOT=.
fi

if [ -z $TOOLS_ROOT ]; then
  TOOLS_ROOT=$DAALA_ROOT
fi

export EXTRA_OPTS=$EXTRA_OPTS

if [ ! -d $DAALA_ROOT ]; then
  echo "Please set DAALA_ROOT to the location of your daala git clone"
  exit 1
fi

case $CODEC in
  daala)
    if [ ! -f $DAALA_ROOT/config.h ]; then
      echo "File not found $DAALA_ROOT/config.h"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    if [ -z "$ENCODER_EXAMPLE" ]; then
      export ENCODER_EXAMPLE=$DAALA_ROOT/examples/encoder_example
    fi

    if [ ! -x "$ENCODER_EXAMPLE" ]; then
      echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    if [ -z "$DUMP_VIDEO" ]; then
      export DUMP_VIDEO=$DAALA_ROOT/examples/dump_video
    fi

    if [ ! -x "$DUMP_VIDEO" ]; then
      echo "Executable not found DUMP_VIDEO=$DUMP_VIDEO"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS="5 7 11 16 25 37 55 81 122 181 270 400"
    fi

    export RD_COLLECT_SUB=$(dirname "$0")/rd_collect_daala.sh
    ;;
  libaom | libaom-rt)
    if [ -z $AOM_ROOT ] || [ ! -d $AOM_ROOT ]; then
      echo "Please set AOM_ROOT to the location of your aom git clone"
      exit 1
    fi

    if [ -z "$AOMENC" ]; then
      export AOMENC=$AOM_ROOT/aomenc
    fi

    if [ -z "$AOMDEC" ]; then
      export AOMDEC=$AOM_ROOT/aomdec
    fi

    if [ ! -x "$AOMENC" ]; then
      echo "Executable not found AOMENC=$AOMENC"
      echo "Do you have the right AOM_ROOT=$AOM_ROOT"
      exit 1
    fi

    if [ ! -x "$AOMDEC" ]; then
      echo "Executable not found AOMDEC=$AOMDEC"
      echo "Do you have the right AOM_ROOT=$AOM_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS="20 32 43 55 63"
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_libaom.sh
    ;;
  vp8)
    if [ -z $LIBVPX_ROOT ] || [ ! -d $LIBVPX_ROOT ]; then
      echo "Please set LIBVPX_ROOT to the location of your libvpx git clone"
      exit 1
    fi

    if [ -z "$VPXENC" ]; then
      export VPXENC=$LIBVPX_ROOT/vpxenc
    fi

    if [ -z "$VPXDEC" ]; then
      export VPXDEC=$LIBVPX_ROOT/vpxdec
    fi

    if [ ! -x "$VPXENC" ]; then
      echo "Executable not found VPXENC=$VPXENC"
      echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
      exit 1
    fi

    if [ ! -x "$VPXDEC" ]; then
      echo "Executable not found VPXDEC=$VPXDEC"
      echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS=$(seq 12 5 64)
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_libvpx.sh
    ;;
  vp9)
    if [ -z $LIBVPX_ROOT ] || [ ! -d $LIBVPX_ROOT ]; then
      echo "Please set LIBVPX_ROOT to the location of your libvpx git clone"
      exit 1
    fi

    if [ -z "$VPXENC" ]; then
      export VPXENC=$LIBVPX_ROOT/vpxenc
    fi

    if [ -z "$VPXDEC" ]; then
      export VPXDEC=$LIBVPX_ROOT/vpxdec
    fi

    if [ ! -x "$VPXENC" ]; then
      echo "Executable not found VPXENC=$VPXENC"
      echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
      exit 1
    fi

    if [ ! -x "$VPXDEC" ]; then
      echo "Executable not found VPXDEC=$VPXDEC"
      echo "Do you have the right LIBVPX_ROOT=$LIBVPX_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS="20 32 43 55 63"
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_libvpx.sh
    ;;
  x264)
    if [ -z $X264_ROOT ] || [ ! -d $X264_ROOT ]; then
      echo "Please set X264_ROOT to the location of your x264 git clone"
      exit 1
    fi

    if [ -z "$X264" ]; then
      export X264=$X264_ROOT/x264
    fi

    if [ ! -x "$X264" ]; then
      echo "Executable not found X264=$X264"
      echo "Do you have the right X264_ROOT=$X264_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS=$(seq 1 51)
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_x264.sh
    ;;
  x265)
    if [ -z $X265_ROOT ] || [ ! -d $X265_ROOT ]; then
      echo "Please set X265_ROOT to the location of your x265 hg checkout"
      exit 1
    fi

    if [ -z "$X265" ]; then
      export X265=$X265_ROOT/build/linux/x265
    fi

    if [ ! -x "$X265" ]; then
      echo "Executable not found X265=$X265"
      echo "Do you have the right X265_ROOT=$X265_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_x265.sh
    ;;
  libjpeg)
    if [ -z "$YUVJPEG" ]; then
      export YUVJPEG=$DAALA_ROOT/tools/yuvjpeg
    fi

    if [ -z "$JPEGYUV" ]; then
      export JPEGYUV=$DAALA_ROOT/tools/jpegyuv
    fi

    if [ ! -x "$YUVJPEG" ]; then
      echo "Executable not found YUVJPEG=$YUVJPEG"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    if [ ! -x "$JPEGYUV" ]; then
      echo "Executable not found JPEGYUV=$JPEGYUV"
      echo "Do you have the right DAALA_ROOT=$DAALA_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_jpeg.sh
    ;;
  mozjpeg)
    if [ -z $MOZJPEG_ROOT ] || [ ! -d $MOZJPEG_ROOT ]; then
      echo "Please set MOZJPEG_ROOT to the location of your mozjpeg git clone"
      exit 1
    fi

    if [ -z "$YUVJPEG" ]; then
      export YUVJPEG=$MOZJPEG_ROOT/yuvjpeg
    fi

    if [ -z "$JPEGYUV" ]; then
      export JPEGYUV=$MOZJPEG_ROOT/jpegyuv
    fi

    if [ ! -x "$YUVJPEG" ]; then
      echo "Executable not found YUVJPEG=$YUVJPEG"
      echo "Do you have the right MOZJPEG_ROOT=$MOZJPEG_ROOT"
      exit 1
    fi

    if [ ! -x "$JPEGYUV" ]; then
      echo "Executable not found JPEGYUV=$JPEGYUV"
      echo "Do you have the right MOZJPEG_ROOT=$MOZJPEG_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_jpeg.sh
    ;;
  theora)
    if [ -z $THEORA_ROOT ] || [ ! -d $THEORA_ROOT ]; then
      echo "Please set THEORA_ROOT to the location of your theora svn checkout"
      exit 1
    fi

    if [ -z "$ENCODER_EXAMPLE" ]; then
      export ENCODER_EXAMPLE=$THEORA_ROOT/examples/encoder_example
    fi

    if [ -z "$DUMP_VIDEO" ]; then
      export DUMP_VIDEO=$THEORA_ROOT/examples/dump_video
    fi

    if [ ! -x "$ENCODER_EXAMPLE" ]; then
      echo "Executable not found ENCODER_EXAMPLE=$ENCODER_EXAMPLE"
      echo "Do you have the right THEORA_ROOT=$THEORA_ROOT"
      exit 1
    fi

    if [ ! -x "$DUMP_VIDEO" ]; then
      echo "Executable not found DUMP_VIDEO=$DUMP_VIDEO"
      echo "Do you have the right THEORA_ROOT=$THEORA_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_theora.sh
    ;;
  webp)
    if [ -z $WEBP_ROOT ] || [ ! -d $WEBP_ROOT ]; then
      echo "Please set WEBP_ROOT to the location of your webp checkout"
      exit 1
    fi

    if [ -z "$CWEBP" ]; then
      export CWEBP=$WEBP_ROOT/examples/cwebp
    fi

    if [ -z "$DWEBP" ]; then
      export DWEBP=$WEBP_ROOT/examples/dwebp
    fi

    if [ ! -x "$CWEBP" ]; then
      echo "Executable not found CWEBP=$CWEBP"
      echo "Do you have the right WEBP_ROOT=$WEBP_ROOT"
      exit 1
    fi

    if [ ! -x "$DWEBP" ]; then
      echo "Executable not found DWEBP=$DWEBP"
      echo "Do you have the right WEBP_ROOT=$WEBP_ROOT"
      exit 1
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_webp.sh
    ;;
  bpg)
    if [ -z $BPG_ROOT ] || [ ! -d $BPG_ROOT ]; then
      echo "Please set BPG_ROOT to the location of your libbpg checkout"
      exit 1
    fi

    if [ -z "$BPGENC" ]; then
      export BPGENC=$BPG_ROOT/bpgenc
    fi

    if [ -z "$BPGDEC" ]; then
      export BPGDEC=$BPG_ROOT/bpgdec
    fi

    if [ ! -x "$BPGENC" ]; then
      echo "Executable not found BPGENC=$BPGENC"
      echo "Do you have the right BPG_ROOT=$BPG_ROOT"
      exit 1
    fi

    if [ ! -x "$BPGDEC" ]; then
      echo "Executable not found BPGDEC=$BPGDEC"
      echo "Do you have the right BPG_ROOT=$BPG_ROOT"
      exit 1
    fi

    export CORES=1
    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_bpg.sh
    ;;
  rav1e)
    if [ -z $RAV1E_ROOT ] || [ ! -d $RAV1E_ROOT ]; then
      echo "Please set RAV1E_ROOT to the location of your rav1e git clone"
      exit 1
    fi

    if [ -z "$RAV1E" ]; then
      export RAV1E="$RAV1E_ROOT/target/release/rav1e"
    fi

    if [ ! -x "$RAV1E" ]; then
      echo "Executable not found RAV1E=$RAV1E"
      echo "Do you have the right RAV1E_ROOT=$RAV1E_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS="80 128 172 220 252"
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_rav1e.sh
    ;;
  svt-av1)
    if [ -z $SVTAV1_ROOT ] || [ ! -d $SVTAV1_ROOT ]; then
      echo "Please set SVTAV1_ROOT to the location of your SVT-AV1 git clone"
      exit 1
    fi

    if [ -z "$SVTAV1" ]; then
      export SVTAV1="$SVTAV1_ROOT/Bin/Release/SvtAv1EncApp"
    fi

    if [ ! -x "$SVTAV1" ]; then
      echo "Executable not found SVTAV1=$SVTAV1"
      echo "Do you have the right SVTAV1_ROOT=$SVTAV1_ROOT"
      exit 1
    fi

    if [ -z "$QPS" ]; then
      QPS="20 32 43 55 63"
    fi

    export RD_COLLECT_SUB=$(dirname $0)/rd_collect_svtav1.sh
    ;;
  *)
    echo "Unknown codec: $CODEC"
    exit 1
esac

if [ -z "$PLANE" ]; then
  export PLANE=0
fi

if [ $PLANE != 0 ] && [ $PLANE != 1 ] && [ $PLANE != 2 ] &&
  [ $PLANE != -1 ]; then
  echo "Invalid plane $PLANE. Must be 0, 1, 2, or -1 (all planes)."
  exit 1
fi

# TODO refactor these out of the daala project into a metrics project

if [ -z "$YUV2YUV4MPEG" ]; then
  export YUV2YUV4MPEG=$TOOLS_ROOT/tools/yuv2yuv4mpeg
fi

if [ -z "$Y4M2PNG" ]; then
  export Y4M2PNG=$TOOLS_ROOT/tools/y4m2png
fi

if [ -z "$Y4M2YUV" ]; then
  export Y4M2YUV=$TOOLS_ROOT/tools/y4m2yuv
fi

if [ -z "$PNG2Y4M" ]; then
  export PNG2Y4M=$TOOLS_ROOT/tools/png2y4m
fi

if [ -z "$DUMP_PSNR" ]; then
  export DUMP_PSNR=$TOOLS_ROOT/tools/dump_psnr
fi

if [ -z "$DUMP_PSNRHVS" ]; then
  export DUMP_PSNRHVS=$TOOLS_ROOT/tools/dump_psnrhvs
fi

if [ -z "$DUMP_SSIM" ]; then
  export DUMP_SSIM=$TOOLS_ROOT/tools/dump_ssim
fi

if [ -z "$DUMP_FASTSSIM" ]; then
  export DUMP_FASTSSIM=$TOOLS_ROOT/tools/dump_fastssim
fi

if [ -z "$DUMP_CIEDE" ]; then
  export DUMP_CIEDE=$(dirname $0)/dump_ciede2000.py
fi

if [ ! -x "$YUV2YUV4MPEG" ]; then
  echo "Executable not found YUV2YUV4MPEG=$YUV2YUV4MPEG"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$Y4M2PNG" ]; then
  echo "Executable not found Y4M2PNG=$Y4M2PNG"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$Y4M2YUV" ]; then
  echo "Executable not found Y4M2YUV=$Y4M2YUV"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$PNG2Y4M" ]; then
  echo "Executable not found PNG2Y4M=$PNG2Y4M"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_PSNR" ]; then
  echo "Executable not found DUMP_PSNR=$DUMP_PSNR"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_PSNRHVS" ]; then
  echo "Executable not found DUMP_PSNRHVS=$DUMP_PSNRHVS"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_SSIM" ]; then
  echo "Executable not found DUMP_FASTSSIM=$DUMP_SSIM"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

if [ ! -x "$DUMP_FASTSSIM" ]; then
  echo "Executable not found DUMP_FASTSSIM=$DUMP_FASTSSIM"
  echo "Do you have the right TOOLS_ROOT=$TOOLS_ROOT"
  exit 1
fi

set +e
temp=$($DUMP_CIEDE)
if [ $? -ne 0 ]; then
  echo "Warning Python dependencies not found. CIEDE2000 will not be computed."
  echo "Required Python dependencies are: numpy, skimage and y4m."
  if [ "$(uname -s)" = "Darwin" ]; then
    DUMP_CIEDE=/usr/bin/true
  else
    DUMP_CIEDE=/bin/true
  fi
fi
set -e

if [ -z "$CORES" ]; then
  if [ "$(uname -s)" = "Darwin" ]; then
    CORES=$(sysctl -n hw.ncpu)
  else
    CORES=$(grep -i processor /proc/cpuinfo | wc -l)
  fi
  #echo "CORES not set, using $CORES"
fi

case $CODEC in
  libaom | libaom-rt | daala | rav1e | svt-av1 | vp8 | vp9 | x264)
    FILES=$(find -L "$@" -type f -name "*.y4m")
    for f in $FILES; do for q in $QPS; do printf "%s\0" $f $q; done; done | xargs -0 -n2 -P$CORES $RD_COLLECT_SUB
    for f in $FILES; do cat $(basename $f)-$CODEC-*.out | sort -n > $(basename $f)-$CODEC.out && rm $(basename $f)-$CODEC-*.out; done
    ;;
  *)
    find -L "$@" -type f -name "*.y4m" -print0 | xargs -0 -n1 -P$CORES $RD_COLLECT_SUB
esac
