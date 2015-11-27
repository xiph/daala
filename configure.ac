# autoconf source script for generating configure

dnl The package_version file will be automatically synced to the git revision
dnl by the update_version script when configured in the repository, but will
dnl remain constant in tarball releases unless it is manually edited.
m4_define([CURRENT_VERSION],
          m4_esyscmd([ if test -e package_version || ./update_version; then
                           . ./package_version
                           printf "$PACKAGE_VERSION"
                       else
                           printf "unknown"
                       fi ]))

AC_INIT([daala],[CURRENT_VERSION],[daala@xiph.org])
AC_CONFIG_SRCDIR([src/decode.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AC_USE_SYSTEM_EXTENSIONS
AC_SYS_LARGEFILE

AC_PROG_CXX

AM_INIT_AUTOMAKE([1.11 foreign no-define subdir-objects])
AM_MAINTAINER_MODE([enable])
LT_INIT

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

dnl Library versioning for libtool.
dnl Please update these for releases.
dnl CURRENT, REVISION, AGE
dnl - library source changed -> increment REVISION
dnl - interfaces added/removed/changed -> increment CURRENT, REVISION = 0
dnl - interfaces added -> increment AGE
dnl - interfaces removed -> AGE = 0

OD_LT_CURRENT=0
OD_LT_REVISION=1
OD_LT_AGE=0

AC_SUBST(OD_LT_CURRENT)
AC_SUBST(OD_LT_REVISION)
AC_SUBST(OD_LT_AGE)

CFLAGS="-std=c89 -pedantic -Wall -Wextra -Wno-long-long \
        -Wshadow -Wno-overlength-strings $CFLAGS $OD_CFLAGS"

# Platform-specific tweaks
case $host in
  *-mingw*)
    # -std=c89 causes some warnings under mingw.
    CC_CHECK_CFLAGS_APPEND([-U__STRICT_ANSI__])
    # We need WINNT>=0x501 (WindowsXP) for getaddrinfo/freeaddrinfo.
    # It's okay to define this even when HTTP support is disabled, as it only
    #  affects header declarations, not linking (unless we actually use some
    #  XP-only functions).
    AC_DEFINE_UNQUOTED(_WIN32_WINNT,0x501,
     [We need at least WindowsXP for getaddrinfo/freaddrinfo])
    host_mingw=true
    ;;
  *-*-darwin*)
    os_darwin=true
    ;;
esac
dnl These two are entirely unused at present.
dnl AM_CONDITIONAL(OD_WIN32, test "$host_mingw" = "true")
dnl AM_CONDITIONAL(OD_OS_DARWIN, test "$os_darwin" = "true")

AC_ARG_ENABLE([assertions],
  AS_HELP_STRING([--enable-assertions], [Enable assertions in code]),,
  enable_assertions=no)
AM_CONDITIONAL([ENABLE_ASSERTIONS], [test "$enable_assertions" = "yes"])
AS_IF([test "$enable_assertions" = "yes"], [
  AC_DEFINE([OD_ENABLE_ASSERTIONS], [1], [Enable assertions in code])
])

AC_ARG_ENABLE([logging],
  AS_HELP_STRING([--enable-logging], [Enable logging]),,
  enable_logging=no)

AS_IF([test "$enable_logging" = "yes"], [
  AC_DEFINE([OD_LOGGING_ENABLED], [1], [Enable logging])
])
dnl Check for examples
AC_ARG_ENABLE([examples],
  AC_HELP_STRING([--disable-examples], [Disable the basic examples]),,
  [enable_examples=yes])

AS_IF([test "$enable_examples" = "yes"], [
    PKG_CHECK_MODULES([OGG], [ogg >= 1.3])
])

AM_CONDITIONAL(ENABLE_EXAMPLES, [test $enable_examples = "yes"])

AC_ARG_ENABLE([player],
  AS_HELP_STRING([--disable-player], [Disable the example player]),,
  enable_player=yes)

AS_IF([test "$enable_examples" = "no"],
      [enable_player = no])

dnl Set LIBM to the math library if needed.
LT_LIB_M

AC_CHECK_FUNCS([gettimeofday ftime],[break])

AS_IF([test "$enable_player" = "yes"], [
  PKG_CHECK_MODULES([SDL], [sdl2])
  PKG_CHECK_MODULES([OGG], [ogg >= 1.3])
])
AM_CONDITIONAL([ENABLE_PLAYER_EXAMPLE], [test "$enable_player" = "yes"])

#CC_ATTRIBUTE_VISIBILITY([default], [
#  CC_FLAG_VISIBILITY([CFLAGS="${CFLAGS} -fvisibility=hidden"])
#])

dnl Check for tools
AC_ARG_ENABLE([tools],
  AC_HELP_STRING([--disable-tools], [Disable the tools]),,
  [enable_tools=yes])

if test "$enable_tools" = "yes" ; then
   AC_CHECK_LIB([jpeg],[jpeg_start_compress],[JPEG_LIBS="-ljpeg"],[AC_MSG_ERROR([libjpeg needed for building tools])])
   AC_SUBST([JPEG_LIBS])
   PKG_CHECK_MODULES([PNG], [libpng])
   PKG_CHECK_MODULES([OGG], [ogg >= 1.3])

   AC_OPENMP
fi

AM_CONDITIONAL(ENABLE_TOOLS, [test $enable_tools = "yes"])

dnl Check for check
AC_ARG_ENABLE([unit-tests],
  AS_HELP_STRING([--disable-unit-tests], [Do not build unit tests]),,
  [enable_unit_tests=yes]
)

if test "$enable_unit_tests" = "yes" ; then
  PKG_CHECK_MODULES([CHECK], [check >= 0.9.8])

  dnl As of version 0.9.10, check does not include -pthread in its Cflags or
  dnl Libs even though it does depend on it. Assuming that check.pc may one day
  dnl be fixed for that, only add it here if it's not there already.  It really
  dnl belongs in LDFLAGS here (since check is a static lib and doesn't expose
  dnl pthread functions in its header), but it can't be added to CHECK_LIBS,
  dnl since automake 1.13 will barf about adding -pthread to *_LDADD rather
  dnl than *_LDFLAGS.  However libtool does also include CFLAGS when linking.
  case $CHECK_CFLAGS in
    *-pthread*)
      ;;
    *)
      CHECK_CFLAGS="-pthread $CHECK_CFLAGS"
      ;;
  esac
fi
AM_CONDITIONAL(ENABLE_UNIT_TESTS, [test $enable_unit_tests = "yes"])

AC_ARG_ENABLE([doc],
  AS_HELP_STRING([--disable-doc], [Do not build API documentation]),,
  [enable_doc=yes]
)
AS_IF([test "$enable_doc" = "yes"], [
  AC_CHECK_PROG([HAVE_DOXYGEN], [doxygen], [yes], [no])
  AC_CHECK_PROG([HAVE_FIG2DEV], [fig2dev], [yes], [no])
],[
  HAVE_DOXYGEN=no
  HAVE_FIG2DEV=no
])
AM_CONDITIONAL([HAVE_DOXYGEN], [test "$HAVE_DOXYGEN" = "yes"])
AM_CONDITIONAL([HAVE_FIG2DEV], [test "$HAVE_FIG2DEV" = "yes"])

AC_ARG_ENABLE([asm],
  AS_HELP_STRING([--disable-asm], [Do not compile assembly versions]),,
  [enable_asm=yes]
)

case $host_cpu in
  i[3456]86)
    cpu_x86=true
    ;;
  x86_64)
    cpu_x86=true
    ;;
  arm*)
    cpu_arm=true
    ;;
esac

TMP_CFLAGS="$CFLAGS"

AS_IF([test "$enable_asm" = "yes" -a "$cpu_x86" = "true"], [
  AC_DEFINE([OD_X86ASM], [1], [Enable asm optimisations])
  enable_gcc_inline_assembly=yes
  AC_DEFINE([OD_GCC_INLINE_ASSEMBLY], [1],
   [Enable gcc inline assembly optimisations])
  CFLAGS="$CFLAGS -msse2"
  AC_TRY_LINK([
#include <xmmintrin.h>
   ], [
    return _mm_cvtsi128_si32(_mm_setzero_si128());
   ], [
    enable_sse2_intrinsics=yes
    AC_DEFINE([OD_SSE2_INTRINSICS], [1],
     [Enable SSE2 intrinsics optimisations])
    CFLAGS="$CFLAGS -msse4.1"
    AC_TRY_LINK([
#include <smmintrin.h>
     ], [
      return _mm_cvtsi128_si32(_mm_mullo_epi32(_mm_setzero_si128(),
       _mm_setzero_si128()));
     ], [
       enable_sse41_intrinsics=yes
       AC_DEFINE([OD_SSE41_INTRINSICS], [1],
        [Enable SSE4.1 intrinsics optimisations])
       CFLAGS="$CFLAGS -mavx2"
       AC_TRY_LINK([
#include <immintrin.h>
       ], [
        return _mm_cvtsi128_si32(_mm_broadcastb_epi8(_mm_setzero_si128()));
       ], [
         enable_avx2_intrinsics=yes
         AC_DEFINE([OD_AVX2_INTRINSICS], [1],
          [Enable AVX2 intrinsics optimisations])
       ], [enable_avx2_intrinsics=no]
      )
     ], [enable_sse41_intrinsics=no]
    )
   ], [enable_sse2_intrinsics=no]
  )
])
CFLAGS="$TMP_CFLAGS"
AM_CONDITIONAL([ENABLE_X86ASM],
 [test "$enable_asm" = "yes" -a "$cpu_x86" = "true"])
AM_CONDITIONAL([ENABLE_GCC_INLINE_ASSEMBLY],
 [test "$enable_asm" = "yes" -a "$cpu_x86" = "true" -a "$enable_gcc_inline_assembly" = "yes"])
AM_CONDITIONAL([ENABLE_SSE2_INTRINSICS],
 [test "$enable_asm" = "yes" -a "$cpu_x86" = "true" -a "$enable_sse2_intrinsics" = "yes"])
AM_CONDITIONAL([ENABLE_SSE41_INTRINSICS],
 [test "$enable_asm" = "yes" -a "$cpu_x86" = "true" -a "$enable_sse41_intrinsics" = "yes"])
AM_CONDITIONAL([ENABLE_AVX2_INTRINSICS],
 [test "$enable_asm" = "yes" -a "$cpu_x86" = "true" -a "$enable_avx2_intrinsics" = "yes"])
AM_CONDITIONAL([ENABLE_ARMASM],
 [test "$enable_asm" = "yes" -a "$cpu_arm" = "true"])

AS_IF([test "$enable_asm" = "yes" -a "$cpu_arm" = "true"], [
  AC_DEFINE([OD_ARMASM], [1], [Enable ARM asm optimisations])
  AC_DEFINE([OD_ARM_MAY_HAVE_NEON], [1], [Enable ARM NEON optimisations])
])

AC_ARG_ENABLE([encoder-check],
  AS_HELP_STRING([--enable-encoder-check], [Compare reconstructed frames]),,
  [enable_encoder_check=no])

AS_IF([test "$enable_encoder_check" = "yes"], [
  AC_DEFINE([OD_ENCODER_CHECK], [1], [Enable comparison of reconstructed frames])
])
AM_CONDITIONAL([ENCODER_CHECK], [test "$enable_encoder_check" = "yes"])

AC_ARG_ENABLE([dct-check-overflow],
  AS_HELP_STRING([--enable-dct-check-overflow], [Check for DCT overflow]),,
  [enable_dct_check_overflow=no])

AS_IF([test "$enable_dct_check_overflow" = "yes"], [
  AC_DEFINE([OD_DCT_CHECK_OVERFLOW], [1], [Check for DCT overflow])
])

AC_ARG_ENABLE([check-asm],
  AS_HELP_STRING([--enable-check-asm], [Validate assembly code]),,
  [enable_check_asm=no])

AS_IF([test "$enable_check_asm" = "yes"], [
  AC_DEFINE([OD_CHECKASM], [1], [Validate assembly code])
])

AC_ARG_ENABLE([mc-animate],
  AS_HELP_STRING([--enable-mc-animate],
    [Animate mc. ./configure ANIFRAME=12 --enable-dump-images --enable-mc-animate
     animates the frame pointed by ANIFRAME]),,
  [enable_mc_animate=no])
AS_IF([test "$enable_mc_animate" = "yes"], [
  AC_DEFINE([OD_ANIMATE], [1], [Enable MC frame animation])
  AC_ARG_VAR([ANIFRAME],[the frame to be set to animate])
  if test "x$ANIFRAME" = "x" ; then
    AC_MSG_NOTICE([==> ANI FRAME IS NOT Set, defaulting to 69])
    AC_DEFINE([ANI_FRAME],[69], [Frame to animate])
  else
    AC_MSG_NOTICE([==> ANI FRAME IS Set to $ANIFRAME])
    AC_DEFINE_UNQUOTED([ANI_FRAME], [$ANIFRAME], [Frame to animate])
  fi
])

AC_ARG_ENABLE([dump-images],
  AS_HELP_STRING([--enable-dump-images], [Dump debugging images]),,
  [enable_dump_images=no])
AS_IF([test "$enable_dump_images" = "yes"], [
  AC_DEFINE([OD_DUMP_IMAGES], [1], [Enable image dumping])
  if test "$enable_tools" != "yes" ; then
    PKG_CHECK_MODULES([PNG], [libpng])
  fi
  PC_PNG_REQUIRES="libpng"
  PC_PNG_LIBS="$PNG_LIBS"
])
AM_CONDITIONAL([DUMP_IMAGES], [test "$enable_dump_images" = "yes"])
AC_SUBST([PC_PNG_REQUIRES])
AC_SUBST([PC_PNG_LIBS])

AC_ARG_ENABLE([dump-recons],
  AS_HELP_STRING([--enable-dump-recons], [Dump reconstructed video]),,
  [enable_dump_recons=no])
AS_IF([test "$enable_dump_recons" = "yes"], [
  AC_DEFINE([OD_DUMP_RECONS], [1], [Enable reconstructed video dumping])
])

AC_ARG_ENABLE([dump-bsize-dist],
  AS_HELP_STRING([--enable-dump-bsize-dist],
   [Dump block size RDO distortion]),, [enable_dump_bsize_dist=no])
AS_IF([test "$enable_dump_bsize_dist" = "yes"], [
  AC_DEFINE([OD_DUMP_BSIZE_DIST], [1],
   [Enable block size RDO distortion dumping])
])

AC_ARG_ENABLE([accounting],
  AS_HELP_STRING([--enable-accounting], [Enable bit accounting]),
  [if test "$host_mingw" = "true"; then AC_MSG_ERROR([accounting not supported on Windows]); fi],
  [enable_accounting=yes])
AS_IF([test "$enable_accounting" = "yes"], [
  AC_DEFINE([OD_ACCOUNTING], [1], [Enable bit accounting])])

AC_ARG_ENABLE([ec-accounting],
  AS_HELP_STRING([--enable-ec-accounting], [Enable entropy coder accounting]),,
  [enable_ec_accounting=no])
AS_IF([test "$enable_ec_accounting" = "yes"], [
  AC_DEFINE([OD_EC_ACCOUNTING], [1], [Enable entropy coder accounting])])

AC_ARG_ENABLE([analyzer],
  AS_HELP_STRING([--enable-analyzer], [Enable the stream analyzer]),,
  enable_analyzer=no)

AS_IF([test "$enable_analyzer" = "yes"], [
  PKG_CHECK_MODULES([OGG], [ogg >= 1.3])
  WX_CONFIG_OPTIONS
  WX_CONFIG_CHECK([2.8.0], [wxWin=1])
  if test "$wxWin" != 1 || test "$wx_ver_ok" != "yes"; then
    AC_MSG_ERROR([
      wxWidgets 2.8.0 must be installed on your system to build the stream
      analyzer. Please install or update your wxWidgets.
    ])
  fi
  if test "$enable_accounting" != "yes"; then
    AC_MSG_ERROR([
      The analyzer requires --enable-accounting to function properly.
    ])
  fi
])
AM_CONDITIONAL([ENABLE_ANALYZER], [test "$enable_analyzer" = "yes"])

AC_CONFIG_FILES([
  Makefile
  daalaenc.pc
  daaladec.pc
  daalaenc-uninstalled.pc
  daaladec-uninstalled.pc
  doc/Doxyfile
  doc/Makefile
])
AC_CONFIG_HEADERS([config.h])
AC_OUTPUT

AC_MSG_NOTICE([
------------------------------------------------------------------------
  $PACKAGE_NAME $PACKAGE_VERSION: Automatic configuration OK.

    Assertions ................... ${enable_assertions}
    Logging ...................... ${enable_logging}
    API documentation ............ ${enable_doc}
    Assembly optimizations ....... ${enable_asm}
    Image dumping ................ ${enable_dump_images}
    MC Animation  ................ ${enable_mc_animate}
    Reconstructed video dumping .. ${enable_dump_recons}
    Check encoder ................ ${enable_encoder_check}
    Check DCT overflow............ ${enable_dct_check_overflow}
    Check assembly................ ${enable_check_asm}
    Bit accounting ............... ${enable_accounting}
    Entropy coder accounting ..... ${enable_ec_accounting}
    Tools ........................ ${enable_tools}
    Unit tests ................... ${enable_unit_tests}
    Examples...................... ${enable_examples}
    Example Player................ ${enable_player}
    Stream Analyzer .............. ${enable_analyzer}
------------------------------------------------------------------------

 Use "make tools" to compile the tools.
 Use "make clean && make debug" to enable assertions and logging
 without needing to reconfigure the source tree.
])
