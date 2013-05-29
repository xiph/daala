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

AC_USE_SYSTEM_EXTENSIONS
AC_SYS_LARGEFILE

AM_INIT_AUTOMAKE([1.10 foreign no-define])
AM_MAINTAINER_MODE([enable])
LT_INIT

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])

dnl Library versioning for libtool.
dnl Please update these for releases.
dnl CURRENT, REVISION, AGE
dnl - library source changed -> increment REVISION
dnl - interfaces added/removed/changed -> increment CURRENT, REVISION = 0
dnl - interfaces added -> increment AGE
dnl - interfaces removed -> AGE = 0

OP_LT_CURRENT=0
OP_LT_REVISION=1
OP_LT_AGE=0

AC_SUBST(OP_LT_CURRENT)
AC_SUBST(OP_LT_REVISION)
AC_SUBST(OP_LT_AGE)

CC_CHECK_CFLAGS_APPEND(
  [-std=c89 -pedantic -Wall -Wextra -Wno-parentheses -Wno-long-long -Wno-overlength-strings])

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
esac
AM_CONDITIONAL(OP_WIN32, test "$host_mingw" = "true")

AC_ARG_ENABLE([assertions],
  AS_HELP_STRING([--enable-assertions], [Enable assertions in code]),,
  enable_assertions=no)

AS_IF([test "$enable_assertions" = "yes"], [
  AC_DEFINE([OP_ENABLE_ASSERTIONS], [1], [Enable assertions in code])
])

PKG_CHECK_MODULES([OGG], [ogg >= 1.3])
PKG_CHECK_MODULES([SDL], [sdl])
PKG_CHECK_MODULES([PNG], [libpng])
PKG_CHECK_MODULES([THEORA], [theoradec])

#CC_ATTRIBUTE_VISIBILITY([default], [
#  CC_FLAG_VISIBILITY([CFLAGS="${CFLAGS} -fvisibility=hidden"])
#])

dnl Check for doxygen
AC_ARG_ENABLE([doc],
  AS_HELP_STRING([--disable-doc], [Do not build API documentation]),,
  [enable_doc=yes]
)
AC_CHECK_PROG(HAVE_DOXYGEN, doxygen, yes, no)
if test "$HAVE_DOXYGEN" != "yes" -o "$enable_doc" != "yes" ; then
  HAVE_DOXYGEN="no"
  enable_doc="no"
fi
AM_CONDITIONAL(HAVE_DOXYGEN, [test $HAVE_DOXYGEN = yes])

AC_CHECK_PROG(HAVE_FIG2DEV, fig2dev, yes, no)
if test "$HAVE_FIG2DEV" != "yes" -o "$enable_doc" != "yes" ; then
  HAVE_FIG2DEV="no"
  enable_doc="no"
fi
AM_CONDITIONAL(HAVE_FIG2DEV, [test $HAVE_FIG2DEV = yes])

AM_CONDITIONAL(ENABLE_X86ASM, [test yes = yes])

AM_CONDITIONAL(ENABLE_DOCS, [test $enable_doc = yes])

AC_DEFINE([OD_ENABLE_ASSERTIONS], [1], [Enable assertions in code])
AC_DEFINE([OD_LOGGING_ENABLED], [1], [Enable logging])

AC_ARG_ENABLE([dump-images],
  AS_HELP_STRING([--disable-dump-images], [Do not dump debugging images]),,
  [enable_dump_images=yes])
AS_IF([test x$enable_dump_images = xyes], [
  AC_DEFINE([OD_DUMP_IMAGES], [1], [Enable image dumping])
])
AM_CONDITIONAL([DUMP_IMAGES], [test x$enable_dump_images = xyes])

DAALA_DEC_LIBS="\$(top_builddir)/src/libdaaladec.la \$(top_builddir)/src/libdaalabase.la -lm"
AC_SUBST(DAALA_DEC_LIBS)

DAALA_ENC_LIBS="\$(top_builddir)/src/libdaalaenc.la \$(top_builddir)/src/libdaalabase.la -lm"
AC_SUBST(DAALA_ENC_LIBS)

AC_CONFIG_FILES([
  Makefile
  daala.pc
  daala-uninstalled.pc
  doc/Makefile
  src/Makefile
  src/tests/Makefile
  examples/Makefile
  include/Makefile
  tools/Makefile
])
AC_OUTPUT

AC_MSG_NOTICE([
------------------------------------------------------------------------
  $PACKAGE_NAME $PACKAGE_VERSION: Automatic configuration OK.

    Assertions ................... ${enable_assertions}

    API documentation ............ ${enable_doc}
------------------------------------------------------------------------
])
