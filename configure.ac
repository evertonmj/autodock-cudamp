#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([autodock], [4.0.1], [autodock@scripps.edu])
AC_CONFIG_SRCDIR([analysis.cc])

AM_INIT_AUTOMAKE()
#AM_INIT_AUTOMAKE([foreign])
#AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for programs.
AC_PROG_AWK

# don't let AC_PROG_CXX force "-g -O2"
saved_CXXFLAGS=${CXXFLAGS}
AC_PROG_CXX
CXXFLAGS=${saved_CXXFLAGS}

AC_PROG_INSTALL
AC_PROG_RANLIB

AC_CANONICAL_HOST

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gethostname pow sqrt strncasecmp])

case $host_os in
  *cygwin* ) CYGWIN=yes;;
         * ) CYGWIN=no;;
esac
AM_CONDITIONAL([WANT_BIGSTACK], [test "$CYGWIN" = yes])

#Cuda compilation parameter
AC_ARG_ENABLE(cuda,
 [  --enable-cuda  Enables Cuda Support (needs nvcc)],
 [ cuda=yes
 ])

if test "x$cuda" = "xyes" ; then 
    echo -n "checking for nvcc... "

    ########################
    # Check for CUDASDK
    ########################
    AC_ARG_WITH([cudasdk],
        [  --with-cudasdk=PATH        prefix where cudasdk is installed [default=auto]])
    if test "x$with_cudasdk" != x; then
      CUDASDK="/usr/local/cuda/NVIDIA_CUDA_SDK"
    else
      CUDASDK="$with_cudasdk"
    fi

    ########################
    # Check for CUDA
    ########################
    AC_ARG_WITH([cuda],
        [  --with-cuda=PATH        prefix where cuda is installed [default=auto]])
    if test "x$with_cuda" = x; then
      CUDA_CFLAGS="-I/usr/local/cuda/include -I$CUDASDK/common/inc"
      CUDA_LIBS=" -L/usr/local/cuda/bin -L/usr/local/cuda/lib -L$CUDASDK/lib -L$CUDASDK/common/lib/linux -lcuda -lcudart"
      NVCC="nvcc"
    else
      CUDA_CFLAGS="-I$with_cuda/include -I$CUDASDK/common/inc"
      CUDA_LIBS=" -L$with_cuda/bin -L$with_cuda/lib -L$CUDASDK/lib -L$CUDASDK/common/lib/linux -lcuda -lcudart"
      NVCC="$with_cuda/bin/nvcc"
    fi

    # ADD TEST FOR NVCC HERE
    NVCCTEST=`which "$NVCC" 2> /dev/null`
    if test "x$NVCCTEST" = x; then
        echo "no"
        AC_MSG_ERROR([no nvcc in path ("$PATH")])
    else
        echo "yes"
        AC_DEFINE(CUDA_READY, 1)
    fi
        
    #NVCCFLAGS="-O2 -use_fast_math -I. -I/home/nldudok1/NVIDIA_CUDA_SDK-2.0/common/inc"
    NVCCFLAGS="-D_DEBUG -g -use_fast_math -I. -I$CUDASDK/common/inc -DCUDA_READY"

    AC_SUBST(CUDASDK)
    AC_SUBST(CUDA_CFLAGS)
    AC_SUBST(CUDA_LIBS)
    AC_SUBST(NVCC)
    AC_SUBST(NVCCFLAGS)

    # Check for CUDA headers usability
    cuda_save_CPPFLAGS=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS $CUDA_CFLAGS $CUDA_LIBS"
#    CXXFLAGS+="$CXXFLAGS $CUDA_CFLAGS"

    AC_CHECK_HEADERS([cuda_runtime.h], [],
                [])
    AC_CHECK_HEADERS([cuda.h], [],
                [])
    AC_LANG_PUSH(C++)
    AC_LANG_PUSH(C++)
#    AC_CHECK_HEADERS([cutil.h], [],
#                [],[])
    AC_LANG_POP
#    CPPFLAGS="$cuda_save_CPPFLAGS"
    CXXFLAGS+=" $CUDA_CFLAGS $CUDA_LIBS "
#    CXX=$NVCC

fi

AM_CONDITIONAL([WANT_CUDA], [test "$cuda" = yes])

# See "automake" manual, version 1.10, section 10.5 (Python)
#AM_PATH_PYTHON([2.1],, [:])
#AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
