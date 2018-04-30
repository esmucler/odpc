#include "autoconfig.h"

#ifdef HAVE_INTTYPES_H
#   include <inttypes.h>
#endif

#ifdef HAVE_STDINT_H
#   include <stdint.h>
#endif

#ifndef __cplusplus
#   undef RESTRICT
#   define RESTRICT restrict
#elif !defined(RESTRICT)
#   define RESTRICT
#endif

#ifndef DEBUG
#   define ARMA_NO_DEBUG
#endif

#ifndef HAVE_OPENMP_CXX
#   define ARMA_DONT_USE_OPENMP
#endif
