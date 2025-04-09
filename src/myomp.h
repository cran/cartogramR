// Compatibility define for LLVM builds of libomp.
#ifdef KMP_VERSION_BUILD
#ifndef _OPENMP
#define _OPENMP KMP_VERSION_BUILD
#endif
#endif

#ifdef _OPENMP
  #ifdef R_DISABLEMATCH
    #undef match
  #endif
  #include <omp.h>
#else
  // for machines with compilers void of openmp support
  #define omp_get_num_threads()  1
  #define omp_get_thread_num()   0
  #define omp_get_max_threads()  1
  #define omp_get_thread_limit() 1
  #define omp_get_num_procs()    1
  #define omp_get_wtime()        0
#endif
