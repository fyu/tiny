///////////////////////////////////////////////////////////////////////////////////////////////
//// 
////  Compiler-specific definitions for the sparse Levenberg - Marquardt minimization algorithm
////  Copyright (C) 2005-2010 Manolis Lourakis (lourakis at ics forth gr)
////  Institute of Computer Science, Foundation for Research & Technology - Hellas
////  Heraklion, Crete, Greece.
////
////  This program is free software; you can redistribute it and/or modify
////  it under the terms of the GNU General Public License as published by
////  the Free Software Foundation; either version 2 of the License, or
////  (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _COMPILER_H_
#define _COMPILER_H_

/* note: intel's icc defines both __ICC & __INTEL_COMPILER.
 * Also, some compilers other than gcc define __GNUC__,
 * therefore gcc should be checked last
 */

#ifdef _MSC_VER // MSVC

#define inline __inline
#define SPLM_FINITE _finite
#if _MSC_VER > 1400 // VC8 provides __restrict
#define restrict __restrict
#else
#define restrict // unsupported
#endif /* _MSC_VER >= 1400 */

typedef signed __int32    int32_t;
typedef unsigned __int32  uint32_t;

#elif defined(__ICC) || defined(__INTEL_COMPILER) || defined(__GNUC__) // ICC, GCC

#define SPLM_FINITE finite
#if defined(__SSE__) // || defined(__pentium4__) || defined(__tune_pentiumpro__)
/* non-temporal & temporal prefetches */
#define __PREFETCH_NT(mem) __asm__ __volatile__ ("prefetchnta %0" : : "m" (*((char *)(mem))))
#define __PREFETCH_T(mem)  __asm__ __volatile__ ("prefetcht0 %0" : : "m" (*((char *)(mem))))
#else
#define __PREFETCH_NT(mem) // empty
#define __PREFETCH_T(mem) // empty
#endif

#include <stdint.h> // 32 bit ints
#define restrict __restrict__ // GCC & ICC support both  __restrict__ & __restrict

#else // other than MSVC, ICC, GCC
/* If the compiler lacks stdint.h, a portable version can be found at
 * http://www.azillionmonkeys.com/qed/pstdint.h
 */
#include <pstdint.h>

#define inline // empty
#define SPLM_FINITE finite // let's hope this will work
#define restrict // empty

#endif /* _MSC_VER */

//#define USE_CLOCK // uncomment to force using clock() for time measurements

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) // branch prediction hints for GCC
/* make sure our GCC supports them */
#if __GNUC__ == 2 && __GNUC_MINOR__ < 96
#define __builtin_expect(x, expected_x) (x)
#endif
#define likely(x)       __builtin_expect((x), 1)
#define unlikely(x)     __builtin_expect((x), 0)
#else // not GCC
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif /* __GNUC__ */

#endif /* _COMPILER_H_ */
