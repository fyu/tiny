/////////////////////////////////////////////////////////////////////////////////
//// 
////  Wrappers to various OS-dependent time measurement functions
////  Copyright (C) 2009 Manolis Lourakis (lourakis at ics forth gr)
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
///////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>

#include "compiler.h"
#include "splm.h"

#if !defined(_WIN32) && !defined(USE_CLOCK) // assume we have getrusage()

#include <sys/resource.h>

/* return user time in seconds */
double splm_gettime()
{
struct rusage usage;

  getrusage(RUSAGE_SELF, &usage);

  return usage.ru_utime.tv_sec+usage.ru_utime.tv_usec/1000000.0;
}

#elif defined(_MSC_VER) && !defined(USE_CLOCK)

#include <windows.h>
#include <mmsystem.h>

/* return elapsed (i.e., wall clock) time in seconds */
double splm_gettime()
{
LARGE_INTEGER tpsec;
LARGE_INTEGER time;

  if(!QueryPerformanceFrequency(&tpsec)){
    fprintf(stderr, "QueryPerformanceFrequency() failed (error %lu)\n", GetLastError());
    return 0;
  }

  if(!QueryPerformanceCounter(&time)){
    fprintf(stderr, "QueryPerformanceCounter() failed (error %lu)\n", GetLastError());
    return 0;
  }

  return ((double)(time.QuadPart))/tpsec.QuadPart;
}

#else // !defined(_MSC_VER) || defined(USE_CLOCK)

#include <time.h>

/* return time in seconds (wraps around approx. every 35.8 min).
 * Note that under unix, clock() returns the user time, while
 * under windows it returns the elapsed (i.e., wall clock) time;
 * on a given system, the elapsed time is longer compared to
 * user time.
 */
double splm_gettime()
{
clock_t time;

  time=clock();

  return ((double)time)/CLOCKS_PER_SEC;
}
#endif
