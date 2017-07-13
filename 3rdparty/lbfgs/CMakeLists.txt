set(CMAKE_C_FLAGS "-msse2 -DUSE_SSE -DHAVE_EMMINTRIN_H -O3 -ffast-math  -Wall")
# set(CMAKE_C_FLAGS " -O3 -ffast-math  -Wall")

add_library(furry_3rdparty_lbfgs lbfgs.c)
