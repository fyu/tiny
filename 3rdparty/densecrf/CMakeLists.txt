include_directories( ../lbfgs/ )

add_library(furry_3rdparty_densecrf util.cpp permutohedral.cpp unary.cpp pairwise.cpp objective.cpp densecrf.cpp labelcompatibility.cpp )

add_library(furry_3rdparty_densecrf_opt optimization.cpp )
target_link_libraries(furry_3rdparty_densecrf_opt furry_3rdparty_lbfgs)
