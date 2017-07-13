add_definitions("-DUSE_64_BIT_PTR_CAST")

add_library(furry_3rdparty_mrf
			mrf.cpp
			ICM.cpp
			GCoptimization.cpp
			graph.cpp
			maxflow.cpp
	       	MaxProdBP.cpp
			LinkedBlockList.cpp
			regions-maxprod.cpp
       		TRW-S.cpp
			BP-S.cpp)
