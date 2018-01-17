
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(INSTALL_INCLUDE_PATH include/furry/common)

INCLUDE_DIRECTORIES(${CGAL_CORE_INCLUDE_DIR}
                    /usr/local/include/)

add_library(
		furry_common
		camera.cc
		cast.cc
		clock.cc
		color_map.cc
		cv.cc
		data.cc
		dir.cc
		draw.cc
		gaps.cc
		geometry.cc
		init.cc
		io.cc
		log.cc
		numeric.cc
		patch.cc
		path.cc
		point.cc
		point_set.cc
		pose.cc
		protobuf.cc
		rand.cc
		str.cc
		tbb.cc
		util.cc
)

target_link_libraries(furry_common ${TBB_LIBRARIES} glog ${GFLAGS_LIBRARIES})

install (TARGETS furry_common DESTINATION lib/${PROJECT_NAME})
install (FILES
        algebra.h
		camera.h
		cast.h
		check.h
		common.h
		color_map.h
		cv.h
		dir.h
		draw.h
		eigen.h
		factory.h
		functional.h
		gaps.h
		geometry.h
		io.h
		io-inl.h
		log.h
		memory.h
		numeric.h
		patch.h
		point.h
		point-inl.h
		point_set.h
		pose.h
		rand.h
		range.h
		stl.h
		str.h
		clock.h
		tbb.h
		DESTINATION ${INSTALL_INCLUDE_PATH})
