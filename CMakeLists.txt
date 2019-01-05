cmake_minimum_required( VERSION 3.0 )

project( EPVH )

find_package( OpenCV )
#find_package( Eigen )
find_package( OpenGL )
find_package( VTK )


if( VTK_FOUND )

include( ${VTK_USE_FILE} )

endif()

set( EIGEN_INCLDUE_DIR "" CACHE PATH "set eigen include dir" )

set( EXECUTABLE_OUTPUT_PATH "" CACHE PATH "set eigen include dir" )
set( LIBRARY_OUTPUT_PATH "" CACHE PATH "set library output path" )

set( EPVH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/externalinclues )
set( EPVH_LINK_LIBRARIES ${OpenCV_LIBS} )

add_definitions( -DDATASET_DIR="${CMAKE_SOURCE_DIR}/datasets" )

include_directories(  ${EPVH_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} ${EIGEN_INCLDUE_DIR}  )

add_library( epvh basevh.cpp camerainfo.cpp epvh.cpp segmentclipper.cpp silhouette.cpp contourio.cpp
                  display3droutines.cpp math2d.cpp silhouette.cpp math3d.cpp sortingpredicates.cpp )

target_link_libraries( epvh ${OpenCV_LIBS} ${VTK_LIBRARIES} )

add_executable( epvhdemo epvh_demo.cpp )

target_link_libraries( epvhdemo epvh )


add_subdirectory(externalinclues)
