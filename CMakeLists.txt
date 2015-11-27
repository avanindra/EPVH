cmake_minimum_required( VERSION 3.0 )

project( EPVH )

find_package( OpenCV )
find_package( Eigen )
find_package( OpenGL )
find_package( VTK )

set( EIGEN_INCLDUE_DIR "" CACHE PATH "set eigen include dir" )

set( EPVH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/externalinclues )
set( EPVH_LINK_LIBRARIES ${OpenCV_LIBS} )

include_directories(  ${EPVH_INCLUDE_DIRS}  )

add_library( epvh basevh.cpp camerainfo.cpp epvh.cpp segmentclipper.cpp silhouette.cpp )


add_subdirectory(externalinclues)
