project(Analysis)
cmake_minimum_required(VERSION 2.8)

#ROOT libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

#Generate the dictionary for vector<vector>
ROOT_GENERATE_DICTIONARY( G__Vector LINKDEF ../inc/LinkDef.h )

add_library(Vector SHARED G__Vector.cxx)
target_link_libraries(Vector ${ROOT_LIBRARIES})

#add_library(minuitlib {$ROOTSYS}/lib/libMinuit2)
#target_link_libraries(minuitlib ${ROOT_LIBRARIES})

#Bring the headers
include_directories(inc)
include_directories($ENV{EIGEN_INSTALL})

#Add the CMakeLists.txt on the subdirectory
add_subdirectory(src ${PROJECT_SOURCE_DIR}/build/sourcebuild)

#Compile custom libraries
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/inc/*.h)

#add here the name of your executable
#add_executable(waveform ${PROJECT_SOURCE_DIR}/waveform/waveform.cc ${sources} ${headers})
#target_link_libraries(waveform ${ROOT_LIBRARIES} Vector CustomClasses)

add_executable(charge_ana_1 ${PROJECT_SOURCE_DIR}/macro/charge_ana_1.cc ${sources} ${headers} )
target_link_libraries( charge_ana_1 ${ROOT_LIBRARIES} Vector CustomClasses -lMinuit )

set(ENV{PATH} "${PROJECT_SOURCE_DIR}/build/$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${PROJECT_SOURCE_DIR}/build/$ENV{LD_LIBRARY_PATH}")
