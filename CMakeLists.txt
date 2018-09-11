cmake_minimum_required(VERSION 3.0)
project(MCRT)

set(CMAKE_CXX_STANDARD 11)

find_package(openmp REQUIRED)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(advanced_sources
    "src/core/envmap.h"
    "src/core/envmap.cpp"
    "src/core/lightpath.h"
    "src/core/lightpath.cpp"
    "src/core/volumetric.h"
    "src/core/volumetric.cpp"
)

set(basic_sources
    "src/core/aabb.h"
    "src/core/aabb.cpp"
    "src/core/bsdf.h"
    "src/core/bsdf.cpp"
    "src/core/camera.h"
    "src/core/camera.cpp"
    "src/core/film.h"
    "src/core/film.cpp"
    "src/core/intersection.h"
    "src/core/intersection.cpp"
    "src/core/light.h"
    "src/core/light.cpp"
    "src/core/material.h"
    "src/core/material.cpp"
    "src/core/mesh.h"
    "src/core/mesh.cpp"
    "src/core/ray.h"
    "src/core/ray.cpp"
    "src/core/scene.h"
    "src/core/scene.cpp"
    "src/core/scenebvh.h"
    "src/core/scenebvh.cpp"
    "src/core/texture.h"
    "src/core/texture.cpp"
    "src/core/triangle.h"
    "src/core/triangle.cpp"
    "src/core/vectors.h"
    "src/core/vectors.cpp"
)

set(utility_sources
    "src/core/consoledebug.h"
    "src/core/consoledebug.cpp"
    "src/core/constants.h"
    "src/core/constants.cpp"
    "src/core/geometry.h"
    "src/core/geometry.cpp"
    "src/core/numeric.h"
    "src/core/numeric.cpp"
    "src/core/sample.h"
    "src/core/sample.cpp"
    "src/core/timer.h"
    "src/core/timer.cpp"
)

include_directories(
    "src/thirdparty/nanort"
    "src/thirdparty/glm/include"
)

find_package(embree)
if (EMBREE_LIBRARY)
    add_definitions(-DWITH_EMBREE=1)
    include_directories(${EMBREE_INCLUDE_DIRS})
    message("with embree")
else()
    add_definitions(-DWITH_EMBREE=0)
    message("without embree")
endif()


add_library(rtcore
    ${advanced_sources}
    ${basic_sources}
    ${utility_sources}
)

source_group("core\\advanced" FILES ${advanced_sources})
source_group("core\\basic" FILES ${basic_sources})
source_group("core\\utility" FILES ${utility_sources})

include_directories(
    src/core
    src/thirdparty
)

add_executable(demo_pt
    "src/examples/pt/demo_pt.cpp"
)

target_link_libraries(demo_pt
    rtcore
    ${EMBREE_LIBRARY}
)

add_executable(demo_bdpt
    "src/examples/bdpt/bdpt_impl.h"
    "src/examples/bdpt/demo_bdpt.cpp"
)

target_link_libraries(demo_bdpt
    rtcore
    ${EMBREE_LIBRARY}
)
