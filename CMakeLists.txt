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

if (0)

add_executable(demo_photonmap
    "src/examples/photonmap/photonmap.h"
    "src/examples/photonmap/photonmap.cpp"
    "src/examples/photonmap/demo_photonmap.cpp"
)

target_link_libraries(demo_photonmap
    rtcore
    ${EMBREE_LIBRARY}
)

add_executable(demo_pssmlt
    "src/examples/pssmlt/demo_pssmlt.cpp"
)

target_link_libraries(demo_pssmlt
    rtcore
    ${EMBREE_LIBRARY}
)

add_executable(demo_vcm
    "src/examples/vcm/vcm_material.h"
    "src/examples/vcm/vcm_misc.h"
    "src/examples/vcm/vcm_scene.h"
    "src/examples/vcm/demo_vcm.cpp"
)

target_link_libraries(demo_vcm
    rtcore
    ${EMBREE_LIBRARY}
)

add_executable(demo_bssrdf
    "src/examples/bssrdf/bssrdf.h"
    "src/examples/bssrdf/bssrdf.cpp"
    "src/examples/bssrdf/mesh_sample.h"
    "src/examples/bssrdf/mesh_sample.cpp"
    "src/examples/bssrdf/octree.h"
    "src/examples/bssrdf/octree.cpp"
    "src/examples/bssrdf/demo_bssrdf.cpp"
)

target_link_libraries(demo_bssrdf
    rtcore
    ${EMBREE_LIBRARY}
)

add_executable(demo_consistentnormal
    "src/examples/consistent_normal/demo_consistent_normal.cpp"
)

target_link_libraries(demo_consistentnormal
    rtcore
    ${EMBREE_LIBRARY}
)

endif()