cmake_minimum_required(VERSION 3.0)
project(MCRT)

set(CMAKE_CXX_STANDARD 11)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (MSVC)
    add_compile_options("/Zi")
    add_link_options("/DEBUG")
endif()

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

set(acceleration_sources
    "src/core/acceleration/bvh.h"
    "src/core/acceleration/bvh.cpp"
    "src/core/acceleration/sahbvh.h"
    "src/core/acceleration/sahbvh.cpp"
    "src/core/acceleration/splitbvh.h"
    "src/core/acceleration/splitbvh.cpp"
)

set(basic_sources
    "src/core/aabb.h"
    "src/core/aabb.cpp"
    "src/core/boundedtriangle.h"
    "src/core/boundedtriangle.cpp"
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
    "src/core/logger.h"
    "src/core/logger.cpp"
)

file(GLOB_RECURSE thirdparty_srcs "src/thirdparty/*")

add_library(thirdparty STATIC ${thirdparty_srcs})

target_include_directories(
    thirdparty PUBLIC
    "src/thirdparty/nanort"
    "src/thirdparty/glm/include"
    "src/thirdparty/nanolog"
)

set_target_properties(thirdparty PROPERTIES FOLDER core)

find_package(embree)
if (EMBREE_LIBRARY)
    add_definitions(-DWITH_EMBREE=1)
    include_directories(${EMBREE_INCLUDE_DIRS})
    message("with embree")
else()
    add_definitions(-DWITH_EMBREE=0)
    message("without embree")
endif()

# rtcore
add_library(rtcore
    ${advanced_sources}
    ${basic_sources}
    ${utility_sources}
    ${acceleration_sources}
)

set_target_properties(rtcore PROPERTIES FOLDER core)

source_group("core\\advanced" FILES ${advanced_sources})
source_group("core\\basic" FILES ${basic_sources})
source_group("core\\acceleration" FILES ${acceleration_sources})
source_group("core\\utility" FILES ${utility_sources})

target_include_directories(rtcore PUBLIC "src/core")

target_link_libraries(rtcore PUBLIC thirdparty)

# path tracing
add_executable(demo_pt
    "src/examples/pt/demo_pt.cpp"
)

target_link_libraries(demo_pt
    rtcore
    ${EMBREE_LIBRARY}
)

set_target_properties(demo_pt PROPERTIES FOLDER demos)

# bdpt
add_executable(demo_bdpt
    "src/examples/bdpt/bdpt_impl.h"
    "src/examples/bdpt/demo_bdpt.cpp"
)

target_link_libraries(demo_bdpt
    rtcore
    ${EMBREE_LIBRARY}
)

set_target_properties(demo_bdpt PROPERTIES FOLDER demos)

# photonmap
add_executable(demo_photonmap
    "src/examples/photonmap/photonmap.h"
    "src/examples/photonmap/photonmap.cpp"
    "src/examples/photonmap/demo_photonmap.cpp"
)

target_link_libraries(demo_photonmap
    rtcore
    ${EMBREE_LIBRARY}
)

set_target_properties(demo_photonmap PROPERTIES FOLDER demos)

# pssmlt
add_executable(demo_pssmlt
    "src/examples/pssmlt/demo_pssmlt.cpp"
)

target_link_libraries(demo_pssmlt
    rtcore
    ${EMBREE_LIBRARY}
)

set_target_properties(demo_pssmlt PROPERTIES FOLDER demos)

# vcm
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

set_target_properties(demo_vcm PROPERTIES FOLDER demos)

# bssrdf
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

set_target_properties(demo_bssrdf PROPERTIES FOLDER demos)

# consistent normal
add_executable(demo_consistentnormal
    "src/examples/consistent_normal/demo_consistent_normal.cpp"
)

target_link_libraries(demo_consistentnormal
    rtcore
    ${EMBREE_LIBRARY}
)

set_target_properties(demo_consistentnormal PROPERTIES FOLDER demos)

# split bvh
add_executable(demo_splitbvh
    "src/examples/splitbvh/demo_splitbvh.cpp"
    "src/examples/splitbvh/debug_mesh.h"
    "src/examples/splitbvh/debug_mesh.cpp"
    "src/examples/splitbvh/debugintersection.h"
    "src/examples/splitbvh/debugintersection.cpp"
)

target_link_libraries(demo_splitbvh
    rtcore
    ${EMBREE_LIBRARY}
)

set_target_properties(demo_splitbvh PROPERTIES FOLDER demos)

# --------------------- tests -----------------------------
enable_testing()

add_test(NAME test_pt COMMAND demo_pt)
add_test(NAME test_bdpt COMMAND demo_bdpt)
add_test(NAME test_photonmap COMMAND demo_photonmap)
add_test(NAME test_pssmlt COMMAND demo_pssmlt)
add_test(NAME test_vcm COMMAND demo_vcm)
add_test(NAME test_bssrdf COMMAND demo_bssrdf)
add_test(NAME test_consistentnormal COMMAND demo_consistentnormal)
# add_test(NAME test_splitbvh COMMAND demo_splitbvh)
