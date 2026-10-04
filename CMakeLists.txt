cmake_minimum_required(VERSION 3.20)
project(gjkepa LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(-std=gnu -ffree-line-length-none -Wall)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(-standard-semantics)
endif()

find_package(OpenMP)

add_library(gjkepa STATIC
    src/geom/geom3d.f90
    src/collision/gjkepa.f90
    src/hull/quickhull.f90
)
set_target_properties(gjkepa PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(gjkepa PUBLIC OpenMP::OpenMP_Fortran)
endif()

add_executable(example_two_polyhedra examples/two_polyhedra.f90)
target_link_libraries(example_two_polyhedra PRIVATE gjkepa)
set_target_properties(example_two_polyhedra PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

enable_testing()
add_executable(test_gjkepa tests/test_gjkepa.f90)
target_link_libraries(test_gjkepa PRIVATE gjkepa)
set_target_properties(test_gjkepa PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
add_test(NAME gjkepa COMMAND test_gjkepa)

add_executable(test_quickhull tests/test_quickhull.f90)
target_link_libraries(test_quickhull PRIVATE gjkepa)
set_target_properties(test_quickhull PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
add_test(NAME quickhull COMMAND test_quickhull)

# Shared C ABI for the Python verification bench (ctypes).
add_library(gjkepa_c SHARED src/collision/gjkepa_capi.f90)
target_link_libraries(gjkepa_c PRIVATE gjkepa)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(gjkepa_c PRIVATE OpenMP::OpenMP_Fortran)
endif()
set_target_properties(gjkepa_c PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    POSITION_INDEPENDENT_CODE ON
)
