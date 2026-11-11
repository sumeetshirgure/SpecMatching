# Builds hc_dump, the h_construct table compiler / vector generator, against the repository's own
# library. The root project is pulled in EXCLUDE_FROM_ALL, so only libspecmatching and libstim are
# compiled. Point FETCHCONTENT_SOURCE_DIR_STIM / _GOOGLETEST at an existing checkout (run_crosscheck.sh
# uses obj/_deps when present) to configure without network access.
#
# (This directory is not called `tools` because the root .gitignore ignores every `tools` path.)
cmake_minimum_required(VERSION 3.14)
project(hc_tablegen CXX)

set(SPECMATCHING_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
add_subdirectory("${SPECMATCHING_ROOT}" specmatching EXCLUDE_FROM_ALL)

add_executable(hc_dump hc_dump.cc)
target_include_directories(hc_dump PRIVATE "${SPECMATCHING_ROOT}" "${SPECMATCHING_ROOT}/src")
target_compile_options(hc_dump PRIVATE -O2)
target_link_libraries(hc_dump libspecmatching libstim)
