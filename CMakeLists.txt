cmake_minimum_required(VERSION 3.4...3.18)
project(SpecMatching2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(src)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SIMD_WIDTH 128)
if (NOT(MSVC))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
         set(ARCH_OPT "-O3" "-mno-avx2")
    else ()
         set(ARCH_OPT "-O3")
    endif ()
else ()
    set(ARCH_OPT "-O2")
    # https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
    STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_MakeAvailable(googletest)

# One source of truth for the stim revision: the fetch below pins it, and `sparse_graph_stats`
# stamps the same string into its `run.log` so a campaign artifact records which stim produced it.
set(SPECMATCHING_STIM_VERSION "v1.16.0")
FetchContent_Declare(stim
        GIT_REPOSITORY https://github.com/quantumlib/Stim.git 
        GIT_TAG ${SPECMATCHING_STIM_VERSION}
      )
FetchContent_MakeAvailable(stim)
if (NOT (MSVC))
    target_compile_options(libstim PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libstim PRIVATE -fPIC ${ARCH_OPT})
endif ()

set(SOURCE_FILES_NO_MAIN
        src/specmatching/sparse_blossom/driver/namespaced_main.cc
        src/specmatching/sparse_blossom/driver/mwpm_decoding.cc
        src/specmatching/sparse_blossom/flooder/graph.cc
        src/specmatching/sparse_blossom/flooder/detector_node.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
        src/specmatching/sparse_blossom/flooder/graph_fill_region.cc
        src/specmatching/sparse_blossom/flooder/match.cc
        src/specmatching/sparse_blossom/flooder/graph_flooder.cc
        src/specmatching/sparse_blossom/matcher/alternating_tree.cc
        src/specmatching/sparse_blossom/matcher/mwpm.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
        src/specmatching/sparse_blossom/tracker/flood_check_event.cc
        src/specmatching/sparse_blossom/diagram/animation_main.cc
        src/specmatching/sparse_blossom/diagram/mwpm_diagram.cc
        src/specmatching/sparse_blossom/search/search_graph.cc
        src/specmatching/sparse_blossom/search/search_detector_node.cc
        src/specmatching/sparse_blossom/search/search_flooder.cc
        src/specmatching/sparse_blossom/driver/user_graph.cc
        src/specmatching/rand/rand_gen.cc
        src/specmatching/spec_matching/truncation/truncated_timeline.cc
        src/specmatching/spec_matching/truncation/harvest.cc
        src/specmatching/spec_matching/truncation/exposed_blossom.cc
        src/specmatching/spec_matching/manifold/ball_tables.cc
        src/specmatching/spec_matching/manifold/ball_serialize.cc
        src/specmatching/spec_matching/manifold/ball_graph.cc
        src/specmatching/spec_matching/manifold/ball_mwpm.cc
        src/specmatching/spec_matching/manifold/ball_decoding.cc
        src/specmatching/spec_matching/fusion/fusion_solve.cc
        src/specmatching/spec_matching/certificate/max_dual.cc
        src/specmatching/spec_matching/escalation/escalate.cc
        src/specmatching/spec_matching/driver/spec_matching_decoding.cc
        )

set(TEST_FILES
        src/specmatching/sparse_blossom/driver/namespaced_main.test.cc
        src/specmatching/sparse_blossom/driver/mwpm_decoding.test.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/varying.test.cc
        src/specmatching/sparse_blossom/flooder/graph.test.cc
        src/specmatching/sparse_blossom/flooder/detector_node.test.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/compressed_edge.test.cc
        src/specmatching/sparse_blossom/flooder/graph_fill_region.test.cc
        src/specmatching/sparse_blossom/flooder/match.test.cc
        src/specmatching/sparse_blossom/flooder/graph_flooder.test.cc
        src/specmatching/sparse_blossom/matcher/alternating_tree.test.cc
        src/specmatching/sparse_blossom/matcher/mwpm.test.cc
        src/specmatching/sparse_blossom/tracker/flood_check_event.test.cc
        src/specmatching/sparse_blossom/tracker/radix_heap_queue.test.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/mwpm_event.test.cc
        src/specmatching/sparse_blossom/tracker/queued_event_tracker.test.cc
        src/specmatching/sparse_blossom/tracker/cyclic.test.cc
        src/specmatching/sparse_blossom/diagram/mwpm_diagram.test.cc
        src/specmatching/sparse_blossom/search/search_graph.test.cc
        src/specmatching/sparse_blossom/search/search_flooder.test.cc
        src/specmatching/sparse_blossom/driver/user_graph.test.cc
        tests/spec_matching/truncation_timeline.test.cc
        tests/spec_matching/truncation_harvest.test.cc
        tests/spec_matching/truncation_harvest_staged.test.cc
        tests/spec_matching/truncation_exposed_blossom.test.cc
        tests/spec_matching/truncation_fuzz.test.cc
        tests/spec_matching/manifold_ball_tables.test.cc
        tests/spec_matching/manifold_ball_graph.test.cc
        tests/spec_matching/manifold_ball_identity.test.cc
        tests/spec_matching/escalation_escalate.test.cc
        tests/spec_matching/stock_on_h_certificate.test.cc
        tests/spec_matching/spec_matching_driver.test.cc
        )

set(PERF_FILES
        src/specmatching/perf/main.perf.cc
        src/specmatching/perf/util.perf.cc
        src/specmatching/sparse_blossom/driver/mwpm_decoding.perf.cc
        src/specmatching/sparse_blossom/driver/user_graph.perf.cc
        src/specmatching/sparse_blossom/flooder_matcher_interop/varying.perf.cc
        src/specmatching/sparse_blossom/tracker/radix_heap_queue.perf.cc
        )

set(PYTHON_API_FILES
        src/specmatching/sparse_blossom/driver/user_graph.pybind.cc
        src/specmatching/rand/rand_gen.pybind.cc
        src/specmatching/spec_matching/driver/spec_matching.pybind.cc
        src/specmatching/specmatching.pybind.cc
        )

add_executable(specmatching src/main.cc ${SOURCE_FILES_NO_MAIN})

target_compile_options(specmatching PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(specmatching PRIVATE -pthread ${ARCH_OPT})
endif ()
target_link_libraries(specmatching libstim)
install(TARGETS specmatching RUNTIME DESTINATION bin)

enable_testing()

add_executable(specmatching_tests ${SOURCE_FILES_NO_MAIN} ${TEST_FILES})
# The spec-matching tests share a helper header under tests/, which is outside the `src` include root.
target_include_directories(specmatching_tests PRIVATE ${CMAKE_CURRENT_LIST_DIR})
# -UNDEBUG so the debug-build invariants (asserts, and the M1 horizon gate counters) are live in the
# test binary even when the tree is configured as a Release build.
target_compile_options(specmatching_tests PRIVATE -fsanitize=address -g -fsanitize=undefined -fno-omit-frame-pointer -coverage -UNDEBUG)
if (NOT (MSVC))
    target_link_options(specmatching_tests PRIVATE -pthread -fsanitize=address -fsanitize=undefined -coverage)
else ()
    target_link_options(specmatching_tests PRIVATE -fsanitize=address -fsanitize=undefined -coverage)
endif ()
target_link_libraries(specmatching_tests GTest::gtest_main GTest::gmock_main libstim)

add_executable(specmatching_perf ${SOURCE_FILES_NO_MAIN} ${PERF_FILES})
target_compile_options(specmatching_perf PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(specmatching_perf PRIVATE -pthread -O3)
endif ()
target_link_libraries(specmatching_perf libstim)

# sparse_graph_stats: the §3 H-structure profiler. Measures the structure of the sparsified graph H
# and its matchability at horizon T — the escalation rate q, the connected components, their sizes,
# degrees and diameters, and each component's own COMPLETE/TRUNCATED status. No latency, no timers,
# no decoder output kept; `grep -c chrono benchmarks/spec_matching/sparse_graph_stats.cc` is 0.
# Not in the python wheel.
#
# The three defines below go into `run.log` so a campaign's artifact records what produced it.
# `git rev-parse` is run at configure time and left as "unknown" outside a git checkout, rather than
# failing the build for a provenance line.
execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE SPECMATCHING_GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
if (NOT SPECMATCHING_GIT_HASH)
    set(SPECMATCHING_GIT_HASH "unknown")
endif ()
string(REPLACE ";" " " SPECMATCHING_BUILD_FLAGS "${ARCH_OPT}")
# `SPECMATCHING_STIM_VERSION` is the tag `FetchContent_Declare(stim ...)` above pins, so the two
# cannot drift apart.

add_executable(sparse_graph_stats benchmarks/spec_matching/sparse_graph_stats.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(sparse_graph_stats PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(sparse_graph_stats PRIVATE ${ARCH_OPT})
target_compile_definitions(sparse_graph_stats PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS}"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(sparse_graph_stats PRIVATE -pthread -O3)
endif ()
target_link_libraries(sparse_graph_stats libstim)

# component_profiler: the per-component run-time profiler. Logs the run time of truncated sparse
# blossom on every connected component of H, keyed by component size and by the component's own
# COMPLETE/TRUNCATED status — one row per component, raw uncorrected nanoseconds. It changes no
# decoder code and runs no escalation; it stops at the trigger. Not in the python wheel.
#
# Same three provenance defines as `sparse_graph_stats`, and for the same reason: the campaign
# artifact has to record what produced it.
add_executable(component_profiler benchmarks/spec_matching/component_profiler.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(component_profiler PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(component_profiler PRIVATE ${ARCH_OPT})
target_compile_definitions(component_profiler PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS}"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(component_profiler PRIVATE -pthread -O3)
endif ()
target_link_libraries(component_profiler libstim)

# load_balancer_profiler: the k-core critical-path latency of the decomposed H solve. Per shot it
# unions H into components, spreads them over k solver cores by LPT under a size-only s^alpha cost
# model, and logs uf + balance + max_c(build_c + solve_c + extract_c) + combine against the stock
# decode on G. Single thread — no thread is ever spawned; the k-core figure is the max over serially
# executed loads. It changes no decoder code and the escalating set is the same for every k, which
# `--verify` checks shot by shot against the monolithic per-component path. Not in the python wheel.
#
# Same three provenance defines as its two companions, and for the same reason: the campaign
# artifact has to record what produced it.
add_executable(load_balancer_profiler benchmarks/spec_matching/load_balancer_profiler.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(load_balancer_profiler PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(load_balancer_profiler PRIVATE ${ARCH_OPT})
target_compile_definitions(load_balancer_profiler PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS}"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(load_balancer_profiler PRIVATE -pthread -O3)
endif ()
target_link_libraries(load_balancer_profiler libstim)

# fusion_lb_profiler: the divide-and-conquer load balancer. An oversized component of H is cut into
# pieces of at most ceil(n/k) defects inside the union-find pass that already finds the components;
# each piece is solved as an independent MWPM problem in which the cut is a virtual boundary, and the
# pieces are fused back by continuing blossom from their combined primal-dual state. Per shot it
# logs the makespan of a simulated k-core list schedule over the measured cost of every leaf, fusion
# and extraction, against the stock decode on G. Single thread — no thread is ever spawned.
#
# It is the one binary here that edits the solver: GraphFlooder gains an edge mask and a per-region
# dual cap (see src/specmatching/sparse_blossom/flooder/graph_flooder.h). Both are inert — nullptr
# and NO_HORIZON — on every instance but this binary's, which is what run.log records. Not in the
# python wheel.
#
# Same three provenance defines as its companions, and for the same reason: the campaign artifact has
# to record what produced it.
add_executable(fusion_lb_profiler benchmarks/spec_matching/fusion_lb_profiler.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(fusion_lb_profiler PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(fusion_lb_profiler PRIVATE ${ARCH_OPT})
target_compile_definitions(fusion_lb_profiler PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS}"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(fusion_lb_profiler PRIVATE -pthread -O3)
endif ()
target_link_libraries(fusion_lb_profiler libstim)

# The same binary with the §11 invariants live: `-UNDEBUG` turns the feasibility asserts of
# `fusion_solve.cc` and the profiler's debug-only checks back on even in a Release tree, which is
# what acceptance test 1 is run with.
add_executable(fusion_lb_profiler_debug benchmarks/spec_matching/fusion_lb_profiler.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(fusion_lb_profiler_debug PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(fusion_lb_profiler_debug PRIVATE ${ARCH_OPT} -UNDEBUG -g)
target_compile_definitions(fusion_lb_profiler_debug PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS} -UNDEBUG"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(fusion_lb_profiler_debug PRIVATE -pthread -O3)
endif ()
target_link_libraries(fusion_lb_profiler_debug libstim)

# lattice_surgery_profiler: docs/lattice_surgery_profiler.md. A 1xN lattice-surgery strip cut into N
# spacetime cubes: each cube runs the fusion solve on its own H with the cross-seam edges masked, a
# global manager fuses every seam component in one layer, every cube's fallback core runs phase 1 from
# t = 0, and a truncated cube is resolved locally by the phase-2 hybrid resume. Per shot it logs the
# makespan of a simulated N(k+2)+1-core schedule against stock decoding of the whole strip. Single
# thread — no thread is ever spawned. Not in the python wheel.
#
# It adds two inert fields to GraphFlooder (a masked boundary slot is skipped, and the outer-layer
# (t1) hook); both are off on every instance but this binary's fallback instances. The circuit
# generator is benchmarks/spec_matching/lattice_surgery_circuit.h. Same provenance defines as its
# companions.
add_executable(lattice_surgery_profiler benchmarks/spec_matching/lattice_surgery_profiler.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(lattice_surgery_profiler PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(lattice_surgery_profiler PRIVATE ${ARCH_OPT})
target_compile_definitions(lattice_surgery_profiler PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS}"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(lattice_surgery_profiler PRIVATE -pthread -O3)
endif ()
target_link_libraries(lattice_surgery_profiler libstim)

# The debug twin of §9.1: `-UNDEBUG` turns on the phase-2 invariant checks (H-only defects at least
# T from Omega after every closure, no H-only dual above T, injections only on uncovered nodes, the
# speculative path's dual feasibility after every completed job) and the solver's own asserts.
add_executable(lattice_surgery_profiler_debug benchmarks/spec_matching/lattice_surgery_profiler.cc ${SOURCE_FILES_NO_MAIN})
target_include_directories(lattice_surgery_profiler_debug PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(lattice_surgery_profiler_debug PRIVATE ${ARCH_OPT} -UNDEBUG -g)
target_compile_definitions(lattice_surgery_profiler_debug PRIVATE
        SPECMATCHING_GIT_HASH="${SPECMATCHING_GIT_HASH}"
        SPECMATCHING_BUILD_FLAGS="${SPECMATCHING_BUILD_FLAGS} -UNDEBUG"
        SPECMATCHING_STIM_VERSION="${SPECMATCHING_STIM_VERSION}")
if (NOT (MSVC))
    target_link_options(lattice_surgery_profiler_debug PRIVATE -pthread -O3)
endif ()
target_link_libraries(lattice_surgery_profiler_debug libstim)

add_library(libspecmatching ${SOURCE_FILES_NO_MAIN})
set_target_properties(libspecmatching PROPERTIES PREFIX "")
target_include_directories(libspecmatching PUBLIC src)
target_compile_options(libspecmatching PRIVATE ${ARCH_OPT} -fPIC)
if(NOT(MSVC))
    target_link_options(libspecmatching PRIVATE -pthread -O3)
endif()
target_link_libraries(libspecmatching libstim)
install(TARGETS libspecmatching LIBRARY DESTINATION)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/" DESTINATION "include" FILES_MATCHING PATTERN "*.h" PATTERN "*.inl")

include(GoogleTest)
gtest_discover_tests(specmatching_tests)

IF(EXISTS ${CMAKE_CURRENT_LIST_DIR}/extern/pybind11/CMakeLists.txt)
    add_subdirectory(extern/pybind11)
    pybind11_add_module(_cpp_specmatching ${PYTHON_API_FILES} ${SOURCE_FILES_NO_MAIN})
    target_link_libraries(_cpp_specmatching PRIVATE libstim)
    target_compile_options(_cpp_specmatching PRIVATE ${ARCH_OPT})
else()
    message("WARNING: Skipped the pybind11 module _cpp_specmatching because the `pybind11` git submodule isn't present. To fix, run `git submodule update --init --recursive`")
endif()

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(_cpp_specmatching
#        PRIVATE VERSION_INFO=${SPECMATCHING_VERSION_INFO})
