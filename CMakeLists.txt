######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.11)
project(GFAse VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

# Compiler options.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (dev)
    message(STATUS "---- Building with debug options ----")

    # Memory debug
    set(CMAKE_CXX_FLAGS "-fexceptions -fsanitize=address")
    set(ASAN_OPTIONS=check_initialization_order=1)
    set(ASAN_OPTIONS=detect_leaks=1)

    add_definitions(-ggdb3 -Og -Wall)       # Debugging + No optimization

else()
    message(STATUS "---- Building with optimization ----")

    # Standard compilation
    add_definitions(-O3 -Wall)              # Much optimization
endif()



#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/align.cpp
        src/Bam.cpp
        src/BinaryIO.cpp
        src/BinarySequence.cpp
        src/binomial.cpp
        src/Bipartition.cpp
        src/Bridges.cpp
        src/BubbleGraph.cpp
        src/chain.cpp
        src/Chainer.cpp
        src/ContactGraph.cpp
        src/Color.cpp
        src/edge.cpp
        src/FixedBinarySequence.cpp
        src/GfaReader.cpp
        src/gfa_to_handle.cpp
        src/graph_utility.cpp
        src/HamiltonianPath.cpp
        src/HamiltonianChainer.cpp
        src/HaplotypePathKmer.cpp
        src/Hasher.cpp
        src/Hasher2.cpp
        src/handle_to_gfa.cpp
        src/IncrementalIdMap.cpp
        src/KmerSets.cpp
        src/misc.cpp
        src/MurmurHash2.cpp
        src/MurmurHash3.cpp
        src/MultiContactGraph.cpp
        src/optimize.cpp
	src/Overlaps.cpp
        src/VectorMultiContactGraph.cpp
        ##        src/OverlapMap.cpp
        src/Phase.cpp
        src/PhaseAssign.cpp
        src/Sequence.cpp
        src/Sam.cpp
        src/SubgraphOverlay.cpp
        src/SvgPlot.cpp
        src/Timer.cpp
        )


project(GFAse)
add_library(GFAse STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(GFAse PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(GFAse PROPERTIES DEFINE_SYMBOL "")


# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

##############################################
# ------------------------------------------ #
# --------- LINKING NATIVE LIBRARY --------- #
# ------------------ OMP ------------------- #
# ------------------------------------------ #
##############################################


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    # TODO: can't we though?

    # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
    # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
    if(EXISTS /usr/local/Cellar/libomp/)
        set(OMPDIR /usr/local/Cellar/libomp/)
    elseif(EXISTS /opt/local/include/libomp)
        set(OMPDIR /opt/local/include/libomp)
    endif()

    set (OpenMP_C "${CMAKE_C_COMPILER}")
    set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I${OMPDIR} -I/usr/local/include -L${OMPDIR} -L/usr/local/lib")
    set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I${OMPDIR} -I/usr/local/include -L${OMPDIR} -L/usr/local/lib")
    set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set (OpenMP_libomp_LIBRARY "omp")
    set (OpenMP_libgomp_LIBRARY "gomp")
    set (OpenMP_libiomp5_LIBRARY "iomp5")

    # and now add the OpenMP parameters to the compile flags
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")

    # Mac needs libdl and libomp when linking the library
    set(PLATFORM_EXTRA_LIB_FLAGS -ldl -lomp)


elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    find_package(OpenMP REQUIRED)

    # To link statically, it is necessary to call lower level commands, as opposed to just
    # target_link_library(target OpenMP::OpenMP_CXX), or else make tries to -lgomp which is apparently
    # bad because it brings in dynamic libraries dlsym, dlopen, dlwtf
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()

endif()


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- htslib ------------------ #
# ------------------------------------------ #
##############################################

message(STATUS "Fetching htslib")

find_package(ZLIB REQUIRED)

# Dependencies
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'z' library")
endif()

find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(WARNING "Couldn't find the 'bz2' library")
endif()

find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(WARNING "Couldn't find the 'curl' library")
endif()

include(FetchContent)

FetchContent_Declare(
        project_htslib
        URL      https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
)

FetchContent_MakeAvailable(project_htslib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/htslib/
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/cram/
        ${CMAKE_BINARY_DIR}/_deps/project_htslib-src/os/
        ${CMAKE_SOURCE_DIR}/external/
)

message(STATUS "htslib_SOURCE_DIR: ${project_htslib_SOURCE_DIR}")


set(HTS_DIR _deps/project_htslib-src)

# HTSLIB configure
add_custom_target(
        BUILD_HTS
        ALL
        WORKING_DIRECTORY ${HTS_DIR}
        COMMAND pwd
        COMMAND autoreconf -i
        COMMAND ./configure --disable-lzma --disable-libcurl --disable-bz2 --disable-s3 --disable-gcs --without-libdeflate --disable-plugins
        COMMAND $(MAKE) print-config
        COMMAND $(MAKE) prefix=${CMAKE_SOURCE_DIR}/external/htslib/ install
)


add_library(htslib STATIC IMPORTED)
set_property(TARGET htslib
        PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/external/htslib/lib/libhts.a)


add_dependencies(htslib BUILD_HTS)
add_dependencies(GFAse htslib)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# --------------- minimap2 ----------------- #
# ------------------------------------------ #
##############################################


message(STATUS "Fetching minimap2")


# https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz

include(FetchContent)

FetchContent_Declare(
        project_minimap2
        GIT_REPOSITORY https://github.com/lh3/minimap2.git
        GIT_TAG v2.24
)

FetchContent_MakeAvailable(project_minimap2)

set(MINIMAP2_SOURCES
        ${project_minimap2_SOURCE_DIR}/ksw2_dispatch.c
        ${project_minimap2_SOURCE_DIR}/lchain.c
        ${project_minimap2_SOURCE_DIR}/format.c
        ${project_minimap2_SOURCE_DIR}/example.c
        ${project_minimap2_SOURCE_DIR}/ksw2_extz2_sse.c
        ${project_minimap2_SOURCE_DIR}/ksw2_extd2_sse.c
        ${project_minimap2_SOURCE_DIR}/kalloc.c
        ${project_minimap2_SOURCE_DIR}/kthread.c
        ${project_minimap2_SOURCE_DIR}/hit.c
        ${project_minimap2_SOURCE_DIR}/ksw2_ll_sse.c
        ${project_minimap2_SOURCE_DIR}/index.c
        ${project_minimap2_SOURCE_DIR}/ksw2_exts2_sse.c
        ${project_minimap2_SOURCE_DIR}/align.c
        ${project_minimap2_SOURCE_DIR}/esterr.c
        ${project_minimap2_SOURCE_DIR}/bseq.c
        ${project_minimap2_SOURCE_DIR}/main.c
        ${project_minimap2_SOURCE_DIR}/map.c
        ${project_minimap2_SOURCE_DIR}/options.c
        ${project_minimap2_SOURCE_DIR}/pe.c
        ${project_minimap2_SOURCE_DIR}/sdust.c
        ${project_minimap2_SOURCE_DIR}/sketch.c
        ${project_minimap2_SOURCE_DIR}/splitidx.c
        ${project_minimap2_SOURCE_DIR}/seed.c
        ${project_minimap2_SOURCE_DIR}/misc.c
        )

message(STATUS "minimap2_SOURCE_DIR: ${project_minimap2_SOURCE_DIR}")


set(MINIMAP2_DIR _deps/project_minimap2-src)

add_custom_target(
        BUILD_MINIMAP2
        ALL
        WORKING_DIRECTORY ${MINIMAP2_DIR}
        COMMAND ""
)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/
        #        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/sse2neon/
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/arm/
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/arm/neon/
        ${CMAKE_BINARY_DIR}/_deps/project_minimap2-src/lib/simde/simde/x86/
)

#add_custom_target(HTSLIB_CONFIGURED DEPENDS "${CMAKE_BINARY_DIR}/_deps/project_htslib-src/config.h")

add_library(minimap2 ${MINIMAP2_SOURCES})
add_dependencies(minimap2 BUILD_MINIMAP2)
add_dependencies(GFAse minimap2)

target_link_libraries(minimap2 pthread bz2 z curl)


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- sparsepp ---------------- #
# ------------------------------------------ #
##############################################


ExternalProject_Add(project_sparsepp
        GIT_REPOSITORY https://github.com/greg7mdp/sparsepp.git
        TIMEOUT 5
        PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND "" # Disable install step, is a header only lib!
        )

add_dependencies(GFAse project_sparsepp)

# Specify include dir
ExternalProject_Get_Property(project_sparsepp SOURCE_DIR)
include_directories(${SOURCE_DIR}/sparsepp)

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ bdsg ------------------ #
# ------------------------------------------ #
##############################################

message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

if (dev)
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
            DOWNLOAD_COMMAND ""
            UPDATE_COMMAND ""
            PREFIX ${CMAKE_SOURCE_DIR}/external/bdsg/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make install
            )
else()
    # Download or update library as an external project
    ExternalProject_Add(project_bdsg
            GIT_REPOSITORY https://github.com/vgteam/libbdsg.git
            PREFIX ${CMAKE_SOURCE_DIR}/external/bdsg/
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/bdsg/ -DRUN_DOXYGEN=OFF -DBUILD_PYTHON_BINDINGS=OFF
            BUILD_IN_SOURCE True
            INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/bdsg/
            INSTALL_COMMAND make install
            )
endif()

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_bdsg INSTALL_DIR)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LIB_SUFFIX "dylib")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LIB_SUFFIX "so")
endif()


# Create new library for external project (so it can be linked with main library)
add_library(bdsg SHARED IMPORTED)
set_property(TARGET bdsg
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libbdsg.${LIB_SUFFIX})

add_library(divsufsort SHARED IMPORTED)
set_property(TARGET divsufsort
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libdivsufsort64.${LIB_SUFFIX})

add_library(libhandlegraph SHARED IMPORTED)
set_property(TARGET libhandlegraph
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhandlegraph.${LIB_SUFFIX})

add_library(libsdsl SHARED IMPORTED)
set_property(TARGET libsdsl
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libsdsl.${LIB_SUFFIX})

# Define library as dependent on the downloaded project
add_dependencies(bdsg
        project_bdsg
        libsdsl
        libhandlegraph
        divsufsort)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(GFAse bdsg)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/bdsg/include/
        external/bdsg/include/handlegraph/
        external/bdsg/include/bdsg/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GFAse PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/bdsg")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")


set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# -------- TESTS --------

set(TESTS
        test_alignment_chain
        test_assign_phase
        test_binomial
        test_bfs
        test_binary_sequence
        test_binary_sequence_performance
        test_bridges
        test_bubble_graph
        test_bubblegraph_chaining
        test_bubble_align
        test_connected_component_finder
        test_contact_graph
        test_chainer
        test_fixed_binary_sequence
        test_fixed_binary_sequence_performance_2
        test_fixed_binary_sequence_sparsepp_performance
        test_gfareader
        test_hamiltonian_chainer
        test_hamiltonian_path
        test_htslib
        test_htslib_bam_reader
        test_incremental_id_io
        test_kmer_unordered_set
	test_overlaps
        test_phase_haplotype_paths
        test_minimap2
        test_minimap2_no_io
        test_multi_contact_graph
        test_multi_contact_graph_io
        test_nonbinary_sequence_performance
        test_nonbinary_sequence_sparsepp_performance
        test_rgb_to_hex
        test_rechain
        test_set_intersection
        test_timer
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            GFAse
            Threads::Threads
            ZLIB::ZLIB
            bdsg
            #            ${Boost_LIBRARIES}
            htslib
            divsufsort
            libhandlegraph
            libsdsl
            minimap2
            #            mummer4
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        assign_phases_via_diploid_alignment
        count_kmers
        create_bandage_path_color_table
        compute_minhash2
        extract_haplotype_kmers
        evaluate_contacts
        evaluate_phasing
        find_bubbles
        find_shasta_bubbles
        get_mapq_distribution_from_sam
        get_degree_stats
        get_path_lengths
        generate_contact_map_from_bam
        locate_kmer_matches
        phase_haplotype_paths
        phase_contacts
        phase_contacts_with_monte_carlo
        rephase_tripartition_with_monte_carlo
        remove_empty_nodes
        solve_maxcut
        split_connected_components
        unzip
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            GFAse
            Threads::Threads
            ZLIB::ZLIB
            bdsg
            divsufsort
            htslib
            libhandlegraph
            libsdsl
            minimap2
            )

endforeach()

#set_target_properties(gfase PROPERTIES LINK_FLAGS "-static" )
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")


# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

## The shared library goes to the bin directory.
#install(TARGETS
#        GFAse
#        gfase
#        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")