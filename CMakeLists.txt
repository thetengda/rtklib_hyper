cmake_minimum_required(VERSION 3.0.0)

project(rtklib)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(-O3)
    set(CMAKE_BUILD_POSTFIX "${CMAKE_RELEASE_POSTFIX}")
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/release)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/release/lib)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_options(-g3 -Wall -Wextra)   
    set(CMAKE_BUILD_POSTFIX "${CMAKE_DEBUG_POSTFIX}")
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/debug)
    set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/debug/lib)
endif()

aux_source_directory(src/ SRC_RTKLIB)
aux_source_directory(src/rcv/ SRC_RTKLIB_RCV)

set(SRC_RTKLIB ${SRC_RTKLIB} ${SRC_RTKLIB_RCV})

# find_package(LAPACK REQUIRED)

add_executable(rnx2rtkp app/rnx2rtkp.c ${SRC_RTKLIB} ${SRC_RTKLIB_RCV})
target_compile_definitions(rnx2rtkp PRIVATE
    ENAGLO ENAGAL ENACMP ENAQZS ENAIRN ENALEO 
    NFREQ=3 MAXNPCV=50 
    OUTSTAT_AMB TRACE)
target_include_directories(rnx2rtkp PRIVATE src/)
if(${LAPACK_FOUND})
    target_compile_definitions(rnx2rtkp PRIVATE LAPACK)
endif()
target_link_libraries(rnx2rtkp m pthread ${LAPACK_LIBRARIES})