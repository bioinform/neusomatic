cmake_minimum_required(VERSION 3.1)
set( CMAKE_CXX_STANDARD 14 )
project (neusomatic LANGUAGES CXX)

#uncomment this for gprof
#set(CMAKE_BUILD_TYPE Debug)
#set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -g -pg")
set(CMAKE_BUILD_TYPE Release)

if (APPLE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
else()
    # We're in sane linux world
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -Wno-unused-variable " )

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GXX_VERSION)
    if (GXX_VERSION VERSION_GREATER 5.3.0 OR GXX_VERSION VERSION_EQUAL 5.3.0)
    else()
        message(FATAL_ERROR "current GCC/G++ version is ${GXX_VERSION}Installation abort. Need GCC/G++ >= 5.3.0")
    endif()

    set( CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -pthread")

    if (APPLE)
        message("-- MAC GCC/G++ version ${GXX_VERSION}")
    elseif (NOT APPLE)
        set( CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -L. ")
    endif()
else()
    message(FATAL_ERROR "GCC/G++ is not found on the system. Please install GCC/G++5.")
endif()

include (../third_party/seqlib.cmake)
include (../third_party/seqan.cmake)
set( CPP_LIBS ${CPP_LIBS} seqlib bwa fml hts bz2)
set( CPP_LIBS ${CPP_LIBS} dl )
set( CPP_LIBS ${CPP_LIBS} c)


include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(../third_party/SeqLib/src/SeqLib )
include_directories(../third_party/SeqLib/src/SeqLib/htslib )
include_directories(../third_party/seqan/src/seqan/include)

link_directories(../third_party/lib)
link_directories(../third_party/SeqLib/src/SeqLib/bin/)

add_subdirectory(cpp)
