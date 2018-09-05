cmake_minimum_required(VERSION 3.1)
set( CMAKE_CXX_STANDARD 14 )
project (neusomatic LANGUAGES CXX)

#uncomment this for gprof
#set(CMAKE_BUILD_TYPE Debug)
#set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -g -pg")
set(CMAKE_BUILD_TYPE Release)
set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -march=native ")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fopenmp " )
set( CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -Wno-unused-variable " )
set( CMAKE_EXE_LINKER_FLAGS " -static-libgcc -static-libstdc++ -L. ")

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
