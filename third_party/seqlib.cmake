include(ExternalProject)
find_package(BZip2 REQUIRED)
if (BZIP2_FOUND)
  get_filename_component(BZIP2_LIB_DIR ${BZIP2_LIBRARIES} DIRECTORY) 
  link_directories(${BZIP2_LIB_DIR})
endif (BZIP2_FOUND)

if (NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/SeqLib)
  message("build new seqlib")

ExternalProject_Add(
  SeqLib
  PREFIX ${CMAKE_CURRENT_LIST_DIR}/SeqLib
  GIT_REPOSITORY https://github.com/walaj/SeqLib.git
  GIT_TAG 5941c68c13abca7271931bb8f6892287f3bf6d12
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ./configure --prefix=${CMAKE_CURRENT_LIST_DIR} LDFLAGS=-L${BZIP2_LIB_DIR}
  BUILD_COMMAND sed -i.bak -e "s/#include <lzma.h>/#undef HAVE_LIBLZMA/g" ${CMAKE_CURRENT_LIST_DIR}/SeqLib/src/SeqLib/htslib/cram/cram_io.c  && make CXXFLAGS='-std=c++14' CPPFLAGS='-I${BZIP2_INCLUDE_DIR}' LDFLAGS='-L${BZIP2_LIB_DIR}'
  INSTALL_COMMAND make install
)

else()

add_custom_target(SeqLib)

endif()

MESSAGE( STATUS "new target: seqlib" )
