#https://seqan.readthedocs.io/en/seqan-v2.0.2/BuildManual/IntegrationWithYourOwnBuildSystem.html
if (NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/seqan)

message("build new seqan ")
ExternalProject_Add(
  seqan
  PREFIX ${CMAKE_CURRENT_LIST_DIR}/seqan
  GIT_REPOSITORY https://github.com/seqan/seqan.git
  CONFIGURE_COMMAND "" 
  BUILD_COMMAND "" 
  INSTALL_COMMAND ""
)

else()

add_custom_target(seqan)

endif()
set (CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_HAS_ZLIB=1 " )

MESSAGE( STATUS "new target: seqan" )
