add_executable(scan_alignments scan_alignments.cpp)
add_executable(msa msa.cpp)
target_link_libraries(scan_alignments ${CPP_LIBS} z)
target_link_libraries(msa ${CPP_LIBS} z)
