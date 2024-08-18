add_library(aspen INTERFACE)
set(ASPEN_INCLUDE_DIR "${CPAM_SOURCE_DIR}/examples/graphs")
target_include_directories(aspen INTERFACE
  $<BUILD_INTERFACE:${ASPEN_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_subdirectory(algorithms)
add_subdirectory(run_batch_updates)
