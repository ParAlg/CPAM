add_executable(run_simultaneous_updates_queries-CPAM-CPAM-Diff run_simultaneous_updates_queries.cc)
  target_link_libraries(run_simultaneous_updates_queries-CPAM-CPAM-Diff PRIVATE parlay CPAM aspen)
  target_compile_definitions(run_simultaneous_updates_queries-CPAM-CPAM-Diff PRIVATE -DUSE_DIFF_ENCODING)

