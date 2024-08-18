add_executable(BC BC.cc)
  target_link_libraries(BC PRIVATE parlay CPAM aspen)

add_executable(BFS BFS.cc)
  target_link_libraries(BFS PRIVATE parlay CPAM aspen)

add_executable(MIS MIS.cc)
  target_link_libraries(MIS PRIVATE parlay CPAM aspen)

add_executable(Flatsnap Flatsnap.cc)
  target_link_libraries(Flatsnap PRIVATE parlay CPAM aspen)
