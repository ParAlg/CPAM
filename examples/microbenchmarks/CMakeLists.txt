function(add_benchmark NAME BLOCK_SIZE)
  add_executable(${NAME}-${BLOCK_SIZE} ${NAME}.cpp)
    target_link_libraries(${NAME}-${BLOCK_SIZE} PRIVATE parlay CPAM)
    target_compile_definitions(${NAME}-${BLOCK_SIZE} PRIVATE -DNO_AUG -DBLOCK_SIZE=${BLOCK_SIZE})
endfunction()

add_benchmark(testParallel 128)
