include_directories(${PROJECT_SOURCE_DIR}/source)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(udmstest
  ${dir}/UDoubleTest.h
  ${dir}/demo1.cpp
  )
target_link_libraries(udmstest
  PRIVATE uncertain
  )

add_executable(utest
  ${dir}/UDoubleTest2.h
  ${dir}/demo2.cpp
  )
target_link_libraries(utest
  PRIVATE uncertain
  )
