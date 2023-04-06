include_directories(${PROJECT_SOURCE_DIR}/source)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(ms_demo
  ${dir}/UDoubleTest.hpp
  ${dir}/ms_demo.cpp
  )
target_link_libraries(ms_demo
  PRIVATE uncertain
  )

add_executable(full_demo
  ${dir}/UDoubleTest2.hpp
  ${dir}/full_demo.cpp
  )
target_link_libraries(full_demo
  PRIVATE uncertain
  )

add_custom_target(examples DEPENDS
  ms_demo
  full_demo
  )
