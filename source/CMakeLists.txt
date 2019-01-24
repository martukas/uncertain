include_directories(.)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(udmstest
  ${dir}/UDoubleMS.h
  ${dir}/UDoubleTest.h
  ${dir}/udmstest.cpp
  )

add_executable(utest
  ${dir}/UDoubleMS.h
  ${dir}/UDoubleTest.h
  ${dir}/UDouble.h
  ${dir}/utest.cpp
  )
