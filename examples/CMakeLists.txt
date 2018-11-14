include_directories(.)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(udmstest
  ../include/functions.h
  ../include/UDoubleMS.h
  ${dir}/UDoubleTest.h
  ${dir}/demo1.cpp
  )

add_executable(utest
  ../include/functions.h
  ../include/UDoubleMS.h
  ../include/UDoubleMSC.h
  ../include/UDoubleCT.h
  ../include/UDoubleEnsemble.h
  ../include/SimpleArray.h
  ../include/ScaledArray.h
  ${dir}/UDoubleTest2.h
  ${dir}/demo2.cpp
  )
