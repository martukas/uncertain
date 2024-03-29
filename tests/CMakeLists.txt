#
# We definitely need a better solution for this!
#
#if (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
#  add_definitions(/wd4251 /wd4275)
#endif ()

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(test_sources
    ${dir}/main.cpp
    ${dir}/functions.cpp
    ${dir}/double_ms.cpp
    #  ${dir}/double_msc.cpp
    #  ${dir}/double_ct.cpp
    #  ${dir}/double_ensemble.cpp
)

set(test_headers ${dir}/test_lib/gtest_print.hpp ${dir}/test_lib/color_bash.hpp)
add_executable(unit_tests EXCLUDE_FROM_ALL ${test_sources} ${test_headers})
add_dependencies(unit_tests uncertain)
target_include_directories(
    unit_tests
    PRIVATE ${GTEST_INCLUDE_DIRS}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${PROJECT_SOURCE_DIR}/source
)
target_link_libraries(
    unit_tests
    PRIVATE uncertain
    PRIVATE ${GTest_LIBRARIES}
    PRIVATE ${CMAKE_THREAD_LIBS_INIT}
    PRIVATE ${COVERAGE_LIBRARIES}
)
enable_coverage(unit_tests)

add_custom_target(
    run_tests
    COMMAND unit_tests "--gtest_output=xml:unit_tests_run.xml"
    DEPENDS unit_tests
)
add_test(NAME "RunGoogleTests" COMMAND run_tests)

create_coverage_targets(
  coverage
  run_tests
  ${PROJECT_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/source/uncertain
  ${PROJECT_BINARY_DIR}/coverage
)
