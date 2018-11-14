set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADERS
  ${dir}/functions.hpp
  ${dir}/double_ct.hpp
  ${dir}/double_ensemble.hpp
  ${dir}/double_ms.hpp
  ${dir}/double_msc.hpp
  ${dir}/scaled_array.hpp
  ${dir}/simple_array.hpp
  ${dir}/source_set.hpp
  )
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/uncertain)

set(uncertain_headers ${HEADERS})
set(uncertain_sources ${dir}/functions.cpp)

#add_doxygen_source_deps(${uncertain_headers})

add_library(uncertain SHARED
  ${uncertain_sources}
  ${uncertain_headers}
  )

#add_doxygen_source_deps(${uncertain_headers})

target_compile_definitions(uncertain PRIVATE DLL_BUILD)

#
# remove the absolute path from the library name
#
set(UNCERTAIN_INTERFACE_LIBS)

target_include_directories(uncertain
  PRIVATE ${PROJECT_SOURCE_DIR}/source
  INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

set_target_properties(uncertain
  PROPERTIES
  VERSION ${PROJECT_VERSION}
  SO_VERSION 0
  CXX_STANDARD_REQUIRED ON
  OUTPUT_NAME uncertain
  )

enable_coverage(uncertain)

target_link_libraries(uncertain
  PUBLIC ${UNCERTAIN_LINKS}
  PRIVATE ${COVERAGE_LIBRARIES} ${UNCERTAIN_PRIVATE_LIBS}
  INTERFACE ${UNCERTAIN_INTERFACE_LIBS})

include(GNUInstallDirs) #Required for CMAKE_INSTALL_XXX

install(TARGETS uncertain
  EXPORT uncertain_targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT runtime
  )

install(EXPORT uncertain_targets
  DESTINATION ${CMAKE_INSTALL_PACKAGEDIR}
  COMPONENT development)
