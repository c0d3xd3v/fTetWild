message(STATUS "Try to install Tetwild libs .....")
FILE(GLOB_RECURSE lib_list *.a)
list(LENGTH lib_list lib_list_count)
FOREACH(_LIB ${lib_list})
    get_filename_component(FN ${_LIB} NAME)
    file(INSTALL ${_LIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
ENDFOREACH()

