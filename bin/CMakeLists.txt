set(BIN_SOURCES
    lauritzen_1992_example.cpp
    neapolitan_example_4_1.cpp)
foreach(bin_source_file ${BIN_SOURCES} )
    # I used a simple string replace, to cut off .cpp.
    string( REPLACE ".cpp" "" bin_name ${bin_source_file} )
    add_executable( ${bin_name} ${bin_source_file} )
    # Make sure YourLib is linked to each app
    target_link_libraries( ${bin_name} CLG )
endforeach( bin_source_file ${BIN_SOURCES} )
