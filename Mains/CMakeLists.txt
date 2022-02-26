# List of main sources

set( MAIN_FILES
  Spectrum_test.cc
  MuonRange_test.cc
  BField_test.cc
  MuStops_test.cc
  Tracks_test.cc
  ELossDists_test.cc
  MoyalDist_test.cc
)
# Generate test targets
foreach( mainfile ${MAIN_FILES} )
    # get the name of the test test from the source path
    # (strip _test.cc, and the path preceding the basename)
    string( REPLACE "_test.cc" "" testnamenoext ${mainfile} )
    get_filename_component(testname ${testnamenoext} NAME)

    # prepend Main_ to the target name to avoid target name clashes,
    # create the executable:
    add_executable( Main_${testname} ${mainfile} )

    # add the project root as an include directory
    # link all test tests to KinKal, MatEnv, and ROOT libraries.
    target_include_directories(Main_${testname} PRIVATE ${PROJECT_SOURCE_DIR}/.. ${KINKAL_INCLUDE_DIR})
    target_link_libraries( Main_${testname}  Spectra Detector General Test ${KINKAL_LIBRARIES} ${ROOT_LIBRARIES} )

    # ensure the test test executable filename is just its test name
    set_target_properties( Main_${testname} PROPERTIES OUTPUT_NAME ${testname})

    # register the target as a test
    add_test (NAME ${testname} COMMAND Main_${testname} )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 250)
    set_tests_properties(${testname} PROPERTIES ENVIRONMENT "TRACKTOY_SOURCE_DIR=${CMAKE_SOURCE_DIR}")

    install( TARGETS Main_${testname}
             RUNTIME DESTINATION bin/ )

endforeach( mainfile ${MAIN_FILES} )
