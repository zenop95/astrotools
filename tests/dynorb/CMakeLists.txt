add_executable(dynorb_unit_tests
    ${TEST_RUNNER}
    AidaTest.cpp
    )

target_link_libraries(dynorb_unit_tests
    ${CPPUNIT_LIBRARIES}
    ${DACE_LIBRARIES}
    atmos
    astro
    dynorb
    odeprop
    )

target_include_directories(dynorb_unit_tests
    PUBLIC ${CMAKE_SOURCE_DIR}/modules
)

add_test(NAME dynorbTest COMMAND dynorb_unit_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/dynorb)