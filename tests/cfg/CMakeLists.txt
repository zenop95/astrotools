add_executable(cfg_unit_tests
    ${TEST_RUNNER}
    EnvironmentTest.cpp
    )

target_link_libraries(cfg_unit_tests
    ${CPPUNIT_LIBRARIES}
    ${DACE_LIBRARIES}
    cfg
    geco
    )

target_include_directories(cfg_unit_tests
    PUBLIC ${CMAKE_SOURCE_DIR}/modules
)

add_test(NAME cfgTest COMMAND cfg_unit_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/cfg)