add_executable(atmos_unit_tests
    ${TEST_RUNNER}
    Nrlmsise00Test.cpp
    )

target_link_libraries(atmos_unit_tests
    ${CPPUNIT_LIBRARIES}
    ${DACE_LIBRARIES}
    atmos
    cfg
    geco
    )

target_include_directories(atmos_unit_tests
    PUBLIC ${CMAKE_SOURCE_DIR}/modules
)

add_test(atmos atmos_unit_tests)