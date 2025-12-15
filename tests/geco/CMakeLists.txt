add_executable(geco_unit_tests
    ${TEST_RUNNER}
    FactoryTest.cpp
    JsonParserTest.cpp
    )

target_link_libraries(geco_unit_tests
    ${CPPUNIT_LIBRARIES}
    ${JSON_LIBRARIES}
    geco
    )

target_include_directories(geco_unit_tests
    PUBLIC ${CMAKE_SOURCE_DIR}/modules
)

add_test(NAME gecoTest COMMAND geco_unit_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/geco)