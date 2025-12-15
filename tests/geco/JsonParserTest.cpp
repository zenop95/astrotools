#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "geco/JsonParser.hpp"

namespace test
{
    class JsonParserTest;
    CPPUNIT_TEST_SUITE_REGISTRATION(JsonParserTest);

    class JsonParserTest : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(JsonParserTest);
          CPPUNIT_TEST(ReadTest);
        CPPUNIT_TEST_SUITE_END();

    public:

        void setUp()
        {
        }

        void tearDown()
        {
            // Nothing
        }

    protected:

        void ReadTest()
        {
            // Open desired input file
            std::ifstream inputFile("./data/jsonParserTest.json", std::ifstream::binary);

            // Read the Json value
            Json::Value root;
            inputFile >> root;

            const Json::Value propagator = root["propagator"];
            geco::CJsonParser parser;

            // Parse double
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e-1, parser.Read<double>(propagator, "tolerance"), 1e-15);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e-2, parser.Read<double>(propagator, "minStep"), 1e-15);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(6e2,  parser.Read<double>(propagator, "maxStep"), 1e-15);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, parser.Read<double>(propagator, "initialStep"), 1e-15);
            CPPUNIT_ASSERT_THROW(parser.Read<double>(propagator, "missing"), std::invalid_argument);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0, parser.Read<double>(propagator, "missing", 42), 1e-15);

            // Parse std::string
            CPPUNIT_ASSERT_EQUAL(std::string("rungeKutta78"), parser.Read<std::string>(propagator, "type"));
            CPPUNIT_ASSERT_EQUAL(std::string("default"), parser.Read<std::string>(propagator, "missing", "default"));

            // Parse boolean
            const Json::Value bools = root["boolean"];
            CPPUNIT_ASSERT_EQUAL(true,  parser.Read<bool>(bools, "A"));
            CPPUNIT_ASSERT_EQUAL(false, parser.Read<bool>(bools, "B"));
            CPPUNIT_ASSERT_EQUAL(true,  parser.Read<bool>(bools, "C", true));
            CPPUNIT_ASSERT_EQUAL(false, parser.Read<bool>(bools, "C", false));

            // Parse boolean
            const Json::Value integers = root["integer"];
            CPPUNIT_ASSERT_EQUAL(42,  parser.Read<int>(integers, "A"));
            CPPUNIT_ASSERT_EQUAL(-37, parser.Read<int>(integers, "B"));
            CPPUNIT_ASSERT_EQUAL(51,  parser.Read<int>(integers, "C", 51));
            CPPUNIT_ASSERT_EQUAL(-28, parser.Read<int>(integers, "C", -28));
        }
    };
}
