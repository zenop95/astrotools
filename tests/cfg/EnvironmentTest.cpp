#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <iomanip>

#include "cfg/Config.hpp"

namespace test
{
    class EnvironmentTest;
    CPPUNIT_TEST_SUITE_REGISTRATION(EnvironmentTest);

    class EnvironmentTest : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(EnvironmentTest);
          CPPUNIT_TEST(DefaultTest);
          CPPUNIT_TEST(NonDefaultTest);
          CPPUNIT_TEST(EnvConfigTest);
        CPPUNIT_TEST_SUITE_END();

    public:

        cfg::CConfig testConfig;

        void setUp()
        {
            // Nothing
            CPPUNIT_ASSERT(not cfg::EnvConfig.IsLoaded());
        }

        void tearDown()
        {
            // Nothing
        }

    protected:

        void DefaultTest()
        {
            CPPUNIT_ASSERT(not cfg::EnvConfig.IsLoaded());

            // Load environment configuration
            CPPUNIT_ASSERT(not testConfig.IsLoaded());
            const std::string& defVal = testConfig().Value("someDefaultValue");
            CPPUNIT_ASSERT(testConfig.IsLoaded());
            CPPUNIT_ASSERT_EQUAL(std::string("defaultValue"), defVal);
        }

        void NonDefaultTest()
        {
            // Load environment configuration
            CPPUNIT_ASSERT(not testConfig.IsLoaded());
            if(not testConfig.IsLoaded())
            {
                testConfig.SetEnvironmentFile("./data/env.json");
            }
            CPPUNIT_ASSERT(not testConfig.IsLoaded());
            const std::string& nonDefVal = testConfig().Value("nonDefaultValue");
            CPPUNIT_ASSERT(testConfig.IsLoaded());
            CPPUNIT_ASSERT_EQUAL(std::string("somFilePath"), nonDefVal);
        }

        void EnvConfigTest()
        {
            // Load environment configuration
            CPPUNIT_ASSERT(not cfg::EnvConfig.IsLoaded());
            const std::string& defVal = cfg::EnvConfig().Value("someDefaultValue");
            CPPUNIT_ASSERT(cfg::EnvConfig.IsLoaded());
            CPPUNIT_ASSERT_EQUAL(std::string("defaultValue"), defVal);
        }
    };
}
