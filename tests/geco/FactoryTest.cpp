#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "geco/Factory.hpp"
#include "geco/JsonParser.hpp"

namespace test
{
    class FactoryTest;
    CPPUNIT_TEST_SUITE_REGISTRATION(FactoryTest);

    class FactoryTest : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(FactoryTest);
          CPPUNIT_TEST(DefaultTest);
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

        class CBaseClass
        {
            public:
            CBaseClass() = default;
            virtual ~CBaseClass() = default;


            virtual CBaseClass* Create(const Json::Value& input, geco::CJsonParser& parser) const = 0;
            virtual std::string Print() const = 0;
        };

        class CDerivedA : public CBaseClass
        {
            private:
            double d_value;
            public:
            CDerivedA() : d_value(-1) {};
            CDerivedA(const Json::Value& input, geco::CJsonParser& parser) : d_value(parser.Read<double>(input,"value")) {};
            virtual CDerivedA* Create(const Json::Value& input, geco::CJsonParser& parser) const override
            {
                return new CDerivedA(input, parser);
            }
            virtual std::string Print() const override
            {
                std::ostringstream stream;
                stream << "Value: " << d_value << std::endl;
                return stream.str();
            }
        };

        class CDerivedB : public CBaseClass
        {
            private:
            double d_value;
            std::string s_string;
            public:
            CDerivedB() : d_value(-1) {};
            CDerivedB(const Json::Value& input, geco::CJsonParser& parser)
            : d_value(parser.Read<double>(input, "someValue")),
              s_string(parser.Read<std::string>(input, "label"))
            {
                // Empty
            };
            virtual CDerivedB* Create(const Json::Value& input, geco::CJsonParser& parser) const override
            {
                return new CDerivedB(input, parser);
            }
            virtual std::string Print() const override
            {
                std::ostringstream stream;
                stream << s_string << ": " << d_value << std::endl;
                return stream.str();
            }
        };

        class CBaseFactory : public geco::CFactory<CBaseClass>
        {
            public:
            CBaseFactory()
            {
                this->m_register["derived A"] = std::make_unique<CDerivedA>();
                this->m_register["derived B"] = std::make_unique<CDerivedB>();
            }
        };

        void DefaultTest()
        {
            // Open desired input file
            std::ifstream inputFile("./data/factoryExample.json", std::ifstream::binary);
            geco::CJsonParser parser;

            // Read the Json value
            Json::Value root;
            inputFile >> root;

            // Create the factory object
            CBaseFactory factory;

            // Create object of type a
            CBaseClass* baseObject = factory.Create(root["exampleA"], parser);
            // Test the print function (polymorphism)
            CPPUNIT_ASSERT_EQUAL(std::string("Value: 42\n"), baseObject->Print());

            // Create object of type b
            baseObject = factory.Create(root["exampleB"], parser);
            // Test the print function (polymorphism)
            CPPUNIT_ASSERT_EQUAL(std::string("Other: 51\n"), baseObject->Print());

            // Test that an exception is thrown when the wrong type is created
            CPPUNIT_ASSERT_THROW(factory.Create(root["wrongExample"], parser), std::runtime_error);

            delete baseObject;
        }
    };
}
