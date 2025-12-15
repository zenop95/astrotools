#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <iostream>
#include <iomanip>

#include <cspice/SpiceUsr.h>
#include <dace/dace.h>

#include "dynorb/DynClass.h"
#include "astro/AstroLibrary.h"
#include "dynorb/AIDA.h"
#include "dynorb/AIDAwrappers.h"
#include "odeprop/RK78Propagator.hpp"

#define _USE_MATH_DEFINES
#define WITH_ALGEBRAICMATRIX

namespace test
{
    class AidaTest;
    CPPUNIT_TEST_SUITE_REGISTRATION(AidaTest);

    class AidaTest : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(AidaTest);
            CPPUNIT_TEST(PropTest);
            CPPUNIT_TEST(JsonPropTest);
        CPPUNIT_TEST_SUITE_END();

    public:

        void setUp()
        {
            // Initialize DA environment
            DACE::DA::init(4,6);
        }

        void tearDown()
        {

        }

    protected:

        void PropTest()
        {
            // Load SPICE kernel (taking it from the configuration file)
            cfg::EnvConfig().LoadSpiceKernels();

            const unsigned int gravOrd = 15;
            std::string gravmodel = "./data/egm2008";
            int AIDA_flags[3] = { 2, 3, 2 };

            // parameters
            double mass     = 1000.;
            double A_drag   = 1;
            double Cd       = 2.2;
            double A_srp    = 1;
            double Cr       = 1.3100000000;

            double Bfactor    = Cd*A_drag/mass;
            double SRPC     = Cr*A_srp/mass;

            DACE::AlgebraicVector<double>  x0cons(6), xfcons(6), x0consECEF(6), xfconsECEF(6), consMapD(6), dx(6);
            DACE::AlgebraicMatrix<double> ECI2ECEF(6,6), ECEF2ECI(6,6);
            DACE::AlgebraicVector<DACE::DA> x0(6), x0ECEF(6), xf(6), xfECEF(6);

            // Initial state
            x0cons[0] =  6.6526831668497343e+03;
            x0cons[1] =  0.0000000000000000e+00;
            x0cons[2] =  0.0000000000000000e+00;
            x0cons[3] =  0.0000000000000000e+00;
            x0cons[4] = -8.8474643737438263e-01;
            x0cons[5] =  7.6917814141592569e+00;
            double Td =  5.3999992602326583e+03;
            double n  =  4.1000000000000000e+00;

            double t0 = 536479200; // ephemeris time i.e. seconds after 1th January 2000
            double tf = t0 + Td*n;

            // Reference: propagation with double
            AIDADynamics<double> aidaCartDynDouble(gravmodel, gravOrd, AIDA_flags, Bfactor, SRPC);
            xfcons = RK78( 6, x0cons, t0, tf, aidaCartDynDouble);

            // Check the final state of the propagation with doubles
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[0]/ 5.3815576787098671e+03, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[1]/-4.2357330689082323e+02, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[2]/ 3.8862021991162565e+03, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[3]/-4.5516388931775174e+00, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[4]/-7.3471781832961347e-01, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[5]/ 6.2167070442332859e+00, 1e-8 );

            // DA initialization of initial state
            x0[0] = x0cons[0] + 1e-4*DACE::DA(1);
            x0[1] = x0cons[1] + 1e-4*DACE::DA(2);
            x0[2] = x0cons[2] + 1e-4*DACE::DA(3);
            x0[3] = x0cons[3] + 1e-7*DACE::DA(4);
            x0[4] = x0cons[4] + 1e-7*DACE::DA(5);
            x0[5] = x0cons[5] + 1e-7*DACE::DA(6);

            // Check of initial state
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[0]), x0cons[0], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[1]), x0cons[1], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[2]), x0cons[2], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[3]), x0cons[3], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[4]), x0cons[4], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[5]), x0cons[5], 1e-15 );

            AIDADynamics<DACE::DA> aidaCartDyn(gravmodel, gravOrd, AIDA_flags, Bfactor, SRPC);
            xf = RK78( 6, x0, t0, tf, aidaCartDyn);

            // Check constant part of final state
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[0])/xfcons[0], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[1])/xfcons[1], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[2])/xfcons[2], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[3])/xfcons[3], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[4])/xfcons[4], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[5])/xfcons[5], 1e-8 );
        }

        void JsonPropTest()
        {
            // Load SPICE kernel (taking it from the configuration file)
            cfg::EnvConfig().LoadSpiceKernels();

            // Open desired input file
            std::ifstream inputFile("./data/input.json", std::ifstream::binary);

            // Read the Json value
            Json::Value root;
            inputFile >> root;

            // Get dynamics
            const Json::Value dynamics = root["dynamics"];

            // parameters
            double mass     = 1000.;
            double A_drag   = 1;
            double Cd       = 2.2;
            double A_srp    = 1;
            double Cr       = 1.3100000000;

            double Bfactor    = Cd*A_drag/mass;
            double SRPC     = Cr*A_srp/mass;

            DACE::AlgebraicVector<double>  x0cons(6), xfcons(6), x0consECEF(6), xfconsECEF(6), consMapD(6), dx(6);
            DACE::AlgebraicMatrix<double> ECI2ECEF(6,6), ECEF2ECI(6,6);
            DACE::AlgebraicVector<DACE::DA> x0(6), x0ECEF(6), xf(6), xfECEF(6);

            // Initial state
            x0cons[0] =  6.6526831668497343e+03;
            x0cons[1] =  0.0000000000000000e+00;
            x0cons[2] =  0.0000000000000000e+00;
            x0cons[3] =  0.0000000000000000e+00;
            x0cons[4] = -8.8474643737438263e-01;
            x0cons[5] =  7.6917814141592569e+00;
            double Td =  5.3999992602326583e+03;
            double n  =  4.1000000000000000e+00;

            double t0 = 536479200; // ephemeris time i.e. seconds after 1th January 2000
            double tf = t0 + Td*n;
            odeprop::RK78Propagator<double> rk78;
            rk78.setTolerance(1e-10);
            rk78.setMinStepSize(1e-2);
            rk78.setMaxStepSize(600);
            rk78.setInitialStepSize(60.0);
            rk78.setMaxAttempts(100);

            // Reference: propagation with double
            dynorb::AIDA<double> aidaDyn(dynamics, geco::CJsonParser());
            aidaDyn.setUncParam("Bfactor",Bfactor);
            aidaDyn.setUncParam("SRPC",SRPC);
            xfcons = rk78.propagate(aidaDyn, x0cons, t0, tf);

            // Check the final state of the propagation with doubles
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[0]/ 5.3815576787098671e+03, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[1]/-4.2357330689082323e+02, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[2]/ 3.8862021991162565e+03, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[3]/-4.5516388931775174e+00, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[4]/-7.3471781832961347e-01, 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, xfcons[5]/ 6.2167070442332859e+00, 1e-8 );

            // DA initialization of initial state
            x0[0] = x0cons[0] + 1e-4*DACE::DA(1);
            x0[1] = x0cons[1] + 1e-4*DACE::DA(2);
            x0[2] = x0cons[2] + 1e-4*DACE::DA(3);
            x0[3] = x0cons[3] + 1e-7*DACE::DA(4);
            x0[4] = x0cons[4] + 1e-7*DACE::DA(5);
            x0[5] = x0cons[5] + 1e-7*DACE::DA(6);

            // Check of initial state
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[0]), x0cons[0], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[1]), x0cons[1], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[2]), x0cons[2], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[3]), x0cons[3], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[4]), x0cons[4], 1e-15 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(DACE::cons(x0[5]), x0cons[5], 1e-15 );

            odeprop::RK78Propagator<DACE::DA> rk78DA;
            rk78DA.setTolerance(1e-10);
            rk78DA.setMinStepSize(1e-2);
            rk78DA.setMaxStepSize(600);
            rk78DA.setInitialStepSize(60.0);
            rk78DA.setMaxAttempts(100);

            dynorb::AIDA<DACE::DA> aidaDynDA(dynamics, geco::CJsonParser());
            aidaDynDA.setUncParam("Bfactor",Bfactor);
            aidaDynDA.setUncParam("SRPC",SRPC);
            xf = rk78DA.propagate(aidaDynDA, x0, t0, tf);

            // Check constant part of final state
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[0])/xfcons[0], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[1])/xfcons[1], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[2])/xfcons[2], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[3])/xfcons[3], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[4])/xfcons[4], 1e-8 );
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, DACE::cons(xf[5])/xfcons[5], 1e-8 );
        }
    };
}
