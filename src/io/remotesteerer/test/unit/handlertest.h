#include <cxxtest/TestSuite.h>
#include <libgeodecomp/io/remotesteerer/handler.h>
#include <libgeodecomp/io/remotesteerer/pipe.h>

using namespace LibGeoDecomp;
using namespace LibGeoDecomp::RemoteSteererHelpers;

namespace LibGeoDecomp {

class HandlerTest : public CxxTest::TestSuite
{
public:
    class MockHandler : public Handler<TestCell<2> >
    {
    public:
        using Handler<TestCell<2> >::GridType;

        MockHandler() :
            RemoteSteererHelpers::Handler<TestCell<2> >("mock")
        {}

        virtual bool operator()(const StringOps::StringVec& parameters, Pipe& pipe, GridType *grid, const Region<Topology::DIM>& validRegion, unsigned step)
        {
            grid->at(Coord<2>(1, 1)).testValue = 4711;
            pipe.addSteeringFeedback("MockHandler mocks you! " + parameters[0]);
            return true;
        }
    };

    void testBasic()
    {
        Pipe pipe;
        MockHandler handler;

        TS_ASSERT_EQUALS("mock", handler.key());
        StringOps::StringVec parameters;
        parameters << "arrrr"
                   << "matey";
        Grid<TestCell<2> > grid(Coord<2>(10, 5));
        Region<2> region;
        region << grid.boundingBox();

        grid[Coord<2>(1, 1)].testValue = -1;
        TS_ASSERT_EQUALS(grid[Coord<2>(1, 1)].testValue, -1);

        bool res = handler(parameters, pipe, &grid, region, 123);
        StringOps::StringVec feedback = pipe.retrieveSteeringFeedback();
        TS_ASSERT_EQUALS(feedback.size(), 1);
        TS_ASSERT_EQUALS(feedback[0], "MockHandler mocks you! arrrr");
        TS_ASSERT_EQUALS(grid[Coord<2>(1, 1)].testValue, 4711.0);
        TS_ASSERT(res);
    }
};

}
