#ifndef LIBGEODECOMP_MISC_CHRONOMETER_H
#define LIBGEODECOMP_MISC_CHRONOMETER_H

#include <libgeodecomp/misc/scopedtimer.h>
#include <libgeodecomp/storage/fixedarray.h>

#include <iomanip>

namespace LibGeoDecomp {

namespace ChronometerHelpers {

/**
 * This class is a tool for counting the number of events and
 * converting their IDs to strings.
 */
template<int ID>
class EventUtil;

class BasicTimerImplementation
{
public:
    template<typename CHRONOMETER>
    BasicTimerImplementation(CHRONOMETER *chrono) :
        totalTimes(chrono->rawTotalTimes()),
        t(ScopedTimer::time())
    {}

    template<typename CHRONOMETER>
    BasicTimerImplementation(CHRONOMETER *chrono, double t) :
        totalTimes(chrono->rawTotalTimes()),
        t(t)
    {}

protected:
    double *totalTimes;

    double elapsed() const
    {
        return ScopedTimer::time() - t;
    }

    double t;
};

template<int ID>
class EventUtil
{
public:
    static const std::size_t EVENT_COUNTER = EventUtil<ID - 1>::EVENT_COUNTER;
    // NUM_EVENTS is only declared here, but not in the template
    // specializations generated by DEFINE_EVNT() below. This is to
    // prevent users (e.g. class Chronometer) from accidentially
    // reading wrong values as all but the last specialization cannot
    // know the total number of events.
    static const std::size_t NUM_EVENTS = EventUtil<ID - 1>::EVENT_COUNTER;

    std::string operator()()
    {
        return "unknown event";
    }
};

#define DEFINE_EVENT(CLASS_NAME, PARENT, EVENT_NAME, EVENT_ID)      \
    namespace ChronometerHelpers {                                  \
                                                                    \
    template<>                                                      \
    class EventUtil<EVENT_ID>                                       \
    {                                                               \
    public:                                                         \
        static const std::size_t EVENT_COUNTER = EVENT_ID + 1;           \
                                                                    \
        std::string operator()()                                    \
        {                                                           \
            return EVENT_NAME;                                      \
        }                                                           \
    };                                                              \
                                                                    \
    }                                                               \
                                                                    \
    class CLASS_NAME ## Implementation :                            \
        protected PARENT ## Implementation                          \
    {                                                               \
    public:                                                         \
        static const int ID = EVENT_ID;                             \
                                                                    \
        template<typename CHRONOMETER>                              \
        CLASS_NAME ## Implementation(CHRONOMETER *chrono) :         \
            PARENT ## Implementation(chrono)                        \
        {}                                                          \
                                                                    \
        template<typename CHRONOMETER>                              \
        CLASS_NAME ## Implementation(CHRONOMETER *chrono, double time) : \
            PARENT ## Implementation(chrono, time)                  \
        {}                                                          \
                                                                    \
        ~CLASS_NAME ## Implementation()                             \
        {                                                           \
            totalTimes[ID] += t;                                    \
        }                                                           \
    };                                                              \
                                                                    \
    class CLASS_NAME : protected CLASS_NAME ## Implementation       \
    {                                                               \
    public:                                                         \
        typedef CLASS_NAME ## Implementation Implementation;        \
        static const int ID = EVENT_ID;                             \
                                                                    \
        template<typename CHRONOMETER>                              \
            CLASS_NAME(CHRONOMETER *chrono) :                       \
            CLASS_NAME ## Implementation(chrono)                    \
        {}                                                          \
                                                                    \
        ~CLASS_NAME()                                               \
        {                                                           \
            t = elapsed();                                          \
        }                                                           \
    };
}

DEFINE_EVENT(TimeTotal,          ChronometerHelpers::BasicTimer,   "total_time",           0)
DEFINE_EVENT(TimeCompute,        ChronometerHelpers::BasicTimer,   "compute_time",         1)
DEFINE_EVENT(TimeComputeInner,   TimeCompute,                      "compute_time_inner",   2)
DEFINE_EVENT(TimeComputeGhost,   TimeCompute,                      "compute_time_ghost",   3)
DEFINE_EVENT(TimePatchAccepters, ChronometerHelpers::BasicTimer,   "patch_accepters_time", 4)
DEFINE_EVENT(TimePatchProviders, ChronometerHelpers::BasicTimer,   "patch_providers_time", 5)
DEFINE_EVENT(TimeCommunication,  ChronometerHelpers::BasicTimer,   "communication_time",   6)
DEFINE_EVENT(TimeInput,          ChronometerHelpers::BasicTimer,   "input_time",           7)
DEFINE_EVENT(TimeOutput,         ChronometerHelpers::BasicTimer,   "output_time",          8)

namespace ChronometerHelpers {

class EventToString
{
public:
    std::string operator()(int id)
    {
        switch(id) {
        case 0:
            return EventUtil<0>()();
        case 1:
            return EventUtil<1>()();
        case 2:
            return EventUtil<2>()();
        case 3:
            return EventUtil<3>()();
        case 4:
            return EventUtil<4>()();
        case 5:
            return EventUtil<5>()();
        case 6:
            return EventUtil<6>()();
        case 7:
            return EventUtil<7>()();
        case 8:
            return EventUtil<8>()();
        case 9:
            return EventUtil<9>()();
        case 10:
            return EventUtil<10>()();
        case 11:
            return EventUtil<11>()();
        case 12:
            return EventUtil<12>()();
        case 13:
            return EventUtil<13>()();
        case 14:
            return EventUtil<14>()();
        case 15:
            return EventUtil<15>()();
        case 16:
            return EventUtil<16>()();
        case 17:
            return EventUtil<17>()();
        case 18:
            return EventUtil<18>()();
        case 19:
            return EventUtil<19>()();
        default:
            throw std::invalid_argument("event not listed in EventToString");
        }
    }
};

}

/**
 * This class can be used to measure execution time of different parts
 * of our code. This is useful to determine the relative load of a
 * node or to find out which part of the algorithm the most time.
 */
class Chronometer
{
public:
    friend class ChronometerTest;
    friend class Serialization;
    friend class Typemaps;

    // measure one time interval per class of events
    static const std::size_t NUM_INTERVALS = ChronometerHelpers::EventUtil<100>::NUM_EVENTS;

    Chronometer() :
        totalTimes(NUM_INTERVALS, 0)
    {
        reset();
    }

    inline const double& operator[](int intervalID) const
    {
        return totalTimes[intervalID];
    }

    inline double& operator[](int intervalID)
    {
        return totalTimes[intervalID];
    }

    /**
     * Aggregates two chronometers. This is good for accumulating
     * measurements from multiple sources.
     */
    Chronometer& operator+=(const Chronometer& other)
    {
        for (std::size_t i = 0; i < NUM_INTERVALS; ++i) {
            totalTimes[i] += other.totalTimes[i];
        }

        return *this;
    }

    Chronometer operator+(const Chronometer& other) const
    {
        Chronometer ret(*this);
        ret += other;
        return ret;
    }

    /**
     * Flushes all time totals to 0.
     */
    void reset()
    {
        std::fill(totalTimes.begin(), totalTimes.end(), 0);
    }

    void cycle()
    {

    }

    double interval(int i) const
    {
        return totalTimes[i];
    }

    template<typename INTERVAL>
    double interval() const
    {
        return totalTimes[INTERVAL::ID];
    }

    /**
     * Returns the ratio of the accumulated times of INTERVAL1 and
     * INTERVAL2.
     *
     * Caveat: returns 0.5 if the second interval is empty. This is
     * because ratio() is typically used for load balancing. For that
     * purpose 0.5 is easier to digest than NaN.
     */
    template<typename INTERVAL1, typename INTERVAL2>
    double ratio()
    {
        const int i1 = INTERVAL1::ID;
        const int i2 = INTERVAL2::ID;

        if (totalTimes[i2] == 0) {
            return 0.5;
        }
        return totalTimes[i1] / totalTimes[i2];
    }

    double *rawTotalTimes()
    {
        return totalTimes.begin();
    }

    template<typename EVENT>
    void addTime(double elapsedTime)
    {
        typename EVENT::Implementation(this, elapsedTime);
    }

    template<typename EVENT>
    void tock(double startTime)
    {
        addTime<EVENT>(ScopedTimer::time() - startTime);
    }

    std::string report()
    {
        std::stringstream buf;

        for (std::size_t i = 0; i < NUM_INTERVALS; ++i) {
            buf << std::left << std::setw(20) << ChronometerHelpers::EventToString()(i)
                << ": " << totalTimes[i] << "s\n";
        }

        return buf.str();
    }

private:
    FixedArray<double, Chronometer::NUM_INTERVALS> totalTimes;
};

}

#endif
