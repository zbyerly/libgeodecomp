#include <libgeodecomp/loadbalancer/mockbalancer.h>

namespace LibGeoDecomp {

MockBalancer::MockBalancer()
{
    events = "";
}

MockBalancer::~MockBalancer()
{
    events += "deleted\n";
}

MockBalancer::WeightVec MockBalancer::balance(
    const MockBalancer::WeightVec& weights,
    const MockBalancer::LoadVec& relativeLoads)
{
    events += "balance() " + toString(weights) + " " +
        toString(relativeLoads) + "\n";
    return weights;
}

std::string MockBalancer::events;

}
