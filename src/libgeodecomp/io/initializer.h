#ifndef LIBGEODECOMP_IO_INITIALIZER_H
#define LIBGEODECOMP_IO_INITIALIZER_H

#include <libgeodecomp/config.h>
#include <libgeodecomp/geometry/adjacencymanufacturer.h>
#include <libgeodecomp/misc/apitraits.h>
#include <libgeodecomp/misc/random.h>
#include <libgeodecomp/storage/gridbase.h>
#include <libgeodecomp/geometry/regionbasedadjacency.h>
#include <stdexcept>

namespace LibGeoDecomp {

/**
 * The initializer sets up the initial state of the grid. For this a
 * Simulator will invoke Initializer::grid(). Keep in mind that grid()
 * might be called multiple times and that for parallel runs each
 * Initializer will be responsible just for a sub-cuboid of the whole
 * grid.
 */
template<typename CELL>
class Initializer : public AdjacencyManufacturer<APITraits::SelectTopology<CELL>::Value::DIM>
{
public:
    typedef typename APITraits::SelectTopology<CELL>::Value Topology;
    typedef CELL Cell;
    typedef typename SharedPtr<Adjacency>::Type AdjacencyPtr;

    static const unsigned NANO_STEPS = APITraits::SelectNanoSteps<CELL>::VALUE;
    static const int DIM = Topology::DIM;

    virtual ~Initializer()
    {}

    /**
     * initializes all cells of the grid at target
     */
    virtual void grid(GridBase<CELL, DIM> *target) = 0;

    /**
     * Allows a Simulator to discover the extent of the whole
     * simulation. Usually Simulations will use 0 as the origin, but
     * there is no oblication to do so.
     */
    virtual CoordBox<DIM> gridBox()
    {
        return CoordBox<DIM>(Coord<DIM>(), gridDimensions());
    }

    /**
     * returns the size of the gridBox().
     */
    virtual Coord<DIM> gridDimensions() const = 0;

    /**
     * yields the logical time step at which the simulation should start
     */
    virtual unsigned startStep() const = 0;

    /**
     * gives the time step at which the simulation should terminate.
     *
     * Example: if startStep is 0 and maxSteps is 10, then the
     * Simulator should start at t=0, update to t=1, update to t=2,
     * ... until it has updated to t=10.
     *
     * If startStep is 5 and maxSteps is 8, then the Simulator is
     * expected to init at t=5, update to t=6, update to t=7, and
     * finally update to t=8.
     */
    virtual unsigned maxSteps() const = 0;

    Coord<DIM> normalize(const Coord<DIM>& coord) const
    {
        return Topology::normalize(coord, gridDimensions());
    }

    /**
     * This seeds the pseudo random number generator (see class
     * Random) with a seed unique to the coordinate. This is handy if
     * an Initializer needs to randomize parameters, but requires
     * these random values to be identical if multiple processes are
     * initializing the same logical grid cell.
     */
    void seedRNG(const Coord<DIM>& coord) const
    {
        // normalization is mandatory as coordinates may lie outside
        // of the grid's bounding box if a model is using periodic
        // boundary conditions (e.g. the torus topology). example: a
        // node is assigned a rectangular region from (0,0) to
        // (50,50). with periodic boundary conditions, its left outer
        // ghost zone would be on the x-coordinate -1, ranging from -1
        // to 51 on the y-axis. the value -1 would then be wrapped
        // around thanks to the torus topology, to align with the
        // grid's dimension.
        std::size_t index = normalize(coord).toIndex(gridDimensions());
        Random::seed(index);
    }

    /**
     * Returns the list of neighboring nodes for at least all IDs
     * stored in the Region. Given a node n_1 \in region, a node n_2
     * is assumed to be a neighbor of n_1, iff there is a directed
     * edge (n_1, n_2) in the adjacency list of the unstructured grid.
     */
    AdjacencyPtr getAdjacency(const Region<DIM>& /* region */) const
    {
        checkTopologyIfAdjacencyIsNeeded(Topology());
        return AdjacencyPtr();
    }

    std::vector<std::size_t> getWeights(const Region<1>& region) const
    {
        std::vector<std::size_t> weights;
	for (Region<1>::Iterator i = region.begin(); i != region.end(); ++i) {
            std::size_t weight = 0;
            weights.push_back(weight);
        }
        return weights;
    }

private:
    template<typename TOPOLOGY>
    void checkTopologyIfAdjacencyIsNeeded(const TOPOLOGY /* unused */) const
    {
        // intentionally left blank
    }

    void checkTopologyIfAdjacencyIsNeeded(const Topologies::Unstructured::Topology /* unused */) const
    {
        std::string message(
            "You are using an unstructured topology but the Initializer did not re-implement getAdjacency()."
            "You really need to implement getAdjacency() otherwise ghost zone (halo) synchronization will not work.");

        LOG(FATAL, message);
        throw std::logic_error(message);
    }

    /**
     * For all nodes n_1 in the Region, return a list of neighbor
     * nodes n_2 where (n_2, n_1) is in the edge list of the
     * unstructured grid.
     *
     * Attention: the direction of the edges in the returned adjacency
     * is the opposite of the direction in the original unstructured
     * grid. See UnstructuredTestInitializer for an example.
     */
    SharedPtr<Adjacency>::Type getReverseAdjacency(const Region<DIM>& /* region */) const
    {
        return makeShared<Adjacency>(new RegionBasedAdjacency());
    }
};

}

#endif
