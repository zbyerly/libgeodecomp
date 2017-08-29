#ifndef LIBGEODECOMP_GEOMETRY_PARTITIONS_PTSCOTCHUNSTRUCTUREDPARTITION_H
#define LIBGEODECOMP_GEOMETRY_PARTITIONS_PTSCOTCHUNSTRUCTUREDPARTITION_H

#include <libgeodecomp/config.h>
#include <libgeodecomp/geometry/partitions/partition.h>
#include <libgeodecomp/geometry/adjacency.h>

#ifdef LIBGEODECOMP_WITH_CPP14
#ifdef LIBGEODECOMP_WITH_SCOTCH

#include <ptscotch.h>

#include <chrono>

namespace LibGeoDecomp {

/**
 * This class is like DistributedPTScotchUnstructuredPartition, but
 * relies on PT-SCOTCH's serial decomposition. This limits
 * scalability, but decouples the code from MPI.
 */
template<int DIM>
class PTScotchUnstructuredPartition : public Partition<DIM>
{
public:
    using Partition<DIM>::startOffsets;
    using Partition<DIM>::weights;
    using typename Partition<DIM>::AdjacencyPtr;

    PTScotchUnstructuredPartition(
            const Coord<DIM>& origin,
            const Coord<DIM>& dimensions,
            const long offset,
            const std::vector<std::size_t>& weights,
            const AdjacencyPtr& adjacency,
            const std::vector<double>& cellWeights) :
        Partition<DIM>(offset, weights),
        adjacency(adjacency),
        origin(origin),
        dimensions(dimensions),
        numCells(dimensions.prod()),
        cellWeights(cellWeights)
    {
        buildRegions();
    }

    Region<DIM> getRegion(const std::size_t node) const override
    {
        return regions.at(node);
    }

private:
    SharedPtr<Adjacency>::Type adjacency;
    Coord<DIM> origin;
    Coord<DIM> dimensions;
    SCOTCH_Num numCells;
    std::vector<Region<DIM> > regions;
    std::vector<double> cellWeights;

    void buildRegions()
    {
        std::vector<SCOTCH_Num> indices(numCells);
        initIndices(indices);

        regions.resize(weights.size());
        createRegions(indices);
    }

    int getCellWeight(int id)
    {
        /* return an integer between 1 and 100 */
        int cellWeight = 100*cellWeights[id];
        std::cout << "cellWeight = " << cellWeight << std::endl;
        return cellWeight;
    }

    void initIndices(std::vector<SCOTCH_Num>& indices)
    {
        // create 2D grid
        SCOTCH_Graph graph;
        int error = SCOTCH_graphInit(&graph);

        SCOTCH_Num numEdges = this->adjacency->size();

        std::cout << "numCellWeights = " << this->cellWeights.size() << std::endl;
        std::cout << "numCells = " << numCells << std::endl;

        std::vector<SCOTCH_Num> verttabGra;
        std::vector<SCOTCH_Num> edgetabGra;
        std::vector<SCOTCH_Num> velotabGra;

        verttabGra.reserve(numCells + 1);
        edgetabGra.reserve(numEdges);

        int currentEdge = 0;
        for (int i = 0; i < numCells; ++i) {
            verttabGra.push_back(currentEdge);

            std::size_t before = edgetabGra.size();
            this->adjacency->getNeighbors(i, &edgetabGra);
            velotabGra.push_back(this->getCellWeight(i));
            currentEdge += edgetabGra.size() - before;
        }

        verttabGra.push_back(currentEdge);

        error = SCOTCH_graphBuild(
                &graph,               /* grafptr Graph structure to fill             */
                0,                    /* baseval Base value                          */
                numCells,             /* vertnbr Number of vertices                  */
                &verttabGra[0],       /* verttab Vertex array [vertnbr or vertnbr+1] */
                nullptr,              /* vendtab Vertex end array [vertnbr]          */
                nullptr,              /* velotab Vertex load array                   */
                nullptr,              /* vlbltab Vertex label array                  */
                numEdges,             /* edgenbr Number of edges (arcs)              */
                &edgetabGra[0],       /* edgetab Edge array [edgenbr]                */
                nullptr);             /* edlotab Edge load array                     */
        if (error) {
            LOG(FAULT, "SCOTCH_graphBuild error: " << error);
        }

        error = SCOTCH_graphCheck(&graph);
        if (error) {
            LOG(FAULT, "SCOTCH_graphCheck error: " << error);
        }

        SCOTCH_Strat *straptr = SCOTCH_stratAlloc();
        error = SCOTCH_stratInit(straptr);
        if (error) {
            LOG(FAULT, "SCOTCH_stratInit error: " << error);
        }

        error = SCOTCH_graphPart(&graph, weights.size(), straptr,& indices[0]);
        if (error) {
            LOG(FAULT, "SCOTCH_graphMap error: " << error);
        }

        SCOTCH_graphExit(&graph);
        SCOTCH_stratExit(straptr);
    }

    void createRegions(const std::vector<SCOTCH_Num>& indices)
    {
        for (int i = 0; i < numCells; ++i) {
            regions[indices[i]] << Coord<1>(i);
        }
    }

};

}

#endif
#endif

#endif
