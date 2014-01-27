#ifndef LIBGEODECOMP_GEOMETRY_PARTITIONS_SCOTCHPARTITION_H
#define LIBGEODECOMP_GEOMETRY_PARTITIONS_SCOTCHPARTITION_H

#include <libgeodecomp/config.h>

#ifdef LIBGEODECOMP_WITH_SCOTCH

#include <libgeodecomp/geometry/partitions/partition.h>

#include <scotch.h>

namespace LibGeoDecomp {

template<int DIM>
class ScotchPartition : public Partition<DIM>
{
public:
    using Partition<DIM>::startOffsets;
    using Partition<DIM>::weights;

    inline ScotchPartition(
        const Coord<DIM>& origin = Coord<DIM>(),
        const Coord<DIM>& dimensions = Coord<DIM>(),
        const long& offset = 0,
        const std::vector<std::size_t>& weights = std::vector<std::size_t>(2)) :
        Partition<DIM>(offset, weights),
        origin(origin),
        dimensions(dimensions),
        cellNbr(dimensions[0] * dimensions[1])
        {
            indices = new SCOTCH_Num[cellNbr];
            regions = new Region<DIM>[weights.size()];
            initIndices();
            createRegions();
        }

    Region<DIM> getRegion(const std::size_t node) const
    {
        return regions[node];
    }

private:
    Coord<DIM> origin;
    Coord<DIM> dimensions;
    SCOTCH_Num * indices;
    SCOTCH_Num const cellNbr;
    Region<DIM> * regions;

    std::vector<std::pair <int,int> > * boxes;

    void initIndices(){
        SCOTCH_Arch arch;
        SCOTCH_archInit(&arch);
        SCOTCH_Num * velotabArch;
        SCOTCH_Num const vertnbrArch = weights.size();
        velotabArch = new SCOTCH_Num[weights.size()];
        for(int i = 0;i<vertnbrArch;++i){
            velotabArch[i] = weights[i];
        }
        SCOTCH_archCmpltw (&arch,vertnbrArch,velotabArch);

        SCOTCH_Strat * straptr = SCOTCH_stratAlloc();;
        SCOTCH_stratInit(straptr);
        SCOTCH_stratGraphMapBuild(straptr,SCOTCH_STRATRECURSIVE,vertnbrArch,0);

        SCOTCH_Graph grafdat;
        SCOTCH_graphInit (&grafdat);

        SCOTCH_Num const edgenbrGra = 2 * (dimensions[0] * (dimensions[1] - 1) +
                                           (dimensions[0] - 1) * dimensions[1]);

        SCOTCH_Num * verttabGra;
        SCOTCH_Num * edgetabGra;
        verttabGra = new SCOTCH_Num[cellNbr + 1];
        edgetabGra = new SCOTCH_Num[edgenbrGra];

        int pointer = 0;
        for(int i = 0;i < cellNbr;++i){
            verttabGra[i] = pointer;
            if(i%dimensions[0] != 0){
                edgetabGra[pointer] = i-1;
                pointer++;
            }
            if(i%dimensions[0] != (dimensions[0]-1)){
                edgetabGra[pointer] = i+1;
                pointer++;
            }
            if(!(i < dimensions[0])){
                edgetabGra[pointer] = i-dimensions[0];
                pointer++;
            }
            if(!(i >= dimensions[0] * (dimensions[1]-1))){
                edgetabGra[pointer] = i+dimensions[0];
                pointer++;
            }
        }
        verttabGra[cellNbr] = pointer;


        SCOTCH_graphBuild(
                          &grafdat,
                          0,
                          cellNbr,
                          verttabGra,
                          verttabGra +1,
                          NULL,
                          NULL,
                          edgenbrGra,
                          edgetabGra,
                          NULL
                          );

        SCOTCH_graphMap (&grafdat,&arch,straptr,indices);
    }

    void createRegions(){
        int rank = indices[0];
        int length = 0;
        int start = 0;
        for(int i = 1;i < cellNbr+1;++i){
            if(rank == indices[i] && i < cellNbr && i%dimensions[0]!=0){
                length++;
            } else {
                length++;
                regions[rank] << CoordBox<DIM>(Coord<DIM>(start%dimensions[0],
                                                          start/dimensions[0]),
                                               Coord<DIM>(length,1));
                rank = indices[i];
                start = i;
                length = 0;
            }
        }
    }

 };
}

#endif

#endif
