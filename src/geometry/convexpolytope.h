#ifndef LIBGEODECOMP_GEOMETRY_CONVEXPOLYHEDRON_H
#define LIBGEODECOMP_GEOMETRY_CONVEXPOLYHEDRON_H

#include <libgeodecomp/geometry/floatcoord.h>
#include <libgeodecomp/geometry/plane.h>
#include <libgeodecomp/misc/random.h>

namespace LibGeoDecomp {

/**
 * ConvexPolytope is an intersection of half-spaces. On the 2D plane
 * this is a convex polygone, in 3D space a convex polyhedron.
 */
template<typename COORD, typename ID = int>
class ConvexPolytope
{
public:
    const static std::size_t SAMPLES = 1000;

    typedef Plane<COORD, ID> EquationType;

    ConvexPolytope(
        const COORD& center,
        const COORD& quadrantSize,
        const COORD& simSpaceDim,
        const double minCellDistance,
        ID id) :
        center(center),
        quadrantSize(quadrantSize),
        simSpaceDim(simSpaceDim),
        minCellDistance(minCellDistance),
        id(id),
        area(simSpaceDim.prod()),
        diameter(*std::max_element(&simSpaceDim[0], &simSpaceDim[0] + 2))
    {
        limits << EquationType(COORD(center[0], 0),              COORD( 0,  1))
               << EquationType(COORD(0, center[1]),              COORD( 1,  0))
               << EquationType(COORD(simSpaceDim[0], center[1]), COORD(-1,  0))
               << EquationType(COORD(center[0], simSpaceDim[1]), COORD( 0, -1));
    }

    ConvexPolytope& operator<<(const EquationType& eq)
    {
        // no need to reinsert if limit already present (would only cause trouble)
        for (typename std::vector<EquationType>::iterator i = limits.begin();
             i != limits.end();
             ++i) {
            if (eq == *i) {
                return *this;
            }
        }

        limits << eq;
        std::vector<COORD > cutPoints = generateCutPoints(limits);


        // we need to set up a kill list to avoid jumping beyond the
        // end of limits.
        std::set<int> deleteSet;

        for (std::size_t i = 0; i < limits.size(); ++i) {
            COORD leftDir = turnLeft90(limits[i].dir);
            int dist1 = (cutPoints[2 * i + 0] - limits[i].base) * leftDir;
            int dist2 = (cutPoints[2 * i + 1] - limits[i].base) * leftDir;
            if (dist2 >= dist1) {
                // twisted differences, deleting
                deleteSet.insert(i);
            }

            for (std::size_t j = 0; j < limits.size(); ++j) {
                if (i != j) {
                    // parallel lines, deleting...
                    if (cutPoint(limits[i], limits[j]) == farAway<2>()) {
                        if (limits[i].dir * limits[j].dir > 0) {
                            double dist1 = (center - limits[i].base) *
                                limits[i].dir;
                            double dist2 = (center - limits[j].base) *
                                limits[i].dir;
                            if (dist2 > dist1) {
                                deleteSet.insert(j);
                            }
                        }
                    }
                }
            }
        }

        std::vector<EquationType> newLimits;
        for (std::size_t i = 0; i < limits.size(); ++i) {
            if (!deleteSet.count(i)) {
                newLimits << limits[i];
            }
        }

        limits = newLimits;

        return *this;
    }

    template<typename POINT>
    ConvexPolytope& operator<<(const POINT& c)
    {
        COORD base = (center + c.center) / 2;
        COORD dir = center - c.center;
        *this << EquationType(base, dir, c.id);
        return *this;
    }

    std::vector<COORD > generateCutPoints(const std::vector<EquationType>& equations) const
    {
        std::vector<COORD > buf(2 * equations.size(), farAway<2>());

        for (std::size_t i = 0; i < equations.size(); ++i) {
            for (std::size_t j = 0; j < equations.size(); ++j) {
                if (i != j) {
                    COORD cut = cutPoint(equations[i], equations[j]);
                    int offset = 2 * i;
                    COORD delta = cut - equations[i].base;
                    COORD turnedDir = turnLeft90(equations[i].dir);
                    double distance =
                        1.0 * delta[0] * turnedDir[0] +
                        1.0 * delta[1] * turnedDir[1];


                    bool isLeftCandidate = true;
                    if (equations[j].dir * turnedDir > 0) {
                        isLeftCandidate = false;
                        offset += 1;
                    }

                    delta = buf[offset] - equations[i].base;
                    double referenceDist =
                        1.0 * delta[0] * turnedDir[0] +
                        1.0 * delta[1] * turnedDir[1];
                    bool flag = false;
                    if (buf[offset] == farAway<2>()) {
                        flag = true;
                    }
                    if (isLeftCandidate  && (distance < referenceDist)) {
                        flag = true;
                    }
                    if (!isLeftCandidate && (distance > referenceDist)) {
                        flag = true;
                    }
                    if (cut == farAway<2>()) {
                        flag = false;
                    }
                    if (flag) {
                        buf[offset] = cut;
                    }
                }
            }
        }

        return buf;
    }

    std::vector<COORD > getShape() const
    {
        std::vector<COORD > cutPoints = generateCutPoints(limits);

        for (std::size_t i = 0; i < cutPoints.size(); ++i) {
            if (cutPoints[i] == farAway<2>()) {
                throw std::logic_error("invalid cut point");
            }
        }

        std::map<double, COORD > points;
        for (typename std::vector<COORD >::iterator i = cutPoints.begin();
             i != cutPoints.end();
             ++i) {
            COORD delta = *i - center;
            double angle = relativeCoordToAngle(delta, cutPoints);

            points[angle] = *i;
        }

        std::vector<COORD > res;
        for (typename std::map<double, COORD >::iterator i = points.begin();
             i != points.end(); ++i) {
            res << i->second;
        }

        if (res.size() < 3) {
            throw std::logic_error("cycle too short");
        }

        return res;
    }

    bool includes(const COORD& c)
    {
        for (std::size_t i = 0; i < limits.size(); ++i) {
            if (!limits[i].isOnTop(c)) {
                return false;
            }
        }
        return true;
    }

    void updateGeometryData()
    {
        std::vector<COORD > cutPoints = generateCutPoints(limits);

        for (std::size_t i = 0; i < limits.size(); ++i) {
            COORD delta = cutPoints[2 * i + 0] - cutPoints[2 * i + 1];
            limits[i].length = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
        }

        COORD min = simSpaceDim;
        COORD max(0, 0);
        for (std::size_t i = 0; i < cutPoints.size(); ++i) {
            COORD& c = cutPoints[i];
            max = c.max(max);
            min = c.min(min);
        }
        COORD delta = max - min;

        int hits = 0;
        for (std::size_t i = 0; i < SAMPLES; ++i) {
            COORD p = COORD(Random::gen_d(delta[0]),
                            Random::gen_d(delta[1])) + min;
            if (includes(p)) {
                ++hits;
            }
        }
        area = 1.0 * hits / SAMPLES * delta.prod();

        double radiusSquared = delta * delta;
        double maxRadiusSquared = quadrantSize * quadrantSize * 0.5;
        if (radiusSquared > maxRadiusSquared) {
            std::cerr << "center: " << center << "\n"
                      << "my diameter: " << diameter << "\n"
                      << "maxRadiusSquared: " << maxRadiusSquared << "\n"
                      << "quadrantSize: " << quadrantSize << "\n"
                      << "cutPoints: " << cutPoints << "\n"
                      << "min: " << min << "\n"
                      << "max: " << max << "\n";

            throw std::logic_error("element too large");
        }

        double newDiameter = *std::max_element(&delta[0], &delta[0] + 2);
        if (newDiameter > diameter) {
            throw std::logic_error("diameter should never ever increase!");

        }

        diameter = newDiameter;
    }

    const COORD& getCenter() const
    {
        return center;
    }

    /**
     * The ConvexPolytope's volume is determined via random sampling
     * integration. Few samples are used, so this method is highly
     * inaccurate.
     */
    const double& getVolume() const
    {
        return area;
    }

    const std::vector<EquationType>& getLimits() const
    {
        return limits;
    }

    const double& getDiameter() const
    {
        return diameter;
    }

private:
    COORD center;
    FloatCoord<2> quadrantSize;
    FloatCoord<2> simSpaceDim;
    double minCellDistance;
    ID id;
    double area;
    double diameter;
    std::vector<EquationType> limits;

    template<int DIM>
    static Coord<DIM> farAway()
    {
        return Coord<DIM>::diagonal(-1);
    }

    double relativeCoordToAngle(const COORD& delta, const std::vector<COORD >& cutPoints) const
    {
        double length = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);

        if (length > 0) {
            double dY = delta[1] / length;
            double angle = asin(dY);

            if (delta[0] < 0) {
                angle = M_PI - angle;
            }

            return angle;
        }

        // If lengths is 0, then we can't deduce the angle
        // from the cutPoint's location. But we know that the
        // center is located on the simulation space's
        // boundary. Hence at least one of the following four
        // cases is true, which we can use to assign a fake
        // angle to the point:
        //
        // 1. all x-coordinates are >= than center[0]
        // 2. all x-coordinates are <= than center[0]
        // 3. all y-coordinates are >= than center[1]
        // 4. all y-coordinates are <= than center[1]

        bool case1 = true;
        bool case2 = true;
        bool case3 = true;
        bool case4 = true;

        for (typename std::vector<COORD >::const_iterator i = cutPoints.begin();
             i != cutPoints.end();
             ++i) {
            if ((*i)[0] < center[0]) {
                case1 = false;
            }
            if ((*i)[0] > center[0]) {
                case2 = false;
            }
            if ((*i)[1] < center[1]) {
                case3 = false;
            }
            if ((*i)[1] > center[1]) {
                case4 = false;
            }
        }

        if (case1) {
            return 1.0 * M_PI;
        }
        if (case2) {
            return 0.0 * M_PI;
        }
        if (case3) {
            return 1.5 * M_PI;
        }
        if (case4) {
            return 0.5 * M_PI;
        }

        throw std::logic_error("oops, boundary case in boundary generation should be logically impossible!");
    }

    static COORD turnLeft90(const COORD& c)
    {
        return COORD(-c[1], c[0]);
    }

    COORD cutPoint(EquationType eq1, EquationType eq2) const
    {
        if (eq1.dir[1] == 0) {
            if (eq2.dir[1] == 0) {
                // throw std::invalid_argument("both lines are vertical")
                return farAway<2>();
            }
            std::swap(eq1, eq2);
        }

        COORD dir1 = turnLeft90(eq1.dir);
        double m1 = 1.0 * dir1[1] / dir1[0];
        double d1 = eq1.base[1] - m1 * eq1.base[0];

        if (eq2.dir[1] == 0) {
            return COORD(eq2.base[0], eq2.base[0] * m1 + d1);
        }

        COORD dir2 = turnLeft90(eq2.dir);
        double m2 = 1.0 * dir2[1] / dir2[0];
        double d2 = eq2.base[1] - m2 * eq2.base[0];

        if (m1 == m2) {
            // throw std::invalid_argument("parallel lines")
            return farAway<2>();
        }

        double x = (d2 - d1) / (m1 - m2);
        double y = d1 + x * m1;

        if ((x < (-10 * simSpaceDim[0])) ||
            (x > ( 10 * simSpaceDim[0])) ||
            (y < (-10 * simSpaceDim[1])) ||
            (y > ( 10 * simSpaceDim[1]))) {
            return farAway<2>();
        }

        return COORD(x, y);
    }
};

}

#endif

