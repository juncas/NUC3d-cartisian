#ifndef centraldifference2nd_h
#define centraldifference2nd_h
#include "fieldOperator.h"


namespace nuc3d
{
    class centraldifference2nd : public differential
    {
    public:
        centraldifference2nd();
        ~centraldifference2nd();
        // 6th order central difference functions
        
        void differentialInner(const Field &,
        const double &,
                               Field &,
                               const int tilesize);
        void differentialBoundaryL(const Field &,
        const double &,
                                   const Field &,
                                   Field &,
                                   const int tilesize);
        
        void differentialBoundaryR(const Field &,
        const double &,
                                   const Field &,
                                   Field &,
                                   const int tilesize);
        
    };
}



#endif