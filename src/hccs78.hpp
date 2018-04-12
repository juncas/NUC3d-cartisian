//
//  hccs78.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/24.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef hccs78_hpp
#define hccs78_hpp
#include "fieldOperator.h"

namespace nuc3d
{
    class hccs78 : public interoplation
    {
        double ss;
        int p;
        
    public:
        hccs78();
        ~hccs78();
        //WENO5-JS functions
        
        void interpolationInner(const Field &,
                                const int,
                                Field &,
                                const int);
        
        void interpolationBoundaryL(const Field &,
                                    const Field &,
                                    const int,
                                    Field &,
                                    const int);
        
        void interpolationBoundaryR(const Field &,
                                    const Field &,
                                    const int,
                                    Field &,
                                    const int);
    private:
        void hccs78p(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void hccs78n(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void hccs78pBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccs78nBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccs78pBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccs78nBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}

#endif /* hccs78_hpp */
