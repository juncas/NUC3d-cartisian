//
//  hccs77.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/24.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#ifndef hccs77_hpp
#define hccs77_hpp
#include "fieldOperator.h"

namespace nuc3d
{
    class hccs77 : public interoplation
    {
        double ss;
        int p;
        
    public:
        hccs77();
        ~hccs77();
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
        void hccs77p(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void hccs77n(const Field & fieldIN,
                      Field & fieldOUT,
                      const int tilesize);
        
        void hccs77pBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccs77nBL(const Field & fieldIN,
                        const Field &boundaryL,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccs77pBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        void hccs77nBR(const Field & fieldIN,
                        const Field &boundaryR,
                        Field & fieldOUT,
                        const int tilesize);
        
        
    };
}

#endif /* hccs77_hpp */
