//
//  block.hpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//
#ifndef block_hpp
#define block_hpp
#include <memory>
#include "field.h"
#include "PDEData3d.hpp"
#include "bufferData.hpp"
#include "gradvector.hpp"

namespace nuc3d
{
    class PDEData3d;
    class EulerData3D;
    class bufferData;
    class physicsModel;
    class MPIComunicator3d_nonblocking;
    class IOController;
    class postproc;
    class fft;
    
    class block
    {
        
    protected:
        int nx;
        int ny;
        int nz;

        double xl;
        double yl;
        double zl;

        double x_origin;
        double y_origin;
        double z_origin;

        double dx;
        double dy;
        double dz;


        int nx_global, ny_global, nz_global;

        int npx,npy,npz;

        int bfsize;

        VectorField xyz_center;
        VectorField xyz;
        
        VectorField OutPutValue_prim;
        VectorField OutPutValue_acoust;

        VectorField fft_re;
        VectorField fft_im;
        
        PDEData3d myPDE;

        /*
         this shared pointer EulerData3D could be:
         - EulerData3D;
         - EulerReactiveData3D;(not finished)
         - NaiverStokesData3d;
         - NaiverStokesReactiveData3d;(not finished)
         */
        std::shared_ptr<EulerData3D> myFluxes;

        std::shared_ptr<postproc> myPost;
        std::shared_ptr<fft> myFFT;
        
        VectorBuffer mybuffer;
        
        double time;
        double dt;
        int istep;
        double RES;
        double wall_time;
        
    public:
        block();
        ~block();
        void initial(fieldOperator3d &,
                     physicsModel &,
                     MPIComunicator3d_nonblocking &,
                     IOController &);
        
        void solve(fieldOperator3d &,
                   physicsModel &,
                   MPIComunicator3d_nonblocking &,
                   IOController &);
        
        void printStatus();
        
        void Post(fieldOperator3d &,
                  physicsModel &,
                  MPIComunicator3d_nonblocking &,
                  IOController &);
        
        void Output(fieldOperator3d &,
                    physicsModel &,
                    MPIComunicator3d_nonblocking &,
                    IOController &);
        int getStep(){return istep;};
        double getTime(){return time;};
        
        
    private:
               
        typedef void (nuc3d::block::*pInitial)(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI);
        
        pInitial myInitial[4]={
            &nuc3d::block::initial_default,
            &nuc3d::block::initial_ivc,
            &nuc3d::block::initial_taylorgreen,
            &nuc3d::block::initial_hit
        };
        
        const std::vector<std::string> VarNameList={
            "rho","u","v",
            "w","pressure","temperature",
            "e","sound_speed"
        };

        void initialData(physicsModel &,fieldOperator3d &);
        
        void getXYZ(int coord_x,int coord_y, int coord_z);
        
        void initialQ(IOController &myIO,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI);

        void initial_default(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI);
        
        void initial_ivc(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI);
        
        void initial_taylorgreen(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI);

        void initial_hit(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI);
        
               
        void outputQ_tecplot(int,physicsModel&);
        void outputGEO_tecplot(int myID);
        void outputQ_binary(int,physicsModel&);
        void inputQ_binary(int,int);
        
    };
    
}
#endif /* block_hpp */
