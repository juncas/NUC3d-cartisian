#ifndef Euler3d_h
#define Euler3d_h
#include "fieldOperator.h"
#include "MPICommunicator.h"
#include "bufferData.hpp"

namespace nuc3d
{
    class EulerData3D;
    class EulerFlux;
    class PDEData3d;
    class physicsModel;
    class bufferData;
    class fft;
    
    class EulerFlux
    {
        friend class physicsModel;
        friend class EulerData3D;
    public:
        VectorField FluxL;
        VectorField FluxR;
        VectorField reconstFluxL;
        VectorField reconstFluxR;
        VectorField reconstFlux;
        double maxEigen;
    public:
        EulerFlux(int,int,int,int);
        ~EulerFlux();
        
    public:
        void combineFluxLR();
    };
    
    class EulerData3D
    {
        friend class physicsModel;
    public:

        VectorField W_Euler;//primitive values (rho,u,v,w,p,...)
        
        VectorField W0_Euler;//acoustics values (T,e,alpha)
        
        //Intermediate fields for Reimann problem
        EulerFlux Flux_xi;
        EulerFlux Flux_eta;
        EulerFlux Flux_zeta;
        
        //invicid derivatives
        VectorField dfdxi;
        VectorField dgdeta;
        VectorField dhdzeta;

        VectorField source; // source term e.g. Forcing, Coriolis force etc.
        
        double dt;
    public:
        int nx;
        int ny;
        int nz;
        double dx;
        double dy;
        double dz;
        EulerData3D( int, int , int , int, double, double ,double  );
        ~EulerData3D();
        
        //used by field operators and Riemann solvers
        EulerFlux& getFluxXi();
        EulerFlux& getFluxEta();
        EulerFlux& getFluxZeta();
        
        //used by field operators and PDE
        virtual VectorField& getDrivativeXi();
        virtual VectorField& getDrivativeEta();
        virtual VectorField& getDrivativeZeta();
        
        //used by Riemann solvers
        VectorField& getPrimatives();
        VectorField& getAcoustics();

        //used by source solver
        VectorField& getSource();

        virtual void solve(PDEData3d &,
         fieldOperator3d &,
         VectorBuffer &,
         physicsModel &,
         MPIComunicator3d_nonblocking &,
         fft &myFFT,
         const VectorField &xyz_center);

    protected:
        void solveRiemann(PDEData3d &,
          physicsModel &);
        
        void solveCon2Prim(PDEData3d &myPDE,
         physicsModel &myModel);
        
        void solveInv(fieldOperator3d &,
          VectorBuffer &,
          MPIComunicator3d_nonblocking &);
        
        void setBoundaryCondition(PDEData3d &,
          physicsModel &myModel,
          std::vector<bufferData> &);
        
        void solveSource(PDEData3d &myPDE,
            physicsModel &myModel,
            MPIComunicator3d_nonblocking &,
            fft &myFFT,
            const VectorField &xyz_center);
        
        void getDt();
        
        void solveRHS(PDEData3d &);
    private:

        void solveInvicidFluxL(EulerFlux &,
         fieldOperator3d &myOP,
         VectorBuffer &,
         MPIComunicator3d_nonblocking &,
         int );
        
        void solveInvicidFluxR(EulerFlux &,
         fieldOperator3d &myOP,
         VectorBuffer &,
         MPIComunicator3d_nonblocking &,
         int );
        
        void setDerivativesInv();
        
        void setDerivativesXiInv();
        void setDerivativesEtaInv();
        void setDerivativesZetaInv();
    };
}

#endif