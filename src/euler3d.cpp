#ifndef euler3d_cpp
#define euler3d_cpp
#include <cmath>
#include "euler3d.h"
#include "physicsModel.h"
#include "PDEData3d.hpp"
#include <cmath>

/**************************************************************************************
 Member functions of class: EulerFlux
 **************************************************************************************/
 nuc3d::EulerFlux::EulerFlux(int nx0,int ny0,int nz0,int neqs):
 FluxL(neqs,Field(nx0,ny0,nz0)),
 FluxR(neqs,Field(nx0,ny0,nz0)),
 reconstFluxL(neqs,Field(nx0+1,ny0,nz0)),
 reconstFluxR(neqs,Field(nx0+1,ny0,nz0)),
 reconstFlux(neqs,Field(nx0+1,ny0,nz0)),
 maxEigen(1.0)
 {

 }

 void nuc3d::EulerFlux::combineFluxLR()
 {
    auto beg=reconstFlux.begin();
    auto end=reconstFlux.end();
    
    for(auto iter=beg;iter!=end;iter++)
    {
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        double *prf=iter->getDataPtr();
        double *prfp=reconstFluxL[iter-beg].getDataPtr();
        double *prfn=reconstFluxR[iter-beg].getDataPtr();
        
        for (int k=0; k<nz0; k++)
        {
            int idx_k=nx0*ny0*k;
            for (int j=0; j<ny0; j++)
            {
                int idx_j=idx_k+nx0*j;
                for (int i=0; i<nx0; i++)
                {
                   int idx_i=idx_j+i;
                   double rfp=prfp[idx_i];
                   double rfn=prfn[idx_i];
                   double rf=rfp+rfn;

                   prf[idx_i]=rf;
               }
           }
       }
   }
}

nuc3d::EulerFlux::~EulerFlux()
{}

/**************************************************************************************
 Member functions of class: EulerData3D
 **************************************************************************************/
 nuc3d::EulerData3D::EulerData3D( int nx0, int ny0, int nz0, int neqs, double dx0, double dy0, double dz0):
 nx(nx0),
 ny(ny0),
 nz(nz0),
 W0_Euler(3,Field(nx0,ny0,nz0)),
 W_Euler(neqs,Field(nx0,ny0,nz0)),
 source(neqs,Field(nx0,ny0,nz0)),
 Flux_xi(nx0,ny0,nz0,neqs),
 Flux_eta(ny0,nz0,nx0,neqs),
 Flux_zeta(nz0,nx0,ny0,neqs),
 dfdxi(neqs,Field(nx0,ny0,nz0)),
 dgdeta(neqs,Field(ny0,nz0,nx0)),
 dhdzeta(neqs,Field(nz0,nx0,ny0)),
 dt(0.0),
 dx(dx0),
 dy(dy0),
 dz(dz0)
 {

 }

 nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeXi()
 {
    return this->dfdxi;
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeEta()
{
    return this->dgdeta;
}

nuc3d::VectorField& nuc3d::EulerData3D::getDrivativeZeta()
{
    return this->dhdzeta;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxXi()
{
    return Flux_xi;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxEta()
{
    return Flux_eta;
}

nuc3d::EulerFlux& nuc3d::EulerData3D::getFluxZeta()
{
    return Flux_zeta;
}

nuc3d::VectorField& nuc3d::EulerData3D::getPrimatives()
{
    return W_Euler;
}


nuc3d::VectorField& nuc3d::EulerData3D::getAcoustics()
{
    return W0_Euler;
}

nuc3d::VectorField& nuc3d::EulerData3D::getSource()
{
    return source;
}


void nuc3d::EulerData3D::setDerivativesInv()
{
    Flux_xi.combineFluxLR();
    Flux_eta.combineFluxLR();
    Flux_zeta.combineFluxLR();
    
    EulerData3D::setDerivativesXiInv();
    EulerData3D::setDerivativesEtaInv();
    EulerData3D::setDerivativesZetaInv();
}

void nuc3d::EulerData3D::setDerivativesXiInv()
{
    VectorField &flux=Flux_xi.reconstFlux;
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        double *f=iter->getDataPtr();
        double *df=dfdxi[iter-flux.begin()].getDataPtr();
        
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        int nx1=dfdxi[iter-flux.begin()].getSizeX();
        int ny1=dfdxi[iter-flux.begin()].getSizeY();
        int nz1=dfdxi[iter-flux.begin()].getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            int idx_k=nx0*ny0*k;
            int idx1_k=nx1*ny1*k;
            for (int j=0; j<ny0; j++)
            {
             int idx_j=idx_k+nx0*j;
             int idx1_j=idx1_k+nx1*j;
             for (int i=1; i<nx0; i++)
             {
                int idx_i=idx_j+i;
                int idx1_i=idx1_j+i;
                double df_loc=f[idx_i]-f[idx_i-1];
                
                df[idx1_i-1]=df_loc/dx;
            }
        }
    }
}
}

void nuc3d::EulerData3D::setDerivativesEtaInv()
{
    VectorField &flux=Flux_eta.reconstFlux;
    
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        double *f=iter->getDataPtr();
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        double *dg=dgdeta[iter-flux.begin()].getDataPtr();
        int nx1=dgdeta[iter-flux.begin()].getSizeX();
        int ny1=dgdeta[iter-flux.begin()].getSizeY();
        int nz1=dgdeta[iter-flux.begin()].getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            int idx_k=nx0*ny0*k;
            int idx1_k=nx1*ny1*k;
            for (int j=0; j<ny0; j++)
            {
             int idx_j=idx_k+nx0*j;
             int idx1_j=idx1_k+nx1*j;
             for (int i=1; i<nx0; i++)
             {
                int idx_i=idx_j+i;
                int idx1_i=idx1_j+i;
                double df_loc=f[idx_i]-f[idx_i-1];
                
                dg[idx1_i-1]=df_loc/dy;
            }
        }
    }
}
}

void nuc3d::EulerData3D::setDerivativesZetaInv()
{
    VectorField &flux=Flux_zeta.reconstFlux;
    
    for (auto iter=flux.begin() ; iter!=flux.end(); iter++)
    {
        double *f=iter->getDataPtr();
        int nx0=iter->getSizeX();
        int ny0=iter->getSizeY();
        int nz0=iter->getSizeZ();
        
        double *dh=dhdzeta[iter-flux.begin()].getDataPtr();
        int nx1=dhdzeta[iter-flux.begin()].getSizeX();
        int ny1=dhdzeta[iter-flux.begin()].getSizeY();
        int nz1=dhdzeta[iter-flux.begin()].getSizeZ();
        
        for (int k=0; k<nz0; k++)
        {
            int idx_k=nx0*ny0*k;
            int idx1_k=nx1*ny1*k;
            for (int j=0; j<ny0; j++)
            {
                int idx_j=idx_k+nx0*j;
                int idx1_j=idx1_k+nx1*j;
                for (int i=1; i<nx0; i++)
                {
                    int idx_i=idx_j+i;
                    int idx1_i=idx1_j+i;
                    double df_loc=f[idx_i]-f[idx_i-1];
                    
                    dh[idx1_i-1]=df_loc/dz;
                }
            }
        }
    }
    
}

nuc3d::EulerData3D::~EulerData3D()
{

}

void nuc3d::EulerData3D::solve(PDEData3d &myPDE,
 fieldOperator3d &myOP,
 std::vector<bufferData> &myBf,
 physicsModel &myModel,
 MPIComunicator3d_nonblocking &myMPI,
 fft &myFFT,
 const VectorField &xyz_center)
{
    nuc3d::EulerData3D::solveCon2Prim(myPDE, myModel);
    nuc3d::EulerData3D::solveRiemann(myPDE, myModel);
    nuc3d::EulerData3D::setBoundaryCondition(myPDE,myModel,myBf);
    nuc3d::EulerData3D::solveInv(myOP,myBf,myMPI);
    nuc3d::EulerData3D::solveRHS(myPDE);
    nuc3d::EulerData3D::solveSource(myPDE,myModel,myMPI,myFFT,xyz_center);
    getDt();
    myMPI.allReduceMin(dt);
    myPDE.setDt(dt);
}

void nuc3d::EulerData3D::solveRiemann(PDEData3d &myPDE,
  physicsModel &myModel)
{
    myModel.solveRiemann(myPDE, this);
}

void nuc3d::EulerData3D::solveCon2Prim(PDEData3d &myPDE,
 physicsModel &myModel)
{
    myModel.solve(myPDE, this);
}

void nuc3d::EulerData3D::solveSource(PDEData3d &myPDE,
 physicsModel &myModel,
 MPIComunicator3d_nonblocking &myMPI,
 fft &myFFT,
 const VectorField &xyz_center)
{
    myModel.solveSource(myPDE, this, myMPI, myFFT, xyz_center);
}

void nuc3d::EulerData3D::setBoundaryCondition(PDEData3d &myPDE,
  physicsModel &myModel,
  std::vector<bufferData> &myBf)
{

}

void nuc3d::EulerData3D::solveInv(fieldOperator3d &myOP,
  std::vector<bufferData> &myBf,
  MPIComunicator3d_nonblocking &myMPI)
{  
    // double t0=MPI_Wtime();
    solveInvicidFluxL(this->getFluxXi(), myOP, myBf,myMPI, 0);
    solveInvicidFluxL(this->getFluxEta(), myOP, myBf,myMPI, 1);
    solveInvicidFluxL(this->getFluxZeta(), myOP, myBf,myMPI, 2);

    // double t1=MPI_Wtime();
    solveInvicidFluxR(this->getFluxXi(), myOP, myBf,myMPI, 0);
    solveInvicidFluxR(this->getFluxEta(), myOP, myBf,myMPI, 1);
    solveInvicidFluxR(this->getFluxZeta(), myOP, myBf,myMPI,2);

    // double t2=MPI_Wtime();
    this->setDerivativesInv();
    // double t3=MPI_Wtime();
    // if(0==myMPI.getMyId())
    //     std::cout<<"fluxL time is ="<<t1-t0<<" "
    // <<"FLuxR time is ="<<t2-t1<<" "
    // <<"SetDer time is ="<<t3-t2<<" "
    // <<std::endl;

    
}


void nuc3d::EulerData3D::solveInvicidFluxL(EulerFlux &myFlux,
 fieldOperator3d &myOP,
 std::vector<bufferData> &myBuff,
 MPIComunicator3d_nonblocking &myMPI,
 int dir)
{
    VectorField &pFlux = myFlux.FluxL;
    VectorField &pReconFlux = myFlux.reconstFluxL;    
    
    for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
    {
        bufferData &bf = myBuff[iter - pFlux.begin()];
        Field &rf = pReconFlux[iter - pFlux.begin()];
        
        myMPI.bufferSendRecv(*iter,bf,dir,static_cast<int>(iter - pFlux.begin()));   
        myOP.reconstructionInner(*iter, dir, 1, rf);    

        myMPI.waitSendRecv(bf,dir);
        
        Field &bfField_L=bf.BufferRecv[dir*2];
        Field &bfField_R=bf.BufferRecv[dir*2+1];
        myOP.reconstructionBoundary(*iter, bfField_L, bfField_R, dir, 1, rf);
    }

    myMPI.barrier();
}

void nuc3d::EulerData3D::solveInvicidFluxR(EulerFlux &myFlux,
 fieldOperator3d &myOP,
 std::vector<bufferData> &myBuff,
 MPIComunicator3d_nonblocking &myMPI,
 int dir)

{
    VectorField &pFlux = myFlux.FluxR;
    VectorField &pReconFlux = myFlux.reconstFluxR;
    
    for (auto iter = pFlux.begin(); iter != pFlux.end(); iter++)
    {
        bufferData &bf = myBuff[iter - pFlux.begin()];
        Field &rf = pReconFlux[iter - pFlux.begin()];
        
        myMPI.bufferSendRecv(*iter,bf,dir,static_cast<int>(iter - pFlux.begin()));        
        myOP.reconstructionInner(*iter, dir, -1, rf);        
        myMPI.waitSendRecv(bf,dir);
        
        Field &bfField_L=bf.BufferRecv[dir*2];
        Field &bfField_R=bf.BufferRecv[dir*2+1];
        
        myOP.reconstructionBoundary(*iter, bfField_L, bfField_R, dir, -1, rf);
    }
    
    myMPI.barrier();
}

void nuc3d::EulerData3D::solveRHS(PDEData3d &myPDE)
{
    VectorField &df_v=this->dfdxi;
    VectorField &dg_v=this->dgdeta;
    VectorField &dh_v=this->dhdzeta;
    auto beg=myPDE.getRHS().begin();
    auto end=myPDE.getRHS().end();
    
    for (auto iter=beg ; iter!=end ; iter++)
    {
        double *df=df_v[iter-beg].getDataPtr();
        double *dg=dg_v[iter-beg].getDataPtr();
        double *dh=dh_v[iter-beg].getDataPtr();
        double *rhs=iter->getDataPtr();
        
        for (int k = 0; k < nz; ++k)
        {
            int idx_k=nx*ny*k;
            for (int j = 0; j < ny; ++j)
            {
                int idx_j=idx_k+nx*j;
                for (int i = 0; i < nx; ++i)
                {                
                    int idx_i=idx_j+i;
                    
                    
                    rhs[idx_i]=df[idx_i];
                }
            }
        }
        
        for (int k = 0; k < nz; ++k)
        {
            int idx_k_xi=nx*ny*k;
            int idx_k_eta=ny*k;
            for (int j = 0; j < ny; ++j)
            {
                int idx_j_xi=idx_k_xi+nx*j;
                int idx_j_eta=idx_k_eta+j;
                for (int i = 0; i < nx; ++i)
                {
                    int idx_i_xi=idx_j_xi+i;
                    int idx_i_eta=idx_j_eta+ny*nz*i;
                    
                    rhs[idx_i_xi]+=dg[idx_i_eta];
                }
            }
        }
        
        for (int k = 0; k < nz; ++k)
        {
            int idx_k_xi=nx*ny*k;
            int idx_k_zeta=k;
            for (int j = 0; j < ny; ++j)
            {
                int idx_j_xi=idx_k_xi+nx*j;
                int idx_j_zeta=idx_k_zeta+nz*nx*j;
                for (int i = 0; i < nx; ++i)
                {
                    int idx_i_xi=idx_j_xi+i;
                    int idx_i_zeta=idx_j_zeta+nz*i;
                    
                    rhs[idx_i_xi]+=dh[idx_i_zeta];
                }
            }
        }
    }
    
}

void nuc3d::EulerData3D::getDt()
{
    double dt_xi=dx/Flux_xi.maxEigen;
    double dt_eta=dy/Flux_eta.maxEigen;
    double dt_zeta=dz/Flux_zeta.maxEigen;
    
    dt=1.0/(1.0/dt_xi+1.0/dt_eta+1.0/dt_zeta);
}

#endif
