//
//  NaiverStokesData3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/3.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "NaiverStokes3d.h"
#include "physicsModel.h"
#include "PDEData3d.hpp"
#include "bufferData.hpp"
#include "gradvector.hpp"
#include <cmath>
/*************************************************************************************
 Member functions of class: NaiverStokesData3d
 **************************************************************************************/
nuc3d::NaiverStokesData3d::NaiverStokesData3d(int nx0, int ny0, int nz0, int neqs, double dx0, double dy0, double dz0):
EulerData3D(nx0,ny0,nz0,neqs,dx0,dy0,dz0),
du(nx0,ny0,nz0),
dv(nx0,ny0,nz0),
dw(nx0,ny0,nz0),
dT(nx0,ny0,nz0),
miu(nx0,ny0,nz0),
coeff(nx0,ny0,nz0),
tau(9,Field(nx0,ny0,nz0)),
Flux_xi_vis(neqs,Field(nx0,ny0,nz0)),
Flux_eta_vis(neqs,Field(ny0,nz0,nx0)),
Flux_zeta_vis(neqs,Field(nz0,nx0,ny0)),
dfvdxi(neqs,Field(nx0,ny0,nz0)),
dgvdeta(neqs,Field(ny0,nz0,nx0)),
dhvdzeta(neqs,Field(nz0,nx0,ny0))
{
    
}

void nuc3d::NaiverStokesData3d::solve(PDEData3d &myPDE,
                                      fieldOperator3d &myOP,
                                      std::vector<bufferData> &myBf,
                                      physicsModel &myModel,
                                      MPIComunicator3d_nonblocking &myMPI,
                                      fft &myFFT,
                                      const VectorField &xyz_center)
{
        // double t0=MPI_Wtime();
    EulerData3D::solveCon2Prim(myPDE, myModel);
        // double t1=MPI_Wtime();
    EulerData3D::solveRiemann(myPDE, myModel);
        // double t2=MPI_Wtime();
    EulerData3D::setBoundaryCondition(myPDE,myModel,myBf);
        // double t3=MPI_Wtime();
    EulerData3D::solveInv(myOP,myBf,myMPI);
        // double t4=MPI_Wtime();
    NaiverStokesData3d::solveVis(myPDE,myOP,myModel,myBf,myMPI);
        // double t5=MPI_Wtime();
    NaiverStokesData3d::solveRHS(myPDE);
        // double t6=MPI_Wtime();
    EulerData3D::solveSource(myPDE,myModel,myMPI,myFFT,xyz_center);
        // double t7=MPI_Wtime();

    // if(0==myMPI.getMyId())
    //     std::cout<<"c2p time is ="<<t1-t0<<" "
    // <<"riemann time is ="<<t2-t1<<" "
    // <<"boundary time is ="<<t3-t2<<" "
    // <<"inv time is ="<<t4-t3<<" "
    // <<"vis time is ="<<t5-t4<<" "
    // <<"rhs time is ="<<t6-t5<<" "
    // <<"source time is ="<<t7-t6<<std::endl;
    this->EulerData3D::getDt();

    myMPI.allReduceMin(dt);
    myPDE.setDt(dt);
}

void nuc3d::NaiverStokesData3d::solveVis(PDEData3d &myPDE,
                                         fieldOperator3d &myOP,
                                         physicsModel &myModel,
                                         std::vector<bufferData> &myBf,
                                         MPIComunicator3d_nonblocking &myMPI
                                         )
{
    solveGrads(myPDE, myOP, myBf, myMPI);
    solveViscousFlux(myModel);
    setDerivativesVis(myOP,myBf,myMPI);
}

void nuc3d::NaiverStokesData3d::solveGrads(PDEData3d &myPDE,
                                           fieldOperator3d &myOP,
                                           std::vector<bufferData> &myBf,
                                           MPIComunicator3d_nonblocking &myMPI)
{
    Field &u=this->EulerData3D::W_Euler[1];
    Field &v=this->EulerData3D::W_Euler[2];
    Field &w=this->EulerData3D::W_Euler[3];
    Field &T=this->EulerData3D::W0_Euler[0];
    
    
    solve_grad(u,myOP,myBf[0],myMPI,du,0);
    
    solve_grad(v,myOP,myBf[1],myMPI,dv,1);
    
    solve_grad(w,myOP,myBf[2],myMPI,dw,2);
    
    solve_grad(T,myOP,myBf[3],myMPI,dT,3);

    myMPI.barrier();
}

void nuc3d::NaiverStokesData3d::solve_grad(Field &myField,
                                           fieldOperator3d &myOP,
                                           bufferData &myBf,
                                           MPIComunicator3d_nonblocking &myMPI,
                                           gradvector &myGrad,
                                           int fdID)
{
    Field &dxi=myGrad.getdxi();
    Field &deta=myGrad.getdeta();
    Field &dzeta=myGrad.getdzeta();
    
    Field &fxi=myGrad.getf_xi();
    Field &feta=myGrad.getf_eta();
    Field &fzeta=myGrad.getf_zeta();
    
    
    myGrad.setGrad(myField);
    
    solveGradXi(fxi,myOP,myBf,myMPI,dxi,fdID);    
    solveGradEta(feta,myOP,myBf,myMPI,deta,fdID);    
    solveGradZeta(fzeta,myOP,myBf,myMPI,dzeta,fdID);
}

void nuc3d::NaiverStokesData3d::solveGradXi(Field &myField,
                                            fieldOperator3d &myOP,
                                            bufferData &myBf,
                                            MPIComunicator3d_nonblocking &myMPI,
                                            Field &dxi,
                                            int fdID)
{
    myMPI.bufferSendRecv(myField, myBf, 0, fdID);
    myOP.differenceInner(myField,dx, 0, dxi);
    myMPI.waitSendRecv(myBf, 0);
    myOP.differenceBoundary(myField,dx, myBf.BufferRecv[0], myBf.BufferRecv[1], 0, dxi);

    myMPI.barrier();
}

void nuc3d::NaiverStokesData3d::solveGradEta(Field &myField,
                                             fieldOperator3d &myOP,
                                             bufferData &myBf,
                                             MPIComunicator3d_nonblocking &myMPI,
                                             Field &deta,
                                             int fdID)
{
    myMPI.bufferSendRecv(myField, myBf, 1, fdID);
    myOP.differenceInner(myField,dy, 1, deta);
    myMPI.waitSendRecv(myBf, 1);
    myOP.differenceBoundary(myField,dy, myBf.BufferRecv[2], myBf.BufferRecv[3], 1,deta);

    myMPI.barrier();
}

void nuc3d::NaiverStokesData3d::solveGradZeta(Field &myField,
                                              fieldOperator3d &myOP,
                                              bufferData &myBf,
                                              MPIComunicator3d_nonblocking &myMPI,
                                              Field &dzeta,
                                              int fdID)
{
    myMPI.bufferSendRecv(myField, myBf, 2, fdID);
    myOP.differenceInner(myField,dz, 2, dzeta);
    myMPI.waitSendRecv(myBf, 2);
    myOP.differenceBoundary(myField,dz, myBf.BufferRecv[4], myBf.BufferRecv[5], 2, dzeta);

    myMPI.barrier();
}

void nuc3d::NaiverStokesData3d::solveViscousFlux(physicsModel &myPhyMod)
{
    
    myPhyMod.getMiu(this->EulerData3D::W0_Euler[0], miu, coeff);
    
    double *pu=this->EulerData3D::W_Euler[1].getDataPtr();
    double *pv=this->EulerData3D::W_Euler[2].getDataPtr();
    double *pw=this->EulerData3D::W_Euler[3].getDataPtr();
    
    double *uxi=du.getdxi().getDataPtr();
    double *vxi=dv.getdxi().getDataPtr();
    double *wxi=dw.getdxi().getDataPtr();
    double *txi=dT.getdxi().getDataPtr();
    
    double *ueta=du.getdeta().getDataPtr();
    double *veta=dv.getdeta().getDataPtr();
    double *weta=dw.getdeta().getDataPtr();
    double *teta=dT.getdeta().getDataPtr();
    
    double *uzeta=du.getdzeta().getDataPtr();
    double *vzeta=dv.getdzeta().getDataPtr();
    double *wzeta=dw.getdzeta().getDataPtr();
    double *tzeta=dT.getdzeta().getDataPtr();
        
    double *pMiu=miu.getDataPtr();
    double *pCoeff=coeff.getDataPtr();
    
    double *flux_xi[5];
    double *flux_eta[5];
    double *flux_zeta[5];
    
    double fv[5];
    double gv[5];
    double hv[5];
    
    flux_xi[0]=Flux_xi_vis[0].getDataPtr();
    flux_xi[1]=Flux_xi_vis[1].getDataPtr();
    flux_xi[2]=Flux_xi_vis[2].getDataPtr();
    flux_xi[3]=Flux_xi_vis[3].getDataPtr();
    flux_xi[4]=Flux_xi_vis[4].getDataPtr();
    
    flux_eta[0]=Flux_eta_vis[0].getDataPtr();
    flux_eta[1]=Flux_eta_vis[1].getDataPtr();
    flux_eta[2]=Flux_eta_vis[2].getDataPtr();
    flux_eta[3]=Flux_eta_vis[3].getDataPtr();
    flux_eta[4]=Flux_eta_vis[4].getDataPtr();
    
    flux_zeta[0]=Flux_zeta_vis[0].getDataPtr();
    flux_zeta[1]=Flux_zeta_vis[1].getDataPtr();
    flux_zeta[2]=Flux_zeta_vis[2].getDataPtr();
    flux_zeta[3]=Flux_zeta_vis[3].getDataPtr();
    flux_zeta[4]=Flux_zeta_vis[4].getDataPtr();
    
    for (int k=0; k<nz; k++)
    {
                int idx_k_xi=nx*ny*k;
                int idx_k_eta=ny*k;
                int idx_k_zeta=k;
        for (int j=0; j<ny; j++)
        {
                int idx_j_xi=idx_k_xi+nx*j;
                int idx_j_eta=idx_k_eta+j;
                int idx_j_zeta=idx_k_zeta+nz*nx*j;
            for (int i=0; i<nx; i++)
            {
                int idx_xi=idx_j_xi+i;
                int idx_eta=idx_j_eta+ny*nz*i;
                int idx_zeta=idx_j_zeta+nz*i;
                                
                //d(u,v,w,T)/d(x,y,z)
                double ux=uxi[idx_xi];
                double uy=ueta[idx_eta];
                double uz=uzeta[idx_zeta];
                
                double vx=vxi[idx_xi];
                double vy=veta[idx_eta];
                double vz=vzeta[idx_zeta];
                
                double wx=wxi[idx_xi];
                double wy=weta[idx_eta];
                double wz=wzeta[idx_zeta];
                
                double tx=txi[idx_xi];
                double ty=teta[idx_eta];
                double tz=tzeta[idx_zeta];
                
                double grad=ux+vy+wz;
                double miu0=pMiu[idx_xi];
                
                double tau_xx=miu0*(2.0*ux-2.0/3.0*grad);
                double tau_xy=miu0*(uy+vx-2.0/3.0*grad);
                double tau_xz=miu0*(uz+wx-2.0/3.0*grad);
                double tau_yy=miu0*(2.0*vy-2.0/3.0*grad);
                double tau_yz=miu0*(vz+wy-2.0/3.0*grad);
                double tau_zz=miu0*(2.0*wz-2.0/3.0*grad);
                
                double coeff0=pCoeff[idx_xi];
                double tau_tx=coeff0*tx;
                double tau_ty=coeff0*ty;
                double tau_tz=coeff0*tz;                
                
                double u=pu[idx_xi];
                double v=pv[idx_xi];
                double w=pw[idx_xi];
                
                flux_xi[0][idx_xi]=0.0;
                flux_xi[1][idx_xi]=tau_xx;
                flux_xi[2][idx_xi]=tau_xy;
                flux_xi[3][idx_xi]=tau_xz;
                flux_xi[4][idx_xi]=u*tau_xx+v*tau_xy+w*tau_xz+tau_tx;
                
                flux_eta[0][idx_eta]=0.0;
                flux_eta[1][idx_eta]=tau_xy;
                flux_eta[2][idx_eta]=tau_yy;
                flux_eta[3][idx_eta]=tau_yz;
                flux_eta[4][idx_eta]=u*tau_xy+v*tau_yy+w*tau_yz+tau_ty;
                
                flux_zeta[0][idx_zeta]=0.0;
                flux_zeta[1][idx_zeta]=tau_xz;
                flux_zeta[2][idx_zeta]=tau_yz;
                flux_zeta[3][idx_zeta]=tau_zz;
                flux_zeta[4][idx_zeta]=u*tau_xz+v*tau_yz+w*tau_zz+tau_tz;
                                
            }
        }
    }
}

void nuc3d::NaiverStokesData3d::setDerivativesVis(fieldOperator3d &myOP,
                                                  std::vector<bufferData> &myBf,
                                                  MPIComunicator3d_nonblocking &myMPI)
{
    setDerivativeXi(myOP,myBf,myMPI);
    setDerivativeEta(myOP,myBf,myMPI);
    setDerivativeZeta(myOP,myBf,myMPI);

    myMPI.barrier();
}

void nuc3d::NaiverStokesData3d::setDerivativeXi(fieldOperator3d &myOP,
                                                std::vector<bufferData> &myBf,
                                                MPIComunicator3d_nonblocking &myMPI)
{
    auto beg=Flux_xi_vis.begin();
    auto end=Flux_xi_vis.end();
    for(auto iter=beg;iter!=end;iter++)
    {
        myMPI.bufferSendRecv(*iter, myBf[iter-beg], 0, static_cast<int>(iter-beg));
        myOP.differenceInner(*iter,dx, 0, dfvdxi[iter-beg]);
        myMPI.waitSendRecv(myBf[iter-beg], 0);
        myOP.differenceBoundary(*iter,dx, myBf[iter-beg].BufferRecv[0], myBf[iter-beg].BufferRecv[1], 0, dfvdxi[iter-beg]);
    }
    

    myMPI.barrier();
    
}

void nuc3d::NaiverStokesData3d::setDerivativeEta(fieldOperator3d &myOP,
                                                 std::vector<bufferData> &myBf,
                                                 MPIComunicator3d_nonblocking &myMPI)
{
    auto beg=Flux_eta_vis.begin();
    auto end=Flux_eta_vis.end();
    for(auto iter=beg;iter!=end;iter++)
    {
        myMPI.bufferSendRecv(*iter, myBf[iter-beg], 1, static_cast<int>(iter-beg));
        myOP.differenceInner(*iter,dy, 1, dgvdeta[iter-beg]);
        myMPI.waitSendRecv(myBf[iter-beg], 1);
        myOP.differenceBoundary(*iter,dy, myBf[iter-beg].BufferRecv[2], myBf[iter-beg].BufferRecv[3], 1, dgvdeta[iter-beg]);
    }

    myMPI.barrier();
    
}

void nuc3d::NaiverStokesData3d::setDerivativeZeta(fieldOperator3d &myOP,
                                                  std::vector<bufferData> &myBf,
                                                  MPIComunicator3d_nonblocking &myMPI)
{
    auto beg=Flux_zeta_vis.begin();
    auto end=Flux_zeta_vis.end();
    for(auto iter=beg;iter!=end;iter++)
    {
        myMPI.bufferSendRecv(*iter, myBf[iter-beg], 2, static_cast<int>(iter-beg));
        myOP.differenceInner(*iter,dz, 2, dhvdzeta[iter-beg]);
        myMPI.waitSendRecv(myBf[iter-beg], 2);
        myOP.differenceBoundary(*iter,dz, myBf[iter-beg].BufferRecv[4], myBf[iter-beg].BufferRecv[5], 2, dhvdzeta[iter-beg]);
    }
    

    myMPI.barrier();
}


void nuc3d::NaiverStokesData3d::solveRHS(PDEData3d &myPDE)
{
    VectorField &df_v=this->EulerData3D::dfdxi;
    VectorField &dg_v=this->EulerData3D::dgdeta;
    VectorField &dh_v=this->EulerData3D::dhdzeta;
    
    VectorField &dfv_v=this->dfvdxi;
    VectorField &dgv_v=this->dgvdeta;
    VectorField &dhv_v=this->dhvdzeta;
    
    auto beg=myPDE.getRHS().begin();
    auto end=myPDE.getRHS().end();
    
    for (auto iter=beg ; iter!=end ; iter++)
    {
        double *df=df_v[iter-beg].getDataPtr();
        double *dg=dg_v[iter-beg].getDataPtr();
        double *dh=dh_v[iter-beg].getDataPtr();
        
        double *dfv=dfv_v[iter-beg].getDataPtr();
        double *dgv=dgv_v[iter-beg].getDataPtr();
        double *dhv=dhv_v[iter-beg].getDataPtr();
        
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
          for (int k = 0; k < nz; ++k)
        {
            int idx_k=nx*ny*k;
            for (int j = 0; j < ny; ++j)
            {
                int idx_j=idx_k+nx*j;
                for (int i = 0; i < nx; ++i)
                {                
                    int idx_i=idx_j+i;
                    
                    
                    rhs[idx_i]-=dfv[idx_i];
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
                    
                    rhs[idx_i_xi]-=dgv[idx_i_eta];
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
                    
                    rhs[idx_i_xi]-=dhv[idx_i_zeta];
                }
            }
        }
        
    }
}

nuc3d::NaiverStokesData3d::~NaiverStokesData3d()
{}

