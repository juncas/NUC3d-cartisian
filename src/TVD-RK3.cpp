#ifndef tvdrk3_cpp
#define tvdrk3_cpp
#include "TVD-RK3.h"

nuc3d::tvdrk3rd::tvdrk3rd():
integration(3)
{
}


nuc3d::tvdrk3rd::~tvdrk3rd()
{
    
}

void nuc3d::tvdrk3rd::integrationAll(const VectorField &RHS, // Right-hand-side: l*dt
                                     const VectorField &Q, // u(nstep)
                                     VectorField &Qi, //u_i
                                     double dt,
                                     int nstep)
{
    switch(nstep)
    {
        case 0:
            rk1st(RHS,Q,Qi,dt);
            break;
        case 1:
            rk2nd(RHS,Q,Qi,dt);
            break;
        case 2:
            rk3rd(RHS,Q,Qi,dt);
            break;
    }
    
}

void nuc3d::tvdrk3rd::rk1st(const VectorField &RHS,
                            const VectorField &u_n, // u_n
                            VectorField &u_i,
                            double dt)
{
    int nx=u_n.begin()->getSizeX();
    int ny=u_n.begin()->getSizeY();
    int nz=u_n.begin()->getSizeZ();
    
    auto beg=u_i.begin();
    auto end=u_i.end();
    
    for (auto iter=beg; iter!=end; ++iter)
    {
        double *pRHS=RHS[iter-beg].getDataPtr();
        double *pu_n=u_n[iter-beg].getDataPtr();
        double *pu_i=iter->getDataPtr();
        for(int k=0;k<nz;++k)
        {
            int idx_k=nx*ny*k;
            for(int j=0;j<ny;++j)
            {
                int idx_j=idx_k+nx*j;
                for(int i=0;i<nx;++i)
                {
                    int idx_i=idx_j+i;

                    double rhs=-pRHS[idx_i];
                    double u=pu_i[idx_i];
                    double u1 = coeff_tvdrk3_alpha0[0][0]*u
                    +coeff_tvdrk3_alpha0[0][1]*rhs*dt;
                    
                    pu_i[idx_i]=u1;
                }
            }
        }
    }
}

void nuc3d::tvdrk3rd::rk2nd(const VectorField &RHS,
                            const VectorField &u_n, // u_n
                            VectorField &u_i,
                            double dt)
{
    int nx=u_n.begin()->getSizeX();
    int ny=u_n.begin()->getSizeY();
    int nz=u_n.begin()->getSizeZ();
    
    auto beg=u_i.begin();
    auto end=u_i.end();
    
    
    for (auto iter=beg; iter!=end; ++iter)
    {
        double *pRHS=RHS[iter-beg].getDataPtr();
        double *pu_n=u_n[iter-beg].getDataPtr();
        double *pu_i=iter->getDataPtr();
        for(int k=0;k<nz;++k)
        {
            int idx_k=nx*ny*k;
            for(int j=0;j<ny;++j)
            {
                int idx_j=idx_k+nx*j;
                for(int i=0;i<nx;++i)
                {
                    int idx_i=idx_j+i;
                    double rhs=-pRHS[idx_i];
                    double u=pu_n[idx_i];
                    double u1=pu_i[idx_i];
                    
                    double u2 = coeff_tvdrk3_alpha0[1][0]*u
                    +coeff_tvdrk3_alpha0[1][1]*(rhs*dt+u1);
                    
                    pu_i[idx_i]=u2;
                }
            }
        }
    }
}

void nuc3d::tvdrk3rd::rk3rd(const VectorField &RHS,
                            const VectorField &u_n, // u_n
                            VectorField &u_i,
                            double dt)
{
    int nx=u_n.begin()->getSizeX();
    int ny=u_n.begin()->getSizeY();
    int nz=u_n.begin()->getSizeZ();
    
    auto beg=u_i.begin();
    auto end=u_i.end();
    
    
    for (auto iter=beg; iter!=end; ++iter)
    {
        double *pRHS=RHS[iter-beg].getDataPtr();
        double *pu_n=u_n[iter-beg].getDataPtr();
        double *pu_i=iter->getDataPtr();
        for(int k=0;k<nz;++k)
        {
            int idx_k=nx*ny*k;
            for(int j=0;j<ny;++j)
            {
                int idx_j=idx_k+nx*j;
                for(int i=0;i<nx;++i)
                {
                    int idx_i=idx_j+i;

                    double rhs=-pRHS[idx_i];
                    double u=pu_n[idx_i];
                    double u2=pu_i[idx_i];
                    
                    double u3 = coeff_tvdrk3_alpha0[2][0]*u
                    +coeff_tvdrk3_alpha0[2][1]*(rhs*dt+u2);
                    
                    pu_i[idx_i]=u3;
                }
            }
        }
    }
}

#endif