//
//  PDEData3d.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/11/4.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "PDEData3d.hpp"
#include <cmath>
/**************************************************************************************
 Member functions of class: PDEData3d
 **************************************************************************************/
nuc3d::PDEData3d::PDEData3d()
{
    
}

nuc3d::PDEData3d::PDEData3d(int nx0,int ny0,int nz0,int neqs):
nEquations(neqs),
Q(neqs,Field(nx0,ny0,nz0)),
Q_Euler(neqs,Field(nx0,ny0,nz0)),
Q_work(neqs,Field(nx0,ny0,nz0)),
RHS(neqs,Field(nx0,ny0,nz0)),
dt_local(0.0),
dt_global(0.0),
res_global(0.0),
res_local(0.0)
{
    
}

void nuc3d::PDEData3d::initPDEData3d(int nx0,int ny0,int nz0,int neqs)
{
    nEquations=neqs;
    for(int i=0;i<neqs;i++)
    {
        Q.push_back(Field(nx0,ny0,nz0));
        Q_Euler.push_back(Field(nx0,ny0,nz0));
        Q_work.push_back(Field(nx0,ny0,nz0));
        RHS.push_back(Field(nx0,ny0,nz0));
    }
    
    dt_local=0.0;
    dt_global=0.0;
    res_global=0.0;
    res_local=0.0;
}


nuc3d::VectorField& nuc3d::PDEData3d::getRHS()
{
    return RHS;
}

nuc3d::VectorField& nuc3d::PDEData3d::getQ()
{
    return Q_work;
}
nuc3d::VectorField& nuc3d::PDEData3d::getQ_clean()
{
    return Q;
}

nuc3d::VectorField& nuc3d::PDEData3d::getQcurrent()
{
    return Q_Euler;
}

void nuc3d::PDEData3d::solve(fieldOperator3d &myOP,
                             MPIComunicator3d_nonblocking &myMPI,
                             double cfl,
                             int step)
{

    if(step==0)
    {
        Q_Euler=Q_work;        
        dt_global=dt_local*cfl;
    }

    myOP.timeIntegral(RHS, Q_Euler,Q_work, dt_global, step);
    
    if(step==2)
    {
        setRES(myMPI);
    }
    
    myMPI.barrier();
    
}

void  nuc3d::PDEData3d::setRES(MPIComunicator3d_nonblocking &myMPI)
{
    auto beg=Q_Euler.begin();
    auto end=Q_Euler.end();
    
    int nx=beg->getSizeX();
    int ny=beg->getSizeY();
    int nz=beg->getSizeZ();
    
    double res_this=0.0;
    static double res_last=0.0;
    
    for(auto iter=beg;iter!=end;iter++)
    {
        double *f=Q_work[iter-beg].getDataPtr();
        double *f0=iter->getDataPtr();
        for (int k=0; k<nz; k++)
        {
            for (int j=0; j<ny; j++)
            {
                for (int i=0; i<nx; i++)
                {
                    int idx_xi=nx*ny*k+nx*j+i;
                    
                    double u=f[idx_xi];
                    double u0=f0[idx_xi];
                    res_this+=std::pow((u-u0),2);
                }
            }
        }
    }

    myMPI.allReduceSum(res_this);
    
    res_global=std::sqrt(res_this/(nx*ny*nz*myMPI.getMyDim(2)*myMPI.getMyDim(1)*myMPI.getMyDim(0)*Q_work.size()))/(dt_global);
}

void nuc3d::PDEData3d::initialQ_work()
{
    auto beg=Q.begin();
    auto end=Q.end();
    
    for(auto iter=beg;iter!=end;iter++)
    {
        int nx=iter->getSizeX();
        int ny=iter->getSizeY();
        int nz=iter->getSizeZ();

        for (int k=0;k<nz;k++)
        {
            for (int j=0;j<ny ;j++ )
            {
                for (int i=0;i<nx ;i++ )
                {
                    double q0=iter->getValue(i, j, k);
                    double qj=q0;
                    
                    Q_work[iter-beg].setValue(i, j, k, qj);
                }
            }
        }
    }
    
}

void nuc3d::PDEData3d::setQ_clean()
{    
    auto beg=Q_work.begin();
    auto end=Q_work.end();
    
    for(auto iter=beg;iter!=end;iter++)
    {
        int nx=iter->getSizeX();
        int ny=iter->getSizeY();
        int nz=iter->getSizeZ();

        for (int k=0;k<nz;k++)
        {
            for (int j=0;j<ny ;j++ )
            {
                for (int i=0;i<nx ;i++ )
                {
                    double q0=iter->getValue(i, j, k);
                    double qj=q0;
                    
                    Q[iter-beg].setValue(i, j, k, qj);
                }
            }
        }
    }
}

void nuc3d::PDEData3d::setDt(double dt)
{
	dt_local=dt;
}
nuc3d::PDEData3d::~PDEData3d()
{}
