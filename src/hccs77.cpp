//
//  hccs77.cpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/24.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//

#include "hccs77.hpp"
#include "schemes.hpp"
#include <cmath>

nuc3d::hccs77::hccs77():
ss(1.0e-20),
p(2)
{

}

nuc3d::hccs77::~hccs77()
{

}

void nuc3d::hccs77::interpolationInner(const Field & fieldIN,
                                     const int uw,//uw=(1,-1)
                                     Field & fieldOUT,
                                     const int tilesize)
{
    switch(uw)
    {
        case 1:
        hccs77p(fieldIN, fieldOUT, tilesize);
        break;
        case -1:
        hccs77n(fieldIN, fieldOUT, tilesize);
        break;
        default:
        std::cout<<"weno error: no such direction"<<std::endl;
        exit(-1);
    }
    
}

void nuc3d::hccs77::hccs77p(const Field & fieldIN,
    Field & fieldOUT,
    const int tilesize)
{

    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();

    
    int ibeg=(tilesize-1);
    int iend=nx-tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;

    double *a=new double[iend-ibeg];
    double *b=new double[iend-ibeg];
    double *c=new double[iend-ibeg];
    double *d=new double[iend-ibeg];
    
    for(int k=kbeg;k<kend;k++)
    {
        double *pData_in_k=pIn+nx*ny*k;
        double *pData_out_k=pOut+nx0*ny0*k;
        
        for(int j=jbeg;j<jend;j++)
        {
            double *pData_in=pData_in_k+nx*j;
            double *pData_out=pData_out_k+nx0*j;

            for(int i=ibeg;i<iend;i++)
            {       
                    double is0=pData_in[i-3]*(547.0*pData_in[i-3]-3882.0*pData_in[i-2]+4642.0*pData_in[i-1]-1854.0*pData_in[i])
                    +pData_in[i-2]*(7043.0*pData_in[i-2]-17246.0*pData_in[i-1]+7042.0*pData_in[i])
                    +pData_in[i-1]*(11003.0*pData_in[i-1]-9402.0*pData_in[i])
                    +2107.0*pData_in[i]*pData_in[i];

                    double q40=coeff_weno7_alpha[0][0]*pData_in[i-3]+coeff_weno7_alpha[0][1]*pData_in[i-2]
                    +coeff_weno7_alpha[0][2]*pData_in[i-1]+coeff_weno7_alpha[0][3]*pData_in[i];

                    double is1=pData_in[i-2]*(267.0*pData_in[i-2]-1642.0*pData_in[i-1]+1602.0*pData_in[i]-494.0*pData_in[i+1])
                    +pData_in[i-1]*(2843.0*pData_in[i-1]-5966.0*pData_in[i]+1922.0*pData_in[i+1])
                    +pData_in[i]*(3443.0*pData_in[i]-2522.0*pData_in[i+1])
                    +547.0*pData_in[i+1]*pData_in[i+1];

                    double q41=coeff_weno7_alpha[1][0]*pData_in[i-2]+coeff_weno7_alpha[1][1]*pData_in[i-1]
                    +coeff_weno7_alpha[1][2]*pData_in[i]+coeff_weno7_alpha[1][3]*pData_in[i+1];

                    double q7=(-pData_in[i-2]+19.0*pData_in[i-1]+239.0*pData_in[i]+159.0*pData_in[i+1]+4.0*pData_in[i+2])/420.0;

                    double is2=pData_in[i-1]*(547.0*pData_in[i-1]-2522.0*pData_in[i]+1922.0*pData_in[i+1]-494.0*pData_in[i+2])
                    +pData_in[i]*(3443.0*pData_in[i]-5966.0*pData_in[i+1]+1602.0*pData_in[i+2])
                    +pData_in[i+1]*(2843.0*pData_in[i+1]-1642.0*pData_in[i+2])
                    +267.0*pData_in[i+2]*pData_in[i+2];

                    double q42=coeff_weno7_alpha[2][0]*pData_in[i-1]+coeff_weno7_alpha[2][1]*pData_in[i]
                    +coeff_weno7_alpha[2][2]*pData_in[i+1]+coeff_weno7_alpha[2][3]*pData_in[i+2];

                    double is3=pData_in[i]*(2107.0*pData_in[i]-9402.0*pData_in[i+1]+7042.0*pData_in[i+2]-1854.0*pData_in[i+3])
                    +pData_in[i+1]*(11003.0*pData_in[i+1]-17246.0*pData_in[i+2]+4642.0*pData_in[i+3])
                    +pData_in[i+2]*(7043.0*pData_in[i+2]-3882.0*pData_in[i+3])
                    +547.0*pData_in[i+3]*pData_in[i+3];

                    double q43=coeff_weno7_alpha[3][0]*pData_in[i]+coeff_weno7_alpha[3][1]*pData_in[i+1]
                    +coeff_weno7_alpha[3][2]*pData_in[i+2]+coeff_weno7_alpha[3][3]*pData_in[i+3];

                    double tau7=std::fabs(is0-is3);
                    double alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),p));
                    double alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),p));
                    double alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),p));
                    double alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),p));

                    double alphaSum=alpha0+alpha1+alpha2+alpha3;

                    double omega0=alpha0/alphaSum;
                    double omega1=alpha1/alphaSum;
                    double omega2=alpha2/alphaSum;
                    double omega3=alpha3/alphaSum;                

                    double theta=1.0/(1.0+std::pow(alphaSum-1.0,p));

                    if((i!=ibeg)&&(i!=(iend-1)))
                    {

                        if(theta>0.5)
                        {
                            a[i-ibeg]=2.0/7.0;
                            b[i-ibeg]=4.0/7.0;
                            c[i-ibeg]=1.0/7.0;
                            d[i-ibeg]=q7;                      

                        }
                        else
                        {
                            a[i-ibeg]=coeff_upwindcompact_c[0]*theta;
                            b[i-ibeg]=1.0+(coeff_upwindcompact_c[1]-1.0)*theta;
                            c[i-ibeg]=coeff_upwindcompact_c[2]*theta;
                            d[i-ibeg]=(omega0*q40+omega1*q41+omega2*q42+omega3*q43)*(1.0-theta)+q7*theta;
                        }

                    }
                    else
                    {
                        a[i-ibeg]=0.0;
                        b[i-ibeg]=1.0;
                        c[i-ibeg]=0.0;
                        d[i-ibeg]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
                    }
                
                
            }
            
            
            for(int i=1;i<(iend-ibeg);i++)
            {
                c[i]=c[i]/(b[i]-c[i-1]*a[i]);
                d[i]=(d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i]);
            }
            
            pData_out[iend]=d[iend-ibeg-1];
            
            for(int i=(iend-ibeg-2);i>=0;i--)
            {
                pData_out[i+ibeg+1]=d[i]-c[i]*pData_out[i+ibeg+2];
            }
            
        }
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
}

void nuc3d::hccs77::hccs77n(const Field & fieldIN,
    Field & fieldOUT,
    const int tilesize)
{

    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int ibeg=(tilesize-1);
    int iend=nx-tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double *a=new double[iend-ibeg];
    double *b=new double[iend-ibeg];
    double *c=new double[iend-ibeg];
    double *d=new double[iend-ibeg];

    
    for(int k=kbeg;k<kend;k++)
    {
        double *pData_in_k=pIn+nx*ny*k;
        double *pData_out_k=pOut+nx0*ny0*k;
        
        for(int j=jbeg;j<jend;j++)
        {
            double *pData_in=pData_in_k+nx*j;
            double *pData_out=pData_out_k+nx0*j;

            for(int i=ibeg;i<iend;i++)
            {       
                    double is0=pData_in[i+4]*(547.0*pData_in[i+4]-3882.0*pData_in[i+3]+4642.0*pData_in[i+2]-1854.0*pData_in[i+1])
                    +pData_in[i+3]*(7043.0*pData_in[i+3]-17246.0*pData_in[i+2]+7042.0*pData_in[i+1])
                    +pData_in[i+2]*(11003.0*pData_in[i+2]-9402.0*pData_in[i+1])
                    +2107.0*pData_in[i+1]*pData_in[i+1];
                    
                    double q40=coeff_weno7_alpha[0][0]*pData_in[i+4]+coeff_weno7_alpha[0][1]*pData_in[i+3]
                    +coeff_weno7_alpha[0][2]*pData_in[i+2]+coeff_weno7_alpha[0][3]*pData_in[i+1];                    

                    double is1=pData_in[i+3]*(267.0*pData_in[i+3]-1642.0*pData_in[i+2]+1602.0*pData_in[i+1]-494.0*pData_in[i])
                    +pData_in[i+2]*(2843.0*pData_in[i+2]-5966.0*pData_in[i+1]+1922.0*pData_in[i])
                    +pData_in[i+1]*(3443.0*pData_in[i+1]-2522.0*pData_in[i])
                    +547.0*pData_in[i]*pData_in[i];
                    
                    double q41=coeff_weno7_alpha[1][0]*pData_in[i+3]+coeff_weno7_alpha[1][1]*pData_in[i+2]
                    +coeff_weno7_alpha[1][2]*pData_in[i+1]+coeff_weno7_alpha[1][3]*pData_in[i];   

                    double is2=pData_in[i+2]*(547.0*pData_in[i+2]-2522.0*pData_in[i+1]+1922.0*pData_in[i]-494.0*pData_in[i-1])
                    +pData_in[i+1]*(3443.0*pData_in[i+1]-5966.0*pData_in[i]+1602.0*pData_in[i-1])
                    +pData_in[i]*(2843.0*pData_in[i]-1642.0*pData_in[i-1])
                    +267.0*pData_in[i-1]*pData_in[i-1];
                    
                    double q42=coeff_weno7_alpha[2][0]*pData_in[i+2]+coeff_weno7_alpha[2][1]*pData_in[i+1]
                    +coeff_weno7_alpha[2][2]*pData_in[i]+coeff_weno7_alpha[2][3]*pData_in[i-1];

                    double is3=pData_in[i+1]*(2107.0*pData_in[i+1]-9402.0*pData_in[i]+7042.0*pData_in[i-1]-1854.0*pData_in[i-2])
                    +pData_in[i]*(11003.0*pData_in[i]-17246.0*pData_in[i-1]+4642.0*pData_in[i-2])
                    +pData_in[i-1]*(7043.0*pData_in[i-1]-3882.0*pData_in[i-2])
                    +547.0*pData_in[i-2]*pData_in[i-2];                             

                    double q43=coeff_weno7_alpha[3][0]*pData_in[i+1]+coeff_weno7_alpha[3][1]*pData_in[i]
                    +coeff_weno7_alpha[3][2]*pData_in[i-1]+coeff_weno7_alpha[3][3]*pData_in[i-2];  

                    double q7=(-pData_in[i+3]+19.0*pData_in[i+2]+239.0*pData_in[i+1]+159.0*pData_in[i]+4.0*pData_in[i-1])/420.0;                

                    double tau7=std::fabs(is0-is3);
                    double alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),p));
                    double alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),p));
                    double alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),p));
                    double alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),p));                
                    
                    double alphaSum=alpha0+alpha1+alpha2+alpha3;
                    
                    double omega0=alpha0/alphaSum;
                    double omega1=alpha1/alphaSum;
                    double omega2=alpha2/alphaSum;
                    double omega3=alpha3/alphaSum;
                    double theta=1.0/(1.0+std::pow(alphaSum-1.0,p));
                    
                    if((i!=ibeg)&&(i!=(iend-1)))
                    {
                        if(theta>0.5)
                        {
                            c[i-ibeg]=2.0/7.0;
                            b[i-ibeg]=4.0/7.0;
                            a[i-ibeg]=1.0/7.0;
                            d[i-ibeg]=q7;                        
                        }
                        else
                        {
                            c[i-ibeg]=coeff_upwindcompact_c[0]*theta;
                            b[i-ibeg]=1.0+(coeff_upwindcompact_c[1]-1.0)*theta;
                            a[i-ibeg]=coeff_upwindcompact_c[2]*theta;
                            d[i-ibeg]=(omega0*q40+omega1*q41+omega2*q42+omega3*q43)*(1.0-theta)+q7*theta;
                        }

                    }
                    else
                    {
                        a[i-ibeg]=0.0;
                        b[i-ibeg]=1.0;
                        c[i-ibeg]=0.0;
                        d[i-ibeg]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
                    }
            }

            
            
            for(int i=1;i<(iend-ibeg);i++)
            {
                c[i]=c[i]/(b[i]-c[i-1]*a[i]);
                d[i]=(d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i]);
            }
            
            pData_out[iend]=d[iend-ibeg-1];
            
            for(int i=(iend-ibeg-2);i>=0;i--)
            {
                pData_out[i+ibeg+1]=d[i]-c[i]*pData_out[i+ibeg+2];
            }
        }
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
}

void nuc3d::hccs77::interpolationBoundaryL(const Field & fieldIN,
   const Field & boundaryL,
                                         const int uw,//uw=(1,-1)
                                         Field & fieldOUT,
                                         const int tilesize)
{
    switch(uw)
    {
        case 1:
        hccs77pBL(fieldIN,boundaryL, fieldOUT, tilesize);
        break;
        case -1:
        hccs77nBL(fieldIN,boundaryL, fieldOUT, tilesize);
        break;
        default:
        std::cout<<"weno error: no such direction"<<std::endl;
        exit(-1);
    }
    
}

void nuc3d::hccs77::hccs77pBL(const Field & fieldIN,
  const Field & boundaryL,
  Field & fieldOUT,
  const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryL.getSizeX();
    int nyBND=boundaryL.getSizeY();
    int nzBND=boundaryL.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                for(int z=-3;z<=3;z++)
                {
                    if((i+z-1)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i+z-1+nxBND);
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z-1;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                    is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                    +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                    +f[2]*(11003.0*f[2]-9402.0*f[3])
                    +2107.0*f[3]*f[3];
                    
                    is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                    +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                    +f[3]*(3443.0*f[3]-2522.0*f[4])
                    +547.0*f[4]*f[4];
                    
                    
                    is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                    +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                    +f[4]*(2843.0*f[4]-1642.0*f[5])
                    +267.0*f[5]*f[5];
                    
                    is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                    +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                    +f[5]*(7043.0*f[5]-3882.0*f[6])
                    +547.0*f[6]*f[6];
                    
                    q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                    +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                    
                    q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                    +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                    
                    q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                    +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                    
                    q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                    +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                    
                    tau7=std::fabs(is0-is3);
                    alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),p));
                    alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),p));
                    alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),p));
                    alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),p));                
                    
                    alphaSum=alpha0+alpha1+alpha2+alpha3;
                    
                    omega0=alpha0/alphaSum;
                    omega1=alpha1/alphaSum;
                    omega2=alpha2/alphaSum;
                    omega3=alpha3/alphaSum;
                    
                    pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;   
                
            }
        }
    }
}

void nuc3d::hccs77::hccs77nBL(const Field & fieldIN,
  const Field & boundaryL,
  Field & fieldOUT,
  const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryL.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryL.getSizeX();
    int nyBND=boundaryL.getSizeY();
    int nzBND=boundaryL.getSizeZ();
    
    int ibeg=0;
    int iend=tilesize;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i;
                
                for(int z=-3;z<=3;z++)
                {
                    if((i-z)<0)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+(i-z+nxBND);
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                
                    is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                    +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                    +f[2]*(11003.0*f[2]-9402.0*f[3])
                    +2107.0*f[3]*f[3];
                    
                    is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                    +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                    +f[3]*(3443.0*f[3]-2522.0*f[4])
                    +547.0*f[4]*f[4];
                    
                    
                    is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                    +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                    +f[4]*(2843.0*f[4]-1642.0*f[5])
                    +267.0*f[5]*f[5];
                    
                    is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                    +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                    +f[5]*(7043.0*f[5]-3882.0*f[6])
                    +547.0*f[6]*f[6];
                    
                    q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                    +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                    
                    q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                    +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                    
                    q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                    +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                    
                    q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                    +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                    
                    tau7=std::fabs(is0-is3);
                    alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),p));
                    alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),p));
                    alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),p));
                    alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),p));                
                    
                    alphaSum=alpha0+alpha1+alpha2+alpha3;
                    
                    omega0=alpha0/alphaSum;
                    omega1=alpha1/alphaSum;
                    omega2=alpha2/alphaSum;
                    omega3=alpha3/alphaSum;
                    
                    pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}

void nuc3d::hccs77::interpolationBoundaryR(const Field & fieldIN,
   const Field & boundaryR,
                                         const int uw,//uw=(1,-1)
                                         Field & fieldOUT,
                                         const int tilesize)
{
    switch(uw)
    {
        case 1:
        hccs77pBR(fieldIN,boundaryR, fieldOUT, tilesize);
        break;
        case -1:
        hccs77nBR(fieldIN,boundaryR, fieldOUT, tilesize);
        break;
        default:
        std::cout<<"weno error: no such direction"<<std::endl;
        exit(-1);
    }
    
}


void nuc3d::hccs77::hccs77pBR(const Field & fieldIN,
  const Field & boundaryR,
  Field & fieldOUT,
  const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryR.getSizeX();
    int nyBND=boundaryR.getSizeY();
    int nzBND=boundaryR.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                for(int z=-3;z<=3;z++)
                {
                    if((i+z)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i+z-nx;
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i+z;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                
                
                    is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                    +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                    +f[2]*(11003.0*f[2]-9402.0*f[3])
                    +2107.0*f[3]*f[3];
                    
                    is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                    +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                    +f[3]*(3443.0*f[3]-2522.0*f[4])
                    +547.0*f[4]*f[4];
                    
                    
                    is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                    +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                    +f[4]*(2843.0*f[4]-1642.0*f[5])
                    +267.0*f[5]*f[5];
                    
                    is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                    +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                    +f[5]*(7043.0*f[5]-3882.0*f[6])
                    +547.0*f[6]*f[6];
                    
                    q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                    +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                    
                    q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                    +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                    
                    q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                    +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                    
                    q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                    +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                    
                    tau7=std::fabs(is0-is3);
                    alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),p));
                    alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),p));
                    alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),p));
                    alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),p));
                    
                    
                    alphaSum=alpha0+alpha1+alpha2+alpha3;
                    
                    omega0=alpha0/alphaSum;
                    omega1=alpha1/alphaSum;
                    omega2=alpha2/alphaSum;
                    omega3=alpha3/alphaSum;
                    
                    pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}

void nuc3d::hccs77::hccs77nBR(const Field & fieldIN,
  const Field & boundaryR,
  Field & fieldOUT,
  const int tilesize)
{
    double *pIn=fieldIN.getDataPtr();
    double *pOut=fieldOUT.getDataPtr();
    double *pBND=boundaryR.getDataPtr();
    
    double omega0,omega1,omega2,omega3;
    double alpha0,alpha1,alpha2,alpha3,alphaSum;
    double is0,is1,is2,is3;
    double q40,q41,q42,q43;
    double tau7;
    
    int nx=fieldIN.getSizeX();
    int ny=fieldIN.getSizeY();
    int nz=fieldIN.getSizeZ();
    
    int nx0=fieldOUT.getSizeX();
    int ny0=fieldOUT.getSizeY();
    int nz0=fieldOUT.getSizeZ();
    
    int nxBND=boundaryR.getSizeX();
    int nyBND=boundaryR.getSizeY();
    int nzBND=boundaryR.getSizeZ();
    
    int ibeg=nx-tilesize;
    int iend=nx;
    int jbeg=0;
    int jend=ny;
    int kbeg=0;
    int kend=nz;
    
    double f[7];
    
    for(int k=kbeg;k<kend;k++)
    {
        for(int j=jbeg;j<jend;j++)
        {
            for(int i=ibeg;i<iend;i++)
            {
                int idx_rf=nx0*ny0*k+nx0*j+i+1;
                
                for(int z=-3;z<=3;z++)
                {
                    if((i-z+1)>=nx)
                    {
                        int idx_BND=nxBND*nyBND*k+nxBND*j+i-z+1-nx;
                        
                        f[3+z]=pBND[idx_BND];
                    }
                    else
                    {
                        int idx_f=nx*ny*k+nx*j+i-z+1;
                        f[3+z]=pIn[idx_f];
                    }
                }
                
                    is0=f[0]*(547.0*f[0]-3882.0*f[1]+4642.0*f[2]-1854.0*f[3])
                    +f[1]*(7043.0*f[1]-17246.0*f[2]+7042.0*f[3])
                    +f[2]*(11003.0*f[2]-9402.0*f[3])
                    +2107.0*f[3]*f[3];
                    
                    is1=f[1]*(267.0*f[1]-1642.0*f[2]+1602.0*f[3]-494.0*f[4])
                    +f[2]*(2843.0*f[2]-5966.0*f[3]+1922.0*f[4])
                    +f[3]*(3443.0*f[3]-2522.0*f[4])
                    +547.0*f[4]*f[4];
                    
                    
                    is2=f[2]*(547.0*f[2]-2522.0*f[3]+1922.0*f[4]-494.0*f[5])
                    +f[3]*(3443.0*f[3]-5966.0*f[4]+1602.0*f[5])
                    +f[4]*(2843.0*f[4]-1642.0*f[5])
                    +267.0*f[5]*f[5];
                    
                    is3=f[3]*(2107.0*f[3]-9402.0*f[4]+7042.0*f[5]-1854.0*f[6])
                    +f[4]*(11003.0*f[4]-17246.0*f[5]+4642.0*f[6])
                    +f[5]*(7043.0*f[5]-3882.0*f[6])
                    +547.0*f[6]*f[6];
                    
                    q40=coeff_weno7_alpha[0][0]*f[0]+coeff_weno7_alpha[0][1]*f[1]
                    +coeff_weno7_alpha[0][2]*f[2]+coeff_weno7_alpha[0][3]*f[3];
                    
                    q41=coeff_weno7_alpha[1][0]*f[1]+coeff_weno7_alpha[1][1]*f[2]
                    +coeff_weno7_alpha[1][2]*f[3]+coeff_weno7_alpha[1][3]*f[4];
                    
                    q42=coeff_weno7_alpha[2][0]*f[2]+coeff_weno7_alpha[2][1]*f[3]
                    +coeff_weno7_alpha[2][2]*f[4]+coeff_weno7_alpha[2][3]*f[5];
                    
                    q43=coeff_weno7_alpha[3][0]*f[3]+coeff_weno7_alpha[3][1]*f[4]
                    +coeff_weno7_alpha[3][2]*f[5]+coeff_weno7_alpha[3][3]*f[6];
                    
                    tau7=std::fabs(is0-is3);
                    alpha0=coeff_weno7_c[0]*(1.0+std::pow(tau7/(ss+is0),p));
                    alpha1=coeff_weno7_c[1]*(1.0+std::pow(tau7/(ss+is1),p));
                    alpha2=coeff_weno7_c[2]*(1.0+std::pow(tau7/(ss+is2),p));
                    alpha3=coeff_weno7_c[3]*(1.0+std::pow(tau7/(ss+is3),p));
                    
                    
                    alphaSum=alpha0+alpha1+alpha2+alpha3;
                    
                    omega0=alpha0/alphaSum;
                    omega1=alpha1/alphaSum;
                    omega2=alpha2/alphaSum;
                    omega3=alpha3/alphaSum;
                    
                    pOut[idx_rf]=omega0*q40+omega1*q41+omega2*q42+omega3*q43;
            }
        }
    }
}
