//
//  block.cpp
//  NUC3d
//
//  Created by Jun Peng on 15/10/20.
//  Copyright © 2015年 Jun Peng. All rights reserved.
//

#include "block.h"
#include "euler3d.h"
#include "eulerReactive3d.h"
#include "NaiverStokes3d.h"
#include "NaiverStokesReactive3d.h"
#include "physicsModel.h"
#include "MPICommunicator.h"
#include "IOcontroller.h"
#include "postproc.hpp"
#include "fft.h"
#include "utilities.h"
#include <cstdlib>
#include <cmath>
#include <random>



nuc3d::block::block():
time(0.0),
dt(0.025),
istep(0),
RES(1.0)
{}

nuc3d::block::~block()
{}

void nuc3d::block::initial(fieldOperator3d &myOP,
                           physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI,
                           IOController &myIO)

{
    
    std::string filename_block = ("inp/block.in");
    
    std::ifstream myFile;
    myFile.open(filename_block);
    if (myFile)
    {
        myFile >> npx >> npy >> npz;
        myMPI.setTopo(npx,npy,npz);
        if(myMPI.getMyId()==0) std::cout<<"MPI process topo:"<<npx<<","<<npy<<","<<npz<<std::endl;

        myFile >> nx_global >> ny_global >> nz_global;
        if(myMPI.getMyId()==0) std::cout<<"Global mesh size:"<<nx_global<<","<<ny_global<<","<<nz_global<<std::endl;

        nx=nx_global/npx;
        ny=ny_global/npy;
        nz=nz_global/npz;

        myFile >> xl >> yl >> zl;
        myFile >> x_origin >> y_origin >> z_origin;
        dx=xl/nx_global;
        dy=yl/ny_global;
        dz=zl/nz_global;
        if(myMPI.getMyId()==0) std::cout<<"Geometry size:"<<xl<<","<<yl<<","<<zl<<std::endl;
        if(myMPI.getMyId()==0) std::cout<<"Geometry delta size:"<<dx<<","<<dy<<","<<dz<<std::endl;
        if(myMPI.getMyId()==0) std::cout<<"local mesh size:"<<nx<<","<<ny<<","<<nz<<std::endl;

        for(int i=0;i<3;i++)
        {
            xyz_center.push_back(Field(nx,ny,nz));
            xyz.push_back(Field(nx+1,ny+1,nz+1));
        }        

        if(myMPI.getMyId()==0) std::cout<<"Start calculating mesh data..."<<std::endl;
        getXYZ(myMPI.getMyCoord(2),myMPI.getMyCoord(1),myMPI.getMyCoord(0));

        if(myMPI.getMyId()==0) std::cout<<"Center data has been calculated..."<<std::endl;
    }
    else
    {
        std::cout<<"File \'"<<filename_block<<"\' does not exist!"<<std::endl;
        exit(-1);
    }
    myFile.close();
    
    initialData(myPhyMod,myOP);    
    myFFT->initFFT(myMPI,nx,ny,nz);

    if(myMPI.getMyId()==0) std::cout<<"Flow field has been allocated!"<<std::endl;

    if(0==(myIO.getStep("startStep")))
        initialQ(myIO,myPhyMod,myMPI);
    else
        inputQ_binary(myMPI.getMyId(),myIO.getStep("startStep"));
    
    if(myMPI.getMyId()==0) std::cout<<"Flow field has been initialized!"<<std::endl;
    
    
}


void nuc3d::block::getXYZ(int coord_x,int coord_y, int coord_z)
{
    double *x_center=xyz_center[0].getDataPtr();
    double *y_center=xyz_center[1].getDataPtr();
    double *z_center=xyz_center[2].getDataPtr();

    double *x=xyz[0].getDataPtr();
    double *y=xyz[1].getDataPtr();
    double *z=xyz[2].getDataPtr();

    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int idx=nx*ny*k+nx*j+i;
                x_center[idx]=xl/npx*coord_x+dx*(0.5+i)-x_origin;
                y_center[idx]=yl/npy*coord_y+dy*(0.5+j)-y_origin;
                z_center[idx]=zl/npz*coord_z+dz*(0.5+k)-z_origin;
            }
        }
    }

    for(int k=0;k<nz+1;k++)
    {
        for(int j=0;j<ny+1;j++)
        {
            for(int i=0;i<nx+1;i++)
            {
                int idx=(nx+1)*(ny+1)*k+(nx+1)*j+i;
                x[idx]=xl/npx*coord_x+dx*i-x_origin;
                y[idx]=yl/npy*coord_y+dy*j-y_origin;
                z[idx]=zl/npz*coord_z+dz*k-z_origin;
            }
        }
    }
}


void nuc3d::block::initialData(physicsModel &myPhy,fieldOperator3d &myOP)
{
    myPDE.initPDEData3d(nx, ny, nz, myPhy.getEqNum());
    if("Euler3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<EulerData3D>(nx,ny,nz,myPhy.getEqNum(),dx,dy,dz);
    }
    else if ("EulerReactive3d"==myPhy.getMyModelName())
    {
        //        myFluxes=std::make_shared<EulerReactiveData3D>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else if ("NavierStokes3d"==myPhy.getMyModelName())
    {
        myFluxes=std::make_shared<NaiverStokesData3d>(nx,ny,nz,myPhy.getEqNum(),dx,dy,dz);
    }
    else if ("NaiverStokesReactive3d"==myPhy.getMyModelName())
    {
        //        myFluxes=std::make_shared<NaiverStokesReactiveData3d>(nx0,ny0,nz0,myPhy.getEqNum());
    }
    else
    {
        std::cout <<"Model Name:"<<myPhy.getMyModelName()
        <<"does not exist!"
        <<std::endl;
        exit(-1);
    }
    
    bfsize=std::fmax(myOP.bufferSize_xi,myOP.bufferSize_eta);
    bfsize=std::fmax(bfsize,myOP.bufferSize_zeta);
    for(int n=0;n<myPhy.getEqNum();n++)
    {
        mybuffer.push_back(bufferData(nx,ny,nz,bfsize));
        OutPutValue_prim.push_back(Field(nx,ny,nz));
        fft_re.push_back(Field(nx,ny,nz));
        fft_im.push_back(Field(nx,ny,nz));
    }
    
    for (int n=0; n<3; n++)
    {
        OutPutValue_acoust.push_back(Field(nx,ny,nz));
    }
    
    myPost=std::make_shared<postproc>(nx,ny,nz,dx,dy,dz);
    myFFT=std::make_shared<fft>();
}



void nuc3d::block::solve(fieldOperator3d &myOP,
                         physicsModel &myPhyMod,
                         MPIComunicator3d_nonblocking &myMPI,
                         IOController &myIO)
{
    int step=0;
    double t0=MPI_Wtime();
    while (myOP.getSteps()!=step)
    {
        myFluxes->solve(myPDE, myOP, mybuffer, myPhyMod, myMPI, *myFFT,xyz_center);
        myPDE.solve(myOP, myMPI,myIO.myTimeController["cfl"],step++);
    }
    double t1=MPI_Wtime();
    
    wall_time=t1-t0;
    istep++;
    
    dt=myPDE.getDt();
    time+=dt;
    
    myIO.myTimeController["dt"]=dt;
    
    
    RES=myPDE.getRes();
}

void nuc3d::block::printStatus()
{
    std::cout<<std::setprecision(6)<<"========step = "<<istep<< "\n time = "<<time<<", dt = "<<dt<<", residual = "<<RES<<", CPU time = "<<wall_time<<"(s) "<<std::endl;
}

void nuc3d::block::Post(fieldOperator3d &myOp,
                        physicsModel &myPhys,
                        MPIComunicator3d_nonblocking &myMPI,
                        IOController &myIO)
{
    myPhys.getPrim(myPDE.getQ(), OutPutValue_prim, OutPutValue_acoust);
    
    myPost->solvePost(OutPutValue_prim,
                      OutPutValue_acoust,
                      xyz,
                      myPhys,
                      myOp, mybuffer, myMPI, myIO, *myFFT,istep,time);

    
}

void nuc3d::block::Output(fieldOperator3d &myOp,
                          physicsModel &myPhys,
                          MPIComunicator3d_nonblocking &myMPI,
                          IOController &myIO)
{
    myPhys.getPrim(myPDE.getQ(), OutPutValue_prim, OutPutValue_acoust);
    
    if(("yes")==(myIO.getType("Binary")))
    {
        outputQ_binary(myMPI.getMyId(),myPhys);
    }
    
    if(("yes")==(myIO.getType("Tecplot")))
    {
        outputQ_tecplot(myMPI.getMyId(),myPhys);
        
    }
    
    if(("yes")==(myIO.getType("PostProc")))
    {
        myPost->OutputPost(OutPutValue_prim,
                           OutPutValue_acoust,
                           xyz,
                           myOp, mybuffer, myMPI, myIO,istep,time);
        
    }

    
}

void nuc3d::block::initialQ(IOController &myIO,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI)
{
    int fp=myIO.getStep("Benchmark");

    VectorField &Q=myPDE.getQ();

    double *prho=Q[0].getDataPtr();
    double *prhou=Q[1].getDataPtr();
    double *prhov=Q[2].getDataPtr();
    double *prhow=Q[3].getDataPtr();
    double *prhoe=Q[4].getDataPtr();
        
    (this->*myInitial[fp])(prho,prhou,prhov,prhow,prhoe,myPhyMod,myMPI);

}

void nuc3d::block::initial_default(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI)
{
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();
    double u_init=myPhyMod.getValue("u_init");
    double v_init=myPhyMod.getValue("v_init");
    double w_init=myPhyMod.getValue("w_init");
    double *px=xyz_center[0].getDataPtr();
    double *py=xyz_center[1].getDataPtr();
    double *pz=xyz_center[2].getDataPtr();

    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int idx=nx*ny*k+nx*j+i;

                rho[idx]=1.0;
                rhou[idx]=u_init;
                rhov[idx]=v_init;
                rhow[idx]=w_init;
                rhoe[idx]=1.0/(mach*mach*gamma)/(gamma-1.0);
            }
        }
    }
    
}

void nuc3d::block::initial_ivc(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI)
{
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();
                
    double pie=4.0*atan(1.0);
    double b=0.5;
    double x_c=5.0;
    double y_c=5.0;    

    double *px=xyz_center[0].getDataPtr();
    double *py=xyz_center[1].getDataPtr();
    double *pz=xyz_center[2].getDataPtr();

    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int idx=nx*ny*k+nx*j+i;
                
                double r=std::sqrt(std::pow(px[idx]-x_c,2)+std::pow(py[idx]-y_c,2));
                double rho0=std::pow((1-(gamma-1.0)*b*b*std::exp(1-r*r)/(8.0*gamma*pie*pie)),2.5);
                double u0=0.5-b/(2.0*pie)*std::exp((1-r*r)/2)*(py[idx]-y_c);
                double v0=b/(2.0*pie)*std::exp((1-r*r)/2)*(px[idx]-x_c);
                double w0=0.0;
                double p0=std::pow(rho0,gamma);

                rho[idx]=rho0;
                rhou[idx]=rho0*u0;
                rhov[idx]=rho0*v0;
                rhow[idx]=rho0*w0;
                rhoe[idx]=p0/(gamma-1.0)+0.5*rho0*(u0*u0+v0*v0+w0*w0);
            }
        }
    }
    
}

void nuc3d::block::initial_taylorgreen(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI)
{
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();
    double *px=xyz_center[0].getDataPtr();
    double *py=xyz_center[1].getDataPtr();
    double *pz=xyz_center[2].getDataPtr();

    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int idx=nx*ny*k+nx*j+i;

                double x0=px[idx];
                double y0=py[idx];
                double z0=pz[idx];

                double rho0=1.0;
                double u0=std::sin(x0)*std::cos(y0)*std::cos(z0);
                double v0=-std::cos(x0)*std::sin(y0)*std::cos(z0);
                double w0=0.0;
                double p0=1.0/(gamma*mach*mach)+((std::cos(2.0*z0)+2.0)*(std::cos(2.0*x0)+std::cos(2.0*y0))-2.0)/16.0;

                rho[idx]=rho0;
                rhou[idx]=rho0*u0;
                rhov[idx]=rho0*v0;
                rhow[idx]=rho0*w0;
                rhoe[idx]=p0/(gamma-1.0)+0.5*rho0*(u0*u0+v0*v0+w0*w0);
            }
        }
    }
}

void nuc3d::block::initial_hit(double *rho,double *rhou,double *rhov,double *rhow,
            double *rhoe,physicsModel &myPhyMod,
                           MPIComunicator3d_nonblocking &myMPI)
{
    /*
        This function generate random phase initial condition
        with energy spectrum:

            E(k)=16*sqrt(2/PI)*(k/k_p)^4*U_0/k_p*exp(-2*(k/k_p)^2)

                            or 

                    E(k)=A*k*U_0^2*exp(-2*(k/k_p)^2)

    */
    double gamma=myPhyMod.getGamma();
    double mach=myPhyMod.getMach();
    double pie=4.0*atan(1.0);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 2.0*pie);

    double k_peak=myPhyMod.getValue("k_peak");
    double mach_t=myPhyMod.getValue("Mach_t");
    double u_rms_init=mach_t/std::pow(mach,2)/std::sqrt(3.0);

    if(0==myMPI.getMyId())
        std::cout<<"Initial u_rms_init = "<< u_rms_init<<" "<<std::endl;

    int npx=myMPI.getMyDim(2);
    int npy=myMPI.getMyDim(1);
    int npz=myMPI.getMyDim(0);

    int px=myMPI.getMyCoord(2);
    int py=myMPI.getMyCoord(1);
    int pz=myMPI.getMyCoord(0);

    VectorField vel_re;
    VectorField vel_im;
    VectorField FFT_vector_re;
    VectorField FFT_vector_im;

    for(int i=0;i<3;i++)
    {
        vel_re.push_back(Field(nx,ny,nz));
        vel_im.push_back(Field(nx,ny,nz));
        FFT_vector_re.push_back(Field(nx,ny,nz));
        FFT_vector_im.push_back(Field(nx,ny,nz));
    }

    double *pu_re=FFT_vector_re[0].getDataPtr();
    double *pv_re=FFT_vector_re[1].getDataPtr();
    double *pw_re=FFT_vector_re[2].getDataPtr();

    double *pu_im=FFT_vector_im[0].getDataPtr();
    double *pv_im=FFT_vector_im[1].getDataPtr();
    double *pw_im=FFT_vector_im[2].getDataPtr();

    double kx=(double)nx*(double)px;
    double ky=(double)ny*(double)py;
    double kz=(double)nz*(double)pz;
    
    int nx_global=npx*nx;
    int ny_global=npy*ny;
    int nz_global=npz*nz;

    for (int k=0; k<nz; k++)
    {
        for (int j=0; j<ny; j++)
        {
            for (int i=0; i<nx; i++)
            {
                int idx=nx*ny*k+nx*j+i;

                double kx0=((kx+(double)i)>=(nx_global/2.0))?kx+(double)i-nx_global:kx+(double)i;
                double ky0=((ky+(double)j)>=(ny_global/2.0))?ky+(double)j-ny_global:ky+(double)j;
                double kz0=((kz+(double)k)>=(nz_global/2.0))?kz+(double)k-nz_global:kz+(double)k;
                double kmag=sqrt(kx0*kx0+ky0*ky0+kz0*kz0);
                double k12=sqrt(kx0*kx0+ky0*ky0);

                if(kmag<1e-6) kmag=1.0;

                double Ek_0=std::pow(u_rms_init,2)*16.0*std::sqrt(2.0/pie)*std::pow(kmag/k_peak,4)/k_peak*std::exp(-2.0*kmag*kmag/k_peak/k_peak);
                double phi1=dis(gen);
                double phi2=dis(gen);
                double phi3=dis(gen);

                double a_re=std::sqrt(2.0*Ek_0/kmag/kmag/4.0/pie)*std::cos(phi1)*std::cos(phi3);
                double a_im=std::sqrt(2.0*Ek_0/kmag/kmag/4.0/pie)*std::sin(phi1)*std::cos(phi3);
                double b_re=std::sqrt(2.0*Ek_0/kmag/kmag/4.0/pie)*std::cos(phi2)*std::sin(phi3);
                double b_im=std::sqrt(2.0*Ek_0/kmag/kmag/4.0/pie)*std::sin(phi2)*std::sin(phi3);

                double k212=(k12<1e-6)?1.0:ky0/k12;
                double k112=(k12<1e-6)?0.0:kx0/k12;

                pu_re[idx]=k212*a_re+k112*kz0/kmag*b_re;
                pv_re[idx]=k212*kz0/kmag*b_re-k112*a_re;                 
                pw_re[idx]=-k12/kmag*b_re;

                pu_im[idx]=k212*a_im+k112*kz0/kmag*b_im;
                pv_im[idx]=k212*kz0/kmag*b_im-k112*a_im;                 
                pw_im[idx]=-k12/kmag*b_im;
            }
        }
    }

    myFFT->solveFFT_backward(FFT_vector_re[0],FFT_vector_im[0],vel_re[0],vel_im[0],myMPI); // FFT of velocity u
    myFFT->solveFFT_backward(FFT_vector_re[1],FFT_vector_im[1],vel_re[1],vel_im[1],myMPI); // FFT of velocity v
    myFFT->solveFFT_backward(FFT_vector_re[2],FFT_vector_im[2],vel_re[2],vel_im[2],myMPI); // FFT of velocity w

    double *pu=vel_re[0].getDataPtr();
    double *pv=vel_re[1].getDataPtr();
    double *pw=vel_re[2].getDataPtr();

    double u2_re,v2_re,w2_re,vel;
    double u2_im,v2_im,w2_im;
    myMPI.globalAveragePower2(vel_re[0],u2_re);
    myMPI.globalAveragePower2(vel_re[1],v2_re);
    myMPI.globalAveragePower2(vel_re[2],w2_re);
    // myMPI.globalAveragePower2(vel_im[0],u2_im);
    // myMPI.globalAveragePower2(vel_im[1],v2_im);
    // myMPI.globalAveragePower2(vel_im[2],w2_im);
    
    vel=std::sqrt((u2_re+v2_re+w2_re)/3.0);
    double scale=u_rms_init/vel;

    if(0==myMPI.getMyId())
        std::cout<<"Calculated u_rms_init = "<< vel <<",scaled by  Scal= "<<scale<<std::endl;


    for(int k=0;k<nz;k++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int i=0;i<nx;i++)
            {
                int idx=nx*ny*k+nx*j+i;

                double rho0=1.0;
                double u0=pu[idx]*scale;
                double v0=pv[idx]*scale;
                double w0=pw[idx]*scale;
                double p0=1.0/(gamma*mach*mach);

                rho[idx]=rho0;
                rhou[idx]=rho0*u0;
                rhov[idx]=rho0*v0;
                rhow[idx]=rho0*w0;
                rhoe[idx]=p0/(gamma-1.0)+0.5*rho0*(u0*u0+v0*v0+w0*w0);
            }
        }
    }

}

void nuc3d::block::outputQ_tecplot(int myID,physicsModel &myPhys)
{
    std::string forename_flow = ("flowData/NUC3d_step_");
    std::string step;
    std::string mid("_id_");
    std::string id;
    std::string tailname = (".dat");
    
    std::stringstream ss_step,ss_id;
    ss_step << istep;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow);
    
    std::string TECplotHeader[2]={"title=NUC3d\n",
        "variables=x,y,z"};
    
    myIOfile<<TECplotHeader[0]
    <<TECplotHeader[1];
    
    for(int i=0;i<VarNameList.size();i++)
    {
        myIOfile<<","<<VarNameList[i];
    }
    
    myIOfile<<"\n Zone I = "<<nx+1<<", J= "<<ny+1<<", K="<<nz+1
    <<"\n DATAPACKING=BLOCK, VARLOCATION=(["<<xyz.size()+1<<"-"
    <<xyz.size()+OutPutValue_prim.size()+OutPutValue_acoust.size()
    <<"]=CELLCENTERED)\n";
    
    for(auto iter=xyz.begin();iter!=xyz.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=OutPutValue_prim.begin();iter!=OutPutValue_prim.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    for(auto iter=OutPutValue_acoust.begin();iter!=OutPutValue_acoust.end();iter++)
    {
        writeField(myIOfile, *iter);
    }
    
    myIOfile.close();
    
}

void nuc3d::block::outputQ_binary(int myID,physicsModel &myPhys)
{
    std::string forename_flow = ("flowData/NUC3d_step_");
    std::string step;
    std::string mid("_ID_");
    std::string id;
    std::string tailname = (".bin");
    
    std::stringstream ss_step,ss_id;
    ss_step << istep;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ofstream myIOfile;
    
    myIOfile.open(filename_flow,std::ios::out|std::ios::binary);
    VectorField &Q=myPDE.getQ();
    myIOfile.write(reinterpret_cast<char *>(&time), sizeof(time));
    for(auto iter=Q.begin();iter!=Q.end();iter++)
    {
	    writeField_binary(myIOfile,*iter);
    }
    
    myIOfile.close();
    
}

void nuc3d::block::inputQ_binary(int myID,int step0)
{
    std::string forename_flow = ("flowData/NUC3d_step_");
    std::string step;
    std::string mid("_ID_");
    std::string id;
    std::string tailname = (".bin");
    
    std::stringstream ss_step,ss_id;
    ss_step << step0;
    ss_step >> step;
    
    ss_id<<myID;
    ss_id>>id;
    
    std::string filename_flow = forename_flow + step + mid + id + tailname;
    
    std::ifstream myIOfile;
    
    myIOfile.open(filename_flow,std::ios::in|std::ios::binary);
    VectorField &Q=myPDE.getQ();
    
    myIOfile.read(reinterpret_cast<char *>(&time), sizeof (time));
    istep=step0;
    for(auto iter=Q.begin();iter!=Q.end();iter++)
    {
        readField_binary(myIOfile,*iter);
    }
    
    myIOfile.close();
    
}
