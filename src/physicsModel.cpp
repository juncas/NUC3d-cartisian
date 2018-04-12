#ifndef physicsModel_cpp
#define physicsModel_cpp
#include "physicsModel.h"
#include "euler3d.h"
#include "PDEData3d.hpp"
#include "fft.h"
#include <cmath>
/**************************************************************************************
 Definition of class communicator: base class for communicators
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
 nuc3d::physicsModel::physicsModel() :
 neqs(5),
 myRiemannMap(
 {
   { "AUSMp",&nuc3d::physicsModel::RiemannAUSMp },
   { "AUSM",&nuc3d::physicsModel::RiemannAUSM },
   { "LF",&nuc3d::physicsModel::RiemannLF}
}
),
 myModelParameters(
 {
  {"Reynolds",1000},
  {"Mach",0.5},
  {"Pt",0.72},
  {"Rossby",0.0},
  {"x_center",3.1415926},
  {"y_center",3.1415926},
  {"z_center",0.0},
  {"well_rad",1.0},
  {"Gamma",1.4},
  {"T_ref",110.4},
  {"T_inf",288.0},
  {"T_wall",1.0},
  {"Forcing_amplitude",0.01},
  {"Forcing_amplitude_hel",0.01},
  {"Forcing_wave_k",1.0},
  {"Mach_t",0.3},
  {"k_peak",8.0}
}
),
 myVisModelMap(
 {
    {"Sutherland",&nuc3d::physicsModel::sutherland},
    {"Constant",&nuc3d::physicsModel::constant}
}
),
 myForcingMap(
 {
    {"None",&nuc3d::physicsModel::AddForcing_None},
    {"Const",&nuc3d::physicsModel::AddForcing_Const},
    {"ABC",&nuc3d::physicsModel::AddForcing_ABC},
    {"TG",&nuc3d::physicsModel::AddForcing_TG}
}
),
 mySourceSwitchMap(
 {
  {"Coriolis","no"},
  {"Centrifugal","no"}
                      //add other source switchers..
}
)
 {
    std::string str;
    std::ifstream file("inp/PhysModel.in");
    //std::cout<<"Start reading PhysModel.in ..."<<std::endl;
    if (file)
    {
        // this order should not change
        file >> str >> neqs;
        file >> str >> myEoSName;
        file >> str >> myRiemannName;
        file >> str >> myModelName;
        file >> str >> myVisModelName;
        file >> str >> myForcingName;
        
        if(neqs!=5)
        {
            std::cout<<"WARNNING: Equation number in present version is 5!!"
            << std::endl;
            neqs=5;
        }
        
        
        
        if(myRiemannMap.find(myRiemannName)==myRiemannMap.end())
        {
            std::cout<<"Riemann Solver name "<<myEoSName<<" does not exist ,using default!"
            << std::endl;
            myEoSName="AUSM";
        }
        
        
        if(myVisModelMap.find(myVisModelName)==myVisModelMap.end())
        {
            std::cout<<"Viscosity model name "<<myVisModelName<<" does not exist ,using default!"
            << std::endl;
            myVisModelName="Sutherland";
            
        }

        if(myForcingMap.find(myForcingName)==myForcingMap.end())
        {
            std::cout<<"Forcing model name "<<myForcingName<<" does not exist ,using default: \"None\"!"
            << std::endl;
            myForcingName="None";
            
        }
        
        while (readPhysFile(file));
        //std::cout<<"PhysModel.in has been read!"<<std::endl;
    }
    else
    {
        std::cout << "IO file \'PhysModel.in\' does not exist!"
        << std::endl;
        exit(0);
    }
    file.close();
}

nuc3d::physicsModel::~physicsModel()
{

}
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
 std::istream& nuc3d::physicsModel::readPhysFile(std::istream& ios)
 {
    std::string word0;
    std::string word1;
    ios >> word0 >> word1;

    if (myModelParameters.find(word0) != myModelParameters.end())
    {
        double value;
        std::istringstream value_strm(word1);
        value_strm >> value;
        myModelParameters[word0] = value;
    }
    else if (mySourceSwitchMap.find(word0) != mySourceSwitchMap.end())
    {
        mySourceSwitchMap[word0] = word1;
    }
    else if(word0.size())
    {
        std::cout << "word " << word0 << " does not exist"<<std::endl;
        exit(0);
    }
    
    return ios;
    
}

void nuc3d::physicsModel::solve(PDEData3d &myPDE,EulerData3D *myEuler)
{
    con2prim(myEoSName,             
       myPDE.getQ(),
       myEuler->W_Euler,
       myEuler->W0_Euler);
    
}

void nuc3d::physicsModel::getPrim(VectorField &Q,VectorField &Prim, VectorField &Acoust)
{
    con2prim(myEoSName,
       Q,
       Prim,
       Acoust);
    
}


void nuc3d::physicsModel::getMiu(Field &T,
   Field &miu,
   Field &coeff)

{
    (this->*myVisModelMap[myVisModelName])(T,
     miu,
     coeff,
     myModelParameters["Reynolds"],
     myModelParameters["Mach"],
     myModelParameters["Pt"],
     myModelParameters["Gamma"],
     myModelParameters["T_ref"],
     myModelParameters["T_inf"]);
}

void nuc3d::physicsModel::sutherland(const Field &T,
   Field &miu,
   Field &coeff,
   double Reynolds,
   double Mach,
   double pt,
   double gamma,
   double T_ref,
   double T_inf)
{
    int nx=T.getSizeX();
    int ny=T.getSizeY();
    int nz=T.getSizeZ();
    
    double *pT=T.getDataPtr();
    double *pMiu=miu.getDataPtr();
    double *pCoeff=coeff.getDataPtr();
    
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                int idx=nx*ny*k+nx*j+i;
                double T_local=pT[idx];
                double non_dim_T_ref=T_ref/T_inf;
                double miu0=(1.0+non_dim_T_ref)/(T_local+non_dim_T_ref)*std::pow(T_local,1.5)/Reynolds;
                double coeff0=miu0/((gamma-1.0)*Mach*Mach*pt);
                
                pMiu[idx]=miu0;
                pCoeff[idx]=coeff0;
            }
        }
    }
}

void nuc3d::physicsModel::constant(const Field &T,
 Field &miu,
 Field &coeff,
 double Reynolds,
 double Mach,
 double pt,
 double gamma,
 double T_ref,
 double T_inf)
{
    int nx=T.getSizeX();
    int ny=T.getSizeY();
    int nz=T.getSizeZ();
    
    double *pT=T.getDataPtr();
    double *pMiu=miu.getDataPtr();
    double *pCoeff=coeff.getDataPtr();
    
    for (int k = 0; k < nz; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                int idx=nx*ny*k+nx*j+i;
                double T_local=pT[idx];
                double non_dim_T_ref=T_ref/T_inf;
                double miu0=1.0/Reynolds;
                double coeff0=miu0/((gamma-1.0)*Mach*Mach*pt);
                
                pMiu[idx]=miu0;
                pCoeff[idx]=coeff0;
            }
        }
    }
    
}

void nuc3d::physicsModel::initial(PDEData3d &myPDE,std::shared_ptr<EulerData3D> myEuler)
{
    prim2con(myEoSName,
       myEuler->W_Euler,
       myPDE.getQ());
}


void nuc3d::physicsModel::con2prim(const std::string &EoSName,
 VectorField &Q_vec,
 VectorField &W_vec,
 VectorField &W0_vec
 )
{
    double rho, u, v, w, E, p;
    double T, e, alpha;    

    auto iterRho = Q_vec.begin();
    auto iterRhoU = Q_vec.begin() + 1;
    auto iterRhoV = Q_vec.begin() + 2;
    auto iterRhoW = Q_vec.begin() + 3;
    auto iterRhoE = Q_vec.begin() + 4;
    
    auto iterRho0 = W_vec.begin();
    auto iterU0   = W_vec.begin() + 1;
    auto iterV0   = W_vec.begin() + 2;
    auto iterW0   = W_vec.begin() + 3;
    auto iterP0   = W_vec.begin() + 4;
    
    auto iterT0 = W0_vec.begin();
    auto iterE0   = W0_vec.begin() + 1;
    auto iterAlpha0   = W0_vec.begin() + 2;
    
    double *pRho=iterRho->getDataPtr();
    double *pRhoU=iterRhoU->getDataPtr();
    double *pRhoV=iterRhoV->getDataPtr();
    double *pRhoW=iterRhoW->getDataPtr();
    double *pRhoE=iterRhoE->getDataPtr();
    
    double *pRho0=iterRho0->getDataPtr();
    double *pU0=iterU0->getDataPtr();
    double *pV0=iterV0->getDataPtr();
    double *pW0=iterW0->getDataPtr();
    double *pP0=iterP0->getDataPtr();
    
    double *pT0=iterT0->getDataPtr();
    double *pE0=iterE0->getDataPtr();
    double *pAlpha0=iterAlpha0->getDataPtr();

    double gamma = myModelParameters["Gamma"];
    double Mach = myModelParameters["Mach"];

    int nx = iterRho->getSizeX();
    int ny = iterRho->getSizeY();
    int nz = iterRho->getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho = pRho[idx_i];
                u = pRhoU[idx_i]/ rho;
                v = pRhoV[idx_i] / rho;
                w = pRhoW[idx_i] / rho;
                E = pRhoE[idx_i];

                p=std::fabs((E-0.5*rho*(u*u+v*v+w*w))*(gamma-1.0));
                e=p/rho/(gamma-1.0);
                alpha=std::sqrt(gamma*p/rho);
                T=gamma*Mach*Mach*p / rho;
                
                pRho0[idx_i]=rho;
                pU0[idx_i]=u;
                pV0[idx_i]=v;
                pW0[idx_i]=w;
                pP0[idx_i]=p;
                
                pT0[idx_i]=T;
                pE0[idx_i]=e;
                pAlpha0[idx_i]=alpha;
            }
        }
    }
}

void nuc3d::physicsModel::prim2con(const std::string &EosName,
 const VectorField &W_vec,
 VectorField &Q_vec)
{
    double rho, u, v, w, E, p;
    double T, e, alpha;    
    double rho0,rhou0,rhov0,rhow0,E0;

    double gamma = myModelParameters["Gamma"];
    double Mach = myModelParameters["Mach"];
    auto iterRho = Q_vec.begin();
    auto iterRhoU = Q_vec.begin() + 1;
    auto iterRhoV = Q_vec.begin() + 2;
    auto iterRhoW = Q_vec.begin() + 3;
    auto iterRhoE = Q_vec.begin() + 4;
    
    auto iterRho0 = W_vec.begin();
    auto iterU0   = W_vec.begin() + 1;
    auto iterV0   = W_vec.begin() + 2;
    auto iterW0   = W_vec.begin() + 3;
    auto iterP0   = W_vec.begin() + 4;

    double *pRho=iterRho->getDataPtr();
    double *pRhoU=iterRhoU->getDataPtr();
    double *pRhoV=iterRhoV->getDataPtr();
    double *pRhoW=iterRhoW->getDataPtr();
    double *pRhoE=iterRhoE->getDataPtr();
    
    double *pRho0=iterRho0->getDataPtr();
    double *pU0=iterU0->getDataPtr();
    double *pV0=iterV0->getDataPtr();
    double *pW0=iterW0->getDataPtr();
    double *pP0=iterP0->getDataPtr();

    int nx = iterRho->getSizeX();
    int ny = iterRho->getSizeY();
    int nz = iterRho->getSizeZ();
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {

                int idx=nx*ny*k+nx*j+i;
                
                rho = pRho0[idx];
                u = pU0[idx]/ rho;
                v = pV0[idx] / rho;
                w = pW0[idx] / rho;
                p = pP0[idx];
                
                pRho[idx]=rho;
                pRhoU[idx]=rho*u;
                pRhoV[idx]=rho*v;
                pRhoW[idx]=rho*w;
                pRhoE[idx]=p/(gamma-1.0)+0.5*rho*(u*u+v*v+w*w);
            }
        }
    }
    
}

void  nuc3d::physicsModel::solveRiemann(PDEData3d &myPDE,EulerData3D *myEuler)
{

    if (myRiemannMap.find(myRiemannName) != myRiemannMap.end())
        (this->*myRiemannMap[myRiemannName])(myEuler->getPrimatives(),
          myEuler->getAcoustics(),
          myPDE.getQ(),
          myEuler->getFluxXi(),
          myEuler->getFluxEta(),
          myEuler->getFluxZeta());
    else
        std::cout << "Riemann Solver " << myRiemannName << " does not exist!" << std::endl;
}

void nuc3d::physicsModel::RiemannLF(const VectorField &W_vec,
    const VectorField &W0_vec,
    const VectorField &Q_vec,
    EulerFlux& myFlux_xi,
    EulerFlux& myFlux_eta,
    EulerFlux& myFlux_zeta)
{    
    double *Rho=Q_vec[0].getDataPtr();
    double *RhoU=Q_vec[1].getDataPtr();
    double *RhoV=Q_vec[2].getDataPtr();
    double *RhoW=Q_vec[3].getDataPtr();
    double *RhoE=Q_vec[4].getDataPtr();
    
    double *rho=W_vec[0].getDataPtr();
    double *u=W_vec[1].getDataPtr();
    double *v=W_vec[2].getDataPtr();
    double *w=W_vec[3].getDataPtr();
    double *p=W_vec[4].getDataPtr();
    
    double *T=W0_vec[0].getDataPtr();
    double *e=W0_vec[1].getDataPtr();
    double *alpha=W0_vec[2].getDataPtr();    
    
    double U0,V0,W0;
    
    double Rho0,Rhou0,Rhov0,Rhow0,Rhoe0;
    double rho0,u0,v0,w0,p0,T0,e0,alpha0;
    
    double MaxEigen_xi=0.0;
    double MaxEigen_eta=0.0;
    double MaxEigen_zeta=0.0;
    double fluxp[5], fluxn[5];
    
    double *flux_xi_l[5];
    double *flux_xi_r[5];
    
    double *flux_eta_l[5];
    double *flux_eta_r[5];
    
    double *flux_zeta_l[5];
    double *flux_zeta_r[5];
    
    flux_xi_l[0]=myFlux_xi.FluxL[0].getDataPtr();
    flux_xi_l[1]=myFlux_xi.FluxL[1].getDataPtr();
    flux_xi_l[2]=myFlux_xi.FluxL[2].getDataPtr();
    flux_xi_l[3]=myFlux_xi.FluxL[3].getDataPtr();
    flux_xi_l[4]=myFlux_xi.FluxL[4].getDataPtr();
    
    flux_xi_r[0]=myFlux_xi.FluxR[0].getDataPtr();
    flux_xi_r[1]=myFlux_xi.FluxR[1].getDataPtr();
    flux_xi_r[2]=myFlux_xi.FluxR[2].getDataPtr();
    flux_xi_r[3]=myFlux_xi.FluxR[3].getDataPtr();
    flux_xi_r[4]=myFlux_xi.FluxR[4].getDataPtr();
    
    flux_eta_l[0]=myFlux_eta.FluxL[0].getDataPtr();
    flux_eta_l[1]=myFlux_eta.FluxL[1].getDataPtr();
    flux_eta_l[2]=myFlux_eta.FluxL[2].getDataPtr();
    flux_eta_l[3]=myFlux_eta.FluxL[3].getDataPtr();
    flux_eta_l[4]=myFlux_eta.FluxL[4].getDataPtr();
    
    flux_eta_r[0]=myFlux_eta.FluxR[0].getDataPtr();
    flux_eta_r[1]=myFlux_eta.FluxR[1].getDataPtr();
    flux_eta_r[2]=myFlux_eta.FluxR[2].getDataPtr();
    flux_eta_r[3]=myFlux_eta.FluxR[3].getDataPtr();
    flux_eta_r[4]=myFlux_eta.FluxR[4].getDataPtr();
    
    flux_zeta_l[0]=myFlux_zeta.FluxL[0].getDataPtr();
    flux_zeta_l[1]=myFlux_zeta.FluxL[1].getDataPtr();
    flux_zeta_l[2]=myFlux_zeta.FluxL[2].getDataPtr();
    flux_zeta_l[3]=myFlux_zeta.FluxL[3].getDataPtr();
    flux_zeta_l[4]=myFlux_zeta.FluxL[4].getDataPtr();
    
    flux_zeta_r[0]=myFlux_zeta.FluxR[0].getDataPtr();
    flux_zeta_r[1]=myFlux_zeta.FluxR[1].getDataPtr();
    flux_zeta_r[2]=myFlux_zeta.FluxR[2].getDataPtr();
    flux_zeta_r[3]=myFlux_zeta.FluxR[3].getDataPtr();
    flux_zeta_r[4]=myFlux_zeta.FluxR[4].getDataPtr();
    
    double flux[5];
    
    int nx = Q_vec[0].getSizeX();
    int ny = Q_vec[0].getSizeY();
    int nz = Q_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=rho[idx_i];
                u0=u[idx_i];
                v0=v[idx_i];
                w0=w[idx_i];
                p0=p[idx_i];

                T0=T[idx_i];
                e0=e[idx_i];
                alpha0=alpha[idx_i];
                
                Rho0=Rho[idx_i];
                Rhou0=RhoU[idx_i];
                Rhov0=RhoV[idx_i];
                Rhow0=RhoW[idx_i];
                Rhoe0=RhoE[idx_i];

                double MaxEigenLocal=std::fabs(u0)+alpha0;
                MaxEigen_xi=std::fmax(MaxEigenLocal ,MaxEigen_xi);

                flux[0]=u0*Rho0;
                flux[1]=u0*Rhou0+p0;
                flux[2]=u0*Rhov0;
                flux[3]=u0*Rhow0;
                flux[4]=u0*(Rhoe0+p0);

                flux_xi_l[0][idx_i] = 0.5*(flux[0]+MaxEigenLocal*Rho0);
                flux_xi_l[1][idx_i] = 0.5*(flux[1]+MaxEigenLocal*Rhou0);
                flux_xi_l[2][idx_i] = 0.5*(flux[2]+MaxEigenLocal*Rhov0);
                flux_xi_l[3][idx_i] = 0.5*(flux[3]+MaxEigenLocal*Rhow0);
                flux_xi_l[4][idx_i] = 0.5*(flux[4]+MaxEigenLocal*Rhoe0);

                flux_xi_r[0][idx_i] = 0.5*(flux[0]-MaxEigenLocal*Rho0);
                flux_xi_r[1][idx_i] = 0.5*(flux[1]-MaxEigenLocal*Rhou0);
                flux_xi_r[2][idx_i] = 0.5*(flux[2]-MaxEigenLocal*Rhov0);
                flux_xi_r[3][idx_i] = 0.5*(flux[3]-MaxEigenLocal*Rhow0);
                flux_xi_r[4][idx_i] = 0.5*(flux[4]-MaxEigenLocal*Rhoe0);                
                
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
                
                rho0=rho[idx_i_xi];
                u0=u[idx_i_xi];
                v0=v[idx_i_xi];
                w0=w[idx_i_xi];
                p0=p[idx_i_xi];

                T0=T[idx_i_xi];
                e0=e[idx_i_xi];
                alpha0=alpha[idx_i_xi];
                
                Rho0=Rho[idx_i_xi];
                Rhou0=RhoU[idx_i_xi];
                Rhov0=RhoV[idx_i_xi];
                Rhow0=RhoW[idx_i_xi];
                Rhoe0=RhoE[idx_i_xi];

                double MaxEigenLocal=std::fabs(v0)+alpha0;
                MaxEigen_eta=std::fmax(MaxEigenLocal ,MaxEigen_eta);

                flux[0]=v0*Rho0;
                flux[1]=v0*Rhou0;
                flux[2]=v0*Rhov0+p0;
                flux[3]=v0*Rhow0;
                flux[4]=v0*(Rhoe0+p0);

                flux_eta_l[0][idx_i_eta] = 0.5*(flux[0]+MaxEigenLocal*Rho0);
                flux_eta_l[1][idx_i_eta] = 0.5*(flux[1]+MaxEigenLocal*Rhou0);
                flux_eta_l[2][idx_i_eta] = 0.5*(flux[2]+MaxEigenLocal*Rhov0);
                flux_eta_l[3][idx_i_eta] = 0.5*(flux[3]+MaxEigenLocal*Rhow0);
                flux_eta_l[4][idx_i_eta] = 0.5*(flux[4]+MaxEigenLocal*Rhoe0);

                flux_eta_r[0][idx_i_eta] = 0.5*(flux[0]-MaxEigenLocal*Rho0);
                flux_eta_r[1][idx_i_eta] = 0.5*(flux[1]-MaxEigenLocal*Rhou0);
                flux_eta_r[2][idx_i_eta] = 0.5*(flux[2]-MaxEigenLocal*Rhov0);
                flux_eta_r[3][idx_i_eta] = 0.5*(flux[3]-MaxEigenLocal*Rhow0);
                flux_eta_r[4][idx_i_eta] = 0.5*(flux[4]-MaxEigenLocal*Rhoe0);                
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

                rho0=rho[idx_i_xi];
                u0=u[idx_i_xi];
                v0=v[idx_i_xi];
                w0=w[idx_i_xi];
                p0=p[idx_i_xi];

                T0=T[idx_i_xi];
                e0=e[idx_i_xi];
                alpha0=alpha[idx_i_xi];
                
                Rho0=Rho[idx_i_xi];
                Rhou0=RhoU[idx_i_xi];
                Rhov0=RhoV[idx_i_xi];
                Rhow0=RhoW[idx_i_xi];
                Rhoe0=RhoE[idx_i_xi];

                double MaxEigenLocal=std::fabs(w0)+alpha0;
                MaxEigen_eta=std::fmax(MaxEigenLocal ,MaxEigen_eta);

                flux[0]=w0*Rho0;
                flux[1]=w0*Rhou0;
                flux[2]=w0*Rhov0;
                flux[3]=w0*Rhow0+p0;
                flux[4]=w0*(Rhoe0+p0);

                flux_zeta_l[0][idx_i_zeta] = 0.5*(flux[0]+MaxEigenLocal*Rho0);
                flux_zeta_l[1][idx_i_zeta] = 0.5*(flux[1]+MaxEigenLocal*Rhou0);
                flux_zeta_l[2][idx_i_zeta] = 0.5*(flux[2]+MaxEigenLocal*Rhov0);
                flux_zeta_l[3][idx_i_zeta] = 0.5*(flux[3]+MaxEigenLocal*Rhow0);
                flux_zeta_l[4][idx_i_zeta] = 0.5*(flux[4]+MaxEigenLocal*Rhoe0);

                flux_zeta_r[0][idx_i_zeta] = 0.5*(flux[0]-MaxEigenLocal*Rho0);
                flux_zeta_r[1][idx_i_zeta] = 0.5*(flux[1]-MaxEigenLocal*Rhou0);
                flux_zeta_r[2][idx_i_zeta] = 0.5*(flux[2]-MaxEigenLocal*Rhov0);
                flux_zeta_r[3][idx_i_zeta] = 0.5*(flux[3]-MaxEigenLocal*Rhow0);
                flux_zeta_r[4][idx_i_zeta] = 0.5*(flux[4]-MaxEigenLocal*Rhoe0);
            }
        }
    }
    
    myFlux_xi.maxEigen=MaxEigen_xi;
    myFlux_eta.maxEigen=MaxEigen_eta;
    myFlux_zeta.maxEigen=MaxEigen_zeta;
};

void nuc3d::physicsModel::RiemannAUSM(const VectorField &W_vec,
    const VectorField &W0_vec,
    const VectorField &Q_vec,
    EulerFlux& myFlux_xi,
    EulerFlux& myFlux_eta,
    EulerFlux& myFlux_zeta)
{    
    double *Rho=Q_vec[0].getDataPtr();
    double *RhoU=Q_vec[1].getDataPtr();
    double *RhoV=Q_vec[2].getDataPtr();
    double *RhoW=Q_vec[3].getDataPtr();
    double *RhoE=Q_vec[4].getDataPtr();
    
    double *rho=W_vec[0].getDataPtr();
    double *u=W_vec[1].getDataPtr();
    double *v=W_vec[2].getDataPtr();
    double *w=W_vec[3].getDataPtr();
    double *p=W_vec[4].getDataPtr();
    
    double *T=W0_vec[0].getDataPtr();
    double *e=W0_vec[1].getDataPtr();
    double *alpha=W0_vec[2].getDataPtr();    
    
    double U0,V0,W0;
    
    double Rho0,Rhou0,Rhov0,Rhow0,Rhoe0;
    double rho0,u0,v0,w0,p0,T0,e0,alpha0;
    
    double MaxEigen_xi=0.0;
    double MaxEigen_eta=0.0;
    double MaxEigen_zeta=0.0;
    double fluxp[5], fluxn[5];
    
    double *flux_xi_l[5];
    double *flux_xi_r[5];
    
    double *flux_eta_l[5];
    double *flux_eta_r[5];
    
    double *flux_zeta_l[5];
    double *flux_zeta_r[5];
    
    flux_xi_l[0]=myFlux_xi.FluxL[0].getDataPtr();
    flux_xi_l[1]=myFlux_xi.FluxL[1].getDataPtr();
    flux_xi_l[2]=myFlux_xi.FluxL[2].getDataPtr();
    flux_xi_l[3]=myFlux_xi.FluxL[3].getDataPtr();
    flux_xi_l[4]=myFlux_xi.FluxL[4].getDataPtr();
    
    flux_xi_r[0]=myFlux_xi.FluxR[0].getDataPtr();
    flux_xi_r[1]=myFlux_xi.FluxR[1].getDataPtr();
    flux_xi_r[2]=myFlux_xi.FluxR[2].getDataPtr();
    flux_xi_r[3]=myFlux_xi.FluxR[3].getDataPtr();
    flux_xi_r[4]=myFlux_xi.FluxR[4].getDataPtr();
    
    flux_eta_l[0]=myFlux_eta.FluxL[0].getDataPtr();
    flux_eta_l[1]=myFlux_eta.FluxL[1].getDataPtr();
    flux_eta_l[2]=myFlux_eta.FluxL[2].getDataPtr();
    flux_eta_l[3]=myFlux_eta.FluxL[3].getDataPtr();
    flux_eta_l[4]=myFlux_eta.FluxL[4].getDataPtr();
    
    flux_eta_r[0]=myFlux_eta.FluxR[0].getDataPtr();
    flux_eta_r[1]=myFlux_eta.FluxR[1].getDataPtr();
    flux_eta_r[2]=myFlux_eta.FluxR[2].getDataPtr();
    flux_eta_r[3]=myFlux_eta.FluxR[3].getDataPtr();
    flux_eta_r[4]=myFlux_eta.FluxR[4].getDataPtr();
    
    flux_zeta_l[0]=myFlux_zeta.FluxL[0].getDataPtr();
    flux_zeta_l[1]=myFlux_zeta.FluxL[1].getDataPtr();
    flux_zeta_l[2]=myFlux_zeta.FluxL[2].getDataPtr();
    flux_zeta_l[3]=myFlux_zeta.FluxL[3].getDataPtr();
    flux_zeta_l[4]=myFlux_zeta.FluxL[4].getDataPtr();
    
    flux_zeta_r[0]=myFlux_zeta.FluxR[0].getDataPtr();
    flux_zeta_r[1]=myFlux_zeta.FluxR[1].getDataPtr();
    flux_zeta_r[2]=myFlux_zeta.FluxR[2].getDataPtr();
    flux_zeta_r[3]=myFlux_zeta.FluxR[3].getDataPtr();
    flux_zeta_r[4]=myFlux_zeta.FluxR[4].getDataPtr();
    
    double flux[5];
    
    int nx = Q_vec[0].getSizeX();
    int ny = Q_vec[0].getSizeY();
    int nz = Q_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=rho[idx_i];
                u0=u[idx_i];
                v0=v[idx_i];
                w0=w[idx_i];
                p0=p[idx_i];

                T0=T[idx_i];
                e0=e[idx_i];
                alpha0=alpha[idx_i];
                
                Rho0=Rho[idx_i];
                Rhou0=RhoU[idx_i];
                Rhov0=RhoV[idx_i];
                Rhow0=RhoW[idx_i];
                Rhoe0=RhoE[idx_i];

                double mach=u0/alpha0;
                double machp = getMachL(mach);
                double machn = getMachR(mach);

                double p_p = getPressureL(mach, p0);
                double p_n = getPressureR(mach, p0);

                double MaxEigenLocal=std::fabs(u0)+alpha0;
                MaxEigen_xi=std::fmax(MaxEigenLocal ,MaxEigen_xi);

                flux_xi_l[0][idx_i] = machp*alpha0*Rho0;
                flux_xi_l[1][idx_i] = machp*alpha0*Rhou0+p_p;
                flux_xi_l[2][idx_i] = machp*alpha0*Rhov0;
                flux_xi_l[3][idx_i] = machp*alpha0*Rhow0;
                flux_xi_l[4][idx_i] = machp*alpha0*Rhoe0 + machp*alpha0*p0;
                
                flux_xi_r[0][idx_i] = machn*alpha0*Rho0;
                flux_xi_r[1][idx_i] = machn*alpha0*Rhou0+p_n;
                flux_xi_r[2][idx_i] = machn*alpha0*Rhov0;
                flux_xi_r[3][idx_i] = machn*alpha0*Rhow0;
                flux_xi_r[4][idx_i] = machn*alpha0*Rhoe0 + machn*alpha0*p0;
                
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
                
                rho0=rho[idx_i_xi];
                u0=u[idx_i_xi];
                v0=v[idx_i_xi];
                w0=w[idx_i_xi];
                p0=p[idx_i_xi];

                T0=T[idx_i_xi];
                e0=e[idx_i_xi];
                alpha0=alpha[idx_i_xi];
                
                Rho0=Rho[idx_i_xi];
                Rhou0=RhoU[idx_i_xi];
                Rhov0=RhoV[idx_i_xi];
                Rhow0=RhoW[idx_i_xi];
                Rhoe0=RhoE[idx_i_xi];

                double mach=v0/alpha0;
                double machp = getMachL(mach);
                double machn = getMachR(mach);

                double p_p = getPressureL(mach, p0);
                double p_n = getPressureR(mach, p0);

                double MaxEigenLocal=std::fabs(v0)+alpha0;
                MaxEigen_eta=std::fmax(MaxEigenLocal ,MaxEigen_eta);

                flux_eta_l[0][idx_i_eta] = machp*alpha0*Rho0;
                flux_eta_l[1][idx_i_eta] = machp*alpha0*Rhou0;
                flux_eta_l[2][idx_i_eta] = machp*alpha0*Rhov0+p_p;
                flux_eta_l[3][idx_i_eta] = machp*alpha0*Rhow0;
                flux_eta_l[4][idx_i_eta] = machp*alpha0*Rhoe0 + machp*alpha0*p0;
                
                flux_eta_r[0][idx_i_eta] = machn*alpha0*Rho0;
                flux_eta_r[1][idx_i_eta] = machn*alpha0*Rhou0;
                flux_eta_r[2][idx_i_eta] = machn*alpha0*Rhov0+p_n;
                flux_eta_r[3][idx_i_eta] = machn*alpha0*Rhow0;
                flux_eta_r[4][idx_i_eta] = machn*alpha0*Rhoe0 + machn*alpha0*p0;            
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

                rho0=rho[idx_i_xi];
                u0=u[idx_i_xi];
                v0=v[idx_i_xi];
                w0=w[idx_i_xi];
                p0=p[idx_i_xi];

                T0=T[idx_i_xi];
                e0=e[idx_i_xi];
                alpha0=alpha[idx_i_xi];
                
                Rho0=Rho[idx_i_xi];
                Rhou0=RhoU[idx_i_xi];
                Rhov0=RhoV[idx_i_xi];
                Rhow0=RhoW[idx_i_xi];
                Rhoe0=RhoE[idx_i_xi];

                double mach=w0/alpha0;
                double machp = getMachL(mach);
                double machn = getMachR(mach);

                double p_p = getPressureL(mach, p0);
                double p_n = getPressureR(mach, p0);

                double MaxEigenLocal=std::fabs(w0)+alpha0;
                MaxEigen_zeta=std::fmax(MaxEigenLocal ,MaxEigen_zeta);

                flux_zeta_l[0][idx_i_zeta] = machp*alpha0*Rho0;
                flux_zeta_l[1][idx_i_zeta] = machp*alpha0*Rhou0;
                flux_zeta_l[2][idx_i_zeta] = machp*alpha0*Rhov0;
                flux_zeta_l[3][idx_i_zeta] = machp*alpha0*Rhow0+p_p;
                flux_zeta_l[4][idx_i_zeta] = machp*alpha0*Rhoe0 + machp*alpha0*p0;
                
                flux_zeta_r[0][idx_i_zeta] = machn*alpha0*Rho0;
                flux_zeta_r[1][idx_i_zeta] = machn*alpha0*Rhou0;
                flux_zeta_r[2][idx_i_zeta] = machn*alpha0*Rhov0;
                flux_zeta_r[3][idx_i_zeta] = machn*alpha0*Rhow0+p_n;
                flux_zeta_r[4][idx_i_zeta] = machn*alpha0*Rhoe0 + machn*alpha0*p0;    
            }
        }
    }
    
    myFlux_xi.maxEigen=MaxEigen_xi;
    myFlux_eta.maxEigen=MaxEigen_eta;
    myFlux_zeta.maxEigen=MaxEigen_zeta;
};


void nuc3d::physicsModel::RiemannAUSMp(const VectorField &W_vec,
    const VectorField &W0_vec,
    const VectorField &Q_vec,
    EulerFlux& myFlux_xi,
    EulerFlux& myFlux_eta,
    EulerFlux& myFlux_zeta)
{     
    double *Rho=Q_vec[0].getDataPtr();
    double *RhoU=Q_vec[1].getDataPtr();
    double *RhoV=Q_vec[2].getDataPtr();
    double *RhoW=Q_vec[3].getDataPtr();
    double *RhoE=Q_vec[4].getDataPtr();
    
    double *rho=W_vec[0].getDataPtr();
    double *u=W_vec[1].getDataPtr();
    double *v=W_vec[2].getDataPtr();
    double *w=W_vec[3].getDataPtr();
    double *p=W_vec[4].getDataPtr();
    
    double *T=W0_vec[0].getDataPtr();
    double *e=W0_vec[1].getDataPtr();
    double *alpha=W0_vec[2].getDataPtr();    
    
    double U0,V0,W0;
    
    double Rho0,Rhou0,Rhov0,Rhow0,Rhoe0;
    double rho0,u0,v0,w0,p0,T0,e0,alpha0;
    
    double MaxEigen_xi=0.0;
    double MaxEigen_eta=0.0;
    double MaxEigen_zeta=0.0;
    double fluxp[5], fluxn[5];
    
    double *flux_xi_l[5];
    double *flux_xi_r[5];
    
    double *flux_eta_l[5];
    double *flux_eta_r[5];
    
    double *flux_zeta_l[5];
    double *flux_zeta_r[5];
    
    flux_xi_l[0]=myFlux_xi.FluxL[0].getDataPtr();
    flux_xi_l[1]=myFlux_xi.FluxL[1].getDataPtr();
    flux_xi_l[2]=myFlux_xi.FluxL[2].getDataPtr();
    flux_xi_l[3]=myFlux_xi.FluxL[3].getDataPtr();
    flux_xi_l[4]=myFlux_xi.FluxL[4].getDataPtr();
    
    flux_xi_r[0]=myFlux_xi.FluxR[0].getDataPtr();
    flux_xi_r[1]=myFlux_xi.FluxR[1].getDataPtr();
    flux_xi_r[2]=myFlux_xi.FluxR[2].getDataPtr();
    flux_xi_r[3]=myFlux_xi.FluxR[3].getDataPtr();
    flux_xi_r[4]=myFlux_xi.FluxR[4].getDataPtr();
    
    flux_eta_l[0]=myFlux_eta.FluxL[0].getDataPtr();
    flux_eta_l[1]=myFlux_eta.FluxL[1].getDataPtr();
    flux_eta_l[2]=myFlux_eta.FluxL[2].getDataPtr();
    flux_eta_l[3]=myFlux_eta.FluxL[3].getDataPtr();
    flux_eta_l[4]=myFlux_eta.FluxL[4].getDataPtr();
    
    flux_eta_r[0]=myFlux_eta.FluxR[0].getDataPtr();
    flux_eta_r[1]=myFlux_eta.FluxR[1].getDataPtr();
    flux_eta_r[2]=myFlux_eta.FluxR[2].getDataPtr();
    flux_eta_r[3]=myFlux_eta.FluxR[3].getDataPtr();
    flux_eta_r[4]=myFlux_eta.FluxR[4].getDataPtr();
    
    flux_zeta_l[0]=myFlux_zeta.FluxL[0].getDataPtr();
    flux_zeta_l[1]=myFlux_zeta.FluxL[1].getDataPtr();
    flux_zeta_l[2]=myFlux_zeta.FluxL[2].getDataPtr();
    flux_zeta_l[3]=myFlux_zeta.FluxL[3].getDataPtr();
    flux_zeta_l[4]=myFlux_zeta.FluxL[4].getDataPtr();
    
    flux_zeta_r[0]=myFlux_zeta.FluxR[0].getDataPtr();
    flux_zeta_r[1]=myFlux_zeta.FluxR[1].getDataPtr();
    flux_zeta_r[2]=myFlux_zeta.FluxR[2].getDataPtr();
    flux_zeta_r[3]=myFlux_zeta.FluxR[3].getDataPtr();
    flux_zeta_r[4]=myFlux_zeta.FluxR[4].getDataPtr();
    
    double flux[5];
    
    int nx = Q_vec[0].getSizeX();
    int ny = Q_vec[0].getSizeY();
    int nz = Q_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=rho[idx_i];
                u0=u[idx_i];
                v0=v[idx_i];
                w0=w[idx_i];
                p0=p[idx_i];

                T0=T[idx_i];
                e0=e[idx_i];
                alpha0=alpha[idx_i];
                
                Rho0=Rho[idx_i];
                Rhou0=RhoU[idx_i];
                Rhov0=RhoV[idx_i];
                Rhow0=RhoW[idx_i];
                Rhoe0=RhoE[idx_i];

                double mach=u0/alpha0;
                double machp = getMachLp(mach);
                double machn = getMachRp(mach);

                double p_p = getPressureLp(mach, p0);
                double p_n = getPressureRp(mach, p0);

                double MaxEigenLocal=std::fabs(u0)+alpha0;
                MaxEigen_xi=std::fmax(MaxEigenLocal ,MaxEigen_xi);

                flux_xi_l[0][idx_i] = machp*alpha0*Rho0;
                flux_xi_l[1][idx_i] = machp*alpha0*Rhou0+p_p;
                flux_xi_l[2][idx_i] = machp*alpha0*Rhov0;
                flux_xi_l[3][idx_i] = machp*alpha0*Rhow0;
                flux_xi_l[4][idx_i] = machp*alpha0*Rhoe0 + machp*alpha0*p0;
                
                flux_xi_r[0][idx_i] = machn*alpha0*Rho0;
                flux_xi_r[1][idx_i] = machn*alpha0*Rhou0+p_n;
                flux_xi_r[2][idx_i] = machn*alpha0*Rhov0;
                flux_xi_r[3][idx_i] = machn*alpha0*Rhow0;
                flux_xi_r[4][idx_i] = machn*alpha0*Rhoe0 + machn*alpha0*p0;
                
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
                
                rho0=rho[idx_i_xi];
                u0=u[idx_i_xi];
                v0=v[idx_i_xi];
                w0=w[idx_i_xi];
                p0=p[idx_i_xi];

                T0=T[idx_i_xi];
                e0=e[idx_i_xi];
                alpha0=alpha[idx_i_xi];
                
                Rho0=Rho[idx_i_xi];
                Rhou0=RhoU[idx_i_xi];
                Rhov0=RhoV[idx_i_xi];
                Rhow0=RhoW[idx_i_xi];
                Rhoe0=RhoE[idx_i_xi];

                double mach=v0/alpha0;
                double machp = getMachLp(mach);
                double machn = getMachRp(mach);

                double p_p = getPressureLp(mach, p0);
                double p_n = getPressureRp(mach, p0);

                double MaxEigenLocal=std::fabs(v0)+alpha0;
                MaxEigen_eta=std::fmax(MaxEigenLocal ,MaxEigen_eta);

                flux_eta_l[0][idx_i_eta] = machp*alpha0*Rho0;
                flux_eta_l[1][idx_i_eta] = machp*alpha0*Rhou0;
                flux_eta_l[2][idx_i_eta] = machp*alpha0*Rhov0+p_p;
                flux_eta_l[3][idx_i_eta] = machp*alpha0*Rhow0;
                flux_eta_l[4][idx_i_eta] = machp*alpha0*Rhoe0 + machp*alpha0*p0;
                
                flux_eta_r[0][idx_i_eta] = machn*alpha0*Rho0;
                flux_eta_r[1][idx_i_eta] = machn*alpha0*Rhou0;
                flux_eta_r[2][idx_i_eta] = machn*alpha0*Rhov0+p_n;
                flux_eta_r[3][idx_i_eta] = machn*alpha0*Rhow0;
                flux_eta_r[4][idx_i_eta] = machn*alpha0*Rhoe0 + machn*alpha0*p0;            
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

                rho0=rho[idx_i_xi];
                u0=u[idx_i_xi];
                v0=v[idx_i_xi];
                w0=w[idx_i_xi];
                p0=p[idx_i_xi];

                T0=T[idx_i_xi];
                e0=e[idx_i_xi];
                alpha0=alpha[idx_i_xi];
                
                Rho0=Rho[idx_i_xi];
                Rhou0=RhoU[idx_i_xi];
                Rhov0=RhoV[idx_i_xi];
                Rhow0=RhoW[idx_i_xi];
                Rhoe0=RhoE[idx_i_xi];

                double mach=w0/alpha0;
                double machp = getMachLp(mach);
                double machn = getMachRp(mach);

                double p_p = getPressureLp(mach, p0);
                double p_n = getPressureRp(mach, p0);

                double MaxEigenLocal=std::fabs(w0)+alpha0;
                MaxEigen_zeta=std::fmax(MaxEigenLocal ,MaxEigen_zeta);

                flux_zeta_l[0][idx_i_zeta] = machp*alpha0*Rho0;
                flux_zeta_l[1][idx_i_zeta] = machp*alpha0*Rhou0;
                flux_zeta_l[2][idx_i_zeta] = machp*alpha0*Rhov0;
                flux_zeta_l[3][idx_i_zeta] = machp*alpha0*Rhow0+p_p;
                flux_zeta_l[4][idx_i_zeta] = machp*alpha0*Rhoe0 + machp*alpha0*p0;
                
                flux_zeta_r[0][idx_i_zeta] = machn*alpha0*Rho0;
                flux_zeta_r[1][idx_i_zeta] = machn*alpha0*Rhou0;
                flux_zeta_r[2][idx_i_zeta] = machn*alpha0*Rhov0;
                flux_zeta_r[3][idx_i_zeta] = machn*alpha0*Rhow0+p_n;
                flux_zeta_r[4][idx_i_zeta] = machn*alpha0*Rhoe0 + machn*alpha0*p0;    
            }
        }
    }
    
    myFlux_xi.maxEigen=MaxEigen_xi;
    myFlux_eta.maxEigen=MaxEigen_eta;
    myFlux_zeta.maxEigen=MaxEigen_zeta;
};

double nuc3d::physicsModel::getMachL(const double &mach)
{
    double MachL;
    
    if (fabs(mach) < 1.0)
        MachL = 0.25*std::pow((mach + 1.0), 2)+0.125*std::pow((mach*mach-1), 2);
    else
        MachL = 0.50*(mach + std::fabs(mach));
    
    return MachL;
    
}

double nuc3d::physicsModel::getMachR(const double &mach)
{

    double MachR;
    
    if (fabs(mach) < 1.0)
        MachR = -0.25*std::pow((mach - 1.0), 2)-0.125*std::pow((mach*mach-1), 2);
    else
        MachR = 0.50*(mach - std::fabs(mach));
    
    return MachR;
    
}

double nuc3d::physicsModel::getPressureL(const double &mach, const double &p)
{
    double pressureL;
    if (fabs(mach) < 1.0)
        pressureL = p*(mach+1.0)/2.0;
    else
        pressureL = 0.50*p*(mach + std::fabs(mach)) / mach;
    
    return pressureL;
}

double nuc3d::physicsModel::getPressureR(const double &mach, const double &p)
{
    double pressureR;
    if (fabs(mach) < 1.0)
        pressureR = -p*(mach-1.0)/2.0;
    else
        pressureR = 0.5*p*(mach - std::fabs(mach)) / mach;
    
    return pressureR;
}


double nuc3d::physicsModel::getMachLp(const double &mach)
{
    double MachL;
    
    if (fabs(mach) < 1.0)
        MachL = 0.25*std::pow((mach + 1.0), 2.0)+0.125*std::pow((mach*mach-1), 2.0);
    else
        MachL = 0.50*(mach + std::fabs(mach));
    
    return MachL;
    
}

double nuc3d::physicsModel::getMachRp(const double &mach)
{

    double MachR;
    
    if (fabs(mach) < 1.0)
        MachR = -0.25*std::pow((mach - 1.0), 2.0)-0.125*std::pow((mach*mach-1), 2.0);
    else
        MachR = 0.50*(mach - std::fabs(mach));
    
    return MachR;
    
}

double nuc3d::physicsModel::getPressureLp(const double &mach, const double &p)
{
    double pressureL;
    if (fabs(mach) < 1.0)
        pressureL = p*(0.25*std::pow((mach + 1.0), 2.0)*(2.0 - mach)+0.1875*mach*std::pow((mach*mach-1.0),2.0));
    else
        pressureL = 0.50*p*(mach + std::fabs(mach)) / mach;
    
    return pressureL;
}

double nuc3d::physicsModel::getPressureRp(const double &mach, const double &p)
{
    double pressureR;
    if (fabs(mach) < 1.0)
        pressureR = p*(0.25*std::pow(mach - 1.0, 2.0)*(2.0 + mach)-0.1875*mach*std::pow((mach*mach-1.0),2.0));
    else
        pressureR = 0.5*p*(mach - std::fabs(mach)) / mach;
    
    return pressureR;
}

void nuc3d::physicsModel::solveSource(PDEData3d &myPDE,
  EulerData3D *myEuler,
  MPIComunicator3d_nonblocking &myMPI,
  fft &myFFT,
  const VectorField &xyz_center)
{
        AddForcing(myEuler->getPrimatives(),
          myEuler->getAcoustics(),
          myPDE.getQ(),
          myPDE.getRHS(),
          myMPI,
          myFFT,
          xyz_center,
          myPDE.getDt());

    if(("no")!=mySourceSwitchMap["Coriolis"])
    {
        // if(0==myMPI.getMyId())
        // std::cout<<"calculating coriolisForce!"<<std::endl;
        
        AddCoriolisForce(myEuler->getPrimatives(),
          myEuler->getAcoustics(),
          myPDE.getQ(),
          myPDE.getRHS());
    }


    if(("no")!=mySourceSwitchMap["Centrifugal"])
    {
        // if(0==myMPI.getMyId())
        // std::cout<<"calculating AddCentrifugalForce!"<<std::endl;

        AddCentrifugalForce(myEuler->getPrimatives(),
          myEuler->getAcoustics(),
          myPDE.getQ(),
          xyz_center,
          myPDE.getRHS());
    }


}


void nuc3d::physicsModel::AddForcing(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  VectorField &RHS,
  MPIComunicator3d_nonblocking &myMPI,
  fft &myFFT,
  const VectorField &xyz_center,
  double dt)
{
    if (myForcingMap.find(myForcingName) != myForcingMap.end())
    {

        // if(0==myMPI.getMyId())
        // std::cout<<"calculating AddForcing!"<<std::endl;

        (this->*myForcingMap[myForcingName])(W_vec,
          W0_vec,
          Q_vec,
          RHS,
          myMPI,
          myFFT,
          xyz_center,
          dt);
    }
    else
    {

        std::cout << "Forcing function " << myRiemannName << " does not exist!" << std::endl;
        exit(-1);
    }

}

void nuc3d::physicsModel::AddForcing_None(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  VectorField &RHS,
  MPIComunicator3d_nonblocking &myMPI,
  fft &myFFT,
  const VectorField &xyz_center,
  double dt)
{

}

void nuc3d::physicsModel::AddForcing_Const(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  VectorField &RHS,
  MPIComunicator3d_nonblocking &myMPI,
  fft &myFFT,
  const VectorField &xyz_center,
  double dt)
{
    int nx=W_vec[0].getSizeX();
    int ny=W_vec[0].getSizeY();
    int nz=W_vec[0].getSizeZ();

    VectorField forcing_re;
    VectorField forcing_im;
    VectorField FFT_vector_re;
    VectorField FFT_vector_im;
    Field A(nx,ny,nz);

    for(int i=0;i<3;i++)
    {
        forcing_re.push_back(Field(nx,ny,nz));
        forcing_im.push_back(Field(nx,ny,nz));
        FFT_vector_re.push_back(Field(nx,ny,nz));
        FFT_vector_im.push_back(Field(nx,ny,nz));
    }

    double *vel_c_re[3];
    double *vel_c_im[3];

    double *vel_s_re[3];
    double *vel_s_im[3];

    double *vel_l_re[3];
    double *vel_l_im[3];

    double *vel_r_re[3];
    double *vel_r_im[3];

    for(int i=0;i<3;i++)
    {
        vel_c_re[i]=new double[nx*ny*nz];
        vel_c_im[i]=new double[nx*ny*nz];

        vel_s_re[i]=new double[nx*ny*nz];
        vel_s_im[i]=new double[nx*ny*nz];

        vel_l_re[i]=new double[nx*ny*nz];
        vel_l_im[i]=new double[nx*ny*nz];

        vel_r_re[i]=new double[nx*ny*nz];
        vel_r_im[i]=new double[nx*ny*nz];
    }

    double forcing_amp=myModelParameters["Forcing_amplitude"];
    double sigma=myModelParameters["Forcing_amplitude_hel"];
    double forcing_k=myModelParameters["Forcing_wave_k"];

    int npx=myMPI.getMyDim(2);
    int npy=myMPI.getMyDim(1);
    int npz=myMPI.getMyDim(0);

    int px=myMPI.getMyCoord(2);
    int py=myMPI.getMyCoord(1);
    int pz=myMPI.getMyCoord(0);
    int nx_global=npx*nx;
    int ny_global=npy*ny;
    int nz_global=npz*nz;

    myFFT.solveFFT_forward(W_vec[1],FFT_vector_re[0],FFT_vector_im[0],myMPI); // FFT of velocity u
    myFFT.solveFFT_forward(W_vec[2],FFT_vector_re[1],FFT_vector_im[1],myMPI); // FFT of velocity v
    myFFT.solveFFT_forward(W_vec[3],FFT_vector_re[2],FFT_vector_im[2],myMPI); // FFT of velocity w

    double *pu_re=FFT_vector_re[0].getDataPtr();
    double *pv_re=FFT_vector_re[1].getDataPtr();
    double *pw_re=FFT_vector_re[2].getDataPtr();

    double *pu_im=FFT_vector_im[0].getDataPtr();
    double *pv_im=FFT_vector_im[1].getDataPtr();
    double *pw_im=FFT_vector_im[2].getDataPtr();

    double *pforcing_u_re=forcing_re[0].getDataPtr();
    double *pforcing_v_re=forcing_re[1].getDataPtr();
    double *pforcing_w_re=forcing_re[2].getDataPtr();

    double *pforcing_u_im=forcing_im[0].getDataPtr();
    double *pforcing_v_im=forcing_im[1].getDataPtr();
    double *pforcing_w_im=forcing_im[2].getDataPtr();

    double kx=(double)nx*(double)px;
    double ky=(double)ny*(double)py;
    double kz=(double)nz*(double)pz;

    double *prho=W_vec[0].getDataPtr();
    double *pu=W_vec[1].getDataPtr();
    double *pv=W_vec[2].getDataPtr(); 
    double *pw=W_vec[3].getDataPtr(); 

    double *prhs1=RHS[1].getDataPtr();  
    double *prhs2=RHS[2].getDataPtr();  
    double *prhs3=RHS[3].getDataPtr();  
    double *prhs4=RHS[4].getDataPtr();  
    
    double rho0,u0,v0,w0;

    double Ek_s_k=0.0;
    double Ek_c_k=0.0;

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx=idx_j+i;

                double kx0=((kx+(double)i)>=(nx_global/2.0))?kx+(double)i-nx_global:kx+(double)i;
                double ky0=((ky+(double)j)>=(ny_global/2.0))?ky+(double)j-ny_global:ky+(double)j;
                double kz0=((kz+(double)k)>=(nz_global/2.0))?kz+(double)k-nz_global:kz+(double)k;

                double kmag=std::sqrt(kx0*kx0+ky0*ky0+kz0*kz0);
                if(kmag<1e-6) kmag=1.0;
                
                double ku_re=(kx0*pu_re[idx]
                    +ky0*pv_re[idx]
                    +kz0*pw_re[idx])/kmag;

                double ku_im=(kx0*pu_im[idx]
                    +ky0*pv_im[idx]
                    +kz0*pw_im[idx])/kmag;
                
                vel_s_re[0][idx]=pu_re[idx]-ku_re*kx0/kmag;
                vel_s_re[1][idx]=pv_re[idx]-ku_re*ky0/kmag;                 
                vel_s_re[2][idx]=pw_re[idx]-ku_re*kz0/kmag;

                vel_s_im[0][idx]=pu_im[idx]-ku_im*kx0/kmag;
                vel_s_im[1][idx]=pv_im[idx]-ku_im*ky0/kmag;                 
                vel_s_im[2][idx]=pw_im[idx]-ku_im*kz0/kmag;

                vel_c_re[0][idx]=pu_re[idx]-vel_s_re[0][idx];
                vel_c_re[1][idx]=pv_re[idx]-vel_s_re[1][idx];                 
                vel_c_re[2][idx]=pw_re[idx]-vel_s_re[2][idx];

                vel_c_im[0][idx]=pu_im[idx]-vel_s_im[0][idx];
                vel_c_im[1][idx]=pv_im[idx]-vel_s_im[1][idx];                 
                vel_c_im[2][idx]=pw_im[idx]-vel_s_im[2][idx];
                
                double kxu0_re=-(ky0*vel_s_im[2][idx]-kz0*vel_s_im[1][idx])/kmag;
                double kxu1_re=(kx0*vel_s_im[2][idx]-kz0*vel_s_im[0][idx])/kmag;
                double kxu2_re=-(kx0*vel_s_im[1][idx]-ky0*vel_s_im[0][idx])/kmag;


                double kxu0_im=(ky0*vel_s_re[2][idx]-kz0*vel_s_re[1][idx])/kmag;
                double kxu1_im=-(kx0*vel_s_re[2][idx]-kz0*vel_s_re[0][idx])/kmag;
                double kxu2_im=(kx0*vel_s_re[1][idx]-ky0*vel_s_re[0][idx])/kmag;

                vel_s_re[0][idx]=(vel_s_re[0][idx]+sigma*kxu0_re);
                vel_s_re[1][idx]=(vel_s_re[1][idx]+sigma*kxu1_re);                 
                vel_s_re[2][idx]=(vel_s_re[2][idx]+sigma*kxu2_re);

                vel_s_im[0][idx]=(vel_s_im[0][idx]+sigma*kxu0_im);
                vel_s_im[1][idx]=(vel_s_im[1][idx]+sigma*kxu1_im);                 
                vel_s_im[2][idx]=(vel_s_im[2][idx]+sigma*kxu2_im);

                if(fabs(kmag-forcing_k)<0.5)
                {

                    Ek_s_k+=0.5*(vel_s_re[0][idx]*vel_s_re[0][idx]
                        +vel_s_re[1][idx]*vel_s_re[1][idx]
                        +vel_s_re[2][idx]*vel_s_re[2][idx]
                        +vel_s_im[0][idx]*vel_s_im[0][idx]
                        +vel_s_im[1][idx]*vel_s_im[1][idx]
                        +vel_s_im[2][idx]*vel_s_im[2][idx]);

                    Ek_c_k+=0.5*(vel_c_re[0][idx]*vel_c_re[0][idx]
                        +vel_c_re[1][idx]*vel_c_re[1][idx]
                        +vel_c_re[2][idx]*vel_c_re[2][idx]
                        +vel_c_im[0][idx]*vel_c_im[0][idx]
                        +vel_c_im[1][idx]*vel_c_im[1][idx]
                        +vel_c_im[2][idx]*vel_c_im[2][idx]);
                }
            }
        }
    }

    myMPI.allReduceSum(Ek_s_k);
    myMPI.allReduceSum(Ek_c_k);

    double theta=std::sqrt(std::fabs(forcing_amp-Ek_c_k)/Ek_s_k);

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx=idx_j+i;

                double kx0=((kx+(double)i)>=(nx_global/2.0))?kx+(double)i-nx_global:kx+(double)i;
                double ky0=((ky+(double)j)>=(ny_global/2.0))?ky+(double)j-ny_global:ky+(double)j;
                double kz0=((kz+(double)k)>=(nz_global/2.0))?kz+(double)k-nz_global:kz+(double)k;
                double kmag=std::sqrt(kx0*kx0+ky0*ky0+kz0*kz0);

                if(fabs(kmag-forcing_k)<0.5)
                {
                    pu_re[idx]=vel_s_re[0][idx]*(theta-1.0);
                    pv_re[idx]=vel_s_re[1][idx]*(theta-1.0);                 
                    pw_re[idx]=vel_s_re[2][idx]*(theta-1.0);

                    pu_im[idx]=vel_s_im[0][idx]*(theta-1.0);
                    pv_im[idx]=vel_s_im[1][idx]*(theta-1.0);                 
                    pw_im[idx]=vel_s_im[2][idx]*(theta-1.0);
                }
                else
                {
                    pu_re[idx]=0.0;
                    pv_re[idx]=0.0;                 
                    pw_re[idx]=0.0;

                    pu_im[idx]=0.0;
                    pv_im[idx]=0.0;                 
                    pw_im[idx]=0.0;

                }

            }
        }
    }

    myFFT.solveFFT_backward(FFT_vector_re[0],FFT_vector_im[0],forcing_re[0],forcing_im[0],myMPI); // FFT of velocity u
    myFFT.solveFFT_backward(FFT_vector_re[1],FFT_vector_im[1],forcing_re[1],forcing_im[1],myMPI); // FFT of velocity v
    myFFT.solveFFT_backward(FFT_vector_re[2],FFT_vector_im[2],forcing_re[2],forcing_im[2],myMPI); // FFT of velocity w

    double *pa=A.getDataPtr();
    double A_aver=0.0;

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=1.0;
                // rho0=prho[idx_i];
                u0=pu[idx_i];
                v0=pv[idx_i];
                w0=pw[idx_i];

                double fx=rho0*pforcing_u_re[idx_i];
                double fy=rho0*pforcing_v_re[idx_i];
                double fz=rho0*pforcing_w_re[idx_i];
                pa[idx_i] = fx*u0+fy*v0+fz*w0;
            }
        }
    }

    myMPI.globalAverage(A,A_aver);

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=1.0;
                // rho0=prho[idx_i];
                u0=pu[idx_i];
                v0=pv[idx_i];
                w0=pw[idx_i];

                double fx=rho0*pforcing_u_re[idx_i];
                double fy=rho0*pforcing_v_re[idx_i];
                double fz=rho0*pforcing_w_re[idx_i];

                prhs1[idx_i]-= fx;
                prhs2[idx_i]-= fy;
                prhs3[idx_i]-= fz;
                prhs4[idx_i]-= fx*u0+fy*v0+fz*w0-A_aver;
            }
        }
    }

    for(int i=0;i<3;i++)
    {
        delete [] vel_c_re[i];
        delete [] vel_c_im[i];

        delete [] vel_s_re[i];
        delete [] vel_s_im[i];

        delete [] vel_l_re[i];
        delete [] vel_l_im[i];
        
        delete [] vel_r_re[i];
        delete [] vel_r_im[i];
    }

}

void nuc3d::physicsModel::AddForcing_ABC(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  VectorField &RHS,
  MPIComunicator3d_nonblocking &myMPI,
  fft &myFFT,
  const VectorField &xyz_center,
  double dt)
{
    double forcing_amp=myModelParameters["Forcing_amplitude"];
    double forcing_k=myModelParameters["Forcing_wave_k"];


    double *prho=W_vec[0].getDataPtr();
    double *pu=W_vec[1].getDataPtr();
    double *pv=W_vec[2].getDataPtr(); 
    double *pw=W_vec[3].getDataPtr(); 

    double *px=xyz_center[0].getDataPtr();
    double *py=xyz_center[1].getDataPtr();
    double *pz=xyz_center[2].getDataPtr();

    double *prhs1=RHS[1].getDataPtr();  
    double *prhs2=RHS[2].getDataPtr();  
    double *prhs3=RHS[3].getDataPtr();  
    double *prhs4=RHS[4].getDataPtr();  
    
    double rho0,u0,v0,w0;

    int nx=W_vec[0].getSizeX();
    int ny=W_vec[0].getSizeY();
    int nz=W_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=1.0;
                // rho0=prho[idx_i];
                u0=pu[idx_i];
                v0=pv[idx_i];
                w0=pw[idx_i];

                double fx=forcing_amp*rho0*(std::cos(forcing_k*py[idx_i])+1.1*std::sin(forcing_k*pz[idx_i]));
                double fy=forcing_amp*rho0*(1.1*std::cos(forcing_k*pz[idx_i])+0.9*std::sin(forcing_k*px[idx_i]));
                double fz=forcing_amp*rho0*(0.9*std::cos(forcing_k*px[idx_i])+std::sin(forcing_k*py[idx_i]));

                prhs1[idx_i]-= fx;
                prhs2[idx_i]-= fy;
                prhs3[idx_i]-= fz;
                prhs4[idx_i]-= fx*u0+fy*v0+fz*w0;
            }
        }
    }

}

void nuc3d::physicsModel::AddForcing_TG(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  VectorField &RHS,
  MPIComunicator3d_nonblocking &myMPI,
  fft &myFFT,
  const VectorField &xyz_center,
  double dt)
{
    double forcing_amp=myModelParameters["Forcing_amplitude"];
    double forcing_k=myModelParameters["Forcing_wave_k"];


    double *prho=W_vec[0].getDataPtr();
    double *pu=W_vec[1].getDataPtr();
    double *pv=W_vec[2].getDataPtr(); 
    double *pw=W_vec[3].getDataPtr(); 

    double *px=xyz_center[0].getDataPtr();
    double *py=xyz_center[1].getDataPtr();
    double *pz=xyz_center[2].getDataPtr();

    double *prhs1=RHS[1].getDataPtr();  
    double *prhs2=RHS[2].getDataPtr();  
    double *prhs3=RHS[3].getDataPtr();  
    double *prhs4=RHS[4].getDataPtr();  
    
    double rho0,u0,v0;

    int nx=W_vec[0].getSizeX();
    int ny=W_vec[0].getSizeY();
    int nz=W_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=1.0;
                // rho0=prho[idx_i];
                u0=pu[idx_i];
                v0=pv[idx_i];
                double fx=  forcing_amp*2.0*(std::sin(forcing_k*px[idx_i])*std::cos(forcing_k*py[idx_i])*std::sin(forcing_k*pz[idx_i]));
                double fy= -forcing_amp*2.0*(std::cos(forcing_k*px[idx_i])*std::sin(forcing_k*py[idx_i])*std::sin(forcing_k*pz[idx_i]));
                prhs1[idx_i]-=fx;
                prhs2[idx_i]-=fy; 
                //prhs3[idx_i]-=0.0;
                prhs4[idx_i]-= (fx*u0+fy*v0);
            }
        }
    }

}

void nuc3d::physicsModel::AddCoriolisForce(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  VectorField &RHS)
{
    double Rossby=myModelParameters["Rossby"];

    double *prho=W_vec[0].getDataPtr();
    double *pu=W_vec[1].getDataPtr();
    double *pv=W_vec[2].getDataPtr(); 

    double *prhs1=RHS[1].getDataPtr();  
    double *prhs2=RHS[2].getDataPtr();  
    
    double rho0,u0,v0;

    int nx=W_vec[0].getSizeX();
    int ny=W_vec[0].getSizeY();
    int nz=W_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=prho[idx_i];
                u0=pu[idx_i];
                v0=pv[idx_i];

                prhs1[idx_i]-= rho0*Rossby*v0;
                prhs2[idx_i]-=-rho0*Rossby*u0;
            }
        }
    }
}

void nuc3d::physicsModel::AddCentrifugalForce(const VectorField &W_vec,
  const VectorField &W0_vec,
  const VectorField &Q_vec,
  const VectorField &xyz_center,
  VectorField &RHS)
{
    double Rossby=myModelParameters["Rossby"];
    double xc=myModelParameters["x_center"];
    double yc=myModelParameters["y_center"];
    double r0=myModelParameters["well_rad"];

    double *prho=W_vec[0].getDataPtr();
    double *pu=W_vec[1].getDataPtr();
    double *pv=W_vec[2].getDataPtr(); 


    double *prhs1=RHS[1].getDataPtr();  
    double *prhs2=RHS[2].getDataPtr();  
    double *prhs3=RHS[3].getDataPtr();  
    double *prhs4=RHS[4].getDataPtr(); 


    double *px=xyz_center[0].getDataPtr();
    double *py=xyz_center[1].getDataPtr();


    double rho0,u0,v0;
    double x0,y0;


    int nx=W_vec[0].getSizeX();
    int ny=W_vec[0].getSizeY();
    int nz=W_vec[0].getSizeZ();

    for (int k = 0; k < nz; ++k)
    {
        int idx_k=nx*ny*k;
        for (int j = 0; j < ny; ++j)
        {
            int idx_j=idx_k+nx*j;
            for (int i = 0; i < nx; ++i)
            {                
                int idx_i=idx_j+i;

                rho0=prho[idx_i];

                u0=pu[idx_i];
                v0=pv[idx_i];

                x0=px[idx_i];
                y0=py[idx_i];
                /*
                    Let us an centrifugal force well: 

                            1/[1.0+(r/r0)^20]

                    to control the effective range of 
                    the centrifugal force.
                */
                    double r=std::sqrt((x0-xc)*(x0-xc)+(y0-yc)*(y0-yc));
                    double fx=rho0*Rossby*Rossby*(x0-xc)/(1.0+std::pow(r/r0,20))/4.0;
                    double fy=rho0*Rossby*Rossby*(y0-yc)/(1.0+std::pow(r/r0,20))/4.0;

                    prhs1[idx_i]-= fx;
                    prhs2[idx_i]-= fy;

                    prhs2[idx_i]-= fx*u0+fy*v0;
                }
            }
        }

    }

/**************************************************************************************
 End of definition
 
 **************************************************************************************/
#endif
