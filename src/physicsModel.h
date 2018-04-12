#ifndef physicsModel_h
#define physicsModel_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <memory>
#include "field.h"

namespace nuc3d
{
  class EulerData3D;
  class EulerFlux;
  class PDEData3d;
  class MPIComunicator3d_nonblocking;
  class fft;

  class physicsModel
  {                
    typedef void (physicsModel::*pRiemann)(const VectorField &W_vec,
      const VectorField &W0_vec,
      const VectorField &Q_vec,
      EulerFlux& myFlux_xi,
      EulerFlux& myFlux_eta,
      EulerFlux& myFlux_zeta);

    typedef void (physicsModel::*pVisMod)(const Field &T,
      Field &miu,
      Field &coeff,
      double Reynolds,
      double Mach,
      double pt,
      double gamma,
      double T_ref,
      double T_inf);

    typedef void (physicsModel::*pForcing)(const VectorField &W_vec,
      const VectorField &W0_vec,
      const VectorField &Q_vec,
      VectorField &source,
      MPIComunicator3d_nonblocking &myMPI,
      fft &myFFT,
      const VectorField &xyz_center,
      double dt);

    int neqs;
    std::string myEoSName;
    std::string myRiemannName;
    std::string myModelName;
    std::string myVisModelName;
    std::string myForcingName;

    std::map<std::string, pRiemann> myRiemannMap;
    std::map<std::string, pVisMod> myVisModelMap;
    std::map<std::string, pForcing> myForcingMap;

    std::map<std::string, std::string> mySourceSwitchMap;        
        std::map<std::string, double> myModelParameters; //parameters for EoS
      public:
        physicsModel();
        ~physicsModel();
        
        std::string getMyModelName(){return myModelName;};
        
        void solve(PDEData3d &, EulerData3D *);
        
        void solveRiemann(PDEData3d &myPDE,EulerData3D *myEuler);

        void solveSource(PDEData3d &myPDE,
          EulerData3D *myEuler,
          MPIComunicator3d_nonblocking &,
          fft &myFFT,
          const VectorField &xyz_center);

        void getMiu(Field &,
          Field &,
          Field &);
        
        int getEqNum(){return neqs;};
        double getMach(){return myModelParameters["Mach"];};
        double getGamma(){return myModelParameters["Gamma"];};
        double getWallTemp(){return myModelParameters["T_wall"];};
        double getValue(std::string name){return myModelParameters[name];};
        
        void initial(PDEData3d &,std::shared_ptr<EulerData3D> );

        void getPrim(VectorField &Q,VectorField &Prim, VectorField &Acoust);
      private:

        std::istream& readPhysFile(std::istream& ios);
//********************************************************    
        /*
          transformation between primative and conservative variables
        */
        void con2prim(const std::string &,
          VectorField &,
          VectorField &,
          VectorField &);
        
        void prim2con(const std::string &,
          const VectorField &,
          VectorField &);
//******************************************************
        /*
          Vicosity functions 
        */
        void sutherland(const Field &T,
          Field &miu,
          Field &coeff,
          double Reynolds,
          double Mach,
          double pt,
          double gamma,
          double T_ref,
          double T_inf);
        
        void constant(const Field &T,
          Field &miu,
          Field &coeff,
          double Reynolds,
          double Mach,
          double pt,
          double gamma,
          double T_ref,
          double T_inf);
//********************************************************   
        /*
        Riemann solvers and asociated functions
        */
        void RiemannAUSM(const VectorField &W_vec,
         const VectorField &W0_vec,
         const VectorField &Q_vec,
         EulerFlux& myFlux_xi,
         EulerFlux& myFlux_eta,
         EulerFlux& myFlux_zeta);
        
        void RiemannAUSMp(const VectorField &W_vec,
         const VectorField &W0_vec,
         const VectorField &Q_vec,
         EulerFlux& myFlux_xi,
         EulerFlux& myFlux_eta,
         EulerFlux& myFlux_zeta);
        
        void RiemannLF(const VectorField &W_vec,
         const VectorField &W0_vec,
         const VectorField &Q_vec,
         EulerFlux& myFlux_xi,
         EulerFlux& myFlux_eta,
         EulerFlux& myFlux_zeta);

        double getPressureL(const double &, const double &);
        double getPressureR(const double &, const double &);
        
        double getMachL(const double &);
        double getMachR(const double &);
        
        double getPressureLp(const double &, const double &);
        double getPressureRp(const double &, const double &);
        
        double getMachLp(const double &);
        double getMachRp(const double &);          
//*************************************************************
        /*
          forcing functions
        */
        void AddForcing(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          VectorField &RHS,
          MPIComunicator3d_nonblocking &myMPI,
          fft &myFFT,
  const VectorField &xyz_center,
          double dt);

        void AddForcing_None(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          VectorField &RHS,
          MPIComunicator3d_nonblocking &myMPI,
          fft &myFFT,
  const VectorField &xyz_center,
          double dt);

        void AddForcing_Const(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          VectorField &RHS,
          MPIComunicator3d_nonblocking &myMPI,
          fft &myFFT,
  const VectorField &xyz_center,
          double dt);

        void AddForcing_ABC(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          VectorField &RHS,
          MPIComunicator3d_nonblocking &myMPI,
          fft &myFFT,
  const VectorField &xyz_center,
          double dt);

        void AddForcing_TG(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          VectorField &RHS,
          MPIComunicator3d_nonblocking &myMPI,
          fft &myFFT,
  const VectorField &xyz_center,
          double dt);

        //...and other types

//****************************************************************
/*
      In this version, only z-direction angle velocity is considered.

      The Rossby Number is given by : Ro=L \Omega / U
      Thus the coriolis force is given by Ro*(0, 0, 1) X (u, v, w)
      ant the centrifugal force is given by Ro^2*(0, 0, 1) X [(0, 0, 1) X (x, y, z)].
      Note that the cnetrifugal forced in incompatible with periodic boundary
      conditions, however, we adjust boundary conditions in future versions.
*/
        /*
          coriolis force
        */
        void AddCoriolisForce(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          VectorField &RHS);

        /*
          centrifugal force
        */
        void AddCentrifugalForce(const VectorField &W_vec,
          const VectorField &W0_vec,
          const VectorField &Q_vec,
          const VectorField &xyz_center,
          VectorField &RHS);

//****************************************************************

        
      };
    }
#endif
