//
//  postproc.hpp
//  NUC3d
//
//  Created by Jun Peng on 16/1/4.
//  Copyright © 2016年 Jun Peng. All rights reserved.
//
/*
    This class is kind of user defined function for postprocessing,
    All data are available for this class, so you can do any processing during
    the calculation.
*/
#ifndef postproc_hpp
#define postproc_hpp

#include <stdio.h>
#include "field.h"
#include "mpi.h"
#include "bufferData.hpp"


namespace nuc3d
{
  class fieldOperator3d;
  class bufferData;
  class MPIComunicator3d_nonblocking;
  class physicsModel;
  class IOController;  
  class fft;
  
  
  class postproc
  {
    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;
    int postSize;
    
    typedef void (postproc::*pSolveAveraged)(VectorField &prims,
     VectorField &acous,
     physicsModel &myPhys,
     fieldOperator3d &myOP,
     VectorBuffer &myBf,
     MPIComunicator3d_nonblocking &myMPI);
    
    const std::vector<std::string> VarNameList_primary={
      "rho","u","v","w","p","T","e","alpha","miu","coeff"       };
      
      const std::vector<std::string> VarNameList_grads={
        "dudx","dudy","dudz",
        "dvdx","dvdy","dvdz",
        "dwdx","dwdy","dwdz",
        "dpdx","dpdy","dpdz",
        "dTdx","dTdy","dTdz",
        "omega0","omega1","omega2",
        "omega0^2","omega1^2","omega2^2"
      };
      
      
      const std::vector<std::string> VarNameList_derived={
        "rhou","rhov","rhow",
        "rhoe",
        "uu","uv","uw",
        "vv","vw","ww",
        "rhouu","rhouv","rhouw",
        "rhovv","rhovw","rhoww",
        "tau_xx","tau_xy","tau_xz",
        "tau_yy","tau_yz","tau_zz",
        "rhotau_xx","rhotau_xy","rhotau_xz",
        "rhotau_yy","rhotau_yz","rhotau_zz"
      };
      
      const std::vector<std::string> VarNameList_TKE={
        "up","vp","wp",
        "pdudx","pdvdy","pdwdz",
        "dtau_xxdx","dtau_yxdx","dtau_zxdx",
        "dtau_xydy","dtau_yydy","dtau_zydy",
        "dtau_xzdz","dtau_yzdz","dtau_zzdz",
        "utau_xx","vtau_yx","wtau_zx","uktau_kx","duktau_kxdx",
        "utau_xy","vtau_yy","wtau_zy","uktau_ky","duktau_kydy",
        "utau_xz","vtau_yz","wtau_zz","uktau_kz","duktau_kzdz",
        "tau_xxdudx","tau_yxdvdx","tau_zxdwdx",
        "tau_xydudy","tau_yydvdy","tau_zydwdy",
        "tau_xzdudz","tau_yzdvdz","tau_zzdwdz",
        "rhouuu","rhouuv","rhouuw",
        "rhovuv","rhouvw",
        "rhouww","rhovvv","rhovvw",
        "rhovww","rhowww",
      };
      
      
        // uv_t denote u'v'
      const std::vector<std::string> VarNameList_TKEbudget={
        "tke",
        "production",
        "turbulen_transportation",
        "pressure_diffusion",
        "p1",
        "p2",
        "p3",
        "viscous_diffusion",
        "v1",
        "v2",
        "v3",
        "convection",
        "drhoudx","drhoudy","drhoudz",
        "drhovdx","drhovdy","drhovdz",
        "drhowdx","drhowdy","drhowdz",
        "up_t","vp_t","wp_t",
        "dupdx_t","dvpdy_t","dwpdz_t",
        "uktau_kx_t","uktau_ky_t","uktau_kz_t",
        "duktau_kxdx_t","duktau_kydy_t","duktau_kzdz_t",
        "rhouk","rhovk","rhowk",
        "drhoukdx","drhovkdy","drhowkdz",
        "u_tke","v_tke","w_tke",
        "du_tkedx","dv_tkedy","dw_tkedz"
      };
      
      const std::vector<std::string> VarNameList_q={
        "dudx","dudy","dudz",
        "dvdx","dvdy","dvdz",
        "dwdx","dwdy","dwdz",
        "omega0","omega1","omega2",
        "enstrophy",
        "Q",
        "helicity",
        "kinetic"
        };// other kinds of simutaneous variables can be added here
                
        Field temp_xi,temp_eta,temp_zeta;
        Field dtempdxi,dtempdeta,dtempdzeta;
        
        VectorField TemporalPrimaryField;
        VectorField TemporalGradsField;
        VectorField TemporalDerivedField;
        VectorField TemporalTkeField;
        VectorField TemporalQField;
        
        VectorField AveragedPrimaryField;
        VectorField AveragedGradsField;
        VectorField AveragedDerivedField;
        VectorField AveragedTkeField;
        
        VectorField TKEbudgetField;

        VectorField FFT_vector_re;
        VectorField FFT_vector_im;
        
        double enstrophy_glb,enstrophy;
        double kinetic_glb,kinetic;
        
        std::vector<int> variableScalarInt;
        std::vector<double> variableScalarDouble;
        
        int averStep;
        
        pSolveAveraged myAverage[3]={
          &postproc::solveAveraged0,
          &postproc::solveAveraged1,
          &postproc::solveAveraged2
        };
        
      public:
        postproc(int,int,int,double,double,double);
        ~postproc();
        
        void initPost();
        
        void solvePost(VectorField &prims,
         VectorField &acous,
         VectorField &xyz,
         physicsModel &myPhys,
         fieldOperator3d &myOP,
         VectorBuffer &myBf,
         MPIComunicator3d_nonblocking &myMPI,
         IOController &myIO,
         fft &myFFT,
         int istep,
         double time);
        
        void OutputPost(VectorField &prims,
          VectorField &acous,
          VectorField &xyz,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI,
          IOController &myIO,
          int istep,
          double time);
      private:
        void solveGrad(Field &myField,
         fieldOperator3d &myOP,
         bufferData &myBf,
         MPIComunicator3d_nonblocking &myMPI,
         int fdID,
         Field &dfdxi,
         Field &dfdeta,
         Field &dfdzeta);
        
        void solveGrad_xyz(Field &Q_x,
         Field &Q_y,
         Field &Q_z,
         Field &Q_xi,
         Field &Q_eta,
         Field &Q_zeta);
        
        void solveGrad_x(Field &Q_x,
         Field &Q_xi,
         Field &Q_eta,
         Field &Q_zeta);
        
        void solveGrad_y(Field &Q_y,
         Field &Q_xi,
         Field &Q_eta,
         Field &Q_zeta);
        
        void solveGrad_z(Field &Q_z,
         Field &Q_xi,
         Field &Q_eta,
         Field &Q_zeta);
        
        void setField(const Field &f);
        
        void solveGrad_df(Field &myField,
         double dxyz,
         fieldOperator3d &myOP,
         bufferData &myBf,
         MPIComunicator3d_nonblocking &myMPI,
         Field &df,
         int dir,
         int fdID);
        
        void solveTemporal(VectorField &prims,
         VectorField &acous,
         physicsModel &myPhys,
         fieldOperator3d &myOP,
         VectorBuffer &myBf,
         MPIComunicator3d_nonblocking &myMPI,
         fft &myFFT,
         int istep,
         double time);
        
        void solveQ(VectorField &prims,
          VectorField &acous,
          physicsModel &myPhys,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI,
          int istep,
          double time);
        
        void solveFFT(VectorField &prims,
          VectorField &acous,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI,
          fft &myFFT,
          int istep,
          double time);

// FFT functions for shell
        void solveFFT_kinetic(VectorField &prims,
        VectorField &acous,
        fieldOperator3d &myOP,
        VectorBuffer &myBf,
        MPIComunicator3d_nonblocking &myMPI,
        fft &myFFT,
        int istep,
        double time);
        void solveFFT_helicity(VectorField &prims,
        VectorField &acous,
        fieldOperator3d &myOP,
        VectorBuffer &myBf,
        MPIComunicator3d_nonblocking &myMPI,
        fft &myFFT,
        int istep,
        double time);
        void solveFFT_enstrophy(VectorField &prims,
        VectorField &acous,
        fieldOperator3d &myOP,
        VectorBuffer &myBf,
        MPIComunicator3d_nonblocking &myMPI,
        fft &myFFT,
        int istep,
        double time);
        void solveFFT_density(VectorField &prims,
        VectorField &acous,
        fieldOperator3d &myOP,
        VectorBuffer &myBf,
        MPIComunicator3d_nonblocking &myMPI,
        fft &myFFT,
        int istep,
        double time);
        void solveFFT_pressure(VectorField &prims,
        VectorField &acous,
        fieldOperator3d &myOP,
        VectorBuffer &myBf,
        MPIComunicator3d_nonblocking &myMPI,
        fft &myFFT,
        int istep,
        double time);
//************************************************************************
        void solveAveraged(VectorField &prims,
         VectorField &acous,
         physicsModel &myPhys,
         fieldOperator3d &myOP,
         VectorBuffer &myBf,
         MPIComunicator3d_nonblocking &myMPI);
        
        void solveAveraged0(VectorField &prims,
          VectorField &acous,
          physicsModel &myPhys,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI);
        
        void solveAveraged1(VectorField &prims,
          VectorField &acous,
          physicsModel &myPhys,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI);
        
        void solveAveraged2(VectorField &prims,
          VectorField &acous,
          physicsModel &myPhys,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI);
        
        void OutputTemporal(VectorField &prims,
          VectorField &acous,
          VectorField &xyz,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI,
          IOController &myIO,
          int istep,
          double time);
        
        void OutputAveraged_bin(VectorField &prims,
          VectorField &acous,
          VectorField &xyz,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI,
          IOController &myIO,
          int istep,
          double time);
        
        void OutputAveraged_tec(VectorField &prims,
          VectorField &acous,
          VectorField &xyz,
          fieldOperator3d &myOP,
          VectorBuffer &myBf,
          MPIComunicator3d_nonblocking &myMPI,
          IOController &myIO,
          int istep,
          double time);

      };
    }

#endif /* postproc_hpp */
