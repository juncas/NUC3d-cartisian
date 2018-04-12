#ifndef MPICommunicator3d_cpp
#define MPICommunicator3d_cpp
#include "MPICommunicator.h"
#include "bufferData.hpp"
# include <iostream>
#define BoundaryFace (-1)
/**************************************************************************************
	Definition of class MPIComunicator3d_nonblocking: non_blocking coommunicator
 **************************************************************************************/
/**************************************************************************************
 Definition of constructors and destructors
 **************************************************************************************/
 nuc3d::MPIComunicator3d_nonblocking::MPIComunicator3d_nonblocking() :
 MPICommunicator()
 {

 };

 nuc3d::MPIComunicator3d_nonblocking::~MPIComunicator3d_nonblocking()
 {
    MPI_Comm_free(&cart_comm);
 };
/**************************************************************************************
 Definition of member functions
 **************************************************************************************/
 void nuc3d::MPIComunicator3d_nonblocking::setTopo(int npx,int npy,int npz)
 {
  if((npx*npy*npz)==nProc)
  {
    cart_dim[0]=npz;
    cart_dim[1]=npy;
    cart_dim[2]=npx;
    int period[3];
    period[0]=1;
    period[1]=1;
    period[2]=1;

        //creating cartisian process topo
        std::cout << "My MPI Id is "<<myProc<<"!"<< std::endl;
        MPI_Cart_create(MPI_COMM_WORLD, 3, cart_dim, period, 1, &cart_comm);
        MPI_Comm_rank(cart_comm, &myProc); // renew procID for possible rank reorder
        MPI_Cart_coords(cart_comm,myProc, 3, cart_coord);

        //set neighbour process numbers
        MPI_Cart_shift(cart_comm, 2, 1, &neighbours[0][0], &neighbours[0][1]); // x direction
        MPI_Cart_shift(cart_comm, 1, 1, &neighbours[1][0], &neighbours[1][1]); // y direction
        MPI_Cart_shift(cart_comm, 0, 1, &neighbours[2][0], &neighbours[2][1]); // z direction

        std::cout << "My MPI Id (px,py,pz) is ("<<cart_coord[2]<<","<<cart_coord[1]<<","<<cart_coord[0]<<"),"<<myProc<<"!"<< std::endl;
      }
      else
      {
        std::cout << "ERROR:process number ("<<npx<<","<<npy<<","<<npz<<") does not match total process number"<<nProc<<"!"<< std::endl;
        MPI_Abort(cart_comm, 999);
      }
    };

    int nuc3d::MPIComunicator3d_nonblocking::getMyCoord(int dir)
    {
      if(dir>2||dir<0)
        {
          std::cout<<"No such direction dir="<<dir<<std::endl;
          exit(-1);          
        }
      else
        return cart_coord[dir];
    };

    int nuc3d::MPIComunicator3d_nonblocking::getMyDim(int dir)
    {
      if(dir>2||dir<0)
        {
          std::cout<<"No such direction dir="<<dir<<std::endl;
          exit(-1);          
        }
      else
        return cart_dim[dir];
    };


    void nuc3d::MPIComunicator3d_nonblocking::bufferSendRecv(Field &myfield, bufferData& buffer,int dir, int FieldID)
    {
      for(int lr=0;lr<2;lr++)
      {
        if(neighbours[dir][lr]!=MPI_PROC_NULL)
        {
          int sourceFaceID=dir*2+lr;
          int destFaceID=dir*2+1-lr;

          int sendTag = setMPISendRecvTag(this->getMyId(), sourceFaceID, FieldID);
          int RecvTag = setMPISendRecvTag(neighbours[dir][lr], destFaceID, FieldID);

          buffer.setBufferSend(myfield, sourceFaceID);

          MPI_Irecv(buffer.BufferRecv[sourceFaceID].getDataPtr(),
            buffer.bufferSize[sourceFaceID],
            MPI_DOUBLE,
            neighbours[dir][lr],
            RecvTag,
            cart_comm,
            &(buffer.myRequestRecv[sourceFaceID]));

          MPI_Isend(buffer.BufferSend[sourceFaceID].getDataPtr(),
            buffer.bufferSize[sourceFaceID],
            MPI_DOUBLE,
            neighbours[dir][lr],
            sendTag,
            cart_comm,
            &(buffer.myRequestSend[sourceFaceID]));
        } 
      }

    };

    void nuc3d::MPIComunicator3d_nonblocking::waitSendRecv(bufferData& buffer,int dir)
    {
      for(int lr=0;lr<2;lr++)
      {
        if(neighbours[dir][lr]!=MPI_PROC_NULL)
        {
          int sourceFaceID=dir*2+lr;

          if((MPI_Wait(&buffer.myRequestSend[sourceFaceID], &buffer.myStatusSend[sourceFaceID])!= MPI_SUCCESS)
           ||(MPI_Wait(&buffer.myRequestRecv[sourceFaceID], &buffer.myStatusRecv[sourceFaceID])!= MPI_SUCCESS))
          {
            std::cout << "Communicating failed" << std::endl;
            MPI_Abort(cart_comm, 999);
          }
        } 
      }

    };

    void nuc3d::MPIComunicator3d_nonblocking::allReduceMin(double &val)
    {  
      double val_global;
     MPI_Allreduce(&val, &val_global, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
     val=val_global;
    }

    void nuc3d::MPIComunicator3d_nonblocking::allReduceMax(double &val)
    {  
      double val_global;
      MPI_Allreduce(&val, &val_global, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
      val=val_global;
    }

    void nuc3d::MPIComunicator3d_nonblocking::allReduceSum(double &val)
    {  
      double val_global;
      MPI_Allreduce(&val, &val_global, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
      val=val_global;
    }

    void nuc3d::MPIComunicator3d_nonblocking::barrier()
    {  
      MPI_Barrier(cart_comm);
    }


    void nuc3d::MPIComunicator3d_nonblocking::globalAverage(const Field &field,double &val)
    {  
     int nx=field.getSizeX();
     int ny=field.getSizeY();
     int nz=field.getSizeZ();

     double sum_local=0.0;
     double sum_global=0.0;
     double *f=field.getDataPtr();

     for (int k=0; k<nz; k++)
     {
      for (int j=0; j<ny; j++)
      {
        for (int i=0; i<nx; i++)
        {
          int idx_xi=nx*ny*k+nx*j+i;
          sum_local+=f[idx_xi];
        }
      }
    }
    MPI_Allreduce(&sum_local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, cart_comm);

    val=sum_global/(cart_dim[0]*cart_dim[1]*cart_dim[2]*nx*ny*nz);
  }

    void nuc3d::MPIComunicator3d_nonblocking::globalAverage(const double *f,double &val, int nx, int ny, int nz)
    {  
     double sum_local=0.0;
     double sum_global=0.0;
     const double *f0=f;

     for (int k=0; k<nz; k++)
     {
      for (int j=0; j<ny; j++)
      {
        for (int i=0; i<nx; i++)
        {
          sum_local+=*(f0++);
        }
      }
    }
    MPI_Allreduce(&sum_local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, cart_comm);

    val=sum_global/(cart_dim[0]*cart_dim[1]*cart_dim[2]*nx*ny*nz);
  }

    void nuc3d::MPIComunicator3d_nonblocking::globalAveragePower2(const Field &field,double &val)
    {  
     int nx=field.getSizeX();
     int ny=field.getSizeY();
     int nz=field.getSizeZ();

     double sum_local=0.0;
     double sum_global=0.0;
     double *f=field.getDataPtr();

     for (int k=0; k<nz; k++)
     {
      for (int j=0; j<ny; j++)
      {
        for (int i=0; i<nx; i++)
        {
          int idx_xi=nx*ny*k+nx*j+i;
          sum_local+=f[idx_xi]*f[idx_xi];
        }
      }
    }
    MPI_Allreduce(&sum_local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, cart_comm);

    val=sum_global/(cart_dim[0]*cart_dim[1]*cart_dim[2]*nx*ny*nz);
  }


    void nuc3d::MPIComunicator3d_nonblocking::SplitWorld(MPI_Comm *new_comm_x,MPI_Comm *new_comm_y,MPI_Comm *new_comm_z)
    {  
      int color_x,color_y,color_z;

      color_x=cart_dim[1]*cart_coord[0]+cart_coord[1]; //npy*k+j i
      color_y=cart_dim[0]*cart_coord[2]+cart_coord[0]; //npz*i+k j
      color_z=cart_dim[2]*cart_coord[1]+cart_coord[2]; //npx*j+i k

      MPI_Comm_split(cart_comm, color_x, cart_coord[2], new_comm_x);
      MPI_Comm_split(cart_comm, color_y, cart_coord[1], new_comm_y);
      MPI_Comm_split(cart_comm, color_z, cart_coord[0], new_comm_z);

    }

    void nuc3d::MPIComunicator3d_nonblocking::allReduceSum(double *val_local,double *val_global, int n)
    {  
      MPI_Allreduce(val_local, val_global, n, MPI_DOUBLE, MPI_SUM, cart_comm);
    }
/**************************************************************************************
 End of definition
 **************************************************************************************/
#endif
