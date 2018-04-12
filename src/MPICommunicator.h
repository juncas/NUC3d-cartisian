#ifndef MPICommunicator3d_h
#define MPICommunicator3d_h
# include "field.h"
#include "MPICommunicatorBase.h"

namespace nuc3d
{
    class bufferData;
    
    class MPIComunicator3d_nonblocking : public MPICommunicator
    {


        int* sendTags;
        int* recvTags;

        MPI_Comm cart_comm;

        int cart_dim[3]; //Total processes number in each direction
        int cart_coord[3]; //Process ID in cartisian frame
        int neighbours[3][2];
        
    public:
        MPIComunicator3d_nonblocking();
        ~MPIComunicator3d_nonblocking();
        
    public:
        /*
         Different from the curved NUC3d, in this version, the 
         setTopo function is used to set global cartisian coordinate.
        */
        void setTopo( int,  int , int);
        void bufferSendRecv(Field &, bufferData &,int,int);
        void waitSendRecv(bufferData& buffer,int);
        void ReduceData(double );
        int getMyCoord(int );        
        int getMyDim(int );
        void allReduceMin(double &val);
        void allReduceMax(double &val);
        void allReduceSum(double &val);
        void allReduceSum(double *val_local,double *val_global, int n);
        void globalAverage(const Field &field,double &val);       
        void globalAverage(const double *f,double &val, int nx, int ny, int nz); 
        void globalAveragePower2(const Field &field,double &val);
        void SplitWorld(MPI_Comm *,MPI_Comm *,MPI_Comm *);
        void barrier();
    private:
        void bufferSendRecvFace(Field &, bufferData &,int,int);
        void waitSendRecvFace(bufferData& buffer,int);
        
        
        int setMPISendRecvTag(
                          int ProcId,
                          int FaceId,
                          int FieldID)
        {
            return FieldID + FaceId * 100 + ProcId * 1000;
        };
                
    };
}
#endif
