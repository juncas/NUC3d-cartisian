#ifndef MPICommunicator_h
#define MPICommunicator_h

# include "mpi.h"
namespace nuc3d 
{
	class MPICommunicator// resposible for transportation only
	{

		
	public:
		int nProc; // Total number of Processes
		
		int myProc; //Global process ID
		MPICommunicator();
		~MPICommunicator();
	public:
		int getMyId() const {return myProc;};
		int getSize() const {return nProc;};
		void FinializeMPI();
		void AbortMPI();
	};
}
#endif