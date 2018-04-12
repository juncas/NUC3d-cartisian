#ifndef fft_h
#define fft_h
#include "field.h"
#include "MPICommunicator.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <memory>
#include "fftw3.h"
#include "fftw3-mpi.h"
/*
	three dimensional fully parallel FFT for multidimensinal block is 
	very challenging. In this class, we use the FFTW3 library to calulate 
	the fft of a field. this the solve funcion is expected to be very slow.
	It should be improved in future versions.
*/

namespace nuc3d
{
	class fft
	{
		int nx,ny,nz;
		int nx_global,ny_global,nz_global;

		fftw_plan forward_plan_x;
		fftw_plan forward_plan_y;
		fftw_plan forward_plan_z;

		fftw_complex *forward_data_x_in;
		fftw_complex *forward_data_y_in;
		fftw_complex *forward_data_z_in;

		fftw_complex *forward_data_x_out;
		fftw_complex *forward_data_y_out;
		fftw_complex *forward_data_z_out;

		fftw_plan backward_plan_x;
		fftw_plan backward_plan_y;
		fftw_plan backward_plan_z;

		fftw_complex *backward_data_x_in;
		fftw_complex *backward_data_y_in;
		fftw_complex *backward_data_z_in;

		fftw_complex *backward_data_x_out;
		fftw_complex *backward_data_y_out;
		fftw_complex *backward_data_z_out;

		ptrdiff_t forward_alloc_local_x,backward_alloc_local_x, local_nx_in, local_nx_start_in, local_nx_out, local_nx_start_out;
		ptrdiff_t forward_alloc_local_y,backward_alloc_local_y, local_ny_in, local_ny_start_in, local_ny_out, local_ny_start_out;
		ptrdiff_t forward_alloc_local_z,backward_alloc_local_z, local_nz_in, local_nz_start_in, local_nz_out, local_nz_start_out;

		MPI_Comm comm_x;
		MPI_Comm comm_y;
		MPI_Comm comm_z;

		double *temp_re;
		double *temp_im;
	public:
		
		fft();
		~fft();
	public:
		void initFFT(MPIComunicator3d_nonblocking &myMPI,int nx,int ny,int nz);

		void solveFFT_forward(const Field &,
			Field &,
			Field &,
			MPIComunicator3d_nonblocking &);

		void solveFFT_backward(const Field &,
			const Field &,
			Field &,
			Field &,
			MPIComunicator3d_nonblocking &);

		void solveShell(const Field &,
			const Field &,
			MPIComunicator3d_nonblocking &,
			std::string name,
			int istep,
			double time);

		void solveShell(const double *pf_re_in,
			MPIComunicator3d_nonblocking &myMPI,
			std::string name,
			int istep,
			double time);

	private:
		void transpose_i2j(const double *, double *);
		void transpose_j2k(const double *, double *);
		void transpose_k2i(const double *, double *);
		void normalization(double *f,int fac);
	};
}
#endif
