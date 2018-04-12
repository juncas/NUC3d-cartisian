#ifndef fft_cpp
#define fft_cpp
#include "fft.h"
#include <cmath>
nuc3d::fft::fft()
{}

nuc3d::fft::~fft()
{

	fftw_destroy_plan( forward_plan_x);
	fftw_destroy_plan( forward_plan_y);
	fftw_destroy_plan( forward_plan_z);

	fftw_free( forward_data_x_in);
	fftw_free( forward_data_y_in);
	fftw_free( forward_data_z_in);
	fftw_free( forward_data_x_out);
	fftw_free( forward_data_y_out);
	fftw_free( forward_data_z_out);

	fftw_destroy_plan( backward_plan_x);
	fftw_destroy_plan( backward_plan_y);
	fftw_destroy_plan( backward_plan_z);

	fftw_free( backward_data_x_in);
	fftw_free( backward_data_y_in);
	fftw_free( backward_data_z_in);
	fftw_free( backward_data_x_out);
	fftw_free( backward_data_y_out);
	fftw_free( backward_data_z_out);

	delete [] temp_re;
	delete [] temp_im;
	
	MPI_Comm_free(&comm_x);
	MPI_Comm_free(&comm_y);
	MPI_Comm_free(&comm_z);
}

void nuc3d::fft::initFFT(MPIComunicator3d_nonblocking &myMPI,int nx0,int ny0,int nz0)
{
	nx_global=nx0*myMPI.getMyDim(2);
	ny_global=ny0*myMPI.getMyDim(1);
	nz_global=nz0*myMPI.getMyDim(0);
	
	nx=nx0;
	ny=ny0;
	nz=nz0;

	myMPI.SplitWorld(&comm_x,&comm_y,&comm_z);

	fftw_mpi_init();

	temp_re=new double [nx*ny*nz];
	temp_im=new double [nx*ny*nz];
//forward data
	// for x direction 
	forward_alloc_local_x = fftw_mpi_local_size_1d(nx_global, 
		comm_x, 
		FFTW_FORWARD,
		FFTW_MEASURE,
		&local_nx_in, 
		&local_nx_start_in,
		&local_nx_out, 
		&local_nx_start_out);

	forward_data_x_in = fftw_alloc_complex( forward_alloc_local_x);
	forward_data_x_out = fftw_alloc_complex( forward_alloc_local_x);

	forward_plan_x = fftw_mpi_plan_dft_1d(nx_global, 
		forward_data_x_in, 
		forward_data_x_out, 
		comm_x,
		FFTW_FORWARD, 
		FFTW_MEASURE);

	// for y direction 
	forward_alloc_local_y = fftw_mpi_local_size_1d(ny_global, 
		comm_y, 
		FFTW_FORWARD,
		FFTW_MEASURE,
		&local_ny_in, 
		&local_ny_start_in,
		&local_ny_out, 
		&local_ny_start_out);

	forward_data_y_in = fftw_alloc_complex( forward_alloc_local_y);
	forward_data_y_out = fftw_alloc_complex( forward_alloc_local_y);

	forward_plan_y = fftw_mpi_plan_dft_1d(ny_global, 
		forward_data_y_in, 
		forward_data_y_out, 
		comm_y,
		FFTW_FORWARD, 
		FFTW_MEASURE);


	// for z direction 
	forward_alloc_local_z = fftw_mpi_local_size_1d(nz_global, 
		comm_z, 
		FFTW_FORWARD,
		FFTW_MEASURE,
		&local_nz_in, 
		&local_nz_start_in,
		&local_nz_out, 
		&local_nz_start_out);

	forward_data_z_in = fftw_alloc_complex( forward_alloc_local_z);
	forward_data_z_out = fftw_alloc_complex( forward_alloc_local_z);

	forward_plan_z = fftw_mpi_plan_dft_1d(nz_global, 
		forward_data_z_in, 
		forward_data_z_out, 
		comm_z,
		FFTW_FORWARD, 
		FFTW_MEASURE);
//backward data

	// for x direction 
	backward_alloc_local_x = fftw_mpi_local_size_1d(nx_global, 
		comm_x, 
		FFTW_BACKWARD,
		FFTW_MEASURE,
		&local_nx_in, 
		&local_nx_start_in,
		&local_nx_out, 
		&local_nx_start_out);

	backward_data_x_in = fftw_alloc_complex( backward_alloc_local_x);
	backward_data_x_out = fftw_alloc_complex( backward_alloc_local_x);

	backward_plan_x = fftw_mpi_plan_dft_1d(nx_global, 
		backward_data_x_in, 
		backward_data_x_out, 
		comm_x,
		FFTW_BACKWARD, 
		FFTW_MEASURE);

	// for y direction 
	backward_alloc_local_y = fftw_mpi_local_size_1d(ny_global, 
		comm_y, 
		FFTW_BACKWARD,
		FFTW_MEASURE,
		&local_ny_in, 
		&local_ny_start_in,
		&local_ny_out, 
		&local_ny_start_out);

	backward_data_y_in = fftw_alloc_complex( backward_alloc_local_y);
	backward_data_y_out = fftw_alloc_complex( backward_alloc_local_y);

	backward_plan_y = fftw_mpi_plan_dft_1d(ny_global, 
		backward_data_y_in, 
		backward_data_y_out, 
		comm_y,
		FFTW_BACKWARD, 
		FFTW_MEASURE);


	// for z direction 
	backward_alloc_local_z = fftw_mpi_local_size_1d( nz_global, 
		comm_z, 
		FFTW_BACKWARD,
		FFTW_MEASURE,
		&local_nz_in, 
		&local_nz_start_in,
		&local_nz_in, 
		&local_nz_start_in);

	backward_data_z_in = fftw_alloc_complex( backward_alloc_local_z);
	backward_data_z_out = fftw_alloc_complex( backward_alloc_local_z);

	backward_plan_z = fftw_mpi_plan_dft_1d(nz_global, 
		backward_data_z_in, 
		backward_data_z_out, 
		comm_z,
		FFTW_BACKWARD, 
		FFTW_MEASURE);
}

void nuc3d::fft::solveFFT_forward(const Field &field_in,
	Field &field_re_out,
	Field &field_im_out,
	MPIComunicator3d_nonblocking &myMPI)
{
	double *pf_in=field_in.getDataPtr();
	double *pf_re_out=field_re_out.getDataPtr();
	double *pf_im_out=field_im_out.getDataPtr();

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;

				forward_data_x_in[i][0]=pf_in[idx_xi] ;
				forward_data_x_in[i][1]=0.0 ;
			}
	
    		myMPI.barrier();
			fftw_execute(forward_plan_x);
	
    		myMPI.barrier();
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;
                                
				temp_re[idx_xi]=forward_data_x_out[i][0];
				temp_im[idx_xi]=forward_data_x_out[i][1];	
			}
		}
	}

	transpose_i2j(temp_re,pf_re_out);
	transpose_i2j(temp_im,pf_im_out);

	for(int i=0;i<nx;i++)
	{
		for(int k=0;k<nz;k++)
		{
			for (int j = 0; j < ny; ++j)
			{
                int idx_eta=ny*nz*i+ny*k+j;

				forward_data_y_in[j][0]=pf_re_out[idx_eta] ;
				forward_data_y_in[j][1]=pf_im_out[idx_eta] ;		
			}

    		myMPI.barrier();
			fftw_execute(forward_plan_y);

    		myMPI.barrier();
			for (int j = 0; j < ny; ++j)
			{
                int idx_eta=ny*nz*i+ny*k+j;
                                
				temp_re[idx_eta]=forward_data_y_out[j][0];
				temp_im[idx_eta]=forward_data_y_out[j][1];		
			}

		}
	}

	transpose_j2k(temp_re,pf_re_out);
	transpose_j2k(temp_im,pf_im_out);

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
		{
			for (int k = 0; k < nz; ++k)
			{
                int idx_zeta=nz*nx*j+nz*i+k;

				forward_data_z_in[k][0]=pf_re_out[idx_zeta] ;
				forward_data_z_in[k][1]=pf_im_out[idx_zeta] ;		
			}

   			myMPI.barrier();
			fftw_execute(forward_plan_z);

    		myMPI.barrier();
			for (int k = 0; k < nz; ++k)
			{
                int idx_zeta=nz*nx*j+nz*i+k;
                                
				temp_re[idx_zeta]=forward_data_z_out[k][0];
				temp_im[idx_zeta]=forward_data_z_out[k][1];		
			}

		}
	}

	transpose_k2i(temp_re,pf_re_out);
	transpose_k2i(temp_im,pf_im_out);

	normalization(pf_re_out,nx_global*ny_global*nz_global);
	normalization(pf_im_out,nx_global*ny_global*nz_global);

}

void nuc3d::fft::solveFFT_backward(const Field &field_re_in,
	const Field &field_im_in,
	Field &field_re_out,
	Field &field_im_out,
	MPIComunicator3d_nonblocking &myMPI)
{

	double *pf_re_in=field_re_in.getDataPtr();
	double *pf_im_in=field_im_in.getDataPtr();
	double *pf_re_out=field_re_out.getDataPtr();
	double *pf_im_out=field_im_out.getDataPtr();

	for(int k=0;k<nz;k++)
	{

		for(int j=0;j<ny;j++)
		{
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;

				backward_data_x_in[i][0]=pf_re_in[idx_xi];
				backward_data_x_in[i][1]=pf_im_in[idx_xi];		
			}

			fftw_execute(backward_plan_x);

			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;
                                
				temp_re[idx_xi]=backward_data_x_out[i][0];
				temp_im[idx_xi]=backward_data_x_out[i][1];		
			}

		}
	}

	transpose_i2j(temp_re,pf_re_out);
	transpose_i2j(temp_im,pf_im_out);

	for(int i=0;i<nx;i++)
	{
		for(int k=0;k<nz;k++)
		{
			for (int j = 0; j < ny; ++j)
			{
                int idx_eta=ny*nz*i+ny*k+j;

				backward_data_y_in[j][0]=pf_re_out[idx_eta] ;
				backward_data_y_in[j][1]=pf_im_out[idx_eta] ;		
			}

			fftw_execute(backward_plan_y);

			for (int j = 0; j < ny; ++j)
			{
                int idx_eta=ny*nz*i+ny*k+j;
                                
				temp_re[idx_eta]=backward_data_y_out[j][0];
				temp_im[idx_eta]=backward_data_y_out[j][1];		
			}

		}
	}

	transpose_j2k(temp_re,pf_re_out);
	transpose_j2k(temp_im,pf_im_out);

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
		{
			for (int k = 0; k < nz; ++k)
			{
                int idx_zeta=nz*nx*j+nz*i+k;

				backward_data_z_in[k][0]=pf_re_out[idx_zeta] ;
				backward_data_z_in[k][1]=pf_im_out[idx_zeta] ;		
			}

			fftw_execute(backward_plan_z);

			for (int k = 0; k < nz; ++k)
			{
                int idx_zeta=nz*nx*j+nz*i+k;
                                
				temp_re[idx_zeta]=backward_data_z_out[k][0];
				temp_im[idx_zeta]=backward_data_z_out[k][1];		
			}

		}
	}

	transpose_k2i(temp_re,pf_re_out);
	transpose_k2i(temp_im,pf_im_out);

}

void nuc3d::fft::solveShell(const double *pf_re_in,
	MPIComunicator3d_nonblocking &myMPI,
	std::string name,
	int istep,
	double time)
{

	int kmax=std::min(std::min(nx_global,ny_global),nz_global)/2;
	int npx=myMPI.getMyDim(2);
	int npy=myMPI.getMyDim(1);
	int npz=myMPI.getMyDim(0);

	int px=myMPI.getMyCoord(2);
	int py=myMPI.getMyCoord(1);
	int pz=myMPI.getMyCoord(0);

	double *power_shell;
	double *power_shell_sum;



	power_shell= new double[(unsigned)(kmax+1)];
	power_shell_sum= new double[(unsigned)(kmax+1)];

	for (int i = 0; i < (kmax); ++i)
	{
		power_shell[i]=0.0;
	}
	
    double kx=(double)nx*(double)px;
    double ky=(double)ny*(double)py;
    double kz=(double)nz*(double)pz;

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;
                double kx0=((kx+(double)i)>=(nx_global/2.0))?kx+(double)i-nx_global:kx+(double)i;
                double ky0=((ky+(double)j)>=(ny_global/2.0))?ky+(double)j-ny_global:ky+(double)j;
                double kz0=((kz+(double)k)>=(nz_global/2.0))?kz+(double)k-nz_global:kz+(double)k;
                double kmag=sqrt(kx0*kx0+ky0*ky0+kz0*kz0);
                if(kmag<1e-6) kmag=1.0;

                if(kmag<kmax)
                {
                	double k0_floor=floor(kmag);
                	int ki=((kmag-k0_floor)-0.5 <= 0 )?k0_floor:k0_floor+1; 
                	power_shell[ki]+=pf_re_in[idx_xi];
                }

			}
		}	
	}	

	myMPI.allReduceSum(power_shell,power_shell_sum,kmax+1);
    
	double k=0.0;

	for (int i = 0; i < (kmax); ++i)
	{
		k+=power_shell_sum[i];    	
	}
	
    if(0==myMPI.getMyId())
    {
	    std::string forename_flow = ("flowDataTemp/power_shell_");
	    std::string step;
	    std::string tailname = (".dat");
	    
	    std::stringstream ss_step,ss_id;
	    ss_step << istep;
	    ss_step >> step;
	    
	    
	    std::string filename_flow = forename_flow + step + "_" + name + tailname;
	    
	    std::ofstream myIOfile;
	    
	    myIOfile.open(filename_flow);

		for (int i = 0; i < (kmax); ++i)
		{    	
    		myIOfile<<std::setprecision(12)<<i<<" "<<power_shell_sum[i]<<" "<<power_shell_sum[i]/k<<"\n";
		}
        myIOfile<<std::endl;
    	myIOfile.close();
    }
    delete [] power_shell;
    delete [] power_shell_sum;
    myMPI.barrier();
}


void nuc3d::fft::solveShell(const Field &f_re_in,
	const Field &f_im_in,
	MPIComunicator3d_nonblocking &myMPI,
	std::string name,
	int istep,
	double time)
{

	double *pf_re_in=f_re_in.getDataPtr();
	double *pf_im_in=f_im_in.getDataPtr();

	int kmax=std::min(std::min(nx_global,ny_global),nz_global)/2;
	int npx=myMPI.getMyDim(2);
	int npy=myMPI.getMyDim(1);
	int npz=myMPI.getMyDim(0);

	int px=myMPI.getMyCoord(2);
	int py=myMPI.getMyCoord(1);
	int pz=myMPI.getMyCoord(0);

	int kx,ky,kz,k0;
	double *power_shell;
	double *power_shell_sum;


	power_shell= new double[(unsigned)(kmax+1)];
	power_shell_sum= new double[(unsigned)(kmax+1)];

	for (int i = 0; i < (kmax); ++i)
	{
		power_shell[i]=0.0;
	}
	
	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;
                kx=nx*px+i;
                ky=ny*py+j;
                kz=nz*pz+k;

                k0=kx*kx+ky*ky+kz*kz;
                if(k0<pow((kmax),2))
                {
                	double k0_floor=floor(sqrt(k0));
                	int ki=((sqrt(k0)-k0_floor)-0.5 <= 0 )?k0_floor:k0_floor+1; 
                	power_shell[ki]+=std::pow(pf_re_in[idx_xi],2)+std::pow(pf_im_in[idx_xi],2);
                }

			}
		}	
	}	

	myMPI.allReduceSum(power_shell,power_shell_sum,kmax+1);

	double k=0.0;

	for (int i = 0; i < (kmax); ++i)
	{
		k+=power_shell_sum[i];    	
	}
    
    if(0==myMPI.getMyId())
    {
	    std::string forename_flow = ("flowDataTemp/power_shell_");
	    std::string step;
	    std::string tailname = (".dat");
	    
	    std::stringstream ss_step,ss_id;
	    ss_step << istep;
	    ss_step >> step;
	    
	    
	    std::string filename_flow = forename_flow + step + "_" + name + tailname;
	    
	    
	    std::ofstream myIOfile;
	    myIOfile.open(filename_flow);

		for (int i = 0; i < (kmax); ++i)
		{    	
    		myIOfile<<std::setprecision(12)<<i<<" "<<power_shell_sum[i]<<" "<<power_shell_sum[i]/k<<"\n";
		}
        myIOfile<<std::endl;
    	myIOfile.close();
    }
    delete [] power_shell;
    delete [] power_shell_sum;
    myMPI.barrier();
}

void nuc3d::fft::transpose_i2j(const double *f_in, double *f_out)
{

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_eta=ny*nz*i+ny*k+j;

                f_out[idx_eta]=f_in[idx_xi];	
			}
		}
	}
}

void nuc3d::fft::transpose_j2k(const double *f_in, double *f_out)
{

	for(int i=0;i<nx;i++)
	{
		for(int k=0;k<nz;k++)
		{
			for (int j = 0; j < ny; ++j)
			{
                int idx_eta=ny*nz*i+ny*k+j;
                int idx_zeta=nz*nx*j+nz*i+k;	
                f_out[idx_zeta]=f_in[idx_eta];				
			}
		}
	}
}

void nuc3d::fft::transpose_k2i(const double *f_in, double *f_out)
{

	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
		{
			for (int k = 0; k < nz; ++k)
			{
                int idx_xi=nx*ny*k+nx*j+i;
                int idx_zeta=nz*nx*j+nz*i+k;	
                f_out[idx_xi]=f_in[idx_zeta];				
			}
		}
	}
}

void nuc3d::fft::normalization(double *f,int fac)
{

	for(int k=0;k<nz;k++)
	{
		for(int j=0;j<ny;j++)
		{
			for (int i = 0; i < nx; ++i)
			{
                int idx_xi=nx*ny*k+nx*j+i;

                f[idx_xi]/=fac;	
			}
		}
	}
}
#endif
