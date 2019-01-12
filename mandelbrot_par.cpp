/**
 * Mandelbrot project
 * Stefano Ottolenghi
 */
 
#include <complex>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <ctime>
#include <chrono>

using namespace std;

int max_iterations = 20;
double step = 0.1;

complex<double> UR =  0.5 +1i;
complex<double> UL = -1 +1i;
complex<double> LL = -1 -1i;
complex<double> LR =  0.5 -1i;

ofstream output;

/*
 * Test convergence for max_iterations of a single complex input point.
 *
 * @param	complex<double> c Input point
 * @return	int Iterations to divergence
 */ 
int mandelbrot(complex<double> c) {
	int n = 1;
	complex<double> z = c;
	
	while(abs(z) < 2 && n < max_iterations) {
		z = pow(z, 2) + c;
		n++;
	}

	return n;
}

/**
 *
 */ 
void scrivi_riga_output(double* riga, int x_size) {
	complex<double> z = { (double) riga[x_size], (double) riga[x_size+1] };

	for(int i = 0; i < x_size; i++) {
		output << real(z) << " " << imag(z) << " " << riga[i] << "\n";
		z = { real(z) + step, imag(z) };
	}
}

int main(int argc, char** argv) {

	int nproc, rank;

	//Params handling
	if(argc < 7) {
		cout << "Usage: mandelbrot <max_iterations> <step_size> <lower_left_real> <lower_left_imaginary> <upper_right_real> <upper_right_imaginary>";
		return 0;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//cout << "MPI Init, size " << nproc << ", rank " << rank << endl << endl;

	max_iterations = atoi(argv[1]);
	step = atof(argv[2]);
	LL = { atof(argv[3]), atof(argv[4]) };
	UR = { atof(argv[5]), atof(argv[6]) };

	LR = { real(UR), imag(LL) };
	UL = { real(LL), imag(UR) };

	int x_size = (abs(LL-LR)/step)+1;
	chrono::high_resolution_clock::time_point tstart, tend;
	chrono::duration<double> diff;

	//cout << x_size << endl << endl;

	//Master thread
	if(rank == 0) {
		tstart = chrono::high_resolution_clock::now();
	
		MPI_Status mpi_status;
		double* recv_buff = (double*) malloc((x_size+2) * sizeof(double));
	
		output.open("output.dat");

		//Span plane region and call mandlebrot
		complex<double> c_send = UL, c_recv = UL;

		for(int i = 1; i < nproc; i++) { //exclude master
			//cout << "Sending point " << c_send << " to node " << i << endl;
			MPI_Send(&c_send, 1, MPI_DOUBLE_COMPLEX, i, 0,  MPI_COMM_WORLD);
			c_send = { real(UL), imag(c_send) - step };

			if(imag(c_send) < imag(LL)) break;
		}

		while(imag(c_recv) >= imag(LL)) {
			MPI_Recv(recv_buff, x_size+2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
			//cout << "Received point " << recv_buff[x_size] << "," << recv_buff[x_size+1] << " from node " << mpi_status.MPI_SOURCE << endl;
			c_recv = { real(UL), imag(c_recv) - step };

			//cout << "Sending point " << c_send << " to node " << mpi_status.MPI_SOURCE << endl;
			MPI_Send(&c_send, 1, MPI_DOUBLE_COMPLEX, mpi_status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			c_send = { real(UL), imag(c_send) - step };
			
			scrivi_riga_output(recv_buff, x_size);
		}

		output.close();

		tend = chrono::high_resolution_clock::now();
		diff = tend - tstart;

	//Worker threads
	} else {
		MPI_Status mpi_status;
		complex<double> c;
		double* man_results = (double*) malloc((x_size+2) * sizeof(double)); //need to be double to store c coord
		
		for(;;) {
			MPI_Recv(&c, 1, MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &mpi_status);
			//cout << "* " << "[" << rank << "]" << " Received point " << c << ", processing..." << endl;

			if(imag(c) < imag(LL)) {
				//cout << "* " << "[" << rank << "]" << " STOPPING thread " << rank << endl;
				break;
			}

			//Store starting point
			man_results[x_size] = real(c);
			man_results[x_size+1] = imag(c);

			int i = 0;
			while(real(c) <= real(UR)) {
				man_results[i] = mandelbrot(c);
				//cout << "  " << c << "," << man_results[i] << "; ";
				c = { real(c) + step, imag(c) };
				i++;
			}

			//cout << "* " << "[" << rank << "]" << " Sending back results from point " << man_results[x_size] << "," << man_results[x_size+1] << endl;
			MPI_Send(man_results, x_size+2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();

	if(rank == 0)
		cout << "\n~~ Completed in " << diff.count() << " seconds ~~\n";
}