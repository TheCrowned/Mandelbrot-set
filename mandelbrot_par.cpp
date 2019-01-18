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

int rows_per_task = 10;

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
void scrivi_riga_output(double* riga, int buff_size) {
	int x_size = buff_size / rows_per_task;
	complex<double> z = { (double) riga[x_size*rows_per_task], (double) riga[x_size*rows_per_task + 1] };

	for(int row = 0; row < rows_per_task; row++) {
		for(int col = 0; col < x_size; col++) {
			output << real(z) << " " << imag(z) << " " << riga[row*x_size + col] << "\n";
			z = { real(z) + step, imag(z) };
		}
		z = { real(z), imag(z) - step };
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

	int x_size = (abs(LL-LR)/step)+1, y_size = ((abs(LL-UL)/step)+1);
	chrono::high_resolution_clock::time_point tstart, tend;
	chrono::duration<double> diff;

	//cout << x_size << endl << endl;

	//Master thread
	if(rank == 0) {
		tstart = chrono::high_resolution_clock::now();

		cout << "Running Mandelbrot routine on " << x_size*y_size << " points with " << nproc-1 << " workers. Each task consists of " << rows_per_task << " rows, each of " << x_size << " columns." << endl << endl;
	
		MPI_Status mpi_status;
		MPI_Request mpi_request;
		double* recv_buff = (double*) malloc((x_size*rows_per_task + 2) * sizeof(double));
	
		output.open("output.dat");

		//Span plane region and call mandlebrot
		complex<double> c_send = UL, c_recv = UL;

		for(int i = 1; i < nproc*2; i++) { //exclude master

			if(i%nproc == 0) continue;
			
			//cout << "Sending point " << c_send << " to node " << i%nproc << endl;
			MPI_Send(&c_send, 1, MPI_DOUBLE_COMPLEX, i%nproc, 0,  MPI_COMM_WORLD);
			c_send = { real(UL), imag(c_send) - step*rows_per_task };

			if(imag(c_send) < imag(LL)) break;
		}

		//Receive from workers and, for each result, send out a new task
		while(imag(c_recv) >= imag(LL)) {
			MPI_Recv(recv_buff, x_size*rows_per_task + 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status); //cannot use non-blocking here, otherwise can't use mpi_status.MPI_SOURCE before scrivi_riga_output
			//cout << "Received point " << recv_buff[x_size*rows_per_task] << "," << recv_buff[x_size*rows_per_task + 1] << " from node " << mpi_status.MPI_SOURCE << endl;
			c_recv = { real(UL), imag(c_recv) - step*rows_per_task };

			cout << ((imag(c_send) - imag(LL))/step) << endl;
			if(rows_per_task > ((imag(c_send) - imag(LL))/step) ) {
				rows_per_task = ((imag(c_send) - imag(LL))/step)+1;
				cout << ((imag(c_send) - imag(LL))/step)+1;
			}

			cout << "Sending point " << c_send << " to node " << mpi_status.MPI_SOURCE << endl;
			MPI_Send(&c_send, 1, MPI_DOUBLE_COMPLEX, mpi_status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			c_send = { real(UL), imag(c_send) - step*rows_per_task };
			
			scrivi_riga_output(recv_buff, x_size*rows_per_task);
		}

		output.close();

		tend = chrono::high_resolution_clock::now();
		diff = tend - tstart;

	//Worker threads
	} else {
		MPI_Status mpi_status;
		MPI_Request mpi_request;
		complex<double> c;
		
		double* man_results = (double*) malloc((x_size*rows_per_task + 2) * sizeof(double)); //need to be double to store c coord

		int count = 0;
		for(;;) {
			MPI_Recv(&c, 1, MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &mpi_status);
			cout << "* " << "[" << rank << "]" << " Received point " << c << ", processing..." << endl;

			if(imag(c) < imag(LL)) {
				//cout << "* " << "[" << rank << "]" << " STOPPING thread " << rank << endl;
				cout << "* Worker " << rank << " processed " << count << " tasks" << endl;
				break;
			}

			//Ensure non-blocking send of result is finished before overwriting its buffer
			if(count != 0 )
				MPI_Wait(&mpi_request, &mpi_status);

			//Store starting point
			man_results[x_size*rows_per_task] = real(c);
			man_results[x_size*rows_per_task + 1] = imag(c);

			int i = 0;
			while(imag(c) >= imag(LL)) {
				while(real(c) <= real(UR)) {
					man_results[i] = mandelbrot(c);
					//cout << "  " << c << "," << man_results[i] << "; ";
					c = { real(c) + step, imag(c) };
					i++;
				}
				c = { real(c), imag(c) - step };
			}

			//cout << "* " << "[" << rank << "]" << " Sending back results from point " << man_results[x_size] << "," << man_results[x_size+1] << endl;
			MPI_Isend(man_results, x_size*rows_per_task + 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpi_request);

			count++;
		}
	}

	MPI_Finalize();

	if(rank == 0)
		cout << "\n~~ Completed in " << diff.count() << " seconds ~~\n";
}