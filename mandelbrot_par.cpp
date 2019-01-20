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

int rows_per_task = 2,
	x_size = 0, y_size = 0; //really init in main

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
 * Writes a "processed line" in output file.
 * Note that a processed line could be made of several rows, depending on
 * the value of riga[2].
 *
 * @param	double* riga Processed line
 * @return	void
 */ 
void scrivi_riga_output(double* riga) {
	int real_rows_per_task = riga[2] / x_size;
	complex<double> z = { (double) riga[0], (double) riga[1] };

	for(int row = 0; row < real_rows_per_task; row++) {
		for(int col = 0; col < x_size; col++) {
			output << real(z) << " " << imag(z) << " " << (int) riga[3 + row*x_size + col] << "\n";
			z = { real(z) + step, imag(z) };
		}
		z = { real(LL), imag(z) - step };
	}
}

int main(int argc, char** argv) {

	int nproc, rank;

	//Params handling
	if(argc < 7) {
		cout << "Usage: mandelbrot <max_iterations> <step_size> <lower_left_real> <lower_left_imaginary> <upper_right_real> <upper_right_imaginary> [<rows_per_task>]";
		return 0;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	max_iterations = atoi(argv[1]);
	step = atof(argv[2]);
	LL = { atof(argv[3]), atof(argv[4]) };
	UR = { atof(argv[5]), atof(argv[6]) };

	LR = { real(UR), imag(LL) };
	UL = { real(LL), imag(UR) };

	if(argc == 8)
		rows_per_task = atoi(argv[7]);

	x_size = (abs(LL-LR)/step)+1, y_size = ((abs(LL-UL)/step)+1);
	chrono::high_resolution_clock::time_point tstart, tend;
	chrono::duration<double> diff;

	//cout << x_size << y_size << endl;

	//Master thread
	if(rank == 0) {
		tstart = chrono::high_resolution_clock::now();

		cout << "Running Mandelbrot routine on " << x_size*y_size << " points with " << nproc-1 << " workers. Each task consists of " << rows_per_task << " rows, each of " << x_size << " columns." << endl << endl;
	
		MPI_Status mpi_status;
		MPI_Request mpi_request;
		double* recv_buff = (double*) malloc((x_size*rows_per_task +2 +1) * sizeof(double));
	
		output.open("output.dat");

		//Span plane region and call mandlebrot
		complex<double> c_send = UL, c_recv = UL;

		for(int i = 0; i < nproc*2; i++) {

			if(i % nproc == 0) continue; //exclude master
			
			//cout << "Sending point " << c_send << " to node " << i%nproc << endl;
			MPI_Send(&c_send, 1, MPI_DOUBLE_COMPLEX, i%nproc, 0,  MPI_COMM_WORLD);
			c_send = { real(UL), imag(c_send) - step*rows_per_task };

			if(imag(c_send) < imag(LL)) break;
		}

		//Receive from workers and, for each result, send out a new task
		while(imag(c_recv) >= imag(LL)) {
			MPI_Recv(recv_buff, x_size*rows_per_task +2 +1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status); //cannot use non-blocking here, otherwise can't use mpi_status.MPI_SOURCE before scrivi_riga_output
			cout << "Received point " << recv_buff[0] << "," << recv_buff[1] << " from node " << mpi_status.MPI_SOURCE << ", buffer real size " << recv_buff[2] << endl;
			c_recv = { real(UL), imag(c_recv) - step*rows_per_task };

			//cout << "Sending point " << c_send << " to node " << mpi_status.MPI_SOURCE << endl;
			MPI_Send(&c_send, 1, MPI_DOUBLE_COMPLEX, mpi_status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			c_send = { real(UL), imag(c_send) - step*rows_per_task };
			
			scrivi_riga_output(recv_buff);
		}

		output.close();

		tend = chrono::high_resolution_clock::now();
		diff = tend - tstart;

	//Worker threads
	} else {
		MPI_Status mpi_status;
		MPI_Request mpi_request;
		complex<double> c;
		
		double* man_results = (double*) malloc((x_size*rows_per_task +2 +1) * sizeof(double)); //need to be double to store c coord

		int count = 0;
		for(;;) {
			MPI_Recv(&c, 1, MPI_DOUBLE_COMPLEX, 0, 0, MPI_COMM_WORLD, &mpi_status);
			//cout << "* " << "[" << rank << "]" << " Received point " << c << ", processing..." << endl;

			if(imag(c) <= imag(LL)-step) { //or it would skip point -1 -1
				/*cout << imag(c)*10000 << imag(LL)*10000 << endl;
				cout << (imag(c) == imag(LL)) << endl;
				cout << (imag(c) < imag(LL)) << endl;
				cout << (imag(c) <= imag(LL)) << endl;*/
				cout << "* " << "[" << rank << "]" << " STOPPING thread " << rank << endl;
				cout << "* Worker " << rank << " processed " << count << " tasks" << endl;
				break;
			}

			//Ensure non-blocking send of result is finished before overwriting its buffer
			if(count != 0 )
				MPI_Wait(&mpi_request, &mpi_status);

			//Store starting point
			man_results[0] = real(c);
			man_results[1] = imag(c);

			int i = 0, j = 0;
			while(imag(c) > imag(LL)-step && i*x_size < rows_per_task*x_size) {
				for(j = 0; j < x_size; j++) {
					man_results[3 +i*x_size +j] = mandelbrot(c);
					//cout << "*  " << c << "," << man_results[3 +i*x_size +j] << "; ";
					
					c = { real(c) + step, imag(c) };
				}
				i++;
				c = { real(LL), imag(c) - step }; //reset x coord, step one y coord
			}

			man_results[2] = (i-1)*x_size +j; //Store buffer real size in buffer (final tasks can be smaller than others)

			//cout << "* " << "[" << rank << "]" << " Sending back results from point " << man_results[0] << "," << man_results[1] << ", size " << man_results[2] << endl;
			MPI_Isend(man_results, x_size*rows_per_task +2 +1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpi_request);

			count++;
		}
	}

	MPI_Finalize();

	if(rank == 0)
		cout << "\n~~ Completed in " << diff.count() << " seconds ~~\n";
}