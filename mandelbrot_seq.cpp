/**
 * Mandelbrot project
 * Stefano Ottolenghi
 */

#include <complex>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>

using namespace std;

#define TRIALS 2
#define DEBUG 0

int max_iterations = 20;
double step = 0.1;

complex<double> UR = { 0.5, +1 };
complex<double> UL = { -1, +1 };
complex<double> LL = { -1, -1 };
complex<double> LR = { 0.5, -1 };

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

int main(int argc, char** argv) {

	chrono::high_resolution_clock::time_point ttstart, tstart, tend;
    chrono::duration<double> diff;
    double min = 1.0e100;

	//Params handling
	if(argc < 7) {
		cout << "Usage: mandelbrot <lower_left_real> <lower_left_imaginary> <upper_right_real> <upper_right_imaginary> <max_iterations> <step_size>" << endl;
		return 0;
	}

	max_iterations = atoi(argv[5]);
	step = atof(argv[6]);
	LL = { atof(argv[1]), atof(argv[2]) };
	UR = { atof(argv[3]), atof(argv[4]) };

	LR = { real(UR), imag(LL) };
	UL = { real(LL), imag(UR) };

	ttstart = chrono::high_resolution_clock::now();
	for(int p = 0; p < TRIALS; p++) {
		tstart = chrono::high_resolution_clock::now();

		//Output
		ofstream output;
		output.open("output.dat");

		//Span plane region and call mandlebrot
		complex<double> c = UL + complex<double> { 0, 0 };
		while(imag(c) >= imag(LL)) {

			while(real(c) <= real(UR)) {
				int man_n = mandelbrot(c);

				if(DEBUG)
					cout << c << " " << man_n << endl;

				output << real(c) << " " << imag(c) << " " << man_n << "\n";
				c = { real(c) + step, imag(c) };
			}
			//cout << endl;

			c = { real(UL), imag(c) - step };
		}

		tend = chrono::high_resolution_clock::now();
		diff = tend - tstart;

		/* take the best performance result */
		if (diff.count() < min) min = diff.count();

		output.close();
	}

	cout << "\n~~ Completed in " << min << " seconds. ~~";
}