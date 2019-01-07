/**
 * Mandelbrot project
 * Stefano Ottolenghi
 */
 
#include <complex>
#include <iostream>
#include <fstream>
using namespace std;

int max_iterations = 20;
double step = 0.1;

complex<double> UR =  0.5 +1i;
complex<double> UL = -1 +1i;
complex<double> LL = -1 -1i;
complex<double> LR =  0.5 -1i;

/*
 *
 */ 
int mandelbrot(complex<double> c) {
	int n = 0;
	complex<double> z = c;
	//cout << z << " ";
	while(abs(z) < 2 && n < max_iterations) {
		z = pow(z, 2) + c;
		//cout << z << " ";
		n++;
	}

	//cout << n ;
	return n;//(n == max_iterations) ? -1 : n;
}

/*void print_matrix(int* matrix, int rows, int cols) {
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			if(matrix[rows*i + j] == -1)
				cout << "#" << " ";
			else
				cout << matrix[rows*i + j] << " ";
		}
		cout << endl;
	}
}*/

int main(int argc, char** argv) {

	//Params handling
	if(argc < 7) {
		cout << "Usage: mandelbrot <max_iterations> <step_size> <lower_left_real> <lower_left_imaginary> <upper_right_real> <upper_right_imaginary>";
		return 0;
	}

	max_iterations = atoi(argv[1]);
	step = atof(argv[2]);
	LL = { atof(argv[3]), atof(argv[4]) };
	UR = { atof(argv[5]), atof(argv[6]) };

	LR = { real(UR), imag(LL) };
	UL = { real(LL), imag(UR) };

	//Memory init
	int rows = abs(UL-UR)/step + 1, cols = (int)abs(UL-LR)/step + 1;
	int* matrix = (int*) malloc(rows * cols * sizeof(int));

	ofstream output;
	output.open("output.dat");

	complex<double> c = UL + complex<double> { 0, 0 };
	while(imag(c) >= imag(LL)) {

		while(real(c) <= real(UR)) {
			int man_n = mandelbrot(c);
			//cout << c << " " << man_n << endl;
			output << real(c) << " " << imag(c) << " " << man_n << "\n";
			c = { real(c) + step, imag(c) };
		}
		//cout << endl;

		c = { real(UL), imag(c) - step };
	}

	output.close();
}

/*
int max_row = 25, max_column = 85, max_iteration = 30;

int main() {
	for(auto row = 0; row < max_row; ++row) {
		for(auto column = 0; column < max_column; ++column) {
			std::complex<float> z, c = {
				(float)column * 2 / max_column - 1.5f,
				(float)row * 2 / max_row - 1				
			};
			
			int iteration = 0;
			while(abs(z) < 2 && ++iteration < max_iteration)
				z = pow(z, 2) + c;
			std::cout << (iteration == max_iteration ? '#' : '.');
		}
		std::cout << '\n';
	}
}
*/