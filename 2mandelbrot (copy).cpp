/**
 * Mandelbrot project
 * Stefano Ottolenghi
 */
 
#include <complex>
#include <iostream>
using namespace std;

int max_iterations = 20;

int mandelbrot_print(complex<double> c) {
	int n = 0;
	complex<double> z = c;
	cout << z << " ";
	while(abs(z) < 2 && n < max_iterations) {
		z = pow(z, 2) + c;
		cout << z << " ";
		n++;
	}

	cout << n ;
	return (n == max_iterations) ? -1 : n;
}

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
	return (n == max_iterations) ? -1 : n;
}

void print_matrix(int* matrix, int rows, int cols) {
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			if(matrix[rows*i + j] == -1)
				cout << "#" << " ";
			else
				cout << matrix[rows*i + j] << " ";
		}
		cout << endl;
	}
}

int main() {
	complex<double> UR =  0.5 +1i;
	complex<double> UL = -1 +1i;
	complex<double> LL = -1 -1i;
	complex<double> LR =  0.5 -1i;

	double step = 0.1;

	//Memory init
	int rows = abs(UL-UR)/step, cols = (int)abs(UL-LR)/step;
	int* matrix = (int*) malloc(rows * cols * sizeof(int));
	//for(int i = 0; i < (int)abs(UL-UR)/step; i++)
	//	output[i] = malloc((int)abs(UL-LR)/step * sizeof(int));
	
	//double x_len = abs(UR-UL), y_len = abs(LR-LL);

	complex<double> c = LL - complex<double> { step, step };
	while(real(c) + step <= real(LR)) {
		c = c+step;
		
		while(imag(c) + step <= imag(UL)) {
			c = { real(c), imag(c) + step };
			int man_n = mandelbrot(c);
			complex<int> aux = {(int)((real(c) - real(LL))*(1/step)), (int)((imag(c) - imag(LL))*(1/step))};
			complex<int> aux2 = {rows - (real(aux)), imag(aux)};
			matrix[(real(aux2)) * cols + imag(aux2)] = man_n;
			cout << c << " " << aux << " " << aux2 << " " << man_n << "; " << endl;
		}
		cout << endl << endl << endl;

		c = { real(c), - imag(c) - step };
	}

	print_matrix(matrix, rows, cols);
	cout << rows << " " << cols << endl;

	cout << mandelbrot_print(complex<double> {0.3,0});
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