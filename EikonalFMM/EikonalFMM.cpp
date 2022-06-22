#include <iostream>
#include "FMM.h"

#define PI 3.14159265359

double func2d(double x, double y) {
	//return 2 * exp(x * x + y * y) * sqrt(x * x + y * y); // Second example
	return 1;
	
	// random
	//return ((double)rand() / RAND_MAX);

	// From article fast methods
	/*if (x > -0.8 && x < -0.75 && y <= 0.9) {
		return INFINITY;
	}
	return 1.;*/

	/*if ((x > -0.95 && x < -0.9 && y <= 0.9) || (x > -0.75 && x < -0.7 && y > -0.9) || (x > -0.55 && x < -0.5 && y <= 0.9) || (x > -0.35 && x < -0.3 && y > -0.9) || (x > -0.15 && x < -0.1 && y <= 0.9)) {
		return INFINITY;
	}
	else {
		return 1.;
	}*/

	//return sqrt(pow(sin(PI + x * PI / 2), 2) + pow(sin(PI + y * PI / 2), 2)); // First example
	
	//return 3 * sqrt((x * x + y * y) / pow(9 + x * x + y * y, 3)); // Third example
}

double func3d(double x, double y, double z) {
	//return 1;
	return sqrt(pow(sin(x), 2) + pow(sin(y), 2) + pow(sin(z), 2));
}

double realSolution2d(double x, double y) {
	return sqrt(x * x + y * y);
}


int main()
{
	int N = 100;
	FMM fmm(-1, 1, N, func2d, func3d);

	// 2D
	Matrix2d origins2d = { {N/2,N/2} };
	fmm.solveFMM("2d_fmm.txt", origins2d);
	fmm.writeError2D(realSolution2d);

	// 3D
	/*Matrix2d origins = { {20,20,20} };
	fmm.solveFMM3("3d_testt.txt", origins, N * N * N);*/

	cout << "End of the program!!!" << endl;
}
