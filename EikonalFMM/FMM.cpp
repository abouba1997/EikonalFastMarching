#include "FMM.h"

FMM::FMM(double _a, double _b, int _N, function2d _f2d, function3d _f3d)
{
	a = _a; b = _b; f2d = _f2d; f3d = _f3d;
	h = (b - a) / (double)_N;
	N = _N + 1;

	x = new double[N];
	y = new double[N];
	z = new double[N];

	space2d = new double[N * N];
	space3d = new double[N * N * N];
	F = new double[N * N];
	F3 = new double[N * N * N];

	int k = 0;
	for (double i = a; i <= b + h; i += h) {
		x[k] = i;
		y[k] = i;
		z[k] = i;
		k++;
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			F[i * N + j] = f2d(x[i], y[j]);
			space2d[i * N + j] = INFINITY;
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				F3[i * N * N + j * N + k] = f3d(x[i], y[j], z[k]);
				space3d[i * N * N + j * N + k] = INFINITY;
			}
		}
	}
}

FMM::~FMM()
{
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] space2d;
	delete[] space3d;
	delete[] F;
	delete[] F3;
}

void FMM::solveFMM(string fname, Matrix2d origins) {
	// Beginning
	// Initialize the heap
	double* far_away = new double[N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			far_away[i * N + j] = 2;
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		vector<int> point = origins[i];
		space2d[(int)point[0] * N + (int)point[1]] = 0;
		far_away[(int)point[0] * N + (int)point[1]] = 0;
	}

	// Get Neighbors of Origins points
	vector<vector<int>> narrowBandIndexOfOrigins = getIndexesNarrowBand(origins);

	// Label Neighbors of Origins points as narrowBand
	for (int i = 0; i < narrowBandIndexOfOrigins.size(); i++) {
		vector<int> row = narrowBandIndexOfOrigins[i];
		far_away[row[0] * N + row[1]] = 1;

		// calculate eikonal for all narrowBand
		calculateEikonal2d(row[0], row[1]);
	}

	// Get all value of narrowBandIndexOfOrigins and save them to the heap

	int iter = 1;

	clock_t tstart = clock();

	while (true) {
		// Get all Indexes of narrowBand
		vector<vector<int>> indexesOnesNarrowBand = getIndexesNarrowBand(far_away);

		// When not finished all element inside narrowband
		if (!indexesOnesNarrowBand.empty()) {
			// Get the minimum element index
			vector<int> minIndex = getIndexesMinNarrowBand(space2d, indexesOnesNarrowBand);
			// Add the minimum to the frozen elements
			far_away[minIndex[0] * N + minIndex[1]] = 0;

			// Get the neighbors of the minvalue
			vector<vector<int>> newNeighborsIndexes;
			for (int j = 0; j < neighbors.size(); j++) {
				vector<int> ind;
				ind.push_back(minIndex[0] + neighbors[j][0]);
				ind.push_back(minIndex[1] + neighbors[j][1]);
				newNeighborsIndexes.push_back(ind);
			}

			// Put 1 where 2 in far_away (label far away as narrowBand)
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				vector<int> current = newNeighborsIndexes[i];
				if (far_away[current[0] * N + current[1]] == 2) {
					far_away[current[0] * N + current[1]] = 1;
				}
			}
			// We got new narrowBand... must get all indexes and recalculate distance field value for new added...

			// Calculate U tentative of that
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				vector<int> currentIndex = newNeighborsIndexes[i];

				calculateEikonal2d(currentIndex[0], currentIndex[1]);
			}
		}

		if (iter == N * N) {
			break;
		}
		iter++;
	}

	delete[] far_away;
	clock_t tfinish = clock();
	std::cout << "Sequential FMM - Time: \tt = " << (tfinish - tstart) / (double)CLOCKS_PER_SEC << endl;
	std::cout << "  FMM - N: \tN = " << N << endl;
	std::cout << "  FMM solving" << endl;
	/*std::cout << "Writing file..." << endl;
	writeFile2d(fname);*/
}

void FMM::writeError2D(function2d realSolution)
{
	double result = 0, sum = 0;
	for (int i = 1; i < N - 1; i++) {
		for (int j = 1; j < N - 1; j++)
		{
			result = abs(space2d[i * N + j] - realSolution(x[i], y[j]));
			if (result > sum) {
				sum = result;
			}
		}
	}
	cout << "N = " << N << "\twith norm = " << sum << std::endl;
}

void FMM::solveFMM3(string filename, Matrix2d origins, int iterations)
{
	// Beginning
	double* far_away = new double[N * N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				far_away[i * N * N + j * N + k] = 2;
			}
		}
	}

	// Get origins fixed data
	for (int i = 0; i < origins.size(); i++) {
		vector<int> point = origins[i];
		space3d[point[0] * N * N + point[1] * N + point[2]] = 0;
		far_away[point[0] * N * N + point[1] * N + point[2]] = 0;
		cout << "Point:\t" << point[0] << "\t" << point[1] << "\t" << point[2] << "\n";
	}

	// Get Neighbors of Origins points
	vector<vector<int>> narrowBandIndexOfOrigins = getIndexesNarrowBand3(origins);

	// Label Neighbors of Origins points as narrowBand
	for (int i = 0; i < narrowBandIndexOfOrigins.size(); i++) {
		vector<int> point = narrowBandIndexOfOrigins[i];
		far_away[point[0] * N * N + point[1] * N + point[2]] = 1;
		//cout << "Point : " << i+1 << " =\t" << point[0] << "\t" << point[1] << "\t" << point[2] << endl;

		// calculate eikonal for all narrowBand
		calculateEikonal3d(point[0], point[1], point[2]);
	}

	int iter = 1;

	while (true) {
		// Get all Indexes of narrowBand
		vector<vector<int>> indexesOnesNarrowBand = getIndexesNarrowBand3(far_away);

		// When not finished all element inside narrowband
		if (!indexesOnesNarrowBand.empty()) {
			// Get the minimum element index
			vector<int> minValueIndex = getIndexesMinNarrowBand3(space3d, indexesOnesNarrowBand);

			// Add the minimum to the frozen elements
			far_away[minValueIndex[0] * N * N + minValueIndex[1] * N + minValueIndex[2]] = 0;

			// Get the neighbors of the minvalue
			vector<vector<int>> newNeighborsIndexes;
			for (int j = 0; j < neighbors3d.size(); j++) {
				vector<int> ind;
				ind.push_back(minValueIndex[0] + neighbors3d[j][0]);
				ind.push_back(minValueIndex[1] + neighbors3d[j][1]);
				ind.push_back(minValueIndex[2] + neighbors3d[j][2]);
				newNeighborsIndexes.push_back(ind);
			}

			// Put 1 where 2 in far_away (label far away as narrowBand)
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				vector<int> current = newNeighborsIndexes[i];
				if (far_away[current[0] * N * N + current[1] * N + current[2]] == 2) {
					far_away[current[0] * N * N + current[1] * N + current[2]] = 1;
				}
			}
			// We got new narrowBand... must get all indexes and recalculate distance field value for new added...

			// Calculate U tentative of that
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				vector<int> currentIndex = newNeighborsIndexes[i];

				calculateEikonal3d(currentIndex[0], currentIndex[1], currentIndex[2]);
			}
		}

		if (iter == iterations) {
			break;
		}
		iter++;
	}
	
	delete[] far_away;

	cout << "Writing file..." << endl;
	writeFile3d(filename);
}

void FMM::writeFile2d(string fname) {
	ofstream file(fname);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			file << x[i] << "\t" << y[j] << "\t" << space2d[i * N + j] << endl;
		}
		file << endl;
	}
	file.close();
}

vector<vector<int>> FMM::getIndexesNarrowBand(double* far) {
	vector<vector<int>> indexes;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (far[i * N + j] == 1) {
				if (i != 0 && i != N - 1 && j != 0 && j != N - 1) {
					vector<int> ind;
					ind.push_back(i);
					ind.push_back(j);
					indexes.push_back(ind);
				}
			}
		}
	}

	return indexes;
}

vector<vector<int>> FMM::getIndexesNarrowBand(vector<vector<int>> origins)
{
	vector<vector<int>> narrowBandIndexOfOrigins;

	for (int i = 0; i < origins.size(); i++) {
		if (origins[i][0] == 0 && origins[i][1] == 0) {
			for (int k = 0; k < top_left_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_left_neighbors[k][0]);
				ind.push_back(origins[i][1] + top_left_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1) {
			for (int k = 0; k < top_right_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_right_neighbors[k][0]);
				ind.push_back(origins[i][1] + top_right_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0) {
			for (int k = 0; k < bottom_left_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_left_neighbors[k][0]);
				ind.push_back(origins[i][1] + bottom_left_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1) {
			for (int k = 0; k < bottom_right_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_right_neighbors[k][0]);
				ind.push_back(origins[i][1] + bottom_right_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0) {
			for (int k = 0; k < top_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary[k][0]);
				ind.push_back(origins[i][1] + top_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1) {
			for (int k = 0; k < bottom_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0) {
			for (int k = 0; k < left_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + left_boundary[k][0]);
				ind.push_back(origins[i][1] + left_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == N - 1) {
			for (int k = 0; k < right_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + right_boundary[k][0]);
				ind.push_back(origins[i][1] + right_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else {
			for (int j = 0; j < neighbors.size(); j++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + neighbors[j][0]);
				ind.push_back(origins[i][1] + neighbors[j][1]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
	}

	return narrowBandIndexOfOrigins;
}

vector<int> FMM::getIndexesMinNarrowBand(double* T, vector<vector<int>> neigh)
{
	vector<int> minIndex = neigh[0];
	double minT = T[minIndex[0] * N + minIndex[1]];

	for (int i = 1; i < neigh.size(); i++) {
		vector<int> currentIndex = neigh[i];

		double currentMin = T[currentIndex[0] * N + currentIndex[1]];

		if (currentMin < minT) {
			minT = currentMin;
			minIndex = currentIndex;
		}
	}

	return minIndex;
}

void FMM::calculateEikonal2d(int i, int j)
{
	// calculate eikonal for all narrowBand
	double minX, minY, calculated;

	if (i == 0) {
		minY = space2d[(i + 1) * N + j];
	}
	else if (i == N - 1) {
		minY = space2d[(i - 1) * N + j];
	}
	else {
		minY = min(space2d[(i + 1) * N + j], space2d[(i - 1) * N + j]);
	}

	if (j == 0) {
		minX = space2d[i * N + j + 1];
	}
	else if (j == N - 1) {
		minX = space2d[i * N + j - 1];
	}
	else {
		minX = min(space2d[i * N + j + 1], space2d[i * N + j - 1]);
	}

	// Result
	if (abs(minX - minY) >= F[i * N + j] * h) {
		calculated = min(minX, minY) + h * F[i * N + j];
	}
	else {
		double sqrtP = 2 * pow(F[i * N + j], 2) * h * h - pow(minX - minY, 2);
		calculated = (minX + minY + sqrt(sqrtP)) / 2.;
	}
	space2d[i * N + j] = min(space2d[i * N + j], calculated);
}

void FMM::writeFile3d(string fname)
{
	ofstream file(fname);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				file << x[i] << "\t" << y[j] << "\t" << z[k] << "\t" << space3d[i * N * N + j * N + k] << endl;
			}
			file << endl;
		}
	}
}

vector<vector<int>> FMM::getIndexesNarrowBand3(double* far)
{
	vector<vector<int>> indexes;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				if (far[i * N * N + j * N + k] == 1) {
					if (i != 0 && i != N - 1 && j != 0 && j != N - 1 && k != 0 && k != N - 1) {
						vector<int> ind;
						ind.push_back(i);
						ind.push_back(j);
						ind.push_back(k);
						indexes.push_back(ind);
					}
					else {
						if (i == 0 && j == 0 && k == 0) {
							vector<int> ind;
							ind.push_back(i);
							ind.push_back(j);
							ind.push_back(k);
							indexes.push_back(ind);
						}
					}
				}
			}
		}
	}

	return indexes;
}

vector<vector<int>> FMM::getIndexesNarrowBand3(vector<vector<int>> origins)
{
	vector<vector<int>> narrowBandIndexOfOrigins;
	for (int i = 0; i < origins.size(); i++) {
		// Les corners
		if (origins[i][0] == 0 && origins[i][1] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_far_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_far_left[k][0]);
				ind.push_back(origins[i][1] + bottom_far_left[k][1]);
				ind.push_back(origins[i][2] + bottom_far_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_near_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_near_left[k][0]);
				ind.push_back(origins[i][1] + bottom_near_left[k][1]);
				ind.push_back(origins[i][2] + bottom_near_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_far_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_far_right[k][0]);
				ind.push_back(origins[i][1] + bottom_far_right[k][1]);
				ind.push_back(origins[i][2] + bottom_far_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_near_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_near_right[k][0]);
				ind.push_back(origins[i][1] + bottom_near_right[k][1]);
				ind.push_back(origins[i][2] + bottom_near_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_far_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_far_left[k][0]);
				ind.push_back(origins[i][1] + top_far_left[k][1]);
				ind.push_back(origins[i][2] + top_far_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_near_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_near_left[k][0]);
				ind.push_back(origins[i][1] + top_near_left[k][1]);
				ind.push_back(origins[i][2] + top_near_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_far_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_far_right[k][0]);
				ind.push_back(origins[i][1] + top_far_right[k][1]);
				ind.push_back(origins[i][2] + top_far_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_near_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_near_right[k][0]);
				ind.push_back(origins[i][1] + top_near_right[k][1]);
				ind.push_back(origins[i][2] + top_near_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		} //top
		else if (origins[i][0] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_far.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_far[k][0]);
				ind.push_back(origins[i][1] + top_boundary_far[k][1]);
				ind.push_back(origins[i][2] + top_boundary_far[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_near.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_near[k][0]);
				ind.push_back(origins[i][1] + top_boundary_near[k][1]);
				ind.push_back(origins[i][2] + top_boundary_near[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_left[k][0]);
				ind.push_back(origins[i][1] + top_boundary_left[k][1]);
				ind.push_back(origins[i][2] + top_boundary_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_right[k][0]);
				ind.push_back(origins[i][1] + top_boundary_right[k][1]);
				ind.push_back(origins[i][2] + top_boundary_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}// Bottom
		else if (origins[i][0] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_far.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_far[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_far[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_far[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_near.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_near[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_near[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_near[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_left[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_left[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_right[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_right[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		} // Middle
		else if (origins[i][0] == 0 && origins[i][1] == 0) {
			for (int k = 0; k < middle_boundary_far_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_far_left[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_far_left[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_far_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1) {
			for (int k = 0; k < middle_boundary_far_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_far_right[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_far_right[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_far_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0) {
			for (int k = 0; k < middle_boundary_near_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_near_left[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_near_left[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_near_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1) {
			for (int k = 0; k < middle_boundary_near_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_near_right[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_near_right[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_near_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}//faces
		else if (origins[i][0] == 0 || origins[i][0] == N - 1) {
			for (int k = 0; k < face_boundary_top_bottom.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + face_boundary_top_bottom[k][0]);
				ind.push_back(origins[i][1] + face_boundary_top_bottom[k][1]);
				ind.push_back(origins[i][2] + face_boundary_top_bottom[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0 || origins[i][1] == N - 1) {
			for (int k = 0; k < face_boundary_left_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + face_boundary_left_right[k][0]);
				ind.push_back(origins[i][1] + face_boundary_left_right[k][1]);
				ind.push_back(origins[i][2] + face_boundary_left_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][2] == 0 || origins[i][2] == N - 1) {
			for (int k = 0; k < face_boundary_far_near.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + face_boundary_far_near[k][0]);
				ind.push_back(origins[i][1] + face_boundary_far_near[k][1]);
				ind.push_back(origins[i][2] + face_boundary_far_near[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		// other
		else {
			for (int k = 0; k < neighbors3d.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + neighbors3d[k][0]);
				ind.push_back(origins[i][1] + neighbors3d[k][1]);
				ind.push_back(origins[i][2] + neighbors3d[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
	}

	return narrowBandIndexOfOrigins;
}

vector<int> FMM::getIndexesMinNarrowBand3(double* T, vector<vector<int>> neigh)
{
	vector<int> minIndex = neigh[0];
	double minT = T[minIndex[0] * N * N + minIndex[1] * N + minIndex[2]];

	for (int i = 1; i < neigh.size(); i++) {
		vector<int> currentIndex = neigh[i];

		double currentMin = T[currentIndex[0] * N * N + currentIndex[1] * N + currentIndex[2]];

		if (currentMin < minT) {
			minT = currentMin;
			minIndex = currentIndex;
		}
	}

	return minIndex;
}

void FMM::calculateEikonal3d(int i, int j, int k)
{
	double txmin, tymin, tzmin, t;
	double b, c;

	if (!isnan(space3d[i * N * N + j * N + k])) {
		// by i
		if (i == 0) {
			tymin = space3d[(i + 1) * N * N + j * N + k];
		}
		else if (i == N - 1) {
			tymin = space3d[(i - 1) * N * N + j * N + k];
		}
		else {
			tymin = min(space3d[(i - 1) * N * N + j * N + k], space3d[(i + 1) * N * N + j * N + k]);
		}

		// by j
		if (j == 0) {
			txmin = space3d[i * N * N + (j + 1) * N + k];
		}
		else if (j == N - 1) {
			txmin = space3d[i * N * N + (j - 1) * N + k];
		}
		else {
			txmin = min(space3d[i * N * N + (j - 1) * N + k], space3d[i * N * N + (j + 1) * N + k]);
		}

		// by k
		if (k == 0) {
			tzmin = space3d[i * N * N + j * N + k + 1];
		}
		else if (k == N - 1) {
			tzmin = space3d[i * N * N + j * N + k - 1];
		}
		else {
			tzmin = min(space3d[i * N * N + j * N + k - 1], space3d[i * N * N + j * N + k + 1]);
		}

		//result
		double l = txmin * txmin + tymin * tymin + tzmin * tzmin;
		b = pow((txmin + tymin + tzmin), 2);
		c = 3 * (l - h * h * f3d(x[i], y[j], z[k]));

		if ((b - c) > 0) {
			t = ((txmin + tymin + tzmin) + sqrt(b - c)) / 3.;
		}
		else {
			t = min3(txmin, tymin, tzmin) + f3d(x[i], y[j], z[k]) * h;
		}

		// set space data
		space3d[i * N * N + j * N + k] = min(space3d[i * N * N + j * N + k], t);
	}
}
