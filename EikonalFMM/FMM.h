#pragma once

#include "Heap.h"

#define Matrix2d vector<vector<int>>
#define function2d function<double(double, double)>
#define function3d function<double(double, double, double)>

#define min2(a, b) ((a) < (b) ? (a) : (b))
#define min3(a, b, c) (min2(min2((a), (b)), (c)))

class FMM
{
	public:
		FMM(double, double, int, function2d, function3d);
		~FMM();
		void solveFMM(string, Matrix2d);
		void writeError2D(function2d);
		void solveFMM3(string, Matrix2d, int);

	private:
		double* space2d, *space3d;
		double a, b, h;
		int N;
		double *x, *y, *z;
		double* F, *F3;
		function2d f2d;
		function3d f3d;

		// 2D
		vector<vector<int>> neighbors = { {1,0},{-1,0},{0,1},{0,-1} };

		vector<vector<int>> top_left_neighbors = { {1,0},{0,1} };
		vector<vector<int>> top_right_neighbors = { {1,0},{0,-1} };
		vector<vector<int>> bottom_left_neighbors = { {-1,0},{0,1} };
		vector<vector<int>> bottom_right_neighbors = { {-1,0},{0,-1} };

		vector<vector<int>> top_boundary = { {-1,0},{1,0},{0,1} };
		vector<vector<int>> bottom_boundary = { {-1,0},{1,0},{0,-1} };
		vector<vector<int>> left_boundary = { {1,0},{0,1},{0,-1} };
		vector<vector<int>> right_boundary = { {-1,0},{0,1},{0,-1} };

		void writeFile2d(string);
		vector<vector<int>> getIndexesNarrowBand(double*);
		vector<vector<int>> getIndexesNarrowBand(vector<vector<int>>);
		vector<int> getIndexesMinNarrowBand(double*, vector<vector<int>>);

		void calculateEikonal2d(int, int);

		// 3D
		Matrix2d neighbors3d = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1} };

		Matrix2d top_far_left = { {1,0,0},{0,1,0},{0,0,-1} };
		Matrix2d top_far_right = { {1,0,0},{0,-1,0},{0,0,-1} };
		Matrix2d top_near_left = { {-1,0,0},{0,1,0},{0,0,-1} };
		Matrix2d top_near_right = { {-1,0,0},{0,-1,0},{0,0,-1} };

		Matrix2d bottom_far_left = { {1,0,0},{0,1,0},{0,0,1} };
		Matrix2d bottom_far_right = { {1,0,0},{0,-1,0},{0,0,1} };
		Matrix2d bottom_near_left = { {-1,0,0},{0,1,0},{0,0,1} };
		Matrix2d bottom_near_right = { {-1,0,0},{0,-1,0},{0,0,1} };

		Matrix2d top_boundary_far = { {1,0,0},{0,1,0},{0,0,-1},{0,-1,0} };
		Matrix2d top_boundary_near = { {-1,0,0},{0,1,0},{0,0,-1},{0,-1,0} };
		Matrix2d top_boundary_left = { {1,0,0},{0,1,0},{-1,0,0},{0,0,-1} };
		Matrix2d top_boundary_right = { {1,0,0},{0,-1,0},{-1,0,0},{0,0,-1} };

		Matrix2d bottom_boundary_far = { {1,0,0},{0,1,0},{0,0,1},{0,-1,0} };
		Matrix2d bottom_boundary_near = { {-1,0,0},{0,1,0},{0,0,1},{0,-1,0} };
		Matrix2d bottom_boundary_left = { {1,0,0},{-1,0,0},{0,0,1},{0,-1,0} };
		Matrix2d bottom_boundary_right = { {1,0,0},{-1,0,0},{0,0,1},{0,1,0} };

		Matrix2d middle_boundary_near_left = { {-1,0,0},{0,1,0},{0,0,1},{0,0,-1} };
		Matrix2d middle_boundary_near_right = { {-1,0,0},{0,-1,0},{0,0,1},{0,0,-1} };
		Matrix2d middle_boundary_far_left = { {1,0,0},{0,1,0},{0,0,1},{0,0,-1} };
		Matrix2d middle_boundary_far_right = { {1,0,0},{0,-1,0},{0,0,1},{0,0,-1} };

		Matrix2d face_boundary_top_bottom = { {0,1,0},{0,-1,0},{0,0,1},{0,0,-1} };
		Matrix2d face_boundary_left_right = { {1,0,0},{-1,0,0},{0,0,1},{0,0,-1} };	
		Matrix2d face_boundary_far_near = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0} };

		void writeFile3d(string);
		vector<vector<int>> getIndexesNarrowBand3(double*);
		vector<vector<int>> getIndexesNarrowBand3(vector<vector<int>>);
		vector<int> getIndexesMinNarrowBand3(double*, vector<vector<int>>);

		void calculateEikonal3d(int, int, int);
};

