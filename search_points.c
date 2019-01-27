#include "mex.h"

typedef struct
{
	double **img_round;
	double *m;
	double *m_p;
	double *n;
	double *n_p;
}laserline;

typedef struct
{
	int *row_mask;
	int *col_mask;
	int row_mask_end[2];
	int col_mask_end[2];
	int *row_max_idx;
	int *col_max_idx;
	double *row_max_val;
	double *col_max_val;
	int global_max_idx[2];
	double global_max_val;
}img_para;

void** fspace_2d(int row, int col, int length)
{
	int i;
	void **b;
	b = (void **)calloc(sizeof(void *), row);
	if (b == NULL) {
		mexPrintf("error1\n");
		exit(1);
	}
	for (i = 0; i < row; i++) {
		b[i] = (void *)calloc(length, col);
		if (b[i] == NULL) {
			mexPrintf("error2\n");
			exit(1);
		}
	}
	return(b);
}

void ffree_2d(void **a, int row)
{
	int i;
	for (i = 0; i < row; i++) {
		free(a[i]);
		a[i] = NULL;
	}
	free(a);
	a = NULL;
}

int round_nngt(double input)
{
	if (input - (int)input < 0.5)
	{
		return((int)input);
	}
	else
	{
		return((int)(input + 1));
	}
}

void process(laserline *test, double **img, int M, int N)
{
	int *row_mask = (int *)calloc(sizeof(int), M);
	int row_mask_end[] = {M-1, 0};
	int *row_max_idx = (int *)calloc(sizeof(int), M);
	double *row_max_val = (double *)calloc(sizeof(double), M);
	int global_max_idx[2];
	double global_max_val = 0;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (row_max_val[i] < img[i][j])
			{
				row_max_val[i] = img[i][j];
				row_max_idx[i] = j;
			}
		}
		// printf("%f\n", row_max_val[i]);
		if (row_max_val[i] > 0)
		{
			// lowest row number
			if (row_mask_end[0] > i)
			{
				row_mask_end[0] = i;
			}
			// highest row number
			if (row_mask_end[1] < i)
			{
				row_mask_end[1] = i;
			}
			// global max value and index
			if (global_max_val < row_max_val[i])
			{
				global_max_val = row_max_val[i];
				global_max_idx[0] = i;
				global_max_idx[1] = row_max_idx[i];
			}
		}
	}

	double th = 0.1;
	double cur_th;
	// gradient and its vertical line gradient
	double k, k1;
	// temp index variable
	int m, n;
	// graybarycenter
	double sum, weighted_sum;
	// subpixel index
	double subpxl_i, subpxl_j;

	// laserline test;
	test->img_round = (double **)fspace_2d(M, N, sizeof(double));
	test->m = (double *)calloc(sizeof(double), M);
	test->m_p = test->m;
	test->n = (double *)calloc(sizeof(double), M);
	test->n_p = test->n;

	// upward
	for (int i = row_mask_end[1]; i >= row_mask_end[0]; i--)
	{
		if (row_max_val[i] <= 0)
		{
			continue;
		}

		sum = row_max_val[i];
		weighted_sum = row_max_val[i]*row_max_idx[i];

		cur_th = th*row_max_val[i];

		m = i;
		// towards the left
		n = row_max_idx[i] - 1;
		while (img[m][n] >= cur_th)
		{
			sum += img[m][n];
			weighted_sum += n*img[m][n];
			n -= 1;
		}

		// towards the right
		n = row_max_idx[i] + 1;
		while (img[m][n] >= cur_th)
		{
			sum += img[m][n];
			weighted_sum += n*img[m][n];
			n += 1;
		}
		subpxl_j = weighted_sum/sum;
		test->img_round[i][round_nngt(subpxl_j)] = 1;
		*test->m_p++ = i;
		*test->n_p++ = subpxl_j;
	}

	free(row_mask);
	free(row_max_idx);
	free(row_max_val);

	// int points_num;
	// points_num = test->m_p - test->m;
	// printf("%d\n", points_num);

	// return &test;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[])
{
	double *input;
	int M, N;
	input = mxGetPr(prhs[0]);
	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);
	// mexPrintf("M = %d\nN = %d\n", M, N);

	// resize to 2d
	double **input_re = (double **)fspace_2d(M, N, sizeof(double));
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			input_re[i][j] = input[j*M + i];
		}
	}

	// define laserline struct
	laserline test;
	laserline *test_p = &test;
	process(test_p, input_re, M, N);


	double *output_img, *output_m, *output_n;
	int points_num;
	points_num = test.m_p - test.m;
	// printf("%d\n", points_num);

	// get output array
	plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(points_num,1,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(points_num,1,mxREAL);
	output_img = mxGetPr(plhs[0]);
	output_m = mxGetPr(plhs[1]); 
	output_n = mxGetPr(plhs[2]); 
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < N; j++)
		{
			output_img[j*M + i] = test.img_round[i][j];
		}
	}
	for (int i = 0; i <= points_num; i++)
	{
		output_m[i] = test.m[i] + 1;
		output_n[i] = test.n[i] + 1;
	}

	free(test.m);
	free(test.n);
	ffree_2d(input_re, M);
	ffree_2d(test.img_round, M);

}
