#include "mex.h"
#include <math.h>
#include <vector>

using namespace std;
const float pi = 3.1415926;

class matrix2d
{
public:
	float * data;
	int M;
	int N;
	
	float & get (int r, int c)
	{
		return data[c*M+r];
	}
	float & get (int ind)
	{
		return data[ind];
	}
	void create(int MM, int NN)
	{
		M = MM;
		N = NN;
		data = new float[M*N];
		memset(data, 0, M*N*sizeof(float));
	}
	void clear()
	{
		memset(data, 0, M*N*sizeof(float));
	}
	void destroy()
	{
		if (data)
		{
			delete [] data;
		}
	}
};

struct xpair
{
	float val;
	int ind;
};

vector<xpair> bubble_sort(matrix2d & arr, float K)
{
	vector<xpair> neibs(K);
	vector<xpair> data(arr.M);
	xpair tmp;

	for (int m=0; m<arr.M; ++m)
	{
		data[m].val = arr.get(m);
		data[m].ind = m;
	}

	for (int k=0; k<K; ++k)
	{
		for (int i=0; i<data.size()-k-1; ++i)
		{
			if ( data[i].val < data[i+1].val )
			{
				tmp = data[i];
				data[i] = data[i+1];
				data[i+1] = tmp;
			}
		}
		// the smallest value is in the end of the vector
		neibs[k] = data[data.size()-k-1];
	}

	return neibs;
}

float normpdf(const float x, const float mu, const float sigma)
{
	return exp( -0.5*pow(x-mu,2)/pow(sigma,2) ) / (sigma * sqrt(2*pi));
}

void parse_input(int nrhs, const mxArray *prhs[], matrix2d& fea, 
	matrix2d& dict, matrix2d& sigma, float& knn, float& C)
{
	fea.data = (float *)mxGetPr(prhs[0]);
	fea.M = mxGetM(prhs[0]);
	fea.N = mxGetN(prhs[0]);

	dict.data = (float *)mxGetPr(prhs[1]);
	dict.M = mxGetM(prhs[1]);
	dict.N = mxGetN(prhs[1]);

	sigma.data = (float *)mxGetPr(prhs[2]);
	sigma.M = mxGetM(prhs[2]);
	sigma.N = mxGetN(prhs[2]);

	knn = *(float *)mxGetPr(prhs[3]);

	C = *(float *)mxGetPr(prhs[4]);
}

void flush_output(int nlhs, mxArray *plhs[], matrix2d & code)
{
	matrix2d out;
	plhs[0] = mxCreateNumericMatrix(code.M, code.N, mxSINGLE_CLASS, mxREAL);
	out.data = (float * )mxGetPr(plhs[0]);
	out.M = code.M;
	out.N = code.N;
	for(int x=0; x<out.M; ++x)
	{
		for(int y=0; y<out.N; ++y)
		{
			out.get(x,y) = code.get(x,y);
		}
	}
}

/*
	prototype:	[code]=rbf_coding(fea, dict, sigma, knn)
	fea: n*D matrix
	dict: m*D matrix, m is the size of dictionary
	sigma: m*1 matrix
	knn: number of nearest neighbor
	C: a factor by sigma
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs == 0)
	{
		mexPrintf("Input must be float.\n");
		mexPrintf("prototype:	[code]=rbf_coding(fea, dict, sigma, knn)\n"); 
		mexPrintf("fea: n*D matrix\n");
		mexPrintf("dict: m*D matrix, m is the size of dictionary\n");
		mexPrintf("sigma: m*1 matrix\n");
		mexPrintf("knn: number of nearest neighbor\n");
		mexPrintf("C: a factor by sigma\n");
		return;
	}
	matrix2d fea, dict, sigma;
	float knn, C;

	parse_input(nrhs, prhs, fea, dict, sigma, knn, C);

	matrix2d rbf_code;
	rbf_code.create(fea.M, dict.M);
	matrix2d arr_dist;
	arr_dist.create(dict.M, 1);
	for(int m_fea=0; m_fea<fea.M; ++m_fea)
	{
		// calclate the L2 distance between the feature 
		// vector and the codebook
		arr_dist.clear();
		for(int m_dict=0; m_dict<dict.M; ++m_dict)
		{
			for(int n=0; n<dict.N; ++n)
			{
				arr_dist.get(m_dict) += pow( fea.get(m_fea,n)-dict.get(m_dict,n), 2 );
			}
			arr_dist.get(m_dict) = sqrt( arr_dist.get(m_dict) );
		}

		// sort the distance array in a ascend way
		// find both the distcancr and the orginal index
		// using the bubble sort algorithm
		vector<xpair> neibs = bubble_sort(arr_dist, knn);
		for (size_t i=0; i<neibs.size(); ++i)
		{
			rbf_code.get(m_fea, neibs[i].ind) = normpdf(neibs[i].val, 0, C*sigma.get(neibs[i].ind) );
		}

	}

	flush_output(nlhs, plhs, rbf_code);

	rbf_code.destroy();
	arr_dist.destroy();
}