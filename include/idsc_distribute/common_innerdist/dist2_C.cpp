/* This is a C verion of following matlab code, Haibin Ling, 09/04/2004
    
	function n2 = dist2(x, c)
	%DIST2	Calculates squared distance between two sets of points.
	%
	%	Description
	%	D = DIST2(X, C) takes two matrices of vectors and calculates the
	%	squared Euclidean distance between them.  Both matrices must be of
	%	the same column dimension.  If X has M rows and N columns, and C has
	%	L rows and N columns, then the result has M rows and L columns.  The
	%	I, Jth entry is the  squared distance from the Ith row of X to the
	%	Jth row of C.
	%
	%	See also
	%	GMMACTIV, KMEANS, RBFFWD
	%
	%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

	[ndata, dimx] = size(x);
	[ncentres, dimc] = size(c);
	if dimx ~= dimc
		error('Data dimension does not match dimension of centres')
	end

	n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
  			ones(ndata, 1) * sum((c.^2)',1) - ...
  			2.*(x*(c'));
*/

#include "common_matlab.h"
void dist2(double* D, double* X1,int m1, double* X2, int m2, int dim);

/*****************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Input data */
	double	*X1		= mxGetPr(prhs[0]);
	double	*X2		= mxGetPr(prhs[1]);
    int		m1		= mxGetM(prhs[0]);
    int		m2		= mxGetM(prhs[1]);
    int		dim		= mxGetN(prhs[0]);
	

    /* Output data */
    mxArray	*pMxD	= mxCreateDoubleMatrix(m1,m2,mxREAL);
    double	*D		= mxGetPr(pMxD);

	/* call algorithm */
	dist2(D, X1,m1, X2, m2, dim);

    /* Return */
    plhs[0] = pMxD;
}


/*****************************************************************************/
void dist2(double* D, double* X1,int m1, double* X2, int m2, int dim)
{
	int i,j,k;
	double	v,t;
	if(!D||!X1||!X2)	return;

	for(i=0;i<m1;i++){
		for(j=0;j<m2;j++){
			v = 0;
			for(k=0;k<dim;k++)	{	
				t	= MAT_GET(X1,k,i,m1)-MAT_GET(X2,k,j,m2);
				v  += t*t;
			}
			MAT_SET(D,j,i,v,m1);
		}
	}
	return;
}
