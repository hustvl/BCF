/* 
	function [XIs,YIs] = uniform_interp(Xs,Ys,n_samp)
	
	Uniformly resampling point chain
	Haibin Ling, 09/04/2004

*/

/* common functions */
#include "common_matlab.h"

void uniform_interp(double *XIs, double *YIs, int n_samp, 
					double *Xs,  double *Ys,  int n_pt);

/*****************************************************************************
function [XIs,YIs] = uniform_interp(Xs,Ys,n_samp)
*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Input data */
	double	*Xs		= mxGetPr(prhs[0]);
	double	*Ys		= mxGetPr(prhs[1]);
    int		n_pt	= mxGetM(prhs[0]);				// Number of points
	int		n_samp	= ROUND(*mxGetPr(prhs[2]));		// Number of interpolation

    /* Output data */
    mxArray	*pMxXIs	= mxCreateDoubleMatrix(n_samp,1,mxREAL);
    mxArray	*pMxYIs	= mxCreateDoubleMatrix(n_samp,1,mxREAL);
    double	*XIs	= mxGetPr(pMxXIs);
    double	*YIs	= mxGetPr(pMxYIs);

	/* call algorithm */
	uniform_interp(XIs, YIs, n_samp,  Xs, Ys, n_pt);

    /* Return */
    plhs[0] = pMxXIs;
    plhs[1] = pMxYIs;
}


/*****************************************************************************
function [XIs,YIs] = uniform_interp(Xs,Ys,n_samp)
*/
void uniform_interp(double *XIs, double *YIs, int n_samp, 
					double *Xs,  double *Ys,  int n_pt)

{
	double	dx,dy,x,y,r;
	int		ii;

	// compute length from first point to all other points
	double	*seg_lens	= new double[n_pt];
	seg_lens[0]	= 0;
	for(ii=1;ii<n_pt;++ii)
	{
		dx	= Xs[ii]-Xs[ii-1];
		dy	= Ys[ii]-Ys[ii-1];
		seg_lens[ii]	= sqrt(dx*dx+dy*dy) + seg_lens[ii-1];
	}

	double	d_len	= seg_lens[n_pt-1]/(n_samp+1);
	int		i_fill	= 0;
	double	cur_len	= d_len;
	for(ii=1;ii<n_pt && i_fill<n_samp;++ii)
	{
		while(cur_len<=seg_lens[ii] && i_fill<n_samp)
		{
			// interpolate a point
			r	= (cur_len-seg_lens[ii-1]) / (seg_lens[ii]-seg_lens[ii-1]);
			x	= Xs[ii-1]+r*(Xs[ii]-Xs[ii-1]);
			y	= Ys[ii-1]+r*(Ys[ii]-Ys[ii-1]);
			XIs[i_fill] = x;
			YIs[i_fill]	= y;
			
			++i_fill;
			cur_len	+= d_len;
		}
	}

	delete	[]seg_lens;
}
