/* 
	Fast way to compute all pair shortest path for 2D binary shape contour

	Haibin Ling, 08/18/2004

*/

/* common functions */
#include "common_matlab.h"

void bellman_ford_allpair(double *dis_mat, double *ang_mat,
						  double *X, double *Y, int n_V,
						  double *E, int n_E);

/*-----------------------------------------------------------------------
     [dis_mat,ang_mat] = bellman_ford_ex_C(X,Y,E)

 ------------------------------------------------------------------------*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* Input data */
	double	*X	= mxGetPr(prhs[0]);
	double	*Y	= mxGetPr(prhs[1]);
    double	*E	= mxGetPr(prhs[2]);
    int		n_V	= mxGetM(prhs[0]);		/* Number of nodes */
    int		n_E	= mxGetM(prhs[2]);		/* Number of edges */

    /* Output data */
    mxArray	*pMxDisMat	= mxCreateDoubleMatrix(n_V,n_V,mxREAL);
    mxArray	*pMxAngMat	= mxCreateDoubleMatrix(n_V,n_V,mxREAL);
    double	*dis_mat	= mxGetPr(pMxDisMat);
    double	*ang_mat	= mxGetPr(pMxAngMat);

	/* call algorithm */
	bellman_ford_allpair(dis_mat, ang_mat,
						 X, Y, n_V,
						 E, n_E);

    /* Return */
    plhs[0] = pMxDisMat;
    plhs[1] = pMxAngMat;
}



//----------------------------------------------------------------------
void bellman_ford_allpair(double *dis_mat, double *ang_mat,
						  double *X, double *Y, int n_V,
						  double *E, int n_E)
{
	const double INF_DIS	= n_V*n_E;
	const double INF_THRE	= INF_DIS-2;
	const double PI	= 3.14159265;
	int		n_E1,i,s,u,v,e;
	double	dx,dy,w,ang;
	bool	bstop;

	int		*Ep1 = new int[n_E];
	int		*Ep2 = new int[n_E];
	double	*Ew	 = new double[n_E];
	
	// initialize dis_mat and ang_mat
	double	*pD,*pA;
	pD	= dis_mat;
	pA	= ang_mat;
	for(u=0;u<n_V*n_V;u++) {
		*(pD++)	= INF_DIS;
		*(pA++)	= -10;
	}

	// set viewable values
	double *pu	= E;
	double *pv	= E+n_E;
	double *pw	= pv+n_E;
	for(e=0;e<n_E;e++)
	{
		u		= (int)pu[e]-1;
		v		= (int)pv[e]-1;
		w		= pw[e];
		Ep1[e]	= u;
		Ep2[e]	= v;
		Ew[e]	= w;
		
		// set dis_mat
		MAT_SET(dis_mat,u,v,w,n_V);
		MAT_SET(dis_mat,v,u,w,n_V);

		// set ang_mat
		dx	= X[v]-X[u];
		dy	= Y[v]-Y[u];
		ang	= atan2(dy,dx);
		MAT_SET(ang_mat,u,v,ang,n_V);
		MAT_SET(ang_mat,v,u,(ang>0)?(ang-PI):(ang+PI),n_V);
	}


	// shortest pathes from every start point s
	int		*U	= new int[2*n_E];
	int		*V	= new int[2*n_E];
	double	*W	= new double[2*n_E];
	for(s=0;s<n_V;s++)
	{
		pD		= dis_mat+s*n_V;
		pA		= ang_mat+s*n_V;
		pD[s]	= 0;
		//pA[s]	= -100;

		/*/ reduce the size of graph
		memcpy(U,Ep1,n_E*sizeof(int));
		memcpy(V,Ep2,n_E*sizeof(int));
		memcpy(W,Ew, n_E*sizeof(double));	
		n_E1	= n_E;/*/
		n_E1	= 0;
		for(e=0;e<n_E;e++)
		{
			u	= Ep1[e];
			v	= Ep2[e];
			w	= Ew[e];
			
			if(pD[u]>INF_THRE)	// this node is not viewable from s
			{
				U[n_E1]	= v;
				V[n_E1]	= u;
				W[n_E1]	= w;
				n_E1++;
			}

			if(pD[v]>INF_THRE)	// this node is not viewable from s
			{
				U[n_E1]	= u;
				V[n_E1]	= v;
				W[n_E1]	= w;
				n_E1++;
			}
		}


		/* Relaxation using standard bellman-ford*/
		bstop	= false;
		for(i=1; i<n_V-1 && !bstop; i++)
		{
			bstop	= 1;
			for(e=0; e<n_E1; ++e)
			{
				/* Relax for each edge */
				u	= U[e];
				v	= V[e];
				w	= W[e];
				if(pD[v]>pD[u]+w)	{
					pD[v]	= pD[u]+w;
					pA[v]	= pA[u];
					//pP[v] = u+1;
					bstop = 0;
				}

				/*
				if(pD[u]>pD[v]+w)	{
					pD[u]	= pD[v]+w;
					pA[u]	= pA[v];
					//pP[v] = u+1;
					bstop = 0;
				}/**/
			}
			/*printf("i=%d\n",i);*/
		}
	}//of start point s


	delete	[]U;
	delete	[]V;
	delete	[]W;
	delete	[]Ep1;
	delete	[]Ep2;
	delete	[]Ew;
}
