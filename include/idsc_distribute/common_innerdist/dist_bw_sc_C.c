#include "common_matlab.h"

/*-----------------------------------------------------------------------
	 functions 
*/
void compu_dist_matrix( double *pSC1, double *pSC2, 
					    int nSamp1, int nSamp2, int nBin1, int nBin2,
						double *cost_mat
					   );



/*-----------------------------------------------------------------------
[dis,cost_mat] = dist_bw_sc_C( sc1, sc2, n_dist, n_theta, n_dis_type)
	
	sc1,sc2	: input shape context of two object, 
			  each COLUMN of sc1,sc2 is a chape contex feature at a given point.
	
	type	:  0 - using the "Hausdorff" distance as in Belongie's paper
			   1 - find the minimum distance with respect all possible global
					rotations
		  
% Compute cost matrix
sc1			= sc1/mean(sum(sc1,1));
sc2			= sc2/mean(sum(sc2,1));
costmat		= hist_cost(sc1,sc2);		% in C code to avoid memory insufficiency

% calculate shape context cost
[a1,b1]=min(costmat,[],1);
[a2,b2]=min(costmat,[],2);
sc_cost=max(mean(a1),mean(a2));
  
 ------------------------------------------------------------------------*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	double	*pSC1,   *pSC2,   *pType;

	mxArray	*pMxCostMat, *pMxDist;	/* output */
	double	*cost_mat,   *pDist;
	double	*cost_mat_tmp;

	double	*sc1,*sc2,*pCur;
	double	*pRowDist;
	int		nType,nSamp1,nSamp2,nBin1,nBin2;
	int		n_theta,n_dist;
	int		r,c,k,k1,k2,ia,id,ir, n_delta, N;
	double	dis,tmp,minDis,dist1,dist2,res_dis;


    /* Analyse input data and parameters*/
    
	pSC1	= mxGetPr(prhs[0]);
    nSamp1	= mxGetN(prhs[0]);
    nBin1	= mxGetM(prhs[0]);
    
    pSC2	= mxGetPr(prhs[1]);
    nSamp2	= mxGetN(prhs[1]);
    nBin2	= mxGetM(prhs[1]);

	n_theta	= 8;
	n_dist	= 8;
	nType	= 0;

	if(nrhs>4) {
		n_dist	= ROUND(*(mxGetPr(prhs[2])));
		n_theta	= ROUND(*(mxGetPr(prhs[3])));
		nType	= ROUND(*(mxGetPr(prhs[4])));
	}
	else if(nrhs>2) {
		nType	= ROUND(*(mxGetPr(prhs[2])));
	}


	/* Create output matrices and initilize */
	
	pMxDist		= mxCreateDoubleMatrix(1,1,mxREAL);
	pDist		= mxGetPr(pMxDist);
	pMxCostMat	= mxCreateDoubleMatrix(nSamp1,nSamp2,mxREAL);
	cost_mat	= mxGetPr(pMxCostMat);


	/*------------------------------------------------
	type	:  0 - using the "Hausdorff" distance as in Belongie's paper
			   1 - find the minimum distance with respect all possible global
					rotations	/
	printf("nType=%d\n",nType);	/**/


	pRowDist = (double*)malloc(nSamp1*sizeof(double));
	if(pRowDist==0)	printf("error when malloc!");

	switch(nType)
	{
	case 1:
		/* compute the distance between the r-th point in shape 1
									and the c-th point in shape 2 */
		cost_mat_tmp	= (double*)malloc(nSamp1*nSamp2*sizeof(double));
		res_dis	= 10000;
		for(ir=0;ir<n_theta;ir+=2)
		{
			/* printf("ir=%d\n",ir);	*/
			for(r=0;r<nSamp1;r++)		pRowDist[r] = 10000;

			n_delta	= ir*n_dist;
			dist2	= 0;
			for(c=0;c<nSamp2;c++)
			{
				sc2		= pSC2+c*nBin2;
				minDis	= 10000;
				for(r=0;r<nSamp1;r++)
				{
					sc1	= pSC1+r*nBin1;
					dis	= 0;

					for(k1=0;k1<nBin1;k1++)
					{
						k2	= k1+n_delta;
						if(k2>=nBin2)	
							k2 = k2-nBin2;
						tmp	= sc1[k1]+sc2[k2];
						if(tmp>0.000001)
							dis	+= (sc1[k1]-sc2[k2])*(sc1[k1]-sc2[k2])/tmp;
					}
					
					MAT_SET(cost_mat_tmp,c,r,dis,nSamp1);
					if(dis<minDis)					minDis		= dis;
					if(dis<pRowDist[r])				pRowDist[r] = dis;
				}
				dist2	+= minDis;
			}
			dist2 /= nSamp2;

			/* compute the shape context distance */
			dist1	= 0;	
			for(r=0;r<nSamp1;r++)		dist1 += pRowDist[r];
			dist1 /= nSamp1;		

			minDis	= (dist1<dist2)? dist2:dist1;

			if(minDis<res_dis) {
				res_dis	= minDis;
				memcpy(cost_mat,cost_mat_tmp,nSamp1*nSamp2*sizeof(double));
			}
		
		}/* of ir=0..n_theta */
		free(cost_mat_tmp);
	break;
	

	case 0:
	default:
		/*
		N	= nBin1*nSamp1;
		for(r=0; r<N; ++r)
			if(pSC1[r]<0.00001)		pSC1[r]=0.00001;
		N	= nBin2*nSamp1;
		for(c=0; c<N; ++c)
			if(pSC2[c]<0.00001)		pSC2[c]=0.00001;
		//*/

		/* compute the distance between the r-th point in shape 1
									and the c-th point in shape 2 */
		for(r=0;r<nSamp1;r++)		pRowDist[r] = 10000;

		dist2	= 0;
		for(c=0;c<nSamp2;c++)
		{
			sc2		= pSC2+c*nBin2;
			minDis	= 10000;
			for(r=0;r<nSamp1;r++)
			{
				sc1	= pSC1+r*nBin1;
				dis	= 0;
				for(k=0;k<nBin1;k++)
				{
					tmp	= sc1[k]+sc2[k];
					if(tmp>0.000001)
						dis	+= (sc1[k]-sc2[k])*(sc1[k]-sc2[k])/tmp;
				}

				MAT_SET(cost_mat,c,r,dis,nSamp1);
				if(dis<minDis)					minDis		= dis;
				if(dis<pRowDist[r])				pRowDist[r] = dis;
			}
			dist2	+= minDis;
		}
		dist2 /= nSamp2;

		/* compute the shape context distance */
		dist1	= 0;	
		for(r=0;r<nSamp1;r++)		dist1 += pRowDist[r];
		dist1 /= nSamp1;		

		res_dis	= (dist1<dist2)?dist2:dist1;
		break;
	}


	/* return */
	free(pRowDist);
	*pDist	= res_dis;
	plhs[0] = pMxDist;
	plhs[1] = pMxCostMat;
}



/*-------------------------------------------------------------*/
void compu_dist_matrix( double *pSC1, double *pSC2, 
					    int nSamp1, int nSamp2, int nBin1, int nBin2,
						double *cost_mat
					   )
{
	double	*sc1,*sc2,*pCur;
	int		r,c,k;
	double	dis,tmp;

	/* compute the distance between the r-th point in shape 1
								and the c-th point in shape 2 */
	pCur	= cost_mat;
	for(c=0;c<nSamp2;c++)
	{
		sc2	= pSC2+c*nBin2;
		for(r=0;r<nSamp1;r++)
		{
			sc1	= pSC1+r*nBin1;
			dis	= 0;
			for(k=0;k<nBin1;k++)
			{
				tmp	= sc1[k]+sc2[k];
				if(tmp>0.000001)
					dis	+= (sc1[k]-sc2[k])*(sc1[k]-sc2[k])/tmp;
			}

			*pCur	= dis;
			pCur++;
		}
	}
}
