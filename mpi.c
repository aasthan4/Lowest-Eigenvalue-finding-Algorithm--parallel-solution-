#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
//#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include "mpi.h"
#define n 100000
#define P 32
double ele(double i, double j)
{
	if (i==j) 
	{
		return (-1.0)/(2*i+1.0);
	}
	else
	{
		return (-1.0)/(10.0*(i+j+1.0));
	}
}
/* Subtract the `struct timeval' values X and Y,
    storing the result in RESULT.
    Return 1 if the difference is negative, otherwise 0. */
int timeval_subtract (struct timeval * result, struct timeval * x, struct timeval * y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
    
  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}
int start(int id, int nproc, int wid)
{
	if (id ==(nproc-1))
	{
		return n-1;
	}
	else
	{
		return (id+1)*wid;
	}
}
int main()
{
// Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Get the rank of the process
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);


    // Print off a hello world message
//    printf("Hello world from processor , rank %d"
  //         " out of %d processors\n", id, nproc);
	

	double as[100];
	for(int x=0; x<100; x++) as[x]=0;
	double alpha=0.0, alpha0=ele(0,0);
	int iter=0, tid, co;
	double threshhold=10E-7, delta0i;
	double c[n], alpha_g;
	double Delta0i=0.0, Ap0i, Ap00, Apii;
	struct timeval ta, tb, tresult;	
	
	gettimeofday ( &ta, NULL );

	int wid=(n-1)/nproc;



	while (fabs(alpha-alpha0)>threshhold)
	{
		alpha0=alpha;
		alpha=0.0;

		
		for(int i=start(id, nproc, wid);i>=(id*wid+1);i--)
		{
			Ap00=ele(0,0);
			Ap0i=ele(0,i);
			Apii=ele(i,i);
			
			if (iter!=0)
			{
				Delta0i=0.0;
				for(int j=1; j<=n-1;j++)
				{
					//printf("%d %d %d %lf\n", iter, i, j, c[j]);
					if (j!=i) Delta0i=Delta0i+ele(i,j)*c[j];
				}
				
				Ap00=ele(0,0)-ele(0,i)*c[i]-Delta0i*c[i];
				Ap0i=Ap0i+Delta0i;
			}
			//if (i==1)
			//	printf("element %lf %lf %lf %lf\n", c[1], Delta0i, alpha0, ele(0,i));
			c[i]=Ap0i/(ele(0,0)+alpha0-ele(i,i));
			//if (i==1)
			//	printf("element after  %lf\n", c[1]);
			alpha=alpha+c[i]*ele(0,i);

		}
		//reduce to 0

		double *ctmp;
		ctmp=(double *) malloc((n*nproc)*sizeof(double));
		double alpha_g=0;
		MPI_Reduce(&alpha, &alpha_g, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
		MPI_Gather(&c, n, MPI_DOUBLE, ctmp, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

		//printf("%lf %lf\n",ctmp[0], c[0]);
		if (id==0)
		{
			c[0]=0.0;
			alpha=alpha_g;
			 co=0;
			//printf("iteration %d alpha %lf\n", iter, alpha+ele(0,0));
			for (int y=0; y<nproc; y++)
			for (int z=(y*n+(y*wid+1)); z<=(y*n+start(y, nproc, wid)); z++)
			{
				c[co]=ctmp[z];
				co++;
				//printf("%lf\n", c[co]); 
			}
			
			//printf("iteration %d ev %lf\n", iter, alpha);
		}
		free(ctmp);
		MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast(&c, n, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		//bradcast alpha

		printf("node %d iteration %d alpha %lf\n", id, iter, alpha+ele(0,0));
		iter++;
	}
	MPI_Finalize();
/* get initial time */
	gettimeofday ( &tb, NULL );
	timeval_subtract ( &tresult, &tb, &ta );
	if (id==0)
	printf ("%lu \t  %lu\n",  tresult.tv_sec, tresult.tv_usec );

	// Finalize the MPI environment.

}
