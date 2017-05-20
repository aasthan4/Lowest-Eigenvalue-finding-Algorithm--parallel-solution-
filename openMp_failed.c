#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#define n 10000
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
int main()
{
	int proc= omp_get_max_threads();
	printf("no of threads on the machine = %d\n", proc);
	proc=P;
	omp_set_num_threads(proc);
	double as[100];

	for(int x=0; x<100; x++) as[x]=0;
	double alpha=0.0, alpha0=ele(0,0);
	int iter=0, tid;
	double threshhold=10E-7, delta0i;
	double c[n];
	double Delta0i=0.0, Ap0i, Ap00, Apii;
struct timeval ta, tb, tresult;	
/* get initial time */
  gettimeofday ( &ta, NULL );
	while (fabs(alpha-alpha0)>threshhold)
	{
		alpha0=alpha;
		alpha=0.0;
		for(int i=n-1;i>=1;i--)
		{
			Ap00=ele(0,0);
			Ap0i=ele(0,i);
			Apii=ele(i,i);
			if (iter!=0)
			{
				Delta0i=0.0;
				#pragma omp parallel for default(shared) private(tid) 	
				for(int j=1; j<=n-1;j++)
				{
					tid = omp_get_thread_num();
					if (j!=i) as[tid]=as[tid]+ele(i,j)*c[j];
				}
				for(int z=0; z<proc; z++)
				{
					Delta0i=Delta0i+as[z];
				}
				for(int x=0; x<100; x++) as[x]=0;
				Ap00=ele(0,0)-ele(0,i)*c[i]-Delta0i*c[i];
				Ap0i=Ap0i+Delta0i;
			}
			c[i]=Ap0i/(ele(0,0)+alpha0-ele(i,i));
			alpha=alpha+c[i]*ele(0,i);
		}
		iter++;
		printf("iteration %d alpha %lf\n", iter, alpha+ele(0,0));
	}
/* get initial time */
	  gettimeofday ( &tb, NULL );
	timeval_subtract ( &tresult, &tb, &ta );
	printf ("Program took  %lu seconds and %lu microseconds\n", tresult.tv_sec, tresult.tv_usec );


}
