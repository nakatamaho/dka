#include <stdio.h>
#include <math.h>
#include <float.h>

#define INFINITY		DBL_MAX
#define MINUS_INFINITY	((-)DBL_MAX)
#define MAX				20
#define max(x,y)		( (x) > (y) ? (x) : (y) )
#define min(x,y)		( (x) < (y) ? (x) : (y) )
#define epsilon			1.0e-10
#define epsilon2		0.00001
#define OFF				0
#define ON				1
#define TRUE			1
#define FALSE			0
#define SolvePolynomial	aberth
#include "complex.h"

complex phi ( complex z, complex w[], int nth, double f[], double fprime[] )
{
  complex tmp1, tmp2, tmp3 ;
  int i, j;
  double r;
 
  tmp1 = SubstituteToPolynomial ( z, fprime, nth-1 );
  if ( ComplexAbs( tmp1 ) > epsilon ) {
	tmp3 = ComplexDiv (  SubstituteToPolynomial ( z, f, nth ),
					   tmp1 );
  } else  return tmp1;

  tmp2.real = tmp2.image = 0.0;
  for ( i = 0 ; i < nth ; i ++ ) {
	tmp1 = ComplexSub ( z, w[i] );
	r = ComplexAbs ( tmp1 );
	if ( r < epsilon ) continue;
	tmp1.real  =  tmp1.real  / r;
	tmp1.image = -tmp1.image / r;
	tmp2 = ComplexAdd ( tmp1, tmp2 ); 
  }
  tmp1 = ComplexMul ( tmp2, tmp3 );
  tmp1.real = - 1.0 + tmp1.real;
  return ComplexDiv ( tmp3 , tmp1 );
}

int aberth ( double f[], int nth, double a, double b, double *answer )
{
  int i, j, flag;						/* for loop */
  double g[MAX],gabs[MAX];
  double fprime[MAX],gabsprime[MAX];
  double beta;					/* kai no jyuusin */
  double radius, tmp1, tmp2;
  double max_error;
  complex roots[MAX];
  complex tmp, tmp3;

  *answer = INFINITY; 
  i = 0 ;
  while ( f[i] < epsilon ) i ++;
  nth = nth - i;
  if ( nth < 1 ) return FALSE;
  for ( j = 0 ; j < nth+1 ; j ++ ) f[j] = f[j+i];

  beta = - f[1] / ( nth * f[0] );
  for ( i = 0 ; i < nth+1 ; i ++ ) g[i]=f[i];
  for ( i = 1 ; i < nth+1 ; i ++ ) 
	for ( j = i ; j > 0  ; j -- ) 
	  g[j] = g[j-1] * beta + g[j];

  for ( i = 0 ; i < nth+1 ; i ++ ) {
	gabs[i] = - fabs(g[i]);
	fprime[i] = ( nth - i ) * f[i];
	gabsprime[i] = (nth -i)*(gabs[i]);
  }
  gabs[0] = - gabs [0];
  gabsprime[0] = -gabsprime[0];
  radius = 0.0 ;
  for ( i = 2 ; i <= nth ; i ++ ) {
	tmp1 = -  ((double)nth * gabs[i] / g[0]) /i;     /*	          | n * g[i] |1/i    */
    if ( fabs ( tmp1 ) <= epsilon ) continue;	     /*	    max   |----------|		 */
   	radius = max ( exp ( log ( tmp1 )) , radius );   /*  i=1,..,n |   g[0]   |		 */			
  }

  while ( 1 ) {
	tmp1 = gabs[0];
	tmp2 = gabsprime[0];
	for ( i = 1 ; i < nth+1 ; i ++ ) tmp1 = radius * tmp1 + gabs[i];
	for ( i = 1 ; i < nth   ; i ++ ) tmp2 = radius * tmp2 + gabsprime[i];
	if ( fabs ( tmp2 ) < epsilon ) { radius = epsilon ; break; }
	tmp1 = tmp1 / tmp2;
	if ( fabs ( tmp1 ) < epsilon ) break;
	radius = radius - tmp1;
  }
  for ( i = 0 ; i < nth ; i ++ ) {
	tmp1 = ( 2.0 * PI * ( i - 1.0 ) + 1.5 ) / (double)nth ;
	roots[i].image = radius * sin ( tmp1 ) ;
	roots[i].real  = beta + radius * cos ( tmp1 ) ;
  }
  j  = 0;
  max_error = 1.0;
  while (  max_error > epsilon ) {
	max_error = 0.0;
	for ( i = 0 ; i < nth ; i ++ ) {
	  tmp = phi ( roots[i], roots, nth, f, fprime );
	  roots[i] = ComplexAdd ( roots [i] , tmp );
	  max_error = max ( max_error, ComplexAbs ( tmp ) );
	}
	j ++ ;
#ifdef DEBUG_ABERTH
	wprintf("%d kaime.\n",j);
	for ( i = 0 ; i < nth ; i ++ ){
	  wprintf("z[%d]=",i); PrintComplex(roots[i]);
	}
	wprintf("\n");
#endif
  }
  flag = FALSE;
  for ( i = 0 ; i < nth ; i ++ ) {
	if ( fabs ( roots[i].image ) < epsilon2 && 
		roots[i].real >= a && roots[i].real <= b ) {
		*answer = min ( *answer, roots[i].real );
		flag = TRUE;
	  }
  }
  return flag;
}
