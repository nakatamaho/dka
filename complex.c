#include "complex.h"

complex ComplexAdd ( complex a, complex b )
{
  a.real  = a.real  + b.real  ;
  a.image = a.image + b.image ;
  return a;
}

complex ComplexSub ( complex a, complex b )
{
  a.real  = a.real  - b.real  ;
  a.image = a.image - b.image ;
  return a;
}

complex ComplexMul ( complex a, complex b ) 
{
  complex tmp;
  tmp.real  =  a.real  * b.real -  a.image * b.image ;
  tmp.image =  a.image * b.real +  a.real  * b.image ;
  return tmp;
}

double ComplexAbs ( complex a )
{
  return ( a.real * a.real + a.image * a.image ) ;
}

complex ComplexConj ( complex a )
{
  a.image = - a.image ;
  return a;
}

complex SubstituteToPolynomial ( complex a, double f[], int nth )
{
  int i;
  complex tmp;
  tmp.real  = f[0] ;
  tmp.image = 0.0  ;
  for ( i = 1 ; i < nth + 1 ; i ++ ) {
	tmp = ComplexMul ( a, tmp );
	tmp.real = f[i] + tmp.real ;
  }
  return tmp;
}

complex ComplexDiv ( complex a, complex b ) 
{
  complex tmp;
  double tmp1;
  tmp1 = ComplexAbs ( b ) ;
  if ( tmp1 == 0.0 ) { 
	printf ("divison error!\n");
	exit ();										/* koreha atode tukurinaosu. */
  }
  tmp.real  = ( a.real  * b.real  +  a.image * b.image ) /tmp1 ;
  tmp.image = ( a.image * b.real  -  a.real  * b.image ) /tmp1 ;
  return tmp;
}

void PrintComplex ( complex a )
{
  printf ("%3.6lf+%3.6lf\n",a.real,a.image);
}
