
typedef struct {
        double real, image;
} COMPLEX;

complex ComplexAdd ( complex a, complex b );
complex ComplexSub ( complex a, complex b );
complex ComplexMul ( complex a, complex b );
double ComplexAbs ( complex a );
complex ComplexConj ( complex a );
complex SubstituteToPolynomial ( complex a, double f[], int n_th );
complex ComplexDiv ( complex a, complex b ); 
void PrintComplex ( complex a );

