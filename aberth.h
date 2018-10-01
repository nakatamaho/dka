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

complex phi ( complex , complex [], int , double [], double []);
int aberth ( double [], int , double , double , double * );
