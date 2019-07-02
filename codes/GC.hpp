#include <algorithm>
struct  MatVirt;

struct MatVirt {
public:
    int n,m;
    virtual  double * addmatmul(double *x,double *Ax) const =0;
    virtual ~MatVirt() {}

    double * matmul(double *x,double *Ax) const
    {
        std::fill(Ax,Ax+n, 0.);
        return addmatmul(x,Ax);
    }
    MatVirt(int nn,int mm) : n(nn),m(mm) {}
    MatVirt(int nn) : n(nn),m(nn) {}
};


int GradientConjugue(MatVirt &A, // fonction et pointeur data pour A
		     MatVirt &C, // fonction et pointeur data pour C
		     double * b, // second membre
		     double * x, // solution qui contient une initialisation
		     int nbitermax,double eps,
                     int niveauimpression)
;

inline double * ProduitMatVec(MatVirt *A,double *x, double *Ax) { return A->matmul(x,Ax);}
inline double * ProduitMatVec(MatVirt &A,double *x, double *Ax) { return A.matmul(x,Ax);}
double * myscopy(int n,const double *x,double *y);
double mysdot(int n,double *x,double *y);
double * mysaxpy(int n,double a,double *x,double *y);
double * myscal(int n,double a,double *x);
// Matrice Identite
double * SetArray(double *x,int n,double c=0.);	
struct MatIdentite : public MatVirt
{
    MatIdentite(int nn):MatVirt(nn) {}
    double * addmatmul(double *x,double *Ax) const {
          mysaxpy(n,1.,x,Ax); return Ax;}
};
