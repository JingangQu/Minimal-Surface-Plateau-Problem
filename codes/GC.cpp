/*
   Jingang QU  Gradient Conjugue
*/
#include <cassert>
#include <cmath>
#include <iostream>
#include "GC.hpp"

#ifdef WITHBLAS
// -DWITHBLAS  (nivau compilation)
#ifdef __cplusplus
extern "C" {
#endif
 // BLAS function
 double ddot_(int *n, double *x, int *incx, double *y, int *incy);
 void  dcopy_(int *n, const double *x, int *incx, double *y,int *incy);
 void  dscal_(int *n, double *alpha, double *x, int *incx);
 void  daxpy_(int *n, double *alpha, double *x, int *incx,
                    double *y, int *incy);
 double dnrm2_(int *n, double *x, int *incx);
 void   idmin_(int *n, double *sx, int *incx);
    
#ifdef __cplusplus
}
#endif
//
double * myscopy(int n,const double *x,double *y)
{
    int un=1;
    dcopy_(&n,x,&un,y,&un);
    return y;
}
double mysdot(int n,double *x,double *y)
{
    int un=1;
    return  ddot_(&n,x,&un,y,&un);
}
double * mysaxpy(int n,double a,double *x,double *y)
{
    int un=1;
    daxpy_(&n,&a,x,&un,y,&un);
    return  y;
}

double * myscal(int n,double a,double *x)
{
    int un=1;
    dscal_(&n,&a,x,&un);
    return  x;
}
#else

double * myscopy(int n,const double *x,double *y)   // 将x拷贝给y
{
    
    for(int i=0;i<n;++i)
        y[i]=x[i];
    return y;
}
double mysdot(int n,double *x,double *y)   // 求两个向量的点积
{
    double s=0;
    for(int i=0;i<n;++i)
        s+= x[i]*y[i];
    return  s;
}
double * mysaxpy(int n,double a,double *x,double *y)   // a * x + y
{
    for(int i=0; i< n; ++i)
        y[i] += a* x[i];
    return  y;
}

double * myscal(int n,double a,double *x)   // a * x
{
    for(int i=0; i<n;++i)
        x[i] *= a;
    return  x;
}
#endif



double * SetArray(double *x,int n,double c)
{
    std::fill(x,x+n,c);
    return x;
}

int GradientConjugue( MatVirt &A, MatVirt &C,
                    double *b, // second membre
                    double *x, // solution qui contient une initialisation
                    int nbitermax,double eps,int niveauimpression)
{   // niveauimpression: 0: pas impression .. 10 a chaque iteration
     int n = A.n, nret=0;
    double NaN = nan("");
    double * G = SetArray(new double[n] ,n,NaN);
    double * CG= SetArray(new double[n],n,NaN);
    double *AH=CG;// meme tableau que CG Attention ??????
    double *H=SetArray(new double[n],n,NaN);
    double rho, gamma, gCgp, gCg, eps2;
    assert( A.m==n && C.m==n && C.n==n);
    niveauimpression=std::min(niveauimpression,10);
    int nprint =std::max((int)((nbitermax+1)*(10.-niveauimpression)/10.),1);
    mysaxpy(n,-1.,b,ProduitMatVec(A,x,G));// G = Ax -b
    gCg = mysdot(n,G,ProduitMatVec(C,G,CG)) ;
    myscal(n,-1,myscopy(n,CG,H)); // H =- CG;
    eps2 = eps*eps*gCg;
    assert( !(gCg != gCg) ); //  verif if NaN => bad Matrix..
    if(gCg < 1e-30)
    {   if(niveauimpression)
          std::cout << " GC: on a converge on 0 iteration ||g||_C^2" << gCg  << std::endl;
        nret = 2;}
    else
    for(int iter=1; iter <= nbitermax; ++iter)
    {
        gCgp = gCg;
        rho = -mysdot(n,G,H)/mysdot(n,H,ProduitMatVec(A,H,AH));
        mysaxpy(n,rho,H,x);
        mysaxpy(n,rho,AH,G);
        gCg = mysdot(n,G,ProduitMatVec(C,G,CG)) ;
        gamma = gCg/ gCgp;
        mysaxpy(n,-1.,CG,myscal(n,gamma,H));
        
        if(gCg < eps2) // We have converged ...
        {
            if(niveauimpression)
                std::cout << " GC: on a converge en iteration " <<iter
                     << " g=" << gCg << " rho=" << rho << " gamma=" <<gamma<<std::endl;;
            nret= 1;
            break;
        }
        else
          if ( (iter % nprint) == 0 )
                std::cout <<"  GC:iteration "<< iter << " rho "<< rho << " gamma "
                     <<gamma<< " ||g||_C^2:" << gCg <<std::endl;
    }
	if( nret==0) 
        std::cout <<"  GC:non convergence  "<< nbitermax << " rho "<< rho << " gamma "
             <<gamma<< " ||g||_C^2:" << gCg <<std::endl;
    delete[] H;
    delete[] CG;
    delete[] G;
    return nret;
}
