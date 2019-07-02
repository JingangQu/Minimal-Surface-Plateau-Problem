/*
   Jingang QU  Probleme de Plateau
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include "EF2d-base.hpp"
#include "GC.hpp"


using namespace std;

class R3 {

public:
    R x, y, z;  // declaration de membre   typedef double R
    R3 () :x(0.),y(0.),z(0.){} // rappel : x(0), y(0) z(0.) sont initialiser via le constructeur de double
    R3 (R a,R b,R c):x(a),y(b),z(c)  {}
    R3 (Vertex a,R b):x(a.x),y(a.y),z(b)  {}
    R3 (const R3 & a,const R3 & b):x(b.x-a.x),y(b.y-a.y),z(b.z-a.z)  {}
    // les autre operateur affectations
    R3 &  operator+=(const R3 & P)  {x += P.x;y += P.y;z += P.z; return *this;}
    R3 &  operator-=(const R3 & P) {x -= P.x;y -= P.y;z -= P.z;return *this;}
    // operateur binaire + -
    R3   operator+(const R3 & P)const {return R3(x+P.x,y+P.y,z+P.z);}
    R3   operator-(const R3 & P)const {return R3(x-P.x,y-P.y,z-P.z);}
    R3   operator^(const R3 & P)const {  // produit mixte
        return R3(y*P.z-z*P.y, z*P.x-x*P.z, x*P.y-y*P.x);
    }
    R norme() const { return sqrt(x*x+y*y+z*z);}
};

inline R3 det_R3(const R3 & A,const R3 & B,const R3 & C)  { return R3(A,B)^R3(A,C);}

// creer V
void build_V(double* V, int N){
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++){
            if (j==i)
                V[i * N + j] = 1;
        }
    }
}

// fonction de borde
double borde(double x, double y, int lab){
    // Probleme 4  Catenoide
    //if (lab == 1)
    //    return 1;
    //if (lab == 2)
    //    return 0.5;

    // Probleme 5
    return cos(x) * cos(y);

    // Probleme 6
    // return sin(2*M_PI*x);
    // return x*x + y*y + 5*cos(x*2*M_PI) + sin(y*3*M_PI);
    // return exp(x) + y*y*y + cos(x*10*M_PI) + sin(y*10*M_PI);
}


int BordeGlobIdx(vector<int> &BordeIdx, Mesh2d &Th){
    int j=0;
    for (int i=0; i< Th.nv; i++){
        Vertex Vi = Th.v[i];
        int LabelVi = Vi.lab;
        if (LabelVi != 0) {   // Déterminer si i est la borne
            BordeIdx.push_back(i);
            j++;
        }
    }
    return j;  // j : Nombre de points sur la borde
}

// creer g   g contient la coordonnée z de tous les points de borde
void build_g(double* g, const Mesh2d &Th, const vector<int> & BordeIdx, int Ng){
    for (int i=0; i<Ng; i++) {
        int j = BordeIdx[i];
        g[j] = borde(Th.v[j].x, Th.v[j].y, Th.v[j].lab);
    }
}

// creer W initial
void build_W0(double* W, const vector<int> & BordeIdx, const double* g, int N){
    for (int i=0; i<N; i++) {
        auto iter = find(BordeIdx.begin(), BordeIdx.end(), i);  // Déterminer si i est la borne
        if (iter != BordeIdx.end())
            W[i] = g[i];
        else
            W[i] = 0;    // Sauf 0, nous pouvons définir ce paramètre comme n'importe quelle valeur
    }
}

// calcul de MatA
void  MatrixA(double* MatA, double &aire, const Mesh2d &Th, const double* V,
          const double* W, const vector<int> & BordeIdx) {
    // Th:maillage, V, W, BordeIdx : Index global des points de borde

    for (int i = 0; i < Th.nv; i++) {    // i : Vertex
        for (int k = 0; k < Th.nt; k++) {  // k : Simplex
            Simplex S = Th[k];
            Vertex P0 = S[0], P1 = S[1], P2 = S[2];
            int ig0 = Th(k,0), ig1 = Th(k,1), ig2 = Th(k,2);  // index global

            R3 A0(P0, W[ig0]), A1(P1, W[ig1]), A2(P2, W[ig2]);
            R3 Nkw = det_R3(A0, A1, A2);
            double norme_Nkw = Nkw.norme();

            // Calculer l'aire de la surface
            if (i==0){
                aire += norme_Nkw/2;
            }

            R2 DNkv = V[ig0 + i * Th.nv]*(P2 - P1).perp() + V[ig1 + i * Th.nv] * (P0 - P2).perp()
                    + V[ig2 + i * Th.nv] * (P1 - P0).perp();
            double a0 = (DNkv.x * (P1.y - P2.y) + DNkv.y * (P2.x - P1.x)) / norme_Nkw;  // coefficient du point ig0
            double a1 = (DNkv.x * (P2.y - P0.y) + DNkv.y * (P0.x - P2.x)) / norme_Nkw;  // coefficient du point ig1
            double a2 = (DNkv.x * (P0.y - P1.y) + DNkv.y * (P1.x - P0.x)) / norme_Nkw;  // coefficient du point ig2

            MatA[ig0 + i * Th.nv] += a0;
            MatA[ig1 + i * Th.nv] += a1;
            MatA[ig2 + i * Th.nv] += a2;
        }
    }
}

// A*U = b
void DecompMatA(double* A, double* b, const double* MatA,
                const vector<int> & LocalIdx, const double* g, int N){
    int Idxi = 0, Idxj = 0;
    for (int i=0; i<N; i++){
        auto iter = find(LocalIdx.begin(), LocalIdx.end(), i);  // Déterminer si i est la borne
        if (iter == LocalIdx.end()) {  // i n'est pas la borne
            for (int j = 0; j<N; j++) {
                auto iter2 = find(LocalIdx.begin(), LocalIdx.end(), j);  // Déterminer si j est la borne
                if (iter2 == LocalIdx.end()) {  // j n'est pas la borne
                    A[Idxj++] = MatA[i * N + j];
                } else {
                    b[Idxi] += -MatA[i * N + j] * g[j];
                }
            }
            Idxi += 1;
        }
    }
}

// MatPleine
struct MatPleine : public MatVirt
{
    double *aij;
    MatPleine(int nn,int mm,double *p):MatVirt(nn,mm) , aij(p) {}
    double * addmatmul(double *x,double *Ax) const
    {
        int k=0;
        for(int j=0;j<m;++j)
            for(int i=0;i<n;++i)
            {
                Ax[i]+= aij[k++]*x[j];
            }
        return Ax;
    }
};


int main()  // int argc,const char ** argv
{
    Mesh2d Th("carre_PI.msh");
	{
        ofstream of("carre_PI");
        for(int k=0; k<Th.nt;++k)
        {
            for(int ip=0; ip<=3;++ip)
            {
                int i3=ip%3;
                //int i = Th(k,i3) == Th(Th[k][i3]);
                of << (R2) Th[k][i3] << " " << Th(Th[k][i3]) <<endl;
            }
            of << "\n\n";
			R2 G = ((R2) Th[k][0]+Th[k][1]+Th[k][2])/3;
			of << G << " " << k << "\n\n\n";
        }
    }

    int N = Th.nv, Ng;
    vector<int> BordeIdx;

    Ng = BordeGlobIdx(BordeIdx, Th);

    auto *W = new double[N];
    auto *MatA = new double[N * N];
    auto *V = new double[N * N];
    auto *A = new double[(N - Ng) * (N - Ng)];
    auto *g = new double[N];
    auto *b = new double[N];
    auto *u = new double[N - Ng];
    auto *pre_u = new double[N - Ng];

    fill(W, W + N, 0);
    fill(V, V + N * N, 0);
    fill(g, g + N, 0);
    fill(u, u + N - Ng, 0);

    build_V(V, N);
    build_g(g, Th, BordeIdx, Ng);
    build_W0(W, BordeIdx, g, N);

    // Algorithme de point fixe
    double error = 100.0, tol = 1.e-6;
    int Iter=0, IterMax = 100;
    vector<double> Vec_Aire;
    double aire;
    while (Iter < IterMax && error > tol) {

        error = 0.0;
        aire = 0.0;
        fill(MatA, MatA + N * N, 0);
        fill(A, A + (N - Ng) * (N - Ng), 0);
        fill(b, b + N, 0);

        MatrixA(MatA, aire, Th, V, W, BordeIdx);
        DecompMatA(A, b, MatA, BordeIdx, g, N);
        Vec_Aire.push_back(aire);

        for (int iu=0; iu<N-Ng; iu++)
            pre_u[iu] = u[iu];

        MatPleine GC_A(N - Ng, N - Ng, A);
        MatIdentite Id(N - Ng);
        GradientConjugue(GC_A, Id, b, u, N - Ng, 1e-9, 10);

        int i = 0;
        for (int j = 0; j < N; ++j) {
            auto itertor = find(BordeIdx.begin(), BordeIdx.end(), j);  // Déterminer si i est la borne
            if (itertor == BordeIdx.end()) {  // i n'est pas la borde
                W[j] = u[i];
                error += pow(u[i]-pre_u[i], 2);
                i++;
            }
        }
        cout << "Iteration : " << Iter << "\n";
        cout << "Error (Norme 2) : " << sqrt(error) << "\n";
        cout << "Aire de surface : " << aire << "\n";
        Iter++;
    }

    {
        ofstream of("U.txt");
        for (int i = 0; i < N; ++i)
            of << W[i] << " ";
    }

    {
        ofstream of("aire.txt");
        for (auto s : Vec_Aire)
            of << s << " ";
    }

    delete []W;
    delete []MatA;
    delete []V;
    delete []g;
    delete []A;
    delete []b;

    return 0;
}