/*
   Provided by teacher
*/

#include <cmath>
#include <cassert>
// Definition de la class R2 
//  sans compilation separe toute les fonctions 
// sont definies dans ce R2.hpp avec des inline 
//  
// definition R (les nombres reals)
// remarque la fonction abort est defini dans 
#include <cstdlib> 
#include <iostream> 
typedef double R;

// The class R2
class R2 {

public:  
  R x,y;  // declaration de membre 
  // les 3 constructeurs ---
  R2 () :x(0.),y(0.) {} // rappel : x(0), y(0)  sont initialiser via le constructeur de double 
  R2 (R a,R b):x(a),y(b)  {}
  R2 (const R2 & a,const R2 & b):x(b.x-a.x),y(b.y-a.y)  {}
  // le constucteur par defaut est inutile
  /*
  R2 (const R2 & a) :x(a.x),y(a.y) { cout << "  const par copy " << endl;} 
  R2 ( R2 & a) :x(a.x),y(a.y) { cout << "  const pas const  par copy " << endl; } 
  */
  // rappel les operator definis dans une class on un parametre
  // cache qui est la classe elle meme (*this)

  // les operateurs affectation
  //  operateur affection (*this) = P est inutil par defaut il fait le travail correctement
  /*
  R2 &  operator=(const R2 & P)  {x = P.x;y = P.y;return *this;}
  */
  // les autre operateur affectations
  R2 &  operator+=(const R2 & P)  {x += P.x;y += P.y;return *this;}
  R2 &  operator-=(const R2 & P) {x -= P.x;y -= P.y;return *this;}
  // operateur binaire + - * , ^ /
  R2   operator+(const R2 & P)const   {return R2(x+P.x,y+P.y);}
  R2   operator-(const R2 & P)const   {return R2(x-P.x,y-P.y);}
  R    operator,(const R2 & P)const  {return  x*P.x+y*P.y;} // produit scalaire
  R    operator^(const R2 & P)const {return  x*P.y-y*P.x;} // produit mixte
  R2   operator*(R c)const {return R2(x*c,y*c);}
  R2   operator/(R c)const {return R2(x/c,y/c);}
  // operateur unaire
  R2   operator-()const  {return R2(-x,-y);} 
  const R2 &  operator+()const  {return *this;}
  // un methode
  R2   perp() const {return R2(-y,x);} // la perpendiculaire
  // les operators  tableau
  // version qui peut modifie la class  via l'adresse de x ou y 
  R & operator[](int i) { if(i==0) return x; else if (i==1) return y; else {assert(0);exit(1);} ;}
  // version qui retourne une reference const qui ne modifie pas la class
  const R & operator[](int i) const { if(i==0) return x; else if (i==1) return y; else {assert(0);exit(1);} ;}

  // les operators  tableaux
  // version qui peut modifie la class  via l'adresse de x ou y 
  R & operator()(int i) { if(i==1) return x; else if (i==2) return y; else {assert(0);abort();} ;}
  // version qui retourne une reference const qui ne modifie pas la class
  const R & operator()(int i) const { if(i==1) return x; else if (i==2) return y; else {assert(0);abort();} ;}
  R norme() const { return std::sqrt(x*x+y*y);}
  //  un variable globale dans la classe donc de nom R2::d 
  static const int  d=2;  // utilisation:  R2::d  (n'ajoute pas de memoire)
};
inline R det(const R2 & A,const R2 & B,const R2 & C)
	{ return R2(A,B)^R2(A,C);}
inline std::ostream& operator <<(std::ostream& f, const R2 & P )
       { f << P.x << ' ' << P.y   ; return f; }
inline std::istream& operator >>(std::istream& f,  R2 & P)
       { f >>  P.x >>  P.y  ; return f; }

inline R2 operator*(R c,const R2 & P) {return P*c;} 
inline R2 perp(const R2 & P) { return R2(-P.y,P.x) ; }
//inline R2 Perp(const R2 & P) { return P.perp(); }  // autre ecriture  de la fonction perp
