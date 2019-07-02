/*
   Provided by teacher
*/

#include "EF2d-base.hpp"
#include <fstream>

Mesh2d::Mesh2d(const char * filename)
{
  std::ifstream  f(filename);
  assert( f);
  int unused, I[4];  // I[4] 用来存储三角形的三个顶点编号
  f >> nv >> nt >> unused;
  assert( f.good());
  t = new Simplex[nt];  // 三角形
  v = new Vertex[nv];   // 顶点
  assert( t && v);
  double mes =0; 
  for(int i=0;i<nv;++i)   // 对顶点赋 ：坐标
    { 
      f >> v[i];    // v包含所有点的二维坐标 x 和 y   nv * R2
      assert( f.good());
    }

  for(int k=0;k<nt;++k)  // 对三角形赋值：三个顶点
    { 
      for(int i=0;i< 4; ++i)
        f >> I[i];
      assert( f.good());
      t[k].build(v,I,-1);
      mes += t[k].mes; 
    }
  std::cout<< " End read " << nv << " " << nt << " mes =" << mes << std::endl; 
}
