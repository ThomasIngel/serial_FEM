#include "hpc.h"
#include <math.h>
void mesh_Neumann(double p1[2], double p2[2], 
                  index typ, double (*fc)(double *, index), double m[3])  
{
  int i;
  double x1[2], x2[2], h2 = 0.0, fac;
  
  for (i = 0 ; i < 2 ; i++ ){
     h2 += (p2[i] - p1[i]) * (p2[i] - p1[i]) ;
     x1[i] = (3.0 * p1[i] + p2[i]) / 4.0;
     x2[i] = (3.0 * p2[i] + p1[i]) / 4.0;
  }
  h2 = sqrt(h2)/4.0;
  fac = fc(x1,typ)*h2;
  m[0] = fac; m[2] = fac;
  fac = fc(x2,typ)*h2;
  m[1] = fac; m[2] += fac;
}

void mesh_vol_elem(double p1[2], double p2[2], double p3[2], 
              index typ, double (*fc)(double *, index), double m[3])  
{
  int i;
  double mid[2], d[2][2], fac;
  
  for (i = 0 ; i < 2 ; i++ ){
     d[0][i] = p1[i]-p3[i];
     d[1][i] = p2[i]-p1[i];
     mid[i] = (p1[i] + p2[i] + p3[i])/2.0;
  }
  fac = fc(mid,typ)/6.0*(d[0][0]*d[1][1]-d[1][0]*d[0][1]);
  for ( i = 0 ; i < 3 ; i++ ){
    m[i] = fac; 
  }
}

/* create rhs - vector for uniform refined macro element
 * ordering w.r.t. [ p1, p2, p3, m1=(p1+p2)/2, m2=(p2+p3)/2, m3=(p3+p1)/2] */
void mesh_vol(double p1[2], double p2[2], double p3[2],
                    index  typ, double (*fc)(double *, index), double b[6]) 
{
  int j;
  double m[3][2], c[3], fac, s[3];
  
  for (j = 0 ; j < 2 ; j++ )
  {
     m[0][j] = ( p1[j] + p2[j] ) / 2.0 ;
     m[1][j] = ( p2[j] + p3[j] ) / 2.0 ;
     m[2][j] = ( p3[j] + p1[j] ) / 2.0 ; 
  }
  mesh_vol_elem(p1, m[0], m[2], typ, fc, s);
  b[0]  = s[0]; b[3]  = s[1]; b[5]  = s[2];
  mesh_vol_elem(p2, m[1], m[0], typ, fc, s);
  b[1]  = s[0]; b[4]  = s[1]; b[3] += s[2];
  mesh_vol_elem(p3, m[2], m[1], typ, fc, s);
  b[2]  = s[0]; b[5] += s[1]; b[4] += s[2];
  mesh_vol_elem(m[1], m[2], m[0], typ, fc, s);
  b[4] += s[0]; b[5] += s[1]; b[3] += s[2];
}

/* compute right hand side using midpoint rule */
void mesh_buildRhs(const mesh *M, double *b, double (*fV)(double *, index), 
                   double (*fN)(double *, index))
{  
  index j, k, nC, nT, nB, *Elem, *Bdry, ind[6] ;
  double *Coord, bx[6] ;
  
  nT = M->nelem; nC = M->ncoord; nB = M->nbdry; 
  Coord = M->coord; Elem = M->elem; Bdry = M->bdry; 
  
  for ( k = 0 ; k < nT; k++)
  {
    for (j = 0 ; j < 3 ; j++) ind[j] =      Elem[7*k+j];                                   
    for (j = 3 ; j < 6 ; j++) ind[j] = nC + Elem[7*k+j];  
    mesh_vol(Coord+2*ind[0],Coord+2*ind[1],Coord+2*ind[2],Elem[7*k+6],fV,bx);
    for (j = 0 ; j < 6 ; j++) b[ind[j]] += bx[j] ;  
  }                    

  for ( k = 0 ; k < nB; k++)
  {
    if (Bdry[4*k+3]) 
    {
      for (j = 0 ; j < 2 ; j++) ind[j] =      Bdry[4*k+j];                                   
      ind[2] = nC + Bdry[4*k+2];  
      mesh_Neumann(Coord+2*ind[0],Coord+2*ind[1],Bdry[4*k+3],fN,bx);
      for (j = 0 ; j < 3 ; j++) b[ind[j]] += bx[j];
    }
  }                    
}
