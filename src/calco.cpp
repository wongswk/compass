#include <Rcpp.h>
#include <math.h> 
#include <iostream>
using namespace Rcpp;
using namespace std;

//' Calculate Cartesian coordinates for 1-4 bonded atoms
//' 
//' @description 
//' Given coordinates for the first three atoms forming a dihedral,
//' the dihedral angle, bond angle, and bond length, calculate the coordinates
//' of the fourth atom forming the dihedral.
//'
//' @param prev_atoms a 3x3 matrix of coordinates for atoms 1-2-3 in dihedral, listed by row
//' @param length bond length between atoms 3-4 in dihedral 
//' @param bAngle planar bond angle between atoms 2-3-4
//' @param tAngle dihedral angle formed by atoms 1-2-3-4
//' 
//' @return
//' Returns the vector of coordinates for the fourth atom in the dihedral
//' 
// [[Rcpp::export]]
NumericVector calCo(NumericMatrix prev_atoms, double length, double bAngle, double tAngle)
{
  double pi=3.1415926; 
  NumericVector su(3);
  NumericVector u2(3);
  NumericVector u3(3);
  NumericVector _SvdV(3);
  NumericVector _u(3);
  NumericVector _v(3);
  NumericVector _w(3);
  NumericVector cur_res(3);
  NumericVector last_res(3);
  NumericVector bfl_res(3);
  double d;
  double dis2;
  NumericVector n_res(3);
  
  cur_res  = prev_atoms(2, _);  
  last_res = prev_atoms(1, _);
  bfl_res  = prev_atoms(0, _);
  
  
  bAngle=bAngle*(pi/180);
  tAngle=tAngle*(pi/180);

  
  _SvdV    = bfl_res-last_res; 
  su       = cur_res-last_res; 
  d        = sqrt((su[0]*su[0])+(su[1]*su[1])+(su[2]*su[2])); 
  _u       = su/d; 
  d        =(_SvdV[0]*_u[0])+(_SvdV[1]*_u[1])+(_SvdV[2]*_u[2]); 
  u3       =NumericVector::create(_u[0]*d,_u[1]*d, _u[2]*d); 
  dis2     = sqrt(pow((_SvdV[0]-u3[0]),2)+pow((_SvdV[1]-u3[1]),2)+pow((_SvdV[2]-u3[2]), 2));
  _v       = (_SvdV-u3)/dis2;
  _w       =NumericVector::create(_u[1]*_v[2]-_u[2]*_v[1], _u[2]*_v[0]-_u[0]*_v[2], _u[0]*_v[1]-_u[1]*_v[0]);
    
  n_res    = cur_res + _u*length*cos(PI-bAngle)
  + _v*length*sin(PI-bAngle)*cos(tAngle)
  + _w*length*sin(PI-bAngle)*sin(tAngle);
  
  return n_res;
  
}



 

 


