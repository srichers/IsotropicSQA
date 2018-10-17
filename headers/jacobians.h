#if !defined(_jacobians)
#define _jacobians

//==========//
// JInverse //
//==========//
MATRIX<double,3,4> JInverse(vector<double> y){
  MATRIX<double,3,4> j;
  
  double cPsi1=cos(y[0]), sPsi1=sin(y[0]);
  double cPsi2=cos(y[1]), sPsi2=sin(y[1]);
  double cPsi3=cos(y[2]), sPsi3=sin(y[2]);
  
  j[0][0]=-sPsi1;       j[0][1]=cPsi1*cPsi2;       j[0][2]=cPsi1*sPsi2*cPsi3; j[0][3]=cPsi1*sPsi2*sPsi3;
  j[1][1]=-sPsi2/sPsi1; j[1][2]=cPsi2*cPsi3/sPsi1; j[1][3]=cPsi2*sPsi3/sPsi1;
  sPsi1*=sPsi2;
  j[2][2]=-sPsi3/sPsi1; j[2][3]=cPsi3/sPsi1;
  
  return j;
}

#endif
