/*
//  Copyright (c) 2018, James Kneller and Sherwood Richers
//
//  This file is part of IsotropicSQA.
//
//  IsotropicSQA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  IsotropicSQA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with IsotropicSQA.  If not, see <http://www.gnu.org/licenses/>.
//
*/

//===//
// D //
//===//
double D(const MATRIX<complex<double>,NF,NF>& Hf){
  return 4.*norm(Hf[e][mu]) +norm(Hf[e][e]-Hf[mu][mu]);
}

//=====//
// f,g //
//=====//
inline double f(const double x){
  return 0.5*x*(1.-0.25*x*(1.-0.5*x*(1.-0.625*x)));
}
inline double g(const double x){
  return 0.5*x*(1.-0.75*x*(1.-0.83333333*x*(1.-0.875*x)));
}

//====//
// k1 //
//====//
inline double k1(const double T,const double sqrtD){
  return 0.5*(T+a1*sqrtD);
}
inline double k1bar(const double Tbar,const double sqrtDbar){
  return 0.5*(Tbar+a1*sqrtDbar);
}
inline double asymptotick1(const double Hee,const double Hmm,const double x,const int s){
  return (1+a1*s)/2.*Hee +(1-a1*s)/2.*Hmm +a1*s/2.*(Hee-Hmm)*f(x);
}
inline double asymptotick1bar(const double Hee,const double Hmm,const double x,const int s){
  return (1+a1*s)/2.*Hee +(1-a1*s)/2.*Hmm +a1*s/2.*(Hee-Hmm)*f(x);
}

//====//
// k2 //
//====//
inline double k2(const double T,const double sqrtD){
  return 0.5*(T+a2*sqrtD);
}
inline double k2bar(const double Tbar,const double sqrtDbar){
  return 0.5*(Tbar+a2*sqrtDbar);
}
inline double asymptotick2(const double Hee,const double Hmm,const double x,const int s){
  return (1+a2*s)/2.*Hee +(1-a2*s)/2.*Hmm +a2*s/2.*(Hee-Hmm)*f(x);
}
inline double asymptotick2bar(const double Hee,const double Hmm,const double x,const int s){
  return (1+a2*s)/2.*Hee +(1-a2*s)/2.*Hmm +a2*s/2.*(Hee-Hmm)*f(x);
}

/* // ********************************************************************* */

//===//
// k //
//===//
array<double,NF> k(const MATRIX<complex<double>,NF,NF>& Hf){
  array<double,NF> k;
  double t=real(Trace(Hf)), sqrtd=sqrt(D(Hf));
  
  k[0]=k1(t,sqrtd);
  k[1]=k2(t,sqrtd);
  return k;
}

array<double,NF> k(const double t,const double sqrtd){
  array<double,NF> k;
  k[0]=k1(t,sqrtd);
  k[1]=k2(t,sqrtd);
  return k;
}

//======//
// kbar //
//======//
array<double,NF> kbar(const MATRIX<complex<double>,NF,NF>& Hf){
  array<double,NF> k;
  double tbar=real(Trace(Hf)), sqrtdbar=sqrt(D(Hf));
  
  k[0]=k1bar(tbar,sqrtdbar);
  k[1]=k2bar(tbar,sqrtdbar);
  
  return k;
}

array<double,NF> kbar(const double tbar,const double sqrtdbar){
  array<double,NF> k;
  k[0]=k1bar(tbar,sqrtdbar);
  k[1]=k2bar(tbar,sqrtdbar);
  return k;
}

//========//
// deltak //
//========//
array<double,1> deltak(const MATRIX<complex<double>,NF,NF>& Hf){
  array<double,1> dk;
  double d=D(Hf), sqrtd=sqrt(d);
  dk[0]=a1*sqrtd;
  return dk;
}

array<double,1> deltak(const double sqrtd){
  array<double,1> dk;
  dk[0]=a1*sqrtd;
  return dk;
}

//===========//
// deltakbar //
//===========//
array<double,1> deltakbar(const MATRIX<complex<double>,NF,NF>& Hfbar){
  array<double,1> dk;
  double dbar=D(Hfbar), sqrtdbar=sqrt(dbar);
  dk[0]=a1*sqrtdbar;
  return dk;
}

array<double,1> deltakbar(const double sqrtdbar){
  array<double,1> dk;
  dk[0]=a1*sqrtdbar;
  return dk;
}

