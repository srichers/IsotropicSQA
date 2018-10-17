/*
//  Copyright (c) 2015, California Institute of Technology and the Regents
//  of the University of California, based on research sponsored by the
//  United States Department of Energy. All rights reserved.
//
//  This file is part of Sedonu.
//
//  Sedonu is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Neither the name of the California Institute of Technology (Caltech)
//  nor the University of California nor the names of its contributors 
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  Sedonu is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Sedonu.  If not, see <http://www.gnu.org/licenses/>.
//
*/

//===//
// D //
//===//
double D(MATRIX<complex<double>,NF,NF> Hf){
  return 4.*norm(Hf[e][mu]) +norm(Hf[e][e]-Hf[mu][mu]);
}

//=====//
// f,g //
//=====//
double f(double x){
  return 0.5*x*(1.-0.25*x*(1.-0.5*x*(1.-0.625*x)));
}
double g(double x){
  return 0.5*x*(1.-0.75*x*(1.-0.83333333*x*(1.-0.875*x)));
}

//====//
// k1 //
//====//
double k1(double T,double sqrtD){
  return 0.5*(T+a1*sqrtD);
}
double k1bar(double Tbar,double sqrtDbar){
  return 0.5*(Tbar+a1*sqrtDbar);
}
double asymptotick1(double Hee,double Hmm,double x,int s){
  return (1+a1*s)/2.*Hee +(1-a1*s)/2.*Hmm +a1*s/2.*(Hee-Hmm)*f(x);
}
double asymptotick1bar(double Hee,double Hmm,double x,int s){
  return (1+a1*s)/2.*Hee +(1-a1*s)/2.*Hmm +a1*s/2.*(Hee-Hmm)*f(x);
}

//====//
// k2 //
//====//
double k2(double T,double sqrtD){
  return 0.5*(T+a2*sqrtD);
}
double k2bar(double Tbar,double sqrtDbar){
  return 0.5*(Tbar+a2*sqrtDbar);
}
double asymptotick2(double Hee,double Hmm,double x,int s){
  return (1+a2*s)/2.*Hee +(1-a2*s)/2.*Hmm +a2*s/2.*(Hee-Hmm)*f(x);
}
double asymptotick2bar(double Hee,double Hmm,double x,int s){
  return (1+a2*s)/2.*Hee +(1-a2*s)/2.*Hmm +a2*s/2.*(Hee-Hmm)*f(x);
}

/* // ********************************************************************* */

//===//
// k //
//===//
vector<double> k(MATRIX<complex<double>,NF,NF> Hf){
  //if(4.*norm(Hf[e][mu])/norm(Hf[e][e]-Hf[mu][mu])<0.0745){ return asymptotick(Hf);}

  vector<double> k(NF);
  double t=real(Trace(Hf)), sqrtd=sqrt(D(Hf));
  
  k[0]=k1(t,sqrtd);
  k[1]=k2(t,sqrtd);
  return k;
}

vector<double> k(double t,double sqrtd){
  vector<double> k(NF);
  k[0]=k1(t,sqrtd);
  k[1]=k2(t,sqrtd);
  return k;
}

//======//
// kbar //
//======//
vector<double> kbar(MATRIX<complex<double>,NF,NF> Hf){
  //if(4.*norm(Hf[e][mu])/norm(Hf[e][e]-Hf[mu][mu])<0.0745){ return asymptotickbar(Hf);}
  
  vector<double> k(NF);
  double tbar=real(Trace(Hf)), sqrtdbar=sqrt(D(Hf));
  
  k[0]=k1bar(tbar,sqrtdbar);
  k[1]=k2bar(tbar,sqrtdbar);
  
  return k;
}

vector<double> kbar(double tbar,double sqrtdbar){
  vector<double> k(NF);
  k[0]=k1bar(tbar,sqrtdbar);
  k[1]=k2bar(tbar,sqrtdbar);
  return k;
}

//========//
// deltak //
//========//
vector<double> deltak(MATRIX<complex<double>,NF,NF> Hf){
  vector<double> dk(1);
  double d=D(Hf), sqrtd=sqrt(d);
  dk[0]=a1*sqrtd;
  return dk;
}

vector<double> deltak(double sqrtd){
  vector<double> dk(1);
  dk[0]=a1*sqrtd;
  return dk;
}

//===========//
// deltakbar //
//===========//
vector<double> deltakbar(MATRIX<complex<double>,NF,NF> Hfbar){
  vector<double> dk(1);
  double dbar=D(Hfbar), sqrtdbar=sqrt(dbar);
  dk[0]=a1*sqrtdbar;
  return dk;
}

vector<double> deltakbar(double sqrtdbar){
  vector<double> dk(1);
  dk[0]=a1*sqrtdbar;
  return dk;
}

