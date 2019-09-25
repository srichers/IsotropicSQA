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
// U //
//===//
MATRIX<complex<double>,NF,NF>
  U(const array<double,1>& dk,
    const array<MATRIX<complex<double>,NF,NF>,NF> &C,
    const array<array<double,NF>,NF>& A){

  MATRIX<complex<double>,NF,NF> u;
  double d;
  array<double,NF> r2;
  
  for(int j=0;j<=NF-1;j++){
    if(j==0) d = -dk[0]; // first column
    if(j==1) d =  dk[0]; // second column
    
    r2[e]  = real(C[j][e ][e ])*d;
    r2[mu] = real(C[j][mu][mu])*d;
    
    if(r2[e]>=r2[mu]){
      u[e ][j] = A[j][e] * C[j][e][e ] / sqrt(r2[e]);
      u[mu][j] = A[j][e] * C[j][e][mu] / sqrt(r2[e]);
    }
    if(r2[mu]>=r2[e]){
      u[e ][j] = A[j][mu] * C[j][mu][e ] / sqrt(r2[mu]);
      u[mu][j] = A[j][mu] * C[j][mu][mu] / sqrt(r2[mu]);
    }
  }
  
  return u;
}

//===//
// C //
//===//
template<flavour a,flavour b> complex<double> C(const MATRIX<complex<double>,NF,NF>& H,const double k);
template<> complex<double> C<e,e>(const MATRIX<complex<double>,NF,NF>& H,const double k){
  return H[mu][mu]-k;}
template<> complex<double> C<e,mu>(const MATRIX<complex<double>,NF,NF>& H,const double k){
  return -H[mu][e];}
template<> complex<double> C<mu,e>(const MATRIX<complex<double>,NF,NF>& H,const double k){
  return -H[e][mu];}
template<> complex<double> C<mu,mu>(const MATRIX<complex<double>,NF,NF>& H,const double k){
  return H[e][e]-k;}
  
//======//
// dCdr //
//======//
template<flavour a,flavour b> complex<double> dCdr(const MATRIX<complex<double>,NF,NF>& dHdr,const double dkdr);
template<> complex<double> dCdr<e,e>(const MATRIX<complex<double>,NF,NF>& dHdr,const double dkdr){
  return dHdr[mu][mu]-dkdr;}
template<> complex<double> dCdr<e,mu>(const MATRIX<complex<double>,NF,NF>& dHdr,const double dkdr){
  return -dHdr[mu][e];}
template<> complex<double> dCdr<mu,e>(const MATRIX<complex<double>,NF,NF>& dHdr,const double dkdr){
  return -dHdr[e][mu];}
template<> complex<double> dCdr<mu,mu>(const MATRIX<complex<double>,NF,NF>& dHdr,const double dkdr){
  return dHdr[e][e]-dkdr;}

//=====================//
// MixingMatrixFactors //
//=====================//
array<array<double,NF>,NF>
MixingMatrixFactors(const array<MATRIX<complex<double>,NF,NF>,NF> &C,
		    const array<MATRIX<complex<double>,NF,NF>,NF> &C0,
		    const array<array<double,NF>,NF> &A0){

  array<array<double,NF>,NF> A(A0);
  for(int j=0;j<=NF-1;j++){
    if(real(C[j][e][mu]*C0[j][e][mu])<0.) A[j][e ] *= -1.;
    if(real(C[j][mu][e]*C0[j][mu][e])<0.) A[j][mu] *= -1.;
  }
  
  return A;
}

//===================//
// Vacuum Potentials //
//===================//
array<array<double,NF>,NE> set_kV(const array<double,NE>& E){
  assert(NF==2);
  array<array<double,NF>,NE> kV;
  for(int i=0;i<NE;i++){
    kV[i][0] = m1*m1           * cgs::constants::c4 /2./E[i];
    kV[i][1] = kV[i][0] + dm21 * cgs::constants::c4 /2./E[i];
  }

  // determine eigenvalue ordering
  if(kV[0][1]>kV[0][0])
    cout<<"\n\nNormal hierarchy" << endl;
  else{
    if(kV[0][1]<kV[0][0])
      cout<<"\n\nInverted hierarchy" << endl;
    else{
      cout<<endl<<endl<<"Neither normal or Inverted"<<endl;
      abort();
    }
  }
  
  return kV;
}

//====================//
// VACUUM HAMILTONIAN //
//====================//
array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>
  Evaluate_HfV(const array<array<double,NF>,NE>& kV,
	       const array<MATRIX<complex<double>,NF,NF>,NM>& UV){

  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV;
  MATRIX<complex<double>,NF,NF> KV;
  
  for(int i=0;i<=NE-1;i++){
    for(int j=0;j<=NF-1;j++) KV[j][j]=kV[i][j];
    HfV[matter][i]=UV[matter]*KV*Adjoint(UV[matter]);
    HfV[antimatter][i]=Conjugate(HfV[matter][i]);
  }
  return HfV;
}

//==================//
// CofactorMatrices //
//==================//
array<MATRIX<complex<double>,NF,NF>,NF>
  CofactorMatrices(const MATRIX<complex<double>,NF,NF>& H,
		   const array<double,NF>& k){
  
  array<MATRIX<complex<double>,NF,NF>,NF> CC;
  for(int j=0;j<=NF-1;j++){
    CC[j][e ][e ] = (H[mu][mu]-k[j]);
    CC[j][e ][mu] = -H[mu][e];
    CC[j][mu][e ] = conj(CC[j][e][mu]);
    CC[j][mu][mu] = (H[e][e]-k[j]);
  }
  
  return CC;
}

//=============//
// Evaluate_UV //
//=============//
// vaccum mixing matrices
array<MATRIX<complex<double>,NF,NF>,NM> Evaluate_UV(void){
  array<MATRIX<complex<double>,NF,NF>,NM> UV;
  UV[matter][0][0] = c12V * exp(-I*alphaV[0]                  );
  UV[matter][0][1] = s12V * exp(-I*alphaV[1]                  );
  UV[matter][1][0] =-s12V * exp( I* betaV[0])*exp(-I*alphaV[0]);
  UV[matter][1][1] = c12V * exp( I* betaV[0])*exp(-I*alphaV[1]);
  
  UV[antimatter]=Conjugate(UV[matter]);
  return UV;
}

//=============//
// Evaluate_CV //
//=============//
// cofactor matrices in vacuum
array<array<MATRIX<complex<double>,NF,NF>,NF>,NE>
  Evaluate_CV(const array<array<double,NF>,NE>& kV,
	      const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& HfV){
  array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV;
  for(int i=0;i<=NE-1;i++){
    CV[i][0][e][mu]=C<e,mu>(HfV[matter][i],kV[i][0]); CV[i][1][e][mu]=C<e,mu>(HfV[matter][i],kV[i][1]);
    CV[i][0][mu][e]=C<mu,e>(HfV[matter][i],kV[i][0]); CV[i][1][mu][e]=C<mu,e>(HfV[matter][i],kV[i][1]);
  }
  return CV;
}

//=============//
// Evaluate_AV //
//=============//
// mixing matrix element prefactors in vacuum
array<array<array<double,NF>,NF>,NE>
  Evaluate_AV(const array<array<double,NF>,NE>& kV,
	      const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& HfV,
	      const array<MATRIX<complex<double>,NF,NF>,NM>& UV){

  array<array<array<double,NF>,NF>,NE> AV;
  
  double Delta;
  for(int i=0;i<=NE-1;i++){
    for(int j=0;j<=NF-1;j++){
      if(j==0) Delta=(kV[i][1]-kV[i][0]);
      if(j==1) Delta=(kV[i][0]-kV[i][1]);
      
      double re2=Delta*real(C<e,e>(HfV[matter][i],kV[i][j]));
      double rmu2=Delta*real(C<mu,mu>(HfV[matter][i],kV[i][j]));
      
      if(norm(UV[matter][e][j])>norm(UV[matter][mu][j]) ){
	AV[i][j][e ]=real( UV[matter][e][j] * sqrt(re2 ) / C<e ,e>(HfV[matter][i],kV[i][j]) );
	AV[i][j][mu]=real( UV[matter][e][j] * sqrt(rmu2) / C<mu,e>(HfV[matter][i],kV[i][j]) );
      }
      if(norm(UV[matter][mu][j])>norm(UV[matter][e][j]) ){
	AV[i][j][e ]=real( UV[matter][mu][j]*sqrt(re2 ) / C<e ,mu>(HfV[matter][i],kV[i][j]) );
	AV[i][j][mu]=real( UV[matter][mu][j]*sqrt(rmu2) / C<mu,mu>(HfV[matter][i],kV[i][j]) );
      }
    }
  }
  return AV;
}

