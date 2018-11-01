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
// U //
//===//
MATRIX<complex<double>,NF,NF>
  U(const vector<double>& dk,
    const vector<MATRIX<complex<double>,NF,NF> > &C,
    const vector<vector<double> >& A){

  MATRIX<complex<double>,NF,NF> u;
  double d;
  vector<double> r2(NF);
  
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
vector<vector<double> >
MixingMatrixFactors(const vector<MATRIX<complex<double>,NF,NF> > &C,
		    const vector<MATRIX<complex<double>,NF,NF> > &C0,
		    const vector<vector<double> > &A0){

  vector<vector<double> > A(A0);
  for(int j=0;j<=NF-1;j++){
    if(real(C[j][e][mu]*C0[j][e][mu])<0.) A[j][e ] *= -1.;
    if(real(C[j][mu][e]*C0[j][mu][e])<0.) A[j][mu] *= -1.;
  }
  
  return A;
}

//===================//
// Vacuum Potentials //
//===================//
vector<vector<double> > set_kV(const vector<double>& E){
  assert(NF==2);
  const unsigned NE = E.size();
  vector<vector<double> > kV(NE,vector<double>(NF));
  for(int i=0;i<NE;i++){
    kV[i][0] = m1*m1             * cgs::constants::c4 /2./E[i];
    kV[i][1] = (kV[i][0] + dm21) * cgs::constants::c4 /2./E[i];
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
vector<vector<MATRIX<complex<double>,NF,NF> > >
  Evaluate_HfV(const vector<vector<double> >& kV,
	       const vector<MATRIX<complex<double>,NF,NF> >& UV){

  const unsigned NE = kV.size();
  vector<vector<MATRIX<complex<double>,NF,NF> > > HfV(NM);
  HfV[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
  HfV[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
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
vector<MATRIX<complex<double>,NF,NF> >
  CofactorMatrices(const MATRIX<complex<double>,NF,NF>& H,
		   const vector<double>& k){
  
  vector<MATRIX<complex<double>,NF,NF> > CC(NF);
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
vector<MATRIX<complex<double>,NF,NF> > Evaluate_UV(void){
  vector<MATRIX<complex<double>,NF,NF> > UV(NM);
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
vector<vector<MATRIX<complex<double>,NF,NF> > >
  Evaluate_CV(const vector<vector<double> >& kV,
	      const vector<vector<MATRIX<complex<double>,NF,NF> > >& HfV){
  const unsigned NE = kV.size();
  vector<vector<MATRIX<complex<double>,NF,NF> > > CV =
    vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF));
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
vector<vector<vector<double> > >
Evaluate_AV(const vector<vector<double> >& kV,
	    const vector<vector<MATRIX<complex<double>,NF,NF> > >& HfV,
	    const vector<MATRIX<complex<double>,NF,NF> >& UV){

  const unsigned NE = kV.size();
  vector<vector<vector<double> > > AV(NE,vector<vector<double> >(NF,vector<double>(NF)));
  
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
