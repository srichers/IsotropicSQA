//===//
// U //
//===//
MATRIX<complex<double>,NF,NF>
  U(vector<double> dk,
    vector<MATRIX<complex<double>,NF,NF> > &C,
    vector<vector<double> > A){

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
template<flavour a,flavour b> complex<double> C(MATRIX<complex<double>,NF,NF> H,double k);
template<> complex<double> C<e,e>(MATRIX<complex<double>,NF,NF> H,double k){
  return H[mu][mu]-k;}
template<> complex<double> C<e,mu>(MATRIX<complex<double>,NF,NF> H,double k){
  return -H[mu][e];}
template<> complex<double> C<mu,e>(MATRIX<complex<double>,NF,NF> H,double k){
  return -H[e][mu];}
template<> complex<double> C<mu,mu>(MATRIX<complex<double>,NF,NF> H,double k){
  return H[e][e]-k;}
  
//======//
// dCdr //
//======//
template<flavour a,flavour b> complex<double> dCdr(MATRIX<complex<double>,NF,NF> dHdr,double dkdr);
template<> complex<double> dCdr<e,e>(MATRIX<complex<double>,NF,NF> dHdr,double dkdr){
  return dHdr[mu][mu]-dkdr;}
template<> complex<double> dCdr<e,mu>(MATRIX<complex<double>,NF,NF> dHdr,double dkdr){
  return -dHdr[mu][e];}
template<> complex<double> dCdr<mu,e>(MATRIX<complex<double>,NF,NF> dHdr,double dkdr){
  return -dHdr[e][mu];}
template<> complex<double> dCdr<mu,mu>(MATRIX<complex<double>,NF,NF> dHdr,double dkdr){
  return dHdr[e][e]-dkdr;}

//=====================//
// MixingMatrixFactors //
//=====================//
vector<vector<double> >
MixingMatrixFactors(vector<MATRIX<complex<double>,NF,NF> > &C,
		    vector<MATRIX<complex<double>,NF,NF> > &C0,
		    vector<vector<double> > A0){

  vector<vector<double> > A(A0);
  for(int j=0;j<=NF-1;j++){
    if(real(C[j][e][mu]*C0[j][e][mu])<0.) A[j][e ] *= -1.;
    if(real(C[j][mu][e]*C0[j][mu][e])<0.) A[j][mu] *= -1.;
  }
  
  return A;
}

//==================//
// CofactorMatrices //
//==================//
vector<MATRIX<complex<double>,NF,NF> >
  CofactorMatrices(MATRIX<complex<double>,NF,NF> H,
		   vector<double> k){
  
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
void Evaluate_UV(void){
  UV[matter][0][0] = c12V * exp(-I*alphaV[0]                  );
  UV[matter][0][1] = s12V * exp(-I*alphaV[1]                  );
  UV[matter][1][0] =-s12V * exp( I* betaV[0])*exp(-I*alphaV[0]);
  UV[matter][1][1] = c12V * exp( I* betaV[0])*exp(-I*alphaV[1]);
  
  UV[antimatter]=Conjugate(UV[matter]);
}

//=============//
// Evaluate_CV //
//=============//
void Evaluate_CV(void){
  for(int i=0;i<=NE-1;i++){
    CV[i][0][e][mu]=C<e,mu>(HfV[matter][i],kV[i][0]); CV[i][1][e][mu]=C<e,mu>(HfV[matter][i],kV[i][1]);
    CV[i][0][mu][e]=C<mu,e>(HfV[matter][i],kV[i][0]); CV[i][1][mu][e]=C<mu,e>(HfV[matter][i],kV[i][1]);
  }
}

//=============//
// Evaluate_AV //
//=============//
void Evaluate_AV(void){
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
}
