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

template<typename T>
T get_parameter(ifstream& fin, const char* name){
  stringstream line;
  string linestring;
  T to_set;
  std::getline(fin, linestring);
  line = stringstream(linestring);
  line >> to_set;
  cout << "PARAMETER " << name << " = " << to_set << endl;
  return to_set;
}

inline double Sign(const double input){
  return input>0 ? 1 : -1;
}

void RungeKuttaCashKarpParameters(int &NRK,int &NOrder,const double* &A,const double** &B,const double* &C,const double* &D)
     { NRK=6; NOrder=5;

       static const double a[]={ 0., 1./5., 3./10., 3./5., 1., 7./8. };
       static const double b0[]={};
       static const double b1[]={ 1./5. };
       static const double b2[]={ 3./40.,9./40. };
       static const double b3[]={ 3./10.,-9./10.,6./5. };
       static const double b4[]={ -11./54.,5./2.,-70./27.,35./27. };
       static const double b5[]={ 1631./55296.,175./512.,575./13824.,44275./110592.,253./4096. };
       static const double* b[]={ b0,b1,b2,b3,b4,b5 };
       static const double c[]={ 37./378.,0.,250./621.,125./594.,0.,512./1771. };
       static const double d[]={ 2825./27648.,0.,18575./48384.,13525./55296.,277./14336.,1./4. };

       A=a; B=b; C=c; D=d;
      } 

inline double Ve(const double rho, const double Ye){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;
}
inline double Vmu(const double rho, const double Ye){ return 0.;}

void getP(const double r,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > >& U0,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf0,
	  const vector<double>& nu,
	  const vector<double>& dnu,
	  vector<MATRIX<complex<double>,NF,NF> >& pmatrixm0matter){

  for(int i=0; i<pmatrixm0matter.size(); i++){
    pmatrixm0matter[i] = Adjoint(U0[matter][i])
      * (fmatrixf0[matter][i] - Conjugate(fmatrixf0[antimatter][i]) )
      * U0[matter][i];
    pmatrixm0matter[i] *= M_SQRT2*cgs::constants::GF /* erg cm^3*/
      * 4.*M_PI*nu[i]*nu[i]*dnu[i]/*Hz^3*/ / pow(cgs::constants::c,3)/*cm^3 Hz^3*/;
  }
}

//===//
// B //
//===//
MATRIX<complex<double>,NF,NF> B(const vector<double>& y){
  MATRIX<complex<double>,NF,NF> s;
  double cPsi1=cos(y[0]),sPsi1=sin(y[0]), cPsi2=cos(y[1]),sPsi2=sin(y[1]), cPsi3=cos(y[2]),sPsi3=sin(y[2]);
  
  s[0][1] = cPsi1 + I*sPsi1*cPsi2;
  sPsi1 *= sPsi2;
  s[0][0] = sPsi1 * (cPsi3 + I*sPsi3);

  s[1][0] = -y[3]*conj(s[0][1]);
  s[1][1] =  y[3]*conj(s[0][0]);

  return s;
}

//===//
// W //
//===//
MATRIX<complex<double>,NF,NF> W(const vector<double>& Y){
  MATRIX<complex<double>,NF,NF> w;
  w[0][0]=exp(-I*M_2PI*Y[4]); w[1][1]=exp(-I*M_2PI*Y[5]);
  return w;
}

//===//
// K //
//===//
void K(const double r,
       const double dr,
       const double rho,
       const double Ye,
       const vector<MATRIX<complex<double>,NF,NF> >& pmatrixm0matter,
       const vector<vector<MATRIX<complex<double>,NF,NF> > > HfV,
       const vector<vector<vector<vector<double> > > > &Y,
       const vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > &C0,
       const vector<vector<vector<vector<double> > > > &A0,
       vector<vector<vector<vector<double> > > > &K){

  const unsigned NE = pmatrixm0matter.size();
  vector<MATRIX<complex<double>,NF,NF> > VfSIE(NE); // contribution to self-interaction potential from each energy
  vector<vector<MATRIX<complex<double>,NF,NF> > > 
    Sa(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)),
    Sabar(NE,vector<MATRIX<complex<double>,NF,NF> >(NS));
  vector<MATRIX<complex<double>,NF,NF> > UWBW(NE), UWBWbar(NE);
  MATRIX<complex<double>,NF,NF> VfSI,VfSIbar;  // self-interaction potential

  MATRIX<complex<double>,NF,NF> VfMSW, VfMSWbar;
  VfMSW[e][e]=Ve(rho,Ye);
  VfMSW[mu][mu]=Vmu(rho,Ye);
  VfMSWbar=-Conjugate(VfMSW);

#pragma omp parallel
  {
  #pragma omp for
  for(int i=0; i<NE; i++){
    MATRIX<complex<double>,NF,NF> Hf  = HfV[matter][i]+VfMSW;
    vector<double> kk  = k(Hf);
    vector<double> dkk = deltak(Hf);
    vector<MATRIX<complex<double>,NF,NF> > CC = CofactorMatrices(Hf,kk);
    vector<vector<double> > AA = MixingMatrixFactors(CC,C0[matter][i],A0[matter][i]);
    MATRIX<complex<double>,NF,NF> UU  = U(dkk,CC,AA);
    MATRIX<complex<double>,NF,NF> BB  = B(Y[matter][i][msw]);
    Sa[i][si] = B(Y[matter][i][si]);
    UWBW[i] = UU * W(Y[matter][i][msw]) * BB * W(Y[matter][i][si]);
    
    MATRIX<complex<double>,NF,NF> Hfbar = HfV[antimatter][i] + VfMSWbar;
    vector<double> kkbar = kbar(Hfbar);
    vector<double> dkkbar = deltakbar(Hfbar);
    MATRIX<complex<double>,NF,NF> UUbar = Conjugate(U(dkkbar,CC,AA));
    MATRIX<complex<double>,NF,NF> BBbar = B(Y[antimatter][i][msw]);
    Sabar[i][si] = B(Y[antimatter][i][si]);
    UWBWbar[i] = UUbar * W(Y[antimatter][i][msw]) *BBbar * W(Y[antimatter][i][si]);
    
    // ****************
    // Matter section *
    // ****************
    for(int j=0;j<=3;j++)
      K[matter][i][msw][j] = 0.;
    K[matter][i][msw][4] = kk[0]*dr/M_2PI/cgs::constants::hbarc;
    K[matter][i][msw][5] = kk[1]*dr/M_2PI/cgs::constants::hbarc;

    // ********************
    // Antimatter section *
    // ********************
    for(int j=0;j<=3;j++)
      K[antimatter][i][msw][j] = 0.;
    K[antimatter][i][msw][4] = kkbar[0]*dr/M_2PI/cgs::constants::hbarc;
    K[antimatter][i][msw][5] = kkbar[1]*dr/M_2PI/cgs::constants::hbarc;

    // *****************************************************************
    // contribution to the self-interaction potential from this energy *
    // *****************************************************************
    MATRIX<complex<double>,NF,NF> Sfm    = UWBW[i] * Sa[i][si];
    VfSIE[i] = Sfm * pmatrixm0matter[i]*Adjoint(Sfm);
    #pragma omp critical
    for(flavour f1=e; f1<=mu; f1++)
      for(flavour f2=e; f2<=mu; f2++)
	VfSI[f1][f2]+=VfSIE[i][f1][f2];
  }//end for loop over i

  #pragma omp single
  VfSIbar=-Conjugate(VfSI);

  // *********************
  // SI part of solution *
  // *********************
  #pragma omp for
  for(int i=0; i<NE; i++){
    MATRIX<double,3,4> JI;
    MATRIX<complex<double>,NF,NF> Ha, HB;
    vector<double> dvdr(4);

    //*********
    // Matter *
    //*********    
    Ha = Adjoint(UWBW[i])*VfSI*UWBW[i];
  
    K[matter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[matter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);
    
    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sa[i][si][1][1] );
    
    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);
    
    JI=JInverse(Y[matter][i][si]);
    
    for(int j=0;j<=2;j++){
      K[matter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[matter][i][si][j]+=JI[j][k]*dvdr[k];
      K[matter][i][si][j]*=dr;
    }
    
    K[matter][i][si][3]=0.;
    
    //*************
    // Antimatter *
    //*************
    Ha=Adjoint(UWBWbar[i])*VfSIbar*UWBWbar[i];

    K[antimatter][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
    K[antimatter][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);

    HB[0][0]=-I/cgs::constants::hbarc*( Ha[0][1]*Sabar[i][si][1][0] );
    HB[0][1]=-I/cgs::constants::hbarc*( Ha[0][1]*Sabar[i][si][1][1] );

    dvdr[0]=real(HB[0][1]);
    dvdr[1]=imag(HB[0][1]);
    dvdr[2]=real(HB[0][0]);
    dvdr[3]=imag(HB[0][0]);

    JI = JInverse(Y[antimatter][i][si]);

    for(int j=0;j<=2;j++){
      K[antimatter][i][si][j]=0.;
      for(int k=j;k<=3;k++) K[antimatter][i][si][j]+=JI[j][k]*dvdr[k];
      K[antimatter][i][si][j]*=dr;
    }

    K[antimatter][i][si][3]=0.;
  }
  }// omp parallel block

}// end of K function


//===========//
// Outputvsr //
//===========//
void Outputvsr(ofstream &foutf, const double r, const vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf){
  foutf << r << "\t";
  const unsigned NE = fmatrixf.size();
  for(int i=0; i<NE; i++)
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  foutf << real( fmatrixf[m][i][f1][f2] ) << "\t";
	  foutf << imag( fmatrixf[m][i][f1][f2] ) << "\t";
	}
  foutf << endl;
  foutf.flush();
}

void Hermitize(MATRIX<complex<double>,2,2>& M, const double accuracy){
  double trace = real(M[e][e] + M[mu][mu]);
  assert(trace>0);

  // matching off-diagonals
  double error = abs(M[e][mu] - conj(M[mu][e]));
  if(error/trace >= accuracy){
    //cout << M << endl;
    assert(error/trace < accuracy);
  }
  complex<double> tmp = 0.5 * (M[e][mu] + conj(M[mu][e]));
  M[mu][e] = tmp;
  M[e][mu] = conj(tmp);
  
  // real on-diagonals
  for(flavour f1=e; f1<=mu; f1++){
    error = abs(imag(M[f1][f1]));
    assert(error/trace < accuracy);
    M[f1][f1] = real(M[f1][f1]);
  }

  // density matrix probabilities
  double ad = abs(M[e][e ] * M[mu][mu]) / trace;
  double bc = abs(M[e][mu] * M[mu][e ]) / trace;
  double det = ad - bc;
  if(det < 0 && bc>0){
    assert( abs(det) < accuracy);
    M[e][mu] *= sqrt(ad/bc);
    M[mu][e] *= sqrt(ad/bc);
  }
}

void unitarize(MATRIX<complex<double>,2,2>& M, const double accuracy){
  // M = ( (a, b), (-e^Iphi b*, e^Iphi a*) )
  //   = ( (a, b), (c, d) )
  double a2 = real(M[e ][e ] * conj(M[e ][e ]) );
  double b2 = real(M[e ][mu] * conj(M[e ][mu]) );
  double c2 = real(M[mu][e ] * conj(M[mu][e ]) );
  double d2 = real(M[mu][mu] * conj(M[mu][mu]) );

  // aa* < 1
  if(a2 > 1.){
    assert( abs(a2-1.) < accuracy);
    M[e][e] /= sqrt(a2);
    a2 = real(M[e ][e ] * conj(M[e ][e ]) );
  }
  
  // aa* + bb* = 1
  assert( abs(a2 + b2 - 1.) < accuracy);
  M[e][mu] *= sqrt( max(0., 1.-a2) / b2 );
  b2 = real(M[e ][mu] * conj(M[e ][mu]) );
  
  // aa* = dd*
  assert( abs(a2-d2) < accuracy );
  M[mu][mu] *= sqrt(a2/d2);
  d2 = real(M[mu][mu] * conj(M[mu][mu]) );

  // bb* = cc*
  assert( abs(b2-c2) < accuracy );
  M[mu][e] *= sqrt(b2/c2);
  c2 = real(M[mu][e ] * conj(M[mu][e ]) );
  
  // det(M) = e^Iphi
  complex<double> newval;
  complex<double> eIphi = M[e][e]*M[mu][mu] - M[e][mu]*M[mu][e];
  double eIphiMag = sqrt(real( eIphi*conj(eIphi) ));
  assert( abs(eIphiMag - 1.) < accuracy);
  eIphi /= eIphiMag;
  newval = eIphi * conj(M[e][e ]);
  assert( abs(M[mu][mu]-newval) < accuracy);
  M[mu][mu] = newval;
  newval = -eIphi * conj(M[e][mu]);
  assert( abs(M[mu][e ]-newval) < accuracy);
  M[mu][e ] = newval;

  // crazy sanity check
  MATRIX<complex<double>,2,2> identity;
  identity = M*Adjoint(M);
  assert( abs( real(identity[e ][e ]*conj(identity[e ][e ])) - 1.) < accuracy);
  assert( abs( real(identity[mu][mu]*conj(identity[mu][mu])) - 1.) < accuracy);
  assert( abs( real(identity[e ][mu]*conj(identity[e ][mu]))     ) < accuracy);
  assert( abs( real(identity[mu][e ]*conj(identity[mu][e ]))     ) < accuracy);
}
