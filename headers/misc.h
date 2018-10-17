double Sign(double input){
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

double Ve(double rho, double Ye){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;
}
double Vmu(double rho, double Ye){ return 0.;}

void getP(const double r,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > >& U0,
	  const vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
	  vector<MATRIX<complex<double>,NF,NF> >& pmatrixm0matter){

  for(int i=0;i<=NE-1;i++){
    pmatrixm0matter[i] = Adjoint(U0[matter][i])
      * (fmatrixf[matter][i] - Conjugate(fmatrixf[antimatter][i]) )
      * U0[matter][i];
    pmatrixm0matter[i] *= M_SQRT2*cgs::constants::GF /* erg cm^3*/
      * 4.*M_PI*nu[i]*nu[i]*dnu[i]/*Hz^3*/ / pow(cgs::constants::c,3)/*cm^3 Hz^3*/;
  }
}

//===//
// B //
//===//
MATRIX<complex<double>,NF,NF> B(vector<double> y){
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
// K //
//===//
void K(double r,
       double dr,
       double rho, double Ye,
       const vector<MATRIX<complex<double>,NF,NF> >& pmatrixm0matter,
       vector<vector<vector<vector<double> > > > &Y,
       vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > &C0,
       vector<vector<vector<vector<double> > > > &A0,
       vector<vector<vector<vector<double> > > > &K){

  MATRIX<complex<double>,NF,NF> VfSI,VfSIbar;  // self-interaction potential
  vector<MATRIX<complex<double>,NF,NF> > VfSIE(NE); // contribution to self-interaction potential from each energy
  MATRIX<complex<double>,NF,NF> VfMSW,VfMSWbar;
  MATRIX<complex<double>,NF,NF> Hf,Hfbar, UU,UUbar;
  vector<double> kk,kkbar, dkk,dkkbar;
  vector<MATRIX<complex<double>,NF,NF> > CC;
  vector<vector<double> > AA;
  MATRIX<complex<double>,NF,NF> BB,BBbar;
  MATRIX<complex<double>,NF,NF> Sfm,Sfmbar;
  vector<vector<MATRIX<complex<double>,NF,NF> > > 
    Sa(NE,vector<MATRIX<complex<double>,NF,NF> >(NS)), 
    Sabar(NE,vector<MATRIX<complex<double>,NF,NF> >(NS));
  vector<MATRIX<complex<double>,NF,NF> > UWBW(NE);
  vector<MATRIX<complex<double>,NF,NF> > UWBWbar(NE);
  double drhodr,dYedr;
  MATRIX<double,3,4> JI;
  int i;
  MATRIX<complex<double>,NF,NF> Ha;
  MATRIX<complex<double>,NF,NF> HB;
  vector<double> dvdr(4);
  // *************
  drhodr=0;
  dYedr=0;
  VfMSW[e][e]=Ve(rho,Ye);
  VfMSW[mu][mu]=Vmu(rho,Ye);
  VfMSWbar=-Conjugate(VfMSW);

#pragma omp parallel for schedule(auto) private(Hf,Hfbar,UU,UUbar,kk,kkbar,dkk,dkkbar,AA,CC,BB,BBbar,Sfm,Sfmbar)
  for(i=0;i<=NE-1;i++){
    Hf  = HfV[matter][i]+VfMSW;
    kk  = k(Hf);
    dkk = deltak(Hf);
    CC  = CofactorMatrices(Hf,kk);
    AA  = MixingMatrixFactors(CC,C0[matter][i],A0[matter][i]);
    UU  = U(dkk,CC,AA);
    BB  = B(Y[matter][i][msw]);
    Sa[i][si] = B(Y[matter][i][si]);
    UWBW[i] = UU * W(Y[matter][i][msw]) * BB * W(Y[matter][i][si]);
    
    Hfbar = HfV[antimatter][i] + VfMSWbar;
    kkbar = kbar(Hfbar);
    dkkbar = deltakbar(Hfbar);
    CC = CofactorMatrices(Hfbar,kkbar);
    AA = MixingMatrixFactors(CC,C0[antimatter][i],A0[antimatter][i]);
    UUbar = Conjugate(U(dkkbar,CC,AA));
    BBbar = B(Y[antimatter][i][msw]);
    Sabar[i][si] = B(Y[antimatter][i][si]);
    UWBWbar[i] = UUbar * W(Y[antimatter][i][msw]) *BBbar * W(Y[antimatter][i][si]);
    
    // ****************
    // Matter section *
    // ****************
    for(int j=0;j<=3;j++)
      K[matter][i][msw][j]=0.;
    
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
    Sfm    = UWBW   [i]*Sa   [i][si];
    Sfmbar = UWBWbar[i]*Sabar[i][si];
    VfSIE[i] =     Sfm   *pmatrixm0matter[i]*Adjoint(Sfm   );
  }//end for loop over i

  // ************************************
  // compute self-interaction potential *
  // ************************************
  for(i=0;i<=NE-1;i++){
    VfSI[e ][e ]+=VfSIE[i][e ][e ];
    VfSI[e ][mu]+=VfSIE[i][e ][mu];
    VfSI[mu][e ]+=VfSIE[i][mu][e ];
    VfSI[mu][mu]+=VfSIE[i][mu][mu];
  }
  VfSIbar=-Conjugate(VfSI);

  // *********************
  // SI part of solution *
  // *********************

#pragma omp parallel for schedule(auto) private(JI) firstprivate(Ha,HB,dvdr)
  for(i=0;i<=NE-1;i++){
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

}// end of K function


//===========//
// Outputvsr //
//===========//
void Outputvsr(ofstream &foutf, const double r, const vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf){
  foutf << r << "\t";
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
