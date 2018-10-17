//=========//
// UpdateC //
//=========//
vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > UpdateC(double r,double Ye){
  vector<MATRIX<complex<double>,NF,NF> > VfMSW(NM);
  MATRIX<complex<double>,NF,NF> Hf,Hfbar;
  vector<double> kk,kkbar;

  vector<vector<vector<MATRIX<complex<double>,NF,NF> > > >
    CC(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF)));

  VfMSW[matter][e][e] = Ve(rho,Ye);
  VfMSW[matter][mu][mu] = Vmu(rho,Ye);

  VfMSW[antimatter] = -VfMSW[matter];

  int i;
#pragma omp parallel for schedule(auto) private(Hf,Hfbar,kk,kkbar)
  for(i=0;i<=NE-1;i++){
    Hf = HfV[matter][i] + VfMSW[matter];
    kk = k(Hf);
    CC[matter][i] = CofactorMatrices(Hf,kk);
    
    Hfbar = HfV[antimatter][i]+VfMSW[antimatter];
    kkbar = kbar(Hfbar);
    CC[antimatter][i] = CofactorMatrices(Hfbar,kkbar);
  }

  return CC;
}

//=========//
// UpdateA //
//=========//
vector<vector<vector<vector<double> > > >
UpdateA(vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C,
	vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C0,
	vector<vector<vector<vector<double> > > > A0){

  vector<vector<vector<vector<double> > > >
    A(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF))));
  
  int i;
#pragma omp parallel for schedule(auto)
  for(i=0;i<=NE-1;i++){
    A[matter][i]=MixingMatrixFactors(C[matter][i],C0[matter][i],A0[matter][i]);
    A[antimatter][i]=MixingMatrixFactors(C[antimatter][i],C0[antimatter][i],A0[antimatter][i]);
  }
  
  return A;
}
