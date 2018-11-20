void evolve_oscillations(State& s, const double rmax, const double accuracy, const double increase, const int step_output){
  
  // vectors of energies and vacuum eigenvalues
  const vector<vector<double> > kV = set_kV(s.eas.E);
  const vector<MATRIX<complex<double>,NF,NF> > UV = Evaluate_UV();
  const vector<vector<MATRIX<complex<double>,NF,NF> > > HfV = Evaluate_HfV(kV,UV);
  const vector<vector<MATRIX<complex<double>,NF,NF> > > CV = Evaluate_CV(kV, HfV);
  const vector<array<array<double,NF>,NF> > AV = Evaluate_AV(kV,HfV,UV);
    
  // MSW potential matrix
  MATRIX<complex<double>,NF,NF> VfMSW0, Hf0;
  VfMSW0[e][e]=Ve(s.rho,s.Ye);
  VfMSW0[mu][mu]=Vmu(s.rho,s.Ye);
    
  // other matrices
  vector<vector< array<MATRIX<complex<double>,NF,NF>,NF> > > C0(NM); // cofactor matrices at initial point
  vector<vector< array<array<double,NF>,NF> > > A0(NM); // mixing matrix element prefactors at initial point
  vector<vector<MATRIX<complex<double>,NF,NF> > > U0(NM); // mixing angles to MSW basis at initial point
  for(int m=matter; m<=antimatter; m++){
    C0[m] = vector<array<MATRIX<complex<double>,NF,NF>,NF> >(s.eas.ng);
    A0[m] = vector<array<array<double,NF>,NF> >(s.eas.ng);
    U0[m] = vector<MATRIX<complex<double>,NF,NF> >(s.eas.ng);
  }

  for(int i=0;i<=s.eas.ng-1;i++){
    MATRIX<complex<double>,NF,NF> Hf0;
    array<double,NF> k0;
    array<double,1> deltak0;

    Hf0=HfV[matter][i]+VfMSW0;
    k0=k(Hf0);
    deltak0=deltak(Hf0);
    C0[matter][i]=CofactorMatrices(Hf0,k0);
    for(int j=0;j<=NF-1;j++){
      if(real(C0[matter][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	A0[matter][i][j][e]=-AV[i][j][e];
      else A0[matter][i][j][e]=AV[i][j][e];
      A0[matter][i][j][mu]=AV[i][j][mu];
    }
    U0[matter][i]=U(deltak0,C0[matter][i],A0[matter][i]);
      
    Hf0=HfV[antimatter][i]-VfMSW0;
    k0=kbar(Hf0);
    deltak0=deltakbar(Hf0);
    C0[antimatter][i]=CofactorMatrices(Hf0,k0);
    for(int j=0;j<=NF-1;j++){
      if(real(C0[antimatter][i][j][mu][e]*CV[i][j][mu][e]) < 0.)
	A0[antimatter][i][j][e]=-AV[i][j][e];
      else A0[antimatter][i][j][e]=AV[i][j][e];
      A0[antimatter][i][j][mu]=AV[i][j][mu];
    }
    U0[antimatter][i]=Conjugate(U(deltak0,C0[antimatter][i],A0[antimatter][i]));
  }
    
  //********************
  // evolved variables *
  //********************
  vector<vector<vector<vector<double> > > > 
    Y(NM,vector<vector<vector<double> > >(s.eas.ng,vector<vector<double> >(NS,vector<double>(NY))));
  for(state m=matter;m<=antimatter;m++)
    for(int i=0;i<=s.eas.ng-1;i++)
      Y[m][i] = YIdentity;
    
  // temporaries
  vector<vector<MATRIX<complex<double>,NF,NF> > > SSMSW(s.fmatrixf), SSSI(s.fmatrixf);
  vector<vector<MATRIX<complex<double>,NF,NF> > > pmatrixm0(s.fmatrixf);
  vector<vector<vector<vector<double> > > >  Ytmp(Y);
  vector<vector<vector<vector<vector<double> > > > > Ks(NRK);
  for(int j=0; j<NRK; j++) Ks[j] = Y;
  
  bool finish=false;
  int next_output = step_output>0 ? rand()%step_output+1 : -1;
  
  do{ // loop over radii
    s.counter++;
    getP(U0,s.fmatrixf,s.eas.nu,s.eas.dnu,pmatrixm0);
    
    double dr = s.dr_osc;
    bool repeat=false;
    do{
      double r = s.r;
      if(r+dr>=rmax){
	finish=true;
	dr = rmax-r;
      }
      else finish=false;
      
      // RK integration for oscillation
      for(int k=0;k<=NRK-1;k++){
	Ytmp=Y;
	for(int l=0;l<=k-1;l++)
	  for(int m=0; m<=1; m++) // 0=matter 1=antimatter
	    for(int i=0;i<=s.eas.ng-1;i++)
	      for(int x=0;x<=1;x++) // 0=msw 1=si
		for(int j=0;j<=NY-1;j++)
		  Ytmp[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];
	  
	K(dr,s.rho,s.Ye,pmatrixm0,HfV,Ytmp,C0,A0,Ks[k]);
      }
	  
      // increment all quantities from oscillation
      Ytmp = Y;
      double maxerror=0.;
      #pragma omp parallel for collapse(2) reduction(max:maxerror)
      for(int m=matter; m<=antimatter; m++){
	for(int i=0;i<=s.eas.ng-1;i++){
	  for(int x=msw; x<=si; x++){
	    for(int j=0;j<=NY-1;j++){
	      double Yerror = 0.;
	      for(int k=0;k<=NRK-1;k++){
		Ytmp[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		Yerror += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
	      }
	      maxerror = max( maxerror, fabs(Yerror) );
	    } // j
	  } // x
	} // i
      } // m
      
      //==============//
      // TIMESTEPPING //
      //==============//
      r += dr;
      if(maxerror>accuracy){
	dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	repeat=true;
      }
      else{
	dr *= increase;
	repeat = false;
	if(maxerror>0) dr *= min( 1.0, pow(accuracy/maxerror,1./max(1,NRKOrder))/increase );
	Y = Ytmp;
	s.r = r;
	if(!finish) s.dr_osc = dr;
      }
	
    }while(repeat==true); // end of RK section
  
    //========================//
    // RESETTING/ACCUMULATING //
    //========================//
#pragma omp parallel for collapse(2)
    for(int m=matter;m<=antimatter;m++){
      for(int i=0;i<=s.eas.ng-1;i++){
	
	// get oscillated f
	SSMSW[m][i] = W(Y[m][i][msw])*B(Y[m][i][msw]);
	SSSI [m][i] = W(Y[m][i][si ])*B(Y[m][i][si ]);
	
	// test that the S matrices are close to diagonal
	if(norm(SSMSW[m][i][0][0])+0.1<norm(SSMSW[m][i][0][1]) or
	   norm(SSSI [m][i][0][0])+0.1<norm(SSSI [m][i][0][1]) or
	   finish){
	  MATRIX<complex<double>,NF,NF> SThisStep = U0[m][i] * SSMSW[m][i]*SSSI[m][i] * Adjoint(U0[m][i]);
	  s.fmatrixf[m][i] = SThisStep * s.fmatrixf[m][i] * Adjoint(SThisStep);
	  Y[m][i] = YIdentity;
	}
	else{ // take modulo 2 pi of phase angles
	  Y[m][i][msw][2]=fmod(Y[m][i][msw][2],M_2PI);
	  Y[m][i][msw][4]=fmod(Y[m][i][msw][4],1.0);
	  Y[m][i][msw][5]=fmod(Y[m][i][msw][5],1.0);
	  
	  Y[m][i][si ][2]=fmod(Y[m][i][si ][2],M_2PI);
	  Y[m][i][si ][4]=fmod(Y[m][i][si ][4],1.0);
	  Y[m][i][si ][5]=fmod(Y[m][i][si ][5],1.0);
	}
      }
    }

    //========//
    // OUTPUT //
    //========//
    if(step_output>0 and (s.counter>=next_output or finish)){
      Outputvsr(s, 0);
      next_output = step_output>0 ? s.counter + rand()%step_output + 1 : -1;
      cout << next_output << endl;
    }
    
  } while(finish==false);
}


void evolve_interactions(State& s, const double rmax, const double accuracy, const double increase){

  vector<vector<MATRIX<complex<double>,NF,NF> > > ftmp;
  array<vector<vector<MATRIX<complex<double>,NF,NF> > >,NRK> dfdr;

  bool finish=false;
  s.counter++;
  do{ // loop over radii
    
    double dr = s.dr_int;
    bool repeat=false;
    do{
      double r = s.r;
      if(r+dr>=rmax){
	finish=true;
	dr = rmax-r;
      }
      else finish=false;

      // RK integration
      for(int k=0;k<=NRK-1;k++){
	ftmp=s.fmatrixf;
	for(int l=0;l<=k-1;l++)
	  for(int m=0; m<=1; m++) // 0=matter 1=antimatter
	    for(int i=0;i<=s.eas.ng-1;i++)
	      ftmp[m][i] += dfdr[l][m][i] * BB[k][l] * dr;
	
	dfdr[k] = my_interact(ftmp, s.rho, s.T, s.Ye, s.eas);
      }
      
      // increment all quantities from oscillation
      double maxerror=0;
      ftmp = s.fmatrixf;
      #pragma omp parallel for collapse(2) reduction(max:maxerror)
      for(int m=matter; m<=antimatter; m++){
	for(int i=0;i<=s.eas.ng-1;i++){
	  MATRIX<complex<double>,NF,NF> err, df;
	  for(int k=0;k<=NRK-1;k++){
	    df += dfdr[k][m][i] * CC[k] * dr;
	    err += dfdr[k][m][i] * (CC[k]-DD[k]) * dr;
	  }
	  ftmp[m][i] += df;
	  
	  double trace = norm(Trace(s.fmatrixf[m][i]));
	  if(trace<=0){
	    cout << i << " " << m << " " << trace << endl;
	    exit(1);
	  }
	  for(flavour f1=e; f1<=mu; f1++){
	    for(flavour f2=e; f2<=mu; f2++){
	      maxerror = max(maxerror, norm(err[f1][f2])/norm(s.fmatrixf[m][i][f1][f2]));//IsospinL(s.fmatrixf[m][i]));
	    }
	  }
	} // i
      } // m
      
      //==============//
      // TIMESTEPPING //
      //==============//
      r+=dr;
      if(maxerror>accuracy){
	dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	repeat=true;
      }
      else{
	dr *= increase;
	repeat = false;
	if(maxerror>0) dr *= min( 1.0, pow(accuracy/maxerror,1./max(1,NRKOrder))/increase );
	s.r = r;
	s.fmatrixf = ftmp;
	if(!finish) s.dr_int = dr;
      }

    }while(repeat==true); // end of RK section
    
  }while(finish==false);

}
