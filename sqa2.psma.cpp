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

#include <complex>
using std::complex;
using std::polar;
using std::abs;
using std::arg;
using std::real;
using std::imag;
#include <cstdarg>
using std::va_list;
#include <cstdlib>
#include<iostream>
using std::cout;
#include<ostream>
using std::ostream;
using std::endl;
using std::flush;
#include<fstream>
using std::ifstream;
using std::ofstream;
#include<sstream>
using std::stringstream;
#include<algorithm>
using std::min;
using std::max;
using std::swap;
using std::lower_bound;
#include<string>
using std::string;
#include <utility>
using std::pair;
#include<functional>
#include<limits>
using std::numeric_limits;
#include<vector>
using std::vector;

// headers
#include "headers/matrix.h"
#include "headers/parameters.h"
#include "headers/flavour basis.h"
#include "headers/eigenvalues.h"
#include "headers/mixing angles.h"
#include "headers/adiabatic basis.h"
#include "headers/jacobians.h"
#include "headers/misc.h"
#include "headers/update.h"
#include "headers/interact.h"
#include "headers/nulib_interface.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
    int in=1;
    string inputfilename;
    ofstream foutf;
    string outputfilename,vfilename, spectrapath, nulibfilename, eosfilename;
    string outputfilenamestem;
    string note;
    
    inputfilename=string(argv[in]);
    ifstream fin(inputfilename.c_str());
    
    // load the nulib table
    fin>>nulibfilename;
    cout << nulibfilename << endl;
    fin >> eosfilename;
    cout << eosfilename << endl;
    eas = EAS(nulibfilename, eosfilename);

    double rho, Ye, temperature/*MeV*/; // rho is the mass density
    fin>>rho;
    fin>>Ye;
    fin>>temperature;
    fin>>outputfilename;
    outputfilenamestem = outputfilename+"/";
    
    double rmax;
    fin>>rmax; // s
    rmax *= cgs::constants::c; // s -> cm
    m1=0.;
    fin>>dm21;
    fin>>theta12V;
    
    alphaV[0]=0.;
    alphaV[1]=0.;
    
    betaV[0]=0.;
    
    double accuracy;
    fin>>accuracy;
    fin>>note;

    cout<<"\n\n*********************************************************\n";
    cout<<"\nrho\t"<<rho;
    cout<<"\nYe\t"<<Ye;
    cout<<"\nT\t"<<temperature;
    cout<<"\noutput\t"<<outputfilename;
    cout<<"\ttmax\t"<<rmax/cgs::constants::c;
    cout<<"\n\nm1\t"<<m1<<"\tdm21^2\t"<<dm21;
    cout<<"\ntheta12V\t"<<theta12V;
    cout<<"\nalpha1V\t"<<alphaV[0]<<"\talpha2V\t"<<alphaV[1];
    cout<<"\nbeta1V\t"<<betaV[0];    
    cout<<"\naccuracy\t"<<accuracy<<"\n";
    cout.flush();
    
    // output filestreams: the arrays of ofstreams cannot use the vector container - bug in g++
    foutf.open((outputfilename+"/f.dat").c_str());
    foutf.precision(12);
    foutf << "# 1:r ";
    for(int i=0; i<NE; i++)
      for(state m=matter; m<=antimatter; m++)
	for(flavour f1=e; f1<=mu; f1++)
	  for(flavour f2=e; f2<=mu; f2++) {
	    int istart = 2*( f2 + f1*2 + m*2*2 + i*2*2*2) + 2;
	    foutf << istart   << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"R\t";
	    foutf << istart+1 << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"I\t";
    }
    foutf << endl;
    foutf.flush();
    
    // unit conversion to cgs
    //Emin *= 1.*mega*cgs::units::eV;
    //Emax *= 1.*mega*cgs::units::eV;
    m1   *= 1.*cgs::units::eV/cgs::constants::c2;
    dm21 *= 1.*cgs::units::eV*cgs::units::eV/cgs::constants::c4;
    theta12V *= M_PI/180.;
    c12V = cos(theta12V);
    s12V = sin(theta12V);
    
    // *************************************************
    // set up global variables defined in parameters.h *
    // *************************************************

    // vectors of energies and vacuum eigenvalues
    set_Ebins(E);
    kV = vector<vector<double> >(NE,vector<double>(NF));
    set_kV(kV);

    // determine eigenvalue ordering
    if(kV[0][1]>kV[0][0]){
      a1=-1;
      a2=+1;
      cout<<"\n\nNormal hierarchy" << endl;
    }
    else{ 
      if(kV[0][1]<kV[0][0]){
	a1=+1; 
	a2=-1; 
	cout<<"\n\nInverted hierarchy" << endl;
      }
      else{ 
	cout<<endl<<endl<<"Neither normal or Inverted"<<endl; 
	abort();
      }
    }
    
    // vaccum mixing matrices and Hamiltonians
    Evaluate_UV();
    
    HfV[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    HfV[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    Evaluate_HfV();
    
    // cofactor matrices in vacuum
    CV=vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF));
    Evaluate_CV();
    
    // mixing matrix element prefactors in vacuum
    AV=vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF)));
    Evaluate_AV();
    
    // **************************************
    // quantities evaluated at inital point *
    // **************************************
    
    // MSW potential matrix
    MATRIX<complex<double>,NF,NF> VfMSW0, Hf0;
    vector<double> k0, deltak0;
    
    VfMSW0[e][e]=Ve(rho,Ye);
    VfMSW0[mu][mu]=Vmu(rho,Ye);
    
    // cofactor matrices at initial point - will be recycled as cofactor matrices at beginning of every step
    vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > 
      C0(NM,vector<vector<MATRIX<complex<double>,NF,NF> > >(NE,vector<MATRIX<complex<double>,NF,NF> >(NF)));

    // mixing matrix element prefactors at initial point - will be recycled like C0
    vector<vector<vector<vector<double> > > > 
      A0(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NF,vector<double>(NF))));
    
    // mixing angles to MSW basis at initial point
    U0[matter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    U0[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(NE);
    
    for(int i=0;i<=NE-1;i++){
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
    
    // density matrices at initial point, rhomatrixm0 - but not rhomatrixf0
    // will be updated whenever discontinuities are crossed and/or S is reset
    vector<MATRIX<complex<double>,NF,NF> > pmatrixm0matter(NE);
    vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf(NM);
    fmatrixf[matter]=fmatrixf[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(NE);

    // yzhu14 density/potential matrices art rmin
    double mixing;
    fin>>mixing;

    // ***************************************
    // quantities needed for the calculation *
    // ***************************************
    double r,r0,dr,drmin;
    double maxerror,increase=3.;
    bool repeat, finish, resetflag, output;
    int counterout,step;
    
    // comment out if not following as a function of r
    fin>>step;

    // do we oscillate and interact?
    int do_oscillate, do_interact;
    fin>>do_oscillate;
    fin>>do_interact;
    initialize(fmatrixf,0,rho,temperature,Ye, mixing, do_interact);
    
    // ***************************************
    // variables followed as a function of r *
    // ***************************************
    
    vector<vector<vector<vector<double> > > > 
      Y(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    vector<vector<vector<vector<double> > > > 
      Y0(NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY))));
    
    // cofactor matrices
    vector<vector<vector<MATRIX<complex<double>,NF,NF> > > > C=C0;;
    // mixing matrix prefactors
    vector<vector<vector<vector<double> > > > A=A0;
        
    // ************************
    // Runge-Kutta quantities *
    // ************************
    int NRK,NRKOrder;
    const double *AA=NULL,**BB=NULL,*CC=NULL,*DD=NULL;
    RungeKuttaCashKarpParameters(NRK,NRKOrder,AA,BB,CC,DD);
    
    vector<vector<vector<vector<vector<double> > > > > 
      Ks(NRK,vector<vector<vector<vector<double> > > >
	 (NM,vector<vector<vector<double> > >(NE,vector<vector<double> >(NS,vector<double>(NY)))));
    
    // temporaries
    MATRIX<complex<double>,NF,NF> SSMSW,SSSI,SThisStep;
    vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf0(NM);
    vector<vector<MATRIX<complex<double>,NF,NF> > > ftmp0, dfdr0, dfdr1;

    // **********************
    // start of calculation *
    // **********************
    
      cout << "t(s)  dt(s)  n_nu(1/ccm)  n_nubar(1/ccm) n_nu-n_nubar(1/ccm)" << endl;
      cout.flush();
      
      // *****************************************
      // initialize at beginning of every domain *
      // *****************************************
      r=0;
      dr=1e-3*cgs::units::cm;
      drmin=4.*r*numeric_limits<double>::epsilon();
      
      for(state m=matter;m<=antimatter;m++)
	for(int i=0;i<=NE-1;i++){
	  Y0[m][i][msw][0] = Y0[m][i][si][0] = M_PI/2.;
	  Y0[m][i][msw][1] = Y0[m][i][si][1] = M_PI/2.;
	  Y0[m][i][msw][2] = Y0[m][i][si][2] = 0.;
	  Y0[m][i][msw][3] = Y0[m][i][si][3] = 1.; // The determinant of the S matrix
	  Y0[m][i][msw][4] = Y0[m][i][si][4] = 0.;
	  Y0[m][i][msw][5] = Y0[m][i][si][5] = 0.;
	}
      Y=Y0;
      
      // *************************************************
      // comment out if not following as a function of r *
      // *************************************************
	
      finish=output=false;
      counterout=1;
      Outputvsr(foutf,r, fmatrixf);
	
      for(state m=matter; m<=antimatter; m++)
	for(int i=0; i<NE; i++)
	  for(flavour f1=e; f1<=mu; f1++)
	    for(flavour f2=e; f2<=mu; f2++)
	      assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
      
      // ***********************
      // start the loop over r *
      // ***********************
      double n0,nbar0;
      do{ 

	// output to stdout
	double intkm = int(r/1e5)*1e5;
	if(r - intkm <= dr){
	  double n=0, nbar=0;
	  double coeff = 4.*M_PI / pow(cgs::constants::c,3);
	  for(int i=0; i<NE; i++){
	    for(flavour f1=e; f1<=mu; f1++){
	      n    += real(fmatrixf[    matter][i][f1][f1]) * nu[i]*nu[i]*dnu[i]*coeff;
	      nbar += real(fmatrixf[antimatter][i][f1][f1]) * nu[i]*nu[i]*dnu[i]*coeff;
	    }
	  }
	  if(r==0){
	    n0=n;
	    nbar0=nbar;
	  }
	  cout << r/cgs::constants::c << " ";
	  cout << dr/cgs::constants::c << " ";
	  cout << n/n0 << " " << nbar/nbar0 << " " << (n-nbar)/(n0-nbar0) << endl;
	  cout.flush();
	}
	  
	// save initial values in case of repeat
	r0=r;
	C0=C;
	A0=A;
	fmatrixf0 = fmatrixf;
	getP(r,U0,fmatrixf,pmatrixm0matter);

	// beginning of RK section
	do{ 
	  repeat=false;
	  maxerror=0.;

	  if(do_oscillate){
	    // RK integration for oscillation
	    for(int k=0;k<=NRK-1;k++){
	      r=r0+AA[k]*dr;
	      Y=Y0;
	      for(state m = matter; m <= antimatter; m++)
		for(int i=0;i<=NE-1;i++)
		  for(solution x=msw;x<=si;x++)
		    for(int j=0;j<=NY-1;j++)
		      for(int l=0;l<=k-1;l++)
			Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];

	      K(r,dr,rho,Ye,pmatrixm0matter,Y,C,A,Ks[k]);
	    }
	  
	    // increment all quantities from oscillation
	    getP(r0+dr,U0,fmatrixf,pmatrixm0matter);
	    Y=Y0;
	    for(state m=matter;m<=antimatter;m++){
	      for(int i=0;i<=NE-1;i++){

		// integrate Y
		for(solution x=msw;x<=si;x++){
		  for(int j=0;j<=NY-1;j++){
		    double Yerror = 0.;
		    for(int k=0;k<=NRK-1;k++){
		      assert(CC[k] == CC[k]);
		      assert(Ks[k][m][i][x][j] == Ks[k][m][i][x][j]);
		      Y[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		      Yerror += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
		      assert(Y[m][i][x][j] == Y[m][i][x][j]);
		    }
		    maxerror = max( maxerror, fabs(Yerror) );
		  }
		}
	      
		// convert fmatrix from flavor basis to mass basis, oscillate, convert back
		SSMSW = W(Y[m][i][msw])*B(Y[m][i][msw]);
		SSSI  = W(Y[m][i][si ])*B(Y[m][i][si ]);
		SThisStep = U0[m][i] * SSMSW*SSSI * Adjoint(U0[m][i]);
		fmatrixf[m][i] = MATRIX<complex<double>,NF,NF>
		  (SThisStep * fmatrixf[m][i] * Adjoint(SThisStep));
	      }
	    }
	    Y=Y0;
	    C=UpdateC(r,rho,Ye);
	    A=UpdateA(C,C0,A0);
	  }
	  r=r0+dr;

	  if(do_interact){
	    // interact with the matter
	    // if fluid changes with r, update opacities here, too
	    dfdr0 = my_interact(fmatrixf, rho, temperature, Ye);
	    ftmp0 = fmatrixf;
	    for(state m=matter; m<=antimatter; m++)
	      for(int i=0; i<NE; i++)
		for(flavour f1=e; f1<=mu; f1++)
		  for(flavour f2=e; f2<=mu; f2++)
		    ftmp0[m][i][f1][f2] += dfdr0[m][i][f1][f2] * dr;
	    dfdr1 = my_interact(ftmp0, rho, temperature, Ye);
	    for(state m=matter; m<=antimatter; m++)
	      for(int i=0; i<NE; i++){
		double trace = real(fmatrixf[m][i][e][e]+fmatrixf[m][i][mu][mu]);
		for(flavour f1=e; f1<=mu; f1++)
		  for(flavour f2=e; f2<=mu; f2++){
		    fmatrixf[m][i][f1][f2] += 0.5*dr
		      * (dfdr0[m][i][f1][f2] + dfdr1[m][i][f1][f2]);
		    maxerror = max(fabs(ftmp0[m][i][f1][f2]-fmatrixf[m][i][f1][f2])/trace,
				   maxerror);
		  }
	      }
	  }
	  
	  // decide whether to accept step, if not adjust step size
	  if(maxerror>accuracy){
	    dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	    if(dr > drmin) repeat=true;
	  }
	  if(maxerror>0)
	    dr = min(dr*pow(accuracy/maxerror,1./max(1,NRKOrder)),increase*dr);
	  else dr *= increase;
	  drmin = 4.*r*numeric_limits<double>::epsilon();
	  dr = max(dr,drmin);
	  if(r+dr > rmax) dr = rmax-r;
	  
	  // reset integration variables to those at beginning of step
	  if(repeat==true){
	    r=r0;
	    C=C0;
	    A=A0;
	    fmatrixf = fmatrixf0;
	  }
	  
	}while(repeat==true); // end of RK section


	// comment out if not following as a function of r
	if(r>=rmax){
	  finish=true;
	  output=true;
	}
	if(counterout==step){
	  output=true;
	  counterout=1;
	}
	else counterout++;
	if(output==true || finish==true){
	  Outputvsr(foutf,r,fmatrixf);
	  output=false;
	}

      } while(finish==false);

  cout<<"\nFinished\n\a"; cout.flush();

  return 0;
}


