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
#include<array>
using std::array;

// headers
#include "headers/matrix.h"
#include "headers/parameters.h"
#include "headers/eigenvalues.h"
#include "headers/mixing angles.h"
#include "headers/jacobians.h"
#include "headers/misc.h"
#include "headers/interact.h"
#include "headers/nulib_interface.h"

//======//
// MAIN //
//======//
int main(int argc, char *argv[]){
  string inputfilename=string(argv[1]);
  ifstream fin(inputfilename.c_str());

  // read in all parameters
  const string nulibfilename = get_parameter<string>(fin,"nulibfilename");
  const string eosfilename = get_parameter<string>(fin,"eosfilename");
  const double rho = get_parameter<double>(fin,"rho");
  const double Ye  = get_parameter<double>(fin,"Ye");
  const double temperature = get_parameter<double>(fin,"temperature"); // MeV
  const string outputfilename = get_parameter<string>(fin,"outputfilename");
  const double rmax = get_parameter<double>(fin,"tmax") * cgs::constants::c; // cm
  const double dr0  = get_parameter<double>(fin,"dt0")  * cgs::constants::c; // cm
  const double accuracy = get_parameter<double>(fin,"accuracy");
  const double increase = get_parameter<double>(fin,"increase"); // factor by which timestep increases if small error
  const double mixing = get_parameter<double>(fin,"mixing");
  const int step_output = get_parameter<int>(fin,"step_output");
  const int step_interact = get_parameter<int>(fin,"step_interact");
  const int do_oscillate = get_parameter<int>(fin,"do_oscillate");
  const int do_interact = get_parameter<int>(fin,"do_interact");
  cout.flush();

  // load the nulib table
  EAS eas(nulibfilename, eosfilename);
    
  // set up output file
  ofstream foutf;
  foutf.open((outputfilename+"/f.dat").c_str());
  foutf.precision(12);
  foutf << "# 1:r ";
  for(int i=0; i<eas.ng; i++)
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  int istart = 2*( f2 + f1*2 + m*2*2 + i*2*2*2) + 2;
	  foutf << istart   << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"R\t";
	  foutf << istart+1 << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"I\t";
	}
  foutf << endl;
  foutf.flush();
    
  // vectors of energies and vacuum eigenvalues
  const vector<vector<double> > kV = set_kV(eas.E);
  const vector<MATRIX<complex<double>,NF,NF> > UV = Evaluate_UV();
  const vector<vector<MATRIX<complex<double>,NF,NF> > > HfV = Evaluate_HfV(kV,UV);
  const vector<vector<MATRIX<complex<double>,NF,NF> > > CV = Evaluate_CV(kV, HfV);
  const vector<vector<vector<double> > > AV = Evaluate_AV(kV,HfV,UV);
    
  // **************************************
  // quantities evaluated at inital point *
  // **************************************
    
  // MSW potential matrix
  MATRIX<complex<double>,NF,NF> VfMSW0, Hf0;
  vector<double> k0, deltak0;
    
  VfMSW0[e][e]=Ve(rho,Ye);
  VfMSW0[mu][mu]=Vmu(rho,Ye);
    
  // cofactor matrices at initial point - will be recycled as cofactor matrices at beginning of every step
  vector<vector< array<MATRIX<complex<double>,NF,NF>,NF> > > 
    C0(NM,vector<array<MATRIX<complex<double>,NF,NF>,NF> >(eas.ng));

  // mixing matrix element prefactors at initial point - will be recycled like C0
  vector<vector<vector<vector<double> > > > 
    A0(NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NF,vector<double>(NF))));
    
  // mixing angles to MSW basis at initial point
  vector<vector<MATRIX<complex<double>,NF,NF> > > U0(NM); 
  U0[matter] = vector<MATRIX<complex<double>,NF,NF> >(eas.ng);
  U0[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(eas.ng);
    
  for(int i=0;i<=eas.ng-1;i++){
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
  vector<MATRIX<complex<double>,NF,NF> > pmatrixm0matter(eas.ng);
  vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf0(NM);
  vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf (NM);
  fmatrixf0[matter]=fmatrixf0[antimatter]=vector<MATRIX<complex<double>,NF,NF> >(eas.ng);
  fmatrixf [matter]=fmatrixf [antimatter]=vector<MATRIX<complex<double>,NF,NF> >(eas.ng);
  initialize(fmatrixf0,eas,0,rho,temperature,Ye, mixing, do_interact);
    
  // ***************************************
  // variables followed as a function of r *
  // ***************************************
  vector<vector<vector<vector<double> > > > 
    Y(NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NS,vector<double>(NY))));
  vector<vector<vector<vector<double> > > > 
    Y0(NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NS,vector<double>(NY))));
  
  // ************************
  // Runge-Kutta quantities *
  // ************************
  int NRK,NRKOrder;
  const double *AA=NULL,**BB=NULL,*CC=NULL,*DD=NULL;
  RungeKuttaCashKarpParameters(NRK,NRKOrder,AA,BB,CC,DD);
    
  vector<vector<vector<vector<vector<double> > > > > 
    Ks(NRK,vector<vector<vector<vector<double> > > >
       (NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NS,vector<double>(NY)))));
    
  // temporaries
  vector<vector<MATRIX<complex<double>,NF,NF> > > ftmp0, dfdr0, dfdr1;

  // **********************
  // start of calculation *
  // **********************
    
      
  // *****************************************
  // initialize at beginning of every domain *
  // *****************************************
  double r=0;
  double r_interact_last = 0;
  double dr = dr0;

  for(state m=matter;m<=antimatter;m++)
    for(int i=0;i<=eas.ng-1;i++)
      Y[m][i] = YIdentity;
  bool finish=false;
  int counter=0;
  fmatrixf = fmatrixf0;
  Outputvsr(foutf,r, fmatrixf);
      
  // ***********************
  // start the loop over r *
  // ***********************
  double n0,nbar0;
  do{ 

    // output to stdout
    double intkm = int(r/1e5)*1e5;
    if(counter%step_output==0){
      double n=0, nbar=0;
      double coeff = 4.*M_PI / pow(cgs::constants::c,3);
      for(int i=0; i<eas.ng; i++){
	for(flavour f1=e; f1<=mu; f1++){
	  n    += real(fmatrixf0[    matter][i][f1][f1]) * eas.nu[i]*eas.nu[i]*eas.dnu[i]*coeff;
	  nbar += real(fmatrixf0[antimatter][i][f1][f1]) * eas.nu[i]*eas.nu[i]*eas.dnu[i]*coeff;
	}
      }
      if(counter==0){
	n0=n;
	nbar0=nbar;
	cout << "iter \t t(s) \t dt(s) \t n_nu("<< n0<<"/ccm) \t n_nubar("<<nbar0<<"/ccm) \t n_nu-n_nubar("<<n0-nbar0<<"/ccm)" << endl;
      }
      cout << counter << "\t";
      cout << r/cgs::constants::c << "\t";
      cout << dr/cgs::constants::c << "\t";
      cout << n/n0 << "\t" << nbar/nbar0 << "\t" << (n-nbar)/(n0-nbar0) << endl;
      cout.flush();
    }
	  
    // save initial values in case of repeat
    double r0=r;
    Y0=Y;
    getP(r,U0,fmatrixf0,eas.nu,eas.dnu,pmatrixm0matter);

    // beginning of RK section
    double maxerror;
    bool repeat, do_reset;
    do{
      do_reset = false;
      repeat = false;
      maxerror=0.;

      if(do_oscillate){

	// RK integration for oscillation
	// if potential changes with r, update potential inside rk loop
	for(int k=0;k<=NRK-1;k++){
	  assert(CC[k] == CC[k]);
	  r=r0+AA[k]*dr;
          #pragma omp parallel for collapse(4)
	  for(int m=0; m<=1; m++){ // 0=matter 1=antimatter
	    for(int i=0;i<=eas.ng-1;i++){
	      for(int x=0;x<=1;x++){ // 0=msw 1=si
		for(int j=0;j<=NY-1;j++){

		  Y[m][i][x][j] = Y0[m][i][x][j];
		  for(int l=0;l<=k-1;l++)
		    Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];
		}
	      }
	    }
	  }
	  
	  K(r,dr,rho,Ye,pmatrixm0matter,HfV,Y,C0,A0,Ks[k]);
	}
	  
	// increment all quantities from oscillation
        #pragma omp parallel for collapse(4) reduction(max:maxerror)
	for(int m=0; m<=1; m++)  // 0=matter 1=antimatter
	  for(int i=0;i<=eas.ng-1;i++)
	    for(int x=0; x<=1; x++)  // 0=msw 1=si
	      for(int j=0;j<=NY-1;j++){

		double Yerror = 0.;
		Y[m][i][x][j] = Y0[m][i][x][j];
		for(int k=0;k<=NRK-1;k++){
		  assert(Ks[k][m][i][x][j] == Ks[k][m][i][x][j]);
		  Y[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		  Yerror += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
		  assert(Y[m][i][x][j] == Y[m][i][x][j]);
		}
		maxerror = max( maxerror, fabs(Yerror) );
	      }

      }
      r=r0+dr;

      // convert fmatrix from flavor basis to mass basis, oscillate, convert back
      #pragma omp parallel for collapse(2) reduction(||:do_reset)
      for(int m=matter; m<=antimatter; m++){
	for(int i=0; i<eas.ng; i++){
	  MATRIX<complex<double>,NF,NF> SSMSW = W(Y[m][i][msw])*B(Y[m][i][msw]);
	  MATRIX<complex<double>,NF,NF> SSSI  = W(Y[m][i][si ])*B(Y[m][i][si ]);
	  MATRIX<complex<double>,NF,NF> SThisStep = U0[m][i] * SSMSW*SSSI * Adjoint(U0[m][i]);
	  fmatrixf[m][i] = SThisStep * fmatrixf0[m][i] * Adjoint(SThisStep);

	  // test that the MSW S matrix is close to diagonal
	  // test the SI S matrix is close to diagonal
	  // test amount of interaction error accumulated
	  if(norm(SSMSW[0][0])+0.1<norm(SSMSW[0][1]) or
	     norm(SSSI [0][0])+0.1<norm(SSSI [0][1]))
	    do_reset = true;
	}
      }

      if(do_interact && counter%step_interact==0){
		
	// interact with the matter
	// if fluid changes with r, update opacities here, too
	double dr_interact = (r-r_interact_last);
	ftmp0 = fmatrixf;
	dfdr0 = my_interact(fmatrixf, rho, temperature, Ye, eas);
        #pragma omp parallel for collapse(2)
	for(int m=matter; m<=antimatter; m++)
	  for(int i=0; i<eas.ng; i++)
	    ftmp0[m][i] += dfdr0[m][i] * dr_interact;

	double interact_impact = 0;
	dfdr1 = my_interact(ftmp0, rho, temperature, Ye, eas);
        #pragma omp parallel for collapse(2)
	#pragma omp reduction(max:maxerror) reduction(max:interact_impact)
	for(int m=matter; m<=antimatter; m++){
	  for(int i=0; i<eas.ng; i++){
	    MATRIX<complex<double>,NF,NF> df = (dfdr0[m][i] + dfdr1[m][i])*0.5*dr_interact;
	    fmatrixf[m][i] += df;
	    
	    double trace = real(fmatrixf[m][i][e][e]+fmatrixf[m][i][mu][mu]);
	    for(flavour f1=e; f1<=mu; f1++){
	      for(flavour f2=e; f2<=mu; f2++){
		double error = fabs(ftmp0[m][i][f1][f2]-fmatrixf[m][i][f1][f2])/trace;
		double impact = fabs(df[f1][f2])/trace;
		maxerror = max(maxerror, do_oscillate?impact:error);
		interact_impact = max(interact_impact, impact);
	      }
	    }
	  }
	}

	if(interact_impact >= 0.1*accuracy) do_reset = true;
	if(do_oscillate) assert(interact_impact < accuracy);
      }
	  
      // decide whether to accept step, if not adjust step size and reset variables
      if(maxerror>accuracy){
	dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	repeat=true;
	r=r0;
	Y=Y0;
      }
      else{
	dr *= increase;
	if(maxerror>0) dr *= min( 1.0, pow(accuracy/maxerror,1./max(1,NRKOrder))/increase );
      }
      dr = max(dr, 4.*r*numeric_limits<double>::epsilon());
      dr = min(dr, rmax-r);
    }while(repeat==true); // end of RK section

    // update fmatrixf0 if necessary
    if(do_reset){
      r_interact_last = r;
      #pragma omp parallel for collapse(2)
      for(int m=0;m<=1;m++){ // 0=matter 1=antimatter
	for(int i=0;i<=eas.ng-1;i++){
	  fmatrixf0[m][i] = fmatrixf[m][i];
	  Y[m][i] = YIdentity;
	}
      }
    }
    else{ // take modulo 2 pi of phase angles
#pragma omp simd collapse(2)
      for(int m=0;m<=1;m++){ // 0=matter 1=antimatter
	for(int i=0;i<=eas.ng-1;i++){
	  Y[m][i][msw][2]=fmod(Y[m][i][msw][2],M_2PI);
	  Y[m][i][msw][4]=fmod(Y[m][i][msw][4],1.0);
	  Y[m][i][msw][5]=fmod(Y[m][i][msw][5],1.0);

	  Y[m][i][si ][2]=fmod(Y[m][i][si ][2],M_2PI);
	  Y[m][i][si ][4]=fmod(Y[m][i][si ][4],1.0);
	  Y[m][i][si ][5]=fmod(Y[m][i][si ][5],1.0);
	}
      }
    }
      
    // output to file
    if(r>=rmax) finish=true;
    if(counter%step_output==0 or finish)
      Outputvsr(foutf,r,fmatrixf);
    counter++;
  } while(finish==false);

  cout<<"\nFinished\n\a"; cout.flush();

  return 0;
}
