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
#include<ctime>

// headers
#include "headers/matrix.h"
#include "headers/parameters.h"
#include "headers/eigenvalues.h"
#include "headers/mixing angles.h"
#include "headers/jacobians.h"
#include "headers/interact.h"
#include "headers/nulib_interface.h"
#include "headers/misc.h"

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
  const int do_oscillate = get_parameter<int>(fin,"do_oscillate");
  const int do_interact = get_parameter<int>(fin,"do_interact");
  cout.flush();

  // load the nulib table
  EAS eas(nulibfilename, eosfilename);
    
  // *******************
  // Mixing Quantities *
  // *******************
    
  // vectors of energies and vacuum eigenvalues
  const vector<vector<double> > kV = set_kV(eas.E);
  const vector<MATRIX<complex<double>,NF,NF> > UV = Evaluate_UV();
  const vector<vector<MATRIX<complex<double>,NF,NF> > > HfV = Evaluate_HfV(kV,UV);
  const vector<vector<MATRIX<complex<double>,NF,NF> > > CV = Evaluate_CV(kV, HfV);
  const vector<array<array<double,NF>,NF> > AV = Evaluate_AV(kV,HfV,UV);
    
  // MSW potential matrix
  MATRIX<complex<double>,NF,NF> VfMSW0, Hf0;
  VfMSW0[e][e]=Ve(rho,Ye);
  VfMSW0[mu][mu]=Vmu(rho,Ye);
    
  // other matrices
  vector<vector< array<MATRIX<complex<double>,NF,NF>,NF> > > C0(NM); // cofactor matrices at initial point
  vector<vector< array<array<double,NF>,NF> > > A0(NM); // mixing matrix element prefactors at initial point
  vector<vector<MATRIX<complex<double>,NF,NF> > > U0(NM); // mixing angles to MSW basis at initial point
  for(int m=matter; m<=antimatter; m++){
    C0[m] = vector<array<MATRIX<complex<double>,NF,NF>,NF> >(eas.ng);
    A0[m] = vector<array<array<double,NF>,NF> >(eas.ng);
    U0[m] = vector<MATRIX<complex<double>,NF,NF> >(eas.ng);
  }

  for(int i=0;i<=eas.ng-1;i++){
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
    Y(NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NS,vector<double>(NY))));
  vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf0(NM);
  fmatrixf0[matter]=fmatrixf0[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(eas.ng);
  initialize(fmatrixf0,eas,rho,temperature,Ye, mixing, do_interact);
  vector<vector<MATRIX<complex<double>,NF,NF> > > dfCumulative(fmatrixf0); // trash initialization
  for(state m=matter;m<=antimatter;m++){
    for(int i=0;i<=eas.ng-1;i++){
      Y[m][i] = YIdentity;
      dfCumulative[m][i] *= 0;
    }
  }
    
  // temporaries
  vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf(fmatrixf0);
  vector<vector<MATRIX<complex<double>,NF,NF> > > dfdr (fmatrixf0), df  (fmatrixf0); // trash initialization
  vector<vector<MATRIX<complex<double>,NF,NF> > > SSMSW(fmatrixf0), SSSI(fmatrixf0); // trash initialization
  vector<vector<MATRIX<complex<double>,NF,NF> > > pmatrixm0(fmatrixf0);                 // trash initialization
  vector<vector<vector<vector<double> > > > 
    Y0(NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NS,vector<double>(NY))));
  getP(U0,fmatrixf0,eas.nu,eas.dnu,pmatrixm0);
  
  // ************************
  // Runge-Kutta quantities *
  // ************************
  int NRK,NRKOrder;
  const double *AA=NULL,**BB=NULL,*CC=NULL,*DD=NULL;
  RungeKuttaCashKarpParameters(NRK,NRKOrder,AA,BB,CC,DD);
  vector<vector<vector<vector<vector<double> > > > > 
    Ks(NRK,vector<vector<vector<vector<double> > > >
       (NM,vector<vector<vector<double> > >(eas.ng,vector<vector<double> >(NS,vector<double>(NY)))));

  // random number generator - prevent aliasing in interactions and output
  srand(time(NULL));
  
  // *****************************************
  // initialize at beginning of every domain *
  // *****************************************
  double r_interact_last=0;
  double dr_interact = dr0;
  double r=0;
  double dr = dr0;
  bool finish=false;
  int counter=0;
  int next_output = rand()%step_output;
  double max_impact=0;

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
  cout << "iter \t t(s) \t dt(s) \t n_nu(1/ccm) \t n_nubar(1/ccm) \t n_nu-n_nubar(1/ccm) \t max_impact" << endl;
  Outputvsr(foutf,r,dr,counter,eas, fmatrixf, 0);

  // ***********************
  // start the loop over r *
  // ***********************
  do{ 
 
    // save initial values in case of repeat
    double r0=r;
    Y0=Y;

    // beginning of RK section
    double maxerror;
    bool repeat, calc_interaction;
    do{
      maxerror=0.;
      getP(U0,fmatrixf0,eas.nu,eas.dnu,pmatrixm0);
      if(do_interact and r_interact_last+dr_interact < r+dr){
	dr = r_interact_last + dr_interact - r;
	calc_interaction = true;
      }
      else calc_interaction = false;

      //==============//
      // DO_OSCILLATE //
      //==============//
      if(do_oscillate){
	// RK integration for oscillation
	for(int k=0;k<=NRK-1;k++){
	  assert(CC[k] == CC[k]);
	  r=r0+AA[k]*dr;
	  Y=Y0;
	  for(int l=0;l<=k-1;l++)
	    for(int m=0; m<=1; m++) // 0=matter 1=antimatter
	      for(int i=0;i<=eas.ng-1;i++)
		for(int x=0;x<=1;x++) // 0=msw 1=si
		  for(int j=0;j<=NY-1;j++)
		    Y[m][i][x][j] += BB[k][l] * Ks[l][m][i][x][j];
	  
	  K(dr,rho,Ye,pmatrixm0,HfV,Y,C0,A0,Ks[k]);
	}
	  
	// increment all quantities from oscillation
	Y=Y0;
#pragma omp parallel for collapse(2) reduction(max:maxerror)
	for(int m=matter; m<=antimatter; m++){
	  for(int i=0;i<=eas.ng-1;i++){

	    for(int x=msw; x<=si; x++){
	      for(int j=0;j<=NY-1;j++){
		double Yerror = 0.;
		for(int k=0;k<=NRK-1;k++){
		  Y[m][i][x][j] += CC[k] * Ks[k][m][i][x][j];
		  Yerror += (CC[k]-DD[k]) * Ks[k][m][i][x][j];
		}
		maxerror = max( maxerror, fabs(Yerror) );
	      } // j
	    } // x

	    // get oscillated f
	    SSMSW[m][i] = W(Y[m][i][msw])*B(Y[m][i][msw]);
	    SSSI [m][i] = W(Y[m][i][si ])*B(Y[m][i][si ]);
	    MATRIX<complex<double>,NF,NF> SThisStep = U0[m][i] * SSMSW[m][i]*SSSI[m][i] * Adjoint(U0[m][i]);
	    fmatrixf[m][i] = SThisStep * fmatrixf0[m][i] * Adjoint(SThisStep);

	  } // i
	} // m
      } // do_oscillate
      else fmatrixf = fmatrixf0;

      //==============//
      // TIMESTEPPING //
      //==============//
      r=r0+dr;
      if(maxerror>accuracy){
	dr *= 0.9 * pow(accuracy/maxerror, 1./(NRKOrder-1.));
	repeat=true;
	r=r0;
	Y=Y0;
      }
      else{
	dr *= increase;
	repeat = false;
	if(maxerror>0) dr *= min( 1.0, pow(accuracy/maxerror,1./max(1,NRKOrder))/increase );
      }
      dr = max(dr, 4.*r*numeric_limits<double>::epsilon());
      dr = min(dr, rmax-r);

    }while(repeat==true); // end of RK section

    //=============//
    // DO_INTERACT //
    //=============//
    if(calc_interaction){
      dfdr = my_interact(fmatrixf, rho, temperature, Ye, eas);
      for(int m=matter;m<=antimatter;m++){ // 0=matter 1=antimatter
	for(int i=0;i<=eas.ng-1;i++){
	  dfCumulative[m][i] += dfdr[m][i] * dr_interact;
	  fmatrixf[m][i] += dfCumulative[m][i];
	}
      }
      max_impact = 0;
    }

    //========================//
    // RESETTING/ACCUMULATING //
    //========================//
    //#pragma omp parallel for collapse(2) reduction(max:max_impact)
    for(int m=matter;m<=antimatter;m++){ // 0=matter 1=antimatter
      for(int i=0;i<=eas.ng-1;i++){
	bool do_reset = false;

	if(calc_interaction){
	  // test amount of interaction error accumulated
	  double trace = norm(Trace(fmatrixf[m][i]));
	  for(flavour f1=e; f1<=mu; f1++){
	    for(flavour f2=e; f2<=mu; f2++){
	      double impact = norm(dfCumulative[m][i][f1][f2])/trace;
	      max_impact = max(max_impact,impact);
	      if(impact >= 0.1*accuracy) do_reset = true;
	    }
	  }
	}

	if(do_oscillate){
	  // test that the S matrices are close to diagonal
	  if(norm(SSMSW[m][i][0][0])+0.1<norm(SSMSW[m][i][0][1]) or
	     norm(SSSI [m][i][0][0])+0.1<norm(SSSI [m][i][0][1])){
	    do_reset = true;
	  }
	}

	if(do_reset){
	  fmatrixf0[m][i] = fmatrixf[m][i];
	  dfCumulative[m][i] *= 0;
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

	// sanity checks
	if(real(fmatrixf[m][i][e ][e ]) > 1. or
	   real(fmatrixf[m][i][mu][mu]) > 1. or
	   real(fmatrixf[m][i][e ][e ]) < 0. or
	   real(fmatrixf[m][i][mu][mu]) < 0. or
	   imag(fmatrixf[m][i][e ][e ]) > accuracy or
	   imag(fmatrixf[m][i][mu][mu]) > accuracy){
	  cout << "m"<<m << " i" << i << endl;
	  cout << fmatrixf[m][i] << endl;
	  exit(1);
	}
      }
    }

    //========//
    // OUTPUT //
    //========//
    if(r>=rmax) finish=true;
    if(calc_interaction){
      if(max_impact > accuracy)
	cout << "WARNING: max_impact = "<<max_impact << endl;
      r_interact_last = r;
      dr_interact *= min(.1*accuracy / max_impact, increase);
    }
    if(counter==next_output or finish){
      Outputvsr(foutf,r,dr,counter,eas, fmatrixf,max_impact);
      next_output = counter + rand()%step_output + 1;
      assert(next_output>counter);
    }
    counter++;

  } while(finish==false);

  cout << endl << "Finished" << endl;
  return 0;
}
