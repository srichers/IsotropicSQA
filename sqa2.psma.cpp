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
#include "headers/evolve.h"

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
  const double rmax = get_parameter<double>(fin,"tmax") * cgs::constants::c; // cm
  const double dr0  = get_parameter<double>(fin,"dt0")  * cgs::constants::c; // cm
  const double accuracy = get_parameter<double>(fin,"accuracy");
  const double target_impact = get_parameter<double>(fin,"target_impact");
  const double increase = get_parameter<double>(fin,"increase"); // factor by which timestep increases if small error
  const double mixing = get_parameter<double>(fin,"mixing");
  const int step_output = get_parameter<int>(fin,"step_output");
  const int do_oscillate = get_parameter<int>(fin,"do_oscillate");
  const int do_interact = get_parameter<int>(fin,"do_interact");
  cout.flush();

  State s;
  s.r=0;
  s.rho = rho;
  s.T = temperature;
  s.Ye = Ye;
  s.eas = EAS(nulibfilename, eosfilename);
  s.fmatrixf.resize(NM);
  s.fmatrixf[matter]=s.fmatrixf[antimatter] = vector<MATRIX<complex<double>,NF,NF> >(s.eas.ng);
  initialize(s.fmatrixf,s.eas,s.rho,s.T,s.Ye, mixing, do_interact);
  s.dr_block = dr0;
  s.dr_osc = dr0;
  s.dr_int = dr0;
  s.counter = 0;
  
  // random number generator - prevent aliasing in interactions and output
  srand(time(NULL));

  // temporary variable
  vector<vector<MATRIX<complex<double>,NF,NF> > > fmatrixf0 = s.fmatrixf;
  
  // set up output file
  s.foutf.open("f.dat");
  s.foutf.precision(12);
  s.foutf << "# 1:r ";
  for(int i=0; i<s.eas.ng; i++)
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) {
	  int istart = 2*( f2 + f1*2 + m*2*2 + i*2*2*2) + 2;
	  s.foutf << istart   << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"R\t";
	  s.foutf << istart+1 << ":ie"<<i<<"m"<<m<<"f"<<f1<<f2<<"I\t";
	}
  s.foutf << endl;
  s.foutf.flush();
  cout << "iter \t t(s) \t dt(s) \t n_nu(1/ccm) \t n_nubar(1/ccm) \t n_nu-n_nubar(1/ccm) \t interact_impact" << endl;
  Outputvsr(s, 0);

  // ***********************
  // start the loop over r *
  // ***********************
  bool finish=false;
  do{
    if(s.r+s.dr_block>=rmax){
      finish=true;
      s.dr_block = rmax-s.r;
    }
    else finish=false;
    
    double impact=0;
    double rstep = s.r + s.dr_block;
    double r0 = s.r;
    if(do_oscillate)
      evolve_oscillations(s, rstep, accuracy, increase, step_output);
    if(do_interact){
      s.r = r0;
      fmatrixf0 = s.fmatrixf;
      evolve_interactions(s, rstep, accuracy, increase);

      // evaluate the net impact
      impact = 0;
      for(int m=matter; m<=antimatter; m++){
	for(int i=0;i<=s.eas.ng-1;i++){
	  double trace = norm(Trace(fmatrixf0[m][i]));
	  for(flavour f1=e; f1<=mu; f1++){
	    for(flavour f2=e; f2<=mu; f2++){
	      impact = max(impact, norm(s.fmatrixf[m][i][f1][f2]-fmatrixf0[m][i][f1][f2])/trace);
	    } // f2
	  } // f1
	} // i
      } // m

      // timestepping
      s.r = rstep;
      if(impact > 10*target_impact)
	cout << "WARNING: impact="<<impact<< endl;
      if(impact>target_impact*2. or impact<target_impact/2.)
	s.dr_block *= min(increase, target_impact/impact);
      
      // output
      Outputvsr(s, impact);
    }
    
    
    // sanity checks
    for(int m=matter;m<=antimatter;m++){
      for(int i=0;i<=s.eas.ng-1;i++){
	if(real(s.fmatrixf[m][i][e ][e ]) > 1. or
	   real(s.fmatrixf[m][i][mu][mu]) > 1. or
	   real(s.fmatrixf[m][i][e ][e ]) < 0. or
	   real(s.fmatrixf[m][i][mu][mu]) < 0. or
	   imag(s.fmatrixf[m][i][e ][e ]) > accuracy or
	   imag(s.fmatrixf[m][i][mu][mu]) > accuracy){
	  cout << "m"<<m << " i" << i << endl;
	  cout << s.fmatrixf[m][i] << endl;
	  exit(1);
	}
      }
    }
    
    
  }while(finish==false);

  cout << endl << "Finished" << endl;
  return 0;
}
