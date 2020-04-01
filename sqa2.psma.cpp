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
#include<hdf5.h>

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
  const string outputfilename = get_parameter<string>(fin,"outputfilename");
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

  cout << "iter \t t(s) \t dt(s) \t n_nu(1/ccm) \t n_nubar(1/ccm) \t n_nu-n_nubar(1/ccm) \t interact_impact" << endl;
  cout.flush();

  // initialize the state
  State s(nulibfilename, eosfilename, rho, Ye, temperature, dr0, mixing, do_interact);
  ifstream tmp_ifstream(outputfilename);
  hid_t output_file;
  if(tmp_ifstream){
    cout << "Recovering from " << outputfilename << endl;
    output_file = recover(outputfilename, s);
  }
  else if(step_output>0) {
    cout << "Starting new file " << outputfilename << endl;
    output_file = setup_file(outputfilename, s);
    write_data(output_file, s, 0);
  }

  // random number generator - prevent aliasing in interactions and output
  srand(time(NULL));

  // temporary variable
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> fmatrixf0 = s.fmatrixf;


  // ***********************
  // start the loop over r *
  // ***********************
  bool finish=false;
  int next_output = step_output>0 ? rand()%step_output+1 : -1;
  do{
    double impact=0;
    double r0 = s.r;
    double rstep = s.r + s.dr_block * min(5., exponential_random());
    if(rstep>=rmax){
      finish=true;
      s.dr_block = rmax-r0;
      rstep = rmax;
    }
    else finish=false;

    if(do_oscillate)
      evolve_oscillations(s, rstep, accuracy, increase);
    if(do_interact){
      s.r = r0;
      fmatrixf0 = s.fmatrixf;
      evolve_interactions(s, rstep, accuracy, increase);

      // evaluate the net impact
      #pragma omp parallel for collapse(2) reduction(max:impact)
      for(int m=matter; m<=antimatter; m++){
	for(int i=0;i<=s.eas.ng-1;i++){
	  double l = IsospinL(fmatrixf0[m][i]);
	  MATRIX<complex<double>,NF,NF> df = s.fmatrixf[m][i] - fmatrixf0[m][i];
	  for(flavour f1=e; f1<=mu; f1++){
	    for(flavour f2=e; f2<=mu; f2++){
	      impact = max(impact, abs(df[f1][f2]) / l);
	    } // f2
	  } // f1
	} // i
      } // m
    }

    // output
    s.counter++;
    if(step_output>0 and (s.counter>=next_output or finish)){
      write_data(output_file, s, impact);
      next_output = step_output>0 ? s.counter + rand()%step_output + 1 : -1;
    }

    // timestepping
    if(do_interact){
      if(impact > target_impact)
	cout << "WARNING: impact="<<impact<< endl;
      double corrected_impact = impact / (s.r-r0) * s.dr_block;
      if(corrected_impact<.1*target_impact)
	s.dr_block *= min(increase, .1*target_impact/corrected_impact);
      if(corrected_impact>.1*target_impact)
	s.dr_block *= .1*target_impact/corrected_impact;
    }

    // sanity checks
    // for(int m=matter;m<=antimatter;m++){
    //   for(int i=0;i<=s.eas.ng-1;i++){
    // 	if(real(s.fmatrixf[m][i][e ][e ]) > 1. or
    // 	   real(s.fmatrixf[m][i][mu][mu]) > 1. or
    // 	   real(s.fmatrixf[m][i][e ][e ]) < 0. or
    // 	   real(s.fmatrixf[m][i][mu][mu]) < 0. or
    // 	   imag(s.fmatrixf[m][i][e ][e ]) > accuracy or
    // 	   imag(s.fmatrixf[m][i][mu][mu]) > accuracy){
    // 	  cout << "m"<<m << " i" << i << endl;
    // 	  cout << s.fmatrixf[m][i] << endl;
    // 	  exit(1);
    // 	}
    //   }
    // }
    #pragma omp parallel for collapse(2)
    for(int m=matter;m<=antimatter;m++)
      for(int i=0;i<=s.eas.ng-1;i++)
	Hermitize(s.fmatrixf[m][i], accuracy);

  }while(finish==false);

  cout << endl << "Finished" << endl;
  return 0;
}
