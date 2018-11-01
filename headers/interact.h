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

#include "nulib_interface.h"

//============//
// Initialize //
//============//
void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
		EAS& eas,
		const double r,
		const double rho,
		const double T,
		const double Ye,
		const double mixing,
		const int do_interact){
  // T should be MeV
  cout << "Setting initial data." << endl;
  cout << "rho = " << rho << " g/ccm" << endl;
  cout << "T = " << T << " MeV" << endl;
  cout << "Ye = " << Ye << endl;
  eas.set(rho,T,Ye,do_interact);
  const unsigned NE = fmatrixf[0].size();
  
  for(int i=0; i<NE; i++){
    for(state m=matter; m<=antimatter; m++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++) 
	  fmatrixf[m][i][f1][f2] = 0;
    double fe = eas.fermidirac(0,i);
    double fa = eas.fermidirac(1,i);
    double fx = eas.fermidirac(2,i);

    fmatrixf[    matter][i][e ][e ] = fe; 
    fmatrixf[    matter][i][mu][mu] = fx;
    fmatrixf[    matter][i][e ][mu] = sqrt(fe*fx)*mixing;
    fmatrixf[    matter][i][mu][e ] = sqrt(fe*fx)*mixing;
    fmatrixf[antimatter][i][e ][e ] = fa;
    fmatrixf[antimatter][i][mu][mu] = fx;
    fmatrixf[antimatter][i][e ][mu] = sqrt(fa*fx)*mixing;
    fmatrixf[antimatter][i][mu][e ] = sqrt(fa*fx)*mixing;
      
    cout << "GROUP " << i << endl;
    cout << "eas.emis \t eas.abs \t BB \t eas.scat \t f" << (eas.do_pair ? "\t eas.pair_emis" : "") << endl;
    for(int is=0; is<4; is++){
      cout << eas.emis(is,i) << "\t";
      cout << eas.abs (is,i) << "\t";
      cout << eas.emis(is,i)/eas.abs(is,i) << "\t";
      cout << eas.scat(is,i) << "\t";
      if(is==0) cout << real(fmatrixf[    matter][i][e ][e ]) << "\t";
      if(is==1) cout << real(fmatrixf[antimatter][i][e ][e ]) << "\t";
      if(is==2) cout << real(fmatrixf[    matter][i][mu][mu]) << "\t";
      if(is==3) cout << real(fmatrixf[antimatter][i][mu][mu]) << "\t";
      if(eas.do_pair){
      	double pair_emis = 0;
      	for(int j=0; j<NE; j++){
	  double conv_to_in_rate = exp(-(eas.E[j]+eas.E[i])/(eas.temperature*1e6*cgs::units::eV));
      	  pair_emis += 4.*M_PI*eas.nu[j]*eas.nu[j]*eas.dnu[j]/cgs::constants::c4 * 0.5*eas.Phi0pair(is,i,j) * conv_to_in_rate;
	}
      	cout << pair_emis << "\t";
      }
      cout << endl;
    }
    cout << endl;
  }

  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++){
	  assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
	  assert(abs(fmatrixf[m][i][f1][f2]) < 1.0);
	}
}

//=============================//
// Self-Interaction Potentials //
//=============================//
MATRIX<double,2,2> avg_matrix(const double eval, const double muval){
  MATRIX<double,2,2> result;
  result[e ][e ] = eval;
  result[mu][mu] = muval;
  result[e ][mu] = result[mu][e] = (eval + muval) / 2.;
  return result;
}
MATRIX<double,2,2> tilde_matrix(const double eval, const double muval){
  MATRIX<double,2,2> result;
  result[e ][e ] = 0;
  result[mu][mu] = 0;
  result[e ][mu] = result[mu][e] = (eval - muval) / (4.*sin2thetaW);
  return result;
}
MATRIX<complex<double>,2,2> blocking_term0(const MATRIX<double,2,2>& Phi0matrix,
					   const MATRIX<complex<double>,2,2>& f,
					   const MATRIX<complex<double>,2,2>& fp){
  MATRIX<complex<double>,2,2> result;
  for(flavour fa=e; fa<=mu; fa++)
    for(flavour fb=e; fb<=mu; fb++){
      result[fa][fb] = 0;
      for(flavour fc=e; fc<=mu; fc++){
	result[fa][fb] += 0.25 * (Phi0matrix[fc][fb]*f[fa][fc]*fp[fc][fb] + Phi0matrix[fa][fc]*fp[fa][fc]*f[fc][fb]);
      }
    }
  return result;
}


//=============//
// MY_INTERACT //
//=============//
vector<vector<MATRIX<complex<double>,NF,NF> > > my_interact
  (const vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
   const double rho, const double T, const double Ye, const EAS& eas){

  const unsigned NE = fmatrixf[0].size();
  vector<vector<MATRIX<complex<double>,NF,NF> > > dfdr(2, vector<MATRIX<complex<double>,NF,NF> >(NE));

  // don't do anything if too sparse
  if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    return dfdr;

  for(state m=matter; m<=antimatter; m++){
    // get nulib species indices
    const int se = (m==matter ? 0 : 1);
    const int sx = (m==matter ? 2 : 3);
    const state mbar = (m==matter ? antimatter : matter);

#pragma omp parallel for
    for(int i=0; i<NE; i++){
      // temporary variables
      MATRIX<complex<double>,NF,NF> Pi_plus, Pi_minus, tmp;

      // make sure dfdr = 0
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++)
	  dfdr[m][i][f1][f2] = 0;
      
      // emission
      dfdr[m][i][e ][e ] += eas.emis(se,i);
      dfdr[m][i][mu][mu] += eas.emis(sx,i);

      // absorption. kappa_abs is <kappa> for absorption
      MATRIX<double,NF,NF> kappa_abs    = avg_matrix(eas.abs(se,i), eas.abs(sx,i));
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++)
	  dfdr[m][i][f1][f2] -= kappa_abs[f1][f2] * fmatrixf[m][i][f1][f2];
    
      // scattering and pair annihilation
      for(int j=0; j<NE; j++){
	
	// reusable variables
	double Phi0e, Phi0x, conv_to_in_rate;
	complex<double> unblock_in, unblock_out;
	MATRIX<double,NF,NF> Phi0avg, Phi0tilde, Phi0;
	MATRIX<complex<double>,NF,NF> block;
      
	// scattering rate from group i to group j
	Phi0e = eas.Phi0scat(se,i,j); // cm^3/s/sr
	Phi0x = eas.Phi0scat(sx,i,j);
	assert(Phi0e >= Phi0x);
	Phi0avg   =   avg_matrix(Phi0e, Phi0x);
	Phi0tilde = tilde_matrix(Phi0e, Phi0x);
	Phi0 = Phi0avg - Phi0tilde;
	block = blocking_term0(Phi0, fmatrixf[m][i], fmatrixf[m][j]);
	conv_to_in_rate = exp((eas.E[j]-eas.E[i])/(eas.temperature*1e6*cgs::units::eV));
	for(flavour f1=e; f1<=mu; f1++){
	  for(flavour f2=e; f2<=mu; f2++){
	    unblock_in  = fmatrixf[m][j][f1][f2] * 0.5*Phi0   [f1][f2]*conv_to_in_rate;
	    unblock_out = fmatrixf[m][i][f1][f2] * 0.5*Phi0avg[f1][f2];
	    dfdr[m][i][f1][f2] += 4.*M_PI*eas.nu[j]*eas.nu[j]*eas.dnu[j]/cgs::constants::c4 *
	      (unblock_in - unblock_out - block[f1][f2]*(conv_to_in_rate - 1.));
	  }
	}
      
	// annihilation from group i and anti group j
	if(eas.do_pair){
	  Phi0e = eas.Phi0pair(se,i,j); // cm^3/s/sr
	  Phi0x = eas.Phi0pair(sx,i,j);
	  assert(Phi0e >= Phi0x);
	  Phi0avg   =   avg_matrix(Phi0e, Phi0x);
	  Phi0tilde = tilde_matrix(Phi0e, Phi0x);
	  Phi0 = Phi0avg - Phi0tilde;
	  block = blocking_term0(Phi0, fmatrixf[m][i], fmatrixf[mbar][j]);
	  conv_to_in_rate = exp(-(eas.E[j]+eas.E[i])/(eas.temperature*1e6*cgs::units::eV));
	  for(flavour f1=e; f1<=mu; f1++){
	    for(flavour f2=e; f2<=mu; f2++){
	      unblock_in  = ((f1==f2 ? 1. : 0.) - fmatrixf[mbar][j][f1][f2])
		* 0.5*Phi0   [f1][f2]*conv_to_in_rate;
	      unblock_out = fmatrixf[m][i][f1][f2]
		* 0.5*Phi0avg[f1][f2]*conv_to_in_rate;
	      dfdr[m][i][f1][f2] += 4.*M_PI*eas.nu[j]*eas.nu[j]*eas.dnu[j]/cgs::constants::c4 *
		(unblock_in - unblock_out + block[f1][f2]*(conv_to_in_rate - 1.));
	    }
	  }
	}
	
	// 4-neutrino scattering
	if(eas.do_nu4scat){
	  for(int j3=0; j3<NE; j3++){
	    int j2 = eas.nu4_bin2(i,j,j3);
	    if(j2<0 or j2>=NE) continue;

	    int index = eas.nu4_kernel_index(i,j,j3);
	    double kernel = __nulibtable_MOD_nulibtable_nu4scat[index] * 0.5;

	    tmp = (1.-fmatrixf[m][j]) * fmatrixf[m][j2];
	    Pi_minus += (Trace(tmp) + tmp) * (1.-fmatrixf[m][j3]) * kernel;
	    tmp = (1.-fmatrixf[m][j3]) * fmatrixf[m][j2];
	    Pi_minus += (Trace(tmp) + tmp) * (1.-fmatrixf[m][j]) * kernel;

	    tmp = fmatrixf[m][j]*(1.-fmatrixf[m][j2]);
	    Pi_plus  += (Trace(tmp) + tmp) * fmatrixf[m][j3] * kernel;
	    tmp = fmatrixf[m][j3]*(1.-fmatrixf[m][j2]);
	    Pi_plus  += (Trace(tmp) + tmp) * fmatrixf[m][j] * kernel;
	  }
	}
	
	// 4-neutrino pair processes
	if(eas.do_nu4pair){
	  for(int j3=0; j3<NE; j3++){
	    int j2 = eas.nu4_bin2(i,j,j3);
	    if(j2<0 or j2>=NE) continue;

	    int index = eas.nu4_kernel_index(i,j,j3);
	    double kernel = __nulibtable_MOD_nulibtable_nu4pair[index] * 0.5;
	    
	    tmp = fmatrixf[mbar][j2] * (1.-fmatrixf[mbar][j]);
	    Pi_minus += (Trace(tmp) + tmp) * (1.-fmatrixf[m][j3]) * kernel;
	    tmp = fmatrixf[mbar][j2] * (1.-fmatrixf[mbar][j3]);
	    Pi_minus += (Trace(tmp) + tmp) * (1.-fmatrixf[m][j]) * kernel;

	    tmp = (1.-fmatrixf[m][j3]) * (1.-fmatrixf[mbar][j]);
	    Pi_minus += (Trace(tmp) + tmp) * fmatrixf[mbar][j2] * kernel;
	    tmp = (1.-fmatrixf[m][j]) * (1.-fmatrixf[mbar][j3]);
	    Pi_minus += (Trace(tmp) + tmp) * fmatrixf[mbar][j2] * kernel;
	    
	    tmp = (1.-fmatrixf[mbar][j2])*fmatrixf[mbar][j];
	    Pi_plus += (Trace(tmp) + tmp) * fmatrixf[m][j3] * kernel;
	    tmp = (1.-fmatrixf[mbar][j2])*fmatrixf[mbar][j3];
	    Pi_plus += (Trace(tmp) + tmp) * fmatrixf[m][j] * kernel;

	    tmp = fmatrixf[m][j3]*fmatrixf[mbar][j];
	    Pi_plus += (Trace(tmp) + tmp) * (1.-fmatrixf[mbar][j2]) * kernel;
	    tmp = fmatrixf[m][j]*fmatrixf[mbar][j3];
	    Pi_plus += (Trace(tmp) + tmp) * (1.-fmatrixf[mbar][j2]) * kernel;
	  }
	}
	  
	assert(j<NE);
      }

      dfdr[m][i] += Pi_plus *(1.-fmatrixf[m][i]) + (1.-fmatrixf[m][i])*Pi_plus ;
      dfdr[m][i] -= Pi_minus*    fmatrixf[m][i]  +     fmatrixf[m][i] *Pi_minus;
    } // end loop over i
  }
  return dfdr;
}
