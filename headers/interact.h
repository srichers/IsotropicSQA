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

#include "nulib_interface.h"

//============//
// Initialize //
//============//
void initialize(array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& fmatrixf,
		EAS& eas,
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

  for(int i=0; i<NE; i++){
    for(int m=matter; m<=antimatter; m++){
      double fe = eas.fermidirac(m, i);
      double fx = eas.fermidirac(2, i);
      double Tr = fe+fx;
      double lmax = min(Tr/2., 1.-Tr/2.);
      assert(lmax >= 0);
      double z = (fe-fx)/2.;
      double xmax = sqrt(lmax*lmax - z*z);
      assert(xmax==xmax);

      fmatrixf[m][i][e ][e ] = fe;
      fmatrixf[m][i][mu][mu] = fx;
      fmatrixf[m][i][e ][mu] = xmax*mixing;
      fmatrixf[m][i][mu][e ] = xmax*mixing;
    }

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
array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> my_interact
  (const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& fmatrixf,
   const double rho, const double T, const double Ye, const EAS& eas){

  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> dfdr;

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

	    if(abs(fmatrixf[m][j2][e][mu]) > 0) // j -> j2
	      assert( abs(fmatrixf[m][j2][e][mu] - conj(fmatrixf[m][j2][mu][e])) / abs(fmatrixf[m][j2][e][mu]) < 1e-5);

	    int index = eas.nu4_kernel_index(i,j,j3);
	    double kernel = __nulibtable_MOD_nulibtable_nu4scat[index];
	    index = eas.nu4_kernel_index(i,j3,j);
	    double other_kernel = __nulibtable_MOD_nulibtable_nu4scat[index];

	    tmp = (1.-fmatrixf[m][j]) * fmatrixf[m][j2];
	    Pi_minus += (Trace(tmp) + tmp) * (1.-fmatrixf[m][j3]) * kernel;

	    tmp = fmatrixf[m][j]*(1.-fmatrixf[m][j2]);
	    Pi_plus  += (Trace(tmp) + tmp) * fmatrixf[m][j3] * kernel;
	  }

		if(abs(Pi_plus[e][mu]) > 0) {
cout<<Pi_plus<< endl;
assert( abs(fmatrixf[m][j][e][mu] - conj(fmatrixf[m][j][mu][e])) / abs(fmatrixf[m][j][e][mu]) < 1e-5);
// assert( abs(fmatrixf[m][j2][e][mu] - conj(fmatrixf[m][j2][mu][e])) / abs(fmatrixf[m][j2][e][mu]) < 1e-5);
assert( abs(Pi_plus[e][mu] - conj(Pi_plus[mu][e])) / abs(Pi_plus[e][mu]) < 1e-5);
		}
	}

	// 4-neutrino pair processes
	if(eas.do_nu4pair){
	  for(int j3=0; j3<NE; j3++){
	    int j2 = eas.nu4_bin2(i,j,j3);
	    if(j2<0 or j2>=NE) continue;

	    int index = eas.nu4_kernel_index(i,j,j3);
	    double kernel = __nulibtable_MOD_nulibtable_nu4pair[index];

	    if(abs(fmatrixf[mbar][j][e][mu]) > 0)
	      assert( abs(fmatrixf[mbar][j][e][mu] - conj(fmatrixf[mbar][j][mu][e])) / abs(fmatrixf[mbar][j][e][mu]) < 1e-5);

	    tmp = fmatrixf[mbar][j2] * (1.-fmatrixf[mbar][j]);
	    Pi_minus += (Trace(tmp) + tmp) * (1.-fmatrixf[m][j3]) * kernel;

	    tmp = (1.-fmatrixf[m][j3]) * (1.-fmatrixf[mbar][j]);
	    Pi_minus += (Trace(tmp) + tmp) * fmatrixf[mbar][j2] * kernel;

	    tmp = (1.-fmatrixf[mbar][j2])*fmatrixf[mbar][j];
	    Pi_plus += (Trace(tmp) + tmp) * fmatrixf[m][j3] * kernel;

	    tmp = fmatrixf[m][j3]*fmatrixf[mbar][j];
	    Pi_plus += (Trace(tmp) + tmp) * (1.-fmatrixf[mbar][j2]) * kernel;
	  }

		if(abs(Pi_plus[e][mu]) > 0) {
cout<<Pi_plus<< endl;
assert( abs(Pi_plus[e][mu] - conj(Pi_plus[mu][e])) / abs(Pi_plus[e][mu]) < 1e-5);
		}
	}

	assert(j<NE);
      }

      if(abs(Pi_plus[e][mu]) > 0) {
	cout<<Pi_plus<< endl;
	assert( abs(Pi_plus[e][mu] - conj(Pi_plus[mu][e])) / abs(Pi_plus[e][mu]) < 1e-5);
			}
      if(abs(Pi_minus[e][mu]) > 0)
	assert( abs(Pi_minus[e][mu] - conj(Pi_minus[mu][e])) / abs(Pi_minus[e][mu]) < 1e-5);
      dfdr[m][i] += Pi_plus *(1.-fmatrixf[m][i]) + (1.-fmatrixf[m][i])*Pi_plus ;
      dfdr[m][i] -= Pi_minus*    fmatrixf[m][i]  +     fmatrixf[m][i] *Pi_minus;
    } // end loop over i
  }
  return dfdr;
}
