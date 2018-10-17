#include "nulib_interface.h"
#include "isospin.h"

//=======//
// Ebins //
//=======//
void set_Ebins(vector<double>& E){
  NE = __nulibtable_MOD_nulibtable_number_groups;
  assert(NE > 0);
  E.resize(NE);
  nu.resize(NE);
  dnu.resize(NE);
  cout << endl;
  cout<<"NE="<<NE << endl;
  for(int i=0;i<NE;i++){
    E[i]            = __nulibtable_MOD_nulibtable_energies[i]*1e6*cgs::units::eV; // erg
    nu[i]           = E[i] / (2.*M_PI*cgs::constants::hbar); // Hz
    double nubottom = __nulibtable_MOD_nulibtable_ebottom [i]*1e6*cgs::units::eV / (2.*M_PI*cgs::constants::hbar); // Hz
    double nutop    = __nulibtable_MOD_nulibtable_etop    [i]*1e6*cgs::units::eV / (2.*M_PI*cgs::constants::hbar); // Hz
    dnu[i]          =     nutop    -     nubottom   ;
    cout << E[i]/(1.e6*cgs::units::eV) << " ";
  }
  cout.flush();
}

//============//
// Initialize //
//============//
void initialize(vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
		double r, double rho, double T, double Ye, double mixing, int do_interact){
  // T should be MeV
  cout << "Setting initial data." << endl;
  cout << "rho = " << rho << " g/ccm" << endl;
  cout << "T = " << T << " MeV" << endl;
  cout << "Ye = " << Ye << endl;
  eas.set(rho,T,Ye,do_interact);
  
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
    cout << "\teas.emis = {" << eas.emis(0,i) << ", " << eas.emis(1,i) << ", " << eas.emis(2,i) << "}" << endl;
    cout << "\teas.abs = {" << eas.abs(0,i) << ", " << eas.abs(1,i) << ", " << eas.abs(2,i) << "}" << endl;
    cout << "\tBB = {" << eas.emis(0,i)/eas.abs(0,i) << ", " << eas.emis(1,i)/eas.abs(1,i) << ", " << eas.emis(2,i)/eas.abs(2,i) << "}" << endl;
    cout << "\teas.scat = {" << eas.scat(0,i) << ", " << eas.scat(1,i) << ", " << eas.scat(2,i) << "}" << endl;

    cout << "\tf = {" << real(fmatrixf[matter][i][e][e]) << ", " << real(fmatrixf[antimatter][i][e][e]) << ", " << real(fmatrixf[matter][i][mu][mu]) << ", " << real(fmatrixf[antimatter][i][mu][mu]) << "}" << endl;
  }
  
  for(state m=matter; m<=antimatter; m++)
    for(int i=0; i<NE; i++)
      for(flavour f1=e; f1<=mu; f1++)
	for(flavour f2=e; f2<=mu; f2++){
	  assert(fmatrixf[m][i][f1][f2] == fmatrixf[m][i][f1][f2]);
	  assert(abs(fmatrixf[m][i][f1][f2]) < 1.0);
	}
}

//===================//
// Vacuum Potentials //
//===================//
double deltaV(const double E){ // erg
  return abs(dm21)*cgs::constants::c4 / (2.*E);
}

void set_kV(vector<vector<double> >& kV){
  assert(NF==2);
  for(int i=0;i<NE;i++){
    kV[i][0] = m1*m1 * cgs::constants::c4 /2./E[i];
    kV[i][1] = kV[i][0] + deltaV(E[i]);
  }
}

//=============================//
// Self-Interaction Potentials //
//=============================//
MATRIX<double,2,2> avg_matrix(double eval, double muval){
  MATRIX<double,2,2> result;
  result[e ][e ] = eval;
  result[mu][mu] = muval;
  result[e ][mu] = result[mu][e] = (eval + muval) / 2.;
  return result;
}
MATRIX<double,2,2> tilde_matrix(double eval, double muval){
  const double sin2thetaW = 0.23122;
  MATRIX<double,2,2> result;
  result[e ][e ] = 0;
  result[mu][mu] = 0;
  result[e ][mu] = result[mu][e] = (eval - muval) / (4.*sin2thetaW);
  return result;
}
MATRIX<complex<double>,2,2> blocking_term0(MATRIX<double,2,2> Phi0matrix, MATRIX<complex<double>,2,2> f, MATRIX<complex<double>,2,2> fp){
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

vector<vector<MATRIX<complex<double>,NF,NF> > > my_interact
  (vector<vector<MATRIX<complex<double>,NF,NF> > >& fmatrixf,
   double rho, double T, double Ye){

  vector<vector<MATRIX<complex<double>,NF,NF> > > dfdr(2, vector<MATRIX<complex<double>,NF,NF> >(NE));

  // don't do anything if too sparse
  if(log10(rho) <= __nulibtable_MOD_nulibtable_logrho_min)
    return dfdr;

  for(state m=matter; m<=antimatter; m++){
    // get nulib species indices
    int se = (m==matter ? 0 : 1);
    int sx = (m==matter ? 2 : 3);

#pragma omp parallel for
    for(int i=0; i<NE; i++){
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
	conv_to_in_rate = exp((E[j]-E[i])/(eas.temperature*1e6*cgs::units::eV));
	for(flavour f1=e; f1<=mu; f1++){
	  for(flavour f2=e; f2<=mu; f2++){
	    unblock_in  = fmatrixf[m][j][f1][f2] * 0.5*Phi0   [f1][f2]*conv_to_in_rate;
	    unblock_out = fmatrixf[m][i][f1][f2] * 0.5*Phi0avg[f1][f2];
	    dfdr[m][i][f1][f2] += 4.*M_PI*nu[j]*nu[j]*dnu[j]/cgs::constants::c4 *
	      (unblock_in - unblock_out - block[f1][f2]*(conv_to_in_rate - 1.));
	  }
	}
      
	// annihilation from group i and anti group j
	if(eas.do_pair){
	  state mbar = (m==matter ? antimatter : matter);
	  Phi0e = eas.Phi0pair(se,i,j); // cm^3/s/sr
	  Phi0x = eas.Phi0pair(sx,i,j);
	  assert(Phi0e >= Phi0x);
	  Phi0avg   =   avg_matrix(Phi0e, Phi0x);
	  Phi0tilde = tilde_matrix(Phi0e, Phi0x);
	  Phi0 = Phi0avg - Phi0tilde;
	  block = blocking_term0(Phi0, fmatrixf[m][i], fmatrixf[mbar][j]);
	  conv_to_in_rate = exp(-(E[j]+E[i])/(eas.temperature*1e6*cgs::units::eV));
	  for(flavour f1=e; f1<=mu; f1++){
	    for(flavour f2=e; f2<=mu; f2++){
	      unblock_in  = ((f1==f2 ? 1. : 0.) - fmatrixf[mbar][j][f1][f2])
		* 0.5*Phi0   [f1][f2]*conv_to_in_rate;
	      unblock_out = fmatrixf[m][i][f1][f2]
		* 0.5*Phi0avg[f1][f2]*conv_to_in_rate;
	      dfdr[m][i][f1][f2] += 4.*M_PI*nu[j]*nu[j]*dnu[j]/cgs::constants::c4 *
		(unblock_in - unblock_out + block[f1][f2]*(conv_to_in_rate - 1.));
	    }
	  }
	}
      
      }
    }
  }
  return dfdr;
}
