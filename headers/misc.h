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

inline double Ve(const double rho, const double Ye){
  return (M_SQRT2*cgs::constants::GF/cgs::constants::Mp)*rho*Ye;
}
inline double Vmu(const double rho, const double Ye){ return 0.;}

class State{
 public:
  double rho, T, Ye;
  double r, dr_block, dr_osc, dr_int;
  int counter;
  EAS eas;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> fmatrixf;

  // temporaries
  array<array<double,NF>,NE> kV;
  array<MATRIX<complex<double>,NF,NF>,NM> UV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> HfV;
  array<array<MATRIX<complex<double>,NF,NF>,NF>,NE> CV;
  array<array<array<double,NF>,NF>,NE> AV;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> U0; // mixing angles to MSW basis at initial point
  array<MATRIX<complex<double>,NF,NF>,NM> VfMSW;
  array<array< array<double,NF>,NE>,NM> k0;
  array<array<MATRIX<complex<double>,NF,NF>,NE>,NM> UWBW;
  array<array<array<MATRIX<complex<double>,NF,NF>,NS>,NE>,NM> Sa;

  
  State(string nulibfilename, string eosfilename, double rho_in, double Ye_in, double T_in, double dr0, double mixing, bool do_interact){
    r=0;
    rho = rho_in;
    T = T_in;
    Ye = Ye_in;
    eas = EAS(nulibfilename, eosfilename);
    initialize(fmatrixf,eas,rho,T,Ye, mixing, do_interact);
    dr_block = dr0;
    dr_osc = dr0;
    dr_int = dr0;
    counter = 0;
  
    // vectors of energies and vacuum eigenvalues
    kV  = set_kV(eas.E);
    UV  = Evaluate_UV();
    HfV = Evaluate_HfV(kV,UV);
    CV  = Evaluate_CV(kV, HfV);
    AV  = Evaluate_AV(kV,HfV,UV);
    
    // MSW potential matrix
    VfMSW[matter][e][e]=Ve(rho,Ye);
    VfMSW[matter][mu][mu]=Vmu(rho,Ye);
    VfMSW[antimatter]=-Conjugate(VfMSW[matter]);
    
    // other matrices
    for(int m=matter; m<=antimatter; m++){
      for(int i=0;i<NE;i++){
	MATRIX<complex<double>,NF,NF> Hf0=HfV[m][i]+ VfMSW[m];
	k0[m][i] = (m==matter? k(Hf0) : kbar(Hf0) );
	array<double,1> deltak0 = (m==matter ? deltak(Hf0) : deltakbar(Hf0) );
	array<MATRIX<complex<double>,NF,NF>,NF> C0 = CofactorMatrices(Hf0,k0[m][i]);
	array<array<double,NF>,NF> A0;
	for(int j=0;j<=NF-1;j++){
	  if(real(C0[j][mu][e]*CV[i][j][mu][e]) < 0.)
	    A0[j][e]=-AV[i][j][e];
	  else A0[j][e]=AV[i][j][e];
	  A0[j][mu]=AV[i][j][mu];
	}
	U0[m][i]=U(deltak0,C0,A0);
	if(m==antimatter) U0[m][i] = Conjugate(U0[m][i]);
      }
    }
  }

  // IO temporaries
  hid_t dset_f, dset_r, dset_dr_osc, dset_dr_int, dset_dr_block;
};

void getP(const State& s, array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& pmatrixm0){
  #pragma omp parallel for collapse(2)
  for(int m=matter; m<=antimatter; m++){
    for(int i=0; i<NE; i++){
      pmatrixm0[m][i] = Adjoint(s.U0[m][i]) * s.fmatrixf[m][i] * s.U0[m][i]
	* (M_SQRT2*cgs::constants::GF /* erg cm^3*/
	   * 4.*M_PI*s.eas.nu[i]*s.eas.nu[i]*s.eas.dnu[i]/*Hz^3*/
	   / pow(cgs::constants::c,3)/*cm^3 Hz^3*/);
    }
  }
}

template<typename T>
T get_parameter(ifstream& fin, const char* name){
  stringstream line;
  string linestring;
  T to_set;
  std::getline(fin, linestring);
  line = stringstream(linestring);
  line >> to_set;
  cout << "PARAMETER " << name << " = " << to_set << endl;
  return to_set;
}

inline double Sign(const double input){
  return input>0 ? 1 : -1;
}

// Cash-Karp RK parameters
const int NRK=6;
const int NRKOrder=5;
static const double AA[]={ 0., 1./5., 3./10., 3./5., 1., 7./8. };
static const double b0[]={};
static const double b1[]={ 1./5. };
static const double b2[]={ 3./40.,9./40. };
static const double b3[]={ 3./10.,-9./10.,6./5. };
static const double b4[]={ -11./54.,5./2.,-70./27.,35./27. };
static const double b5[]={ 1631./55296.,175./512.,575./13824.,44275./110592.,253./4096. };
static const double* BB[]={ b0,b1,b2,b3,b4,b5 };
static const double CC[]={ 37./378.,0.,250./621.,125./594.,0.,512./1771. };
static const double DD[]={ 2825./27648.,0.,18575./48384.,13525./55296.,277./14336.,1./4. };

//===//
// B //
//===//
MATRIX<complex<double>,NF,NF> B(const array<double,NY>& y){
  MATRIX<complex<double>,NF,NF> s;
  double cPsi1=cos(y[0]),sPsi1=sin(y[0]), cPsi2=cos(y[1]),sPsi2=sin(y[1]), cPsi3=cos(y[2]),sPsi3=sin(y[2]);
  
  s[0][1] = cPsi1 + I*sPsi1*cPsi2;
  sPsi1 *= sPsi2;
  s[0][0] = sPsi1 * (cPsi3 + I*sPsi3);

  s[1][0] = -y[3]*conj(s[0][1]);
  s[1][1] =  y[3]*conj(s[0][0]);

  return s;
}

//===//
// W //
//===//
MATRIX<complex<double>,NF,NF> W(const array<double,NY>& Y){
  MATRIX<complex<double>,NF,NF> w;
  w[0][0]=exp(-I*M_2PI*Y[4]); w[1][1]=exp(-I*M_2PI*Y[5]);
  return w;
}

//===//
// K //
//===//
void K(const double dr,
       State& s,
       const array<array<MATRIX<complex<double>,NF,NF>,NE>,NM>& pmatrixm0,
       const array<array<array<array<double,NY>,NS>,NE>,NM> &Y,
       array<array<array<array<double,NY>,NS>,NE>,NM> &K){

  array<MATRIX<complex<double>,NF,NF>,NM> VfSI;  // self-interaction potential

#pragma omp parallel
  {
  #pragma omp for collapse(2)
  for(int m=matter; m<=antimatter; m++){
    for(int i=0; i<NE; i++){
      MATRIX<complex<double>,NF,NF> BB  = B(Y[m][i][msw]);
      s.Sa[m][i][si] = B(Y[m][i][si]);
      s.UWBW[m][i] = s.U0[m][i] * W(Y[m][i][msw]) * BB * W(Y[m][i][si]);

      for(int j=0;j<=3;j++)
	K[m][i][msw][j] = 0.;
      K[m][i][msw][4] = s.k0[m][i][0]*dr/M_2PI/cgs::constants::hbarc;
      K[m][i][msw][5] = s.k0[m][i][1]*dr/M_2PI/cgs::constants::hbarc;

      // contribution to the self-interaction potential from this energy *
      MATRIX<complex<double>,NF,NF> Sfm = s.UWBW[m][i] * s.Sa[m][i][si];
      MATRIX<complex<double>,NF,NF> VfSIE = Sfm*pmatrixm0[m][i]*Adjoint(Sfm);
      if(m==antimatter) VfSIE = -Conjugate(VfSIE);
      #pragma omp critical
      VfSI[matter] += VfSIE;
    }
  }

  #pragma omp single
  VfSI[antimatter]=-Conjugate(VfSI[matter]);

  // *********************
  // SI part of solution *
  // *********************
  #pragma omp for collapse(2)
  for(int m=matter; m<=antimatter; m++){
    for(int i=0; i<NE; i++){

      MATRIX<complex<double>,NF,NF> Ha = Adjoint(s.UWBW[m][i]) * VfSI[m] * s.UWBW[m][i];

      K[m][i][si][4]=dr*real(Ha[0][0])/(M_2PI*cgs::constants::hbarc);
      K[m][i][si][5]=dr*real(Ha[1][1])/(M_2PI*cgs::constants::hbarc);

      complex<double> HB0 =-I/cgs::constants::hbarc*( Ha[0][1]*s.Sa[m][i][si][1][0] );
      complex<double> HB1 =-I/cgs::constants::hbarc*( Ha[0][1]*s.Sa[m][i][si][1][1] );

      double dvdr[4];
      dvdr[0]=real(HB1);
      dvdr[1]=imag(HB1);
      dvdr[2]=real(HB0);
      dvdr[3]=imag(HB0);

      MATRIX<double,3,4> JI = JInverse(Y[m][i][si]);

      for(int j=0;j<=2;j++){
	K[m][i][si][j]=0.;
	for(int k=j;k<=3;k++) K[m][i][si][j]+=JI[j][k]*dvdr[k];
	K[m][i][si][j]*=dr;
      }

      K[m][i][si][3]=0.;
    }
  }
  }// omp parallel block

}// end of K function


//============//
// setup_file //
//============//
hid_t setup_file(string filename, State& s){
  hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t file_space;
  hsize_t ndims;
  hsize_t dims[6] = {0, NM, NE, NF, NF, 2};
  hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_layout(plist, H5D_CHUNKED);

  // FMATRIXF //
  ndims = 6;
  file_space = H5Screate_simple(ndims, dims, max_dims);
  H5Pset_chunk(plist, ndims, chunk_dims);
  H5Dcreate(file, "fmatrixf", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  s.dset_f = H5Dopen(file, "fmatrixf", H5P_DEFAULT);

  // RADIUS/TIME //
  ndims = 1;
  file_space = H5Screate_simple(ndims, dims, max_dims);
  H5Pset_chunk(plist, ndims, chunk_dims);
  H5Dcreate(file, "r(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(file, "dr_block(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(file, "dr_osc(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Dcreate(file, "dr_int(cm)", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
  s.dset_r        = H5Dopen(file, "r(cm)", H5P_DEFAULT);
  s.dset_dr_osc   = H5Dopen(file, "dr_osc(cm)", H5P_DEFAULT);
  s.dset_dr_int   = H5Dopen(file, "dr_int(cm)", H5P_DEFAULT);
  s.dset_dr_block = H5Dopen(file, "dr_block(cm)", H5P_DEFAULT);

  // clear resources
  H5Sclose(file_space);
  H5Pclose(plist);
  
  return file;
}

//============//
// write_data //
//============//
void write_data(const hid_t file, State& s, const double impact){
  // output to stdout
  double n=0, nbar=0;
  double coeff = 4.*M_PI / pow(cgs::constants::c,3);
  for(int i=0; i<NE; i++){
    double dnu3 = s.eas.nu[i]*s.eas.nu[i]*s.eas.dnu[i];
    for(flavour f1=e; f1<=mu; f1++){
      n    += real(s.fmatrixf[    matter][i][f1][f1]) * dnu3 * coeff;
      nbar += real(s.fmatrixf[antimatter][i][f1][f1]) * dnu3 * coeff;
    }
  }
  cout << s.counter << "\t";
  cout << s.r/cgs::constants::c << "\t";
  cout << s.dr_osc/cgs::constants::c << "\t";
  cout << s.dr_int/cgs::constants::c << "\t";
  cout << s.dr_block/cgs::constants::c << "\t";
  cout << n << "\t" << nbar << "\t" << (n-nbar) << "\t";
  cout << impact << endl;
  cout.flush();


  
  // output to file //

  hid_t mem_space, file_space;
  hsize_t ndims;

  // create the memory space
  ndims = 6;
  mem_space = H5Screate_simple(ndims, chunk_dims, NULL);
  file_space = H5Dget_space (s.dset_f);
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(file_space, dims, NULL);
  dims[0]++;
  hsize_t start[6] = {dims[0]-1, 0, 0, 0, 0, 0};

  // fmatrixf
  H5Dset_extent(s.dset_f, dims);
  file_space = H5Dget_space(s.dset_f);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dwrite(s.dset_f, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.fmatrixf);

  // 1D stuff
  ndims = 1;
  mem_space = H5Screate_simple(ndims, chunk_dims, NULL);
  
  // r
  H5Dset_extent(s.dset_r, dims);
  file_space = H5Dget_space(s.dset_r);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dwrite(s.dset_r, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.r);

  // dr_osc
  H5Dset_extent(s.dset_dr_osc, dims);
  file_space = H5Dget_space(s.dset_dr_osc);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dwrite(s.dset_dr_osc, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dr_osc);

  // dr_int
  H5Dset_extent(s.dset_dr_int, dims);
  file_space = H5Dget_space(s.dset_dr_int);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dwrite(s.dset_dr_int, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dr_int);

  // dr_block
  H5Dset_extent(s.dset_dr_block, dims);
  file_space = H5Dget_space(s.dset_dr_block);
  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dwrite(s.dset_dr_block, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dr_block);

  // free resources
  H5Sclose(file_space);
  H5Sclose(mem_space);
  H5Fflush(file,H5F_SCOPE_LOCAL);
}

//=========//
// recover //
//=========//
hid_t recover(const string filename, State& s){
  hsize_t start[6] = {0, 0, 0, 0, 0, 0};
  hsize_t dims[6];
  hsize_t ndims;
  hid_t file_space, mem_space;

  hid_t file = H5Fopen (filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  // fmatrixf
  s.dset_f = H5Dopen (file, "fmatrixf", H5P_DEFAULT);
  ndims=6;
  mem_space = H5Screate_simple(ndims, chunk_dims, NULL);
  file_space = H5Dget_space (s.dset_f);
  H5Sget_simple_extent_dims(file_space, dims, NULL);
  start[0] = dims[0]-1;
  assert(dims[1]==NM);
  assert(dims[2]==NE);
  assert(dims[3]==NF);
  assert(dims[4]==NF);
  assert(dims[5]==2);
  H5Sselect_hyperslab (file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dread (s.dset_f, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.fmatrixf);

  // 1D stuff
  ndims=1;
  mem_space = H5Screate_simple(ndims, chunk_dims, NULL);
  
  // r
  s.dset_r = H5Dopen (file, "r(cm)", H5P_DEFAULT);
  file_space = H5Dget_space (s.dset_r);
  H5Sselect_hyperslab (file_space, H5S_SELECT_SET, start, NULL, chunk_dims, NULL);
  H5Dread (s.dset_r, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.r);

  // dr_osc
  s.dset_dr_osc = H5Dopen (file, "dr_osc(cm)", H5P_DEFAULT);
  H5Dread (s.dset_dr_osc, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dr_osc);
  
  // dr_int
  s.dset_dr_int = H5Dopen (file, "dr_int(cm)", H5P_DEFAULT);
  H5Dread (s.dset_dr_int, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dr_int);
  
  // dr_block
  s.dset_dr_block = H5Dopen (file, "dr_block(cm)", H5P_DEFAULT);
  H5Dread (s.dset_dr_block, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &s.dr_block);

  s.counter = dims[0]-1;
  
  // clear resources
  H5Sclose(file_space);
  H5Sclose(mem_space);

  return file;
}


void Hermitize(MATRIX<complex<double>,2,2>& M, const double accuracy){
  double trace = real(M[e][e] + M[mu][mu]);
  assert(trace>0);

  // matching off-diagonals
  double error = abs(M[e][mu] - conj(M[mu][e]));
  if(error/trace >= accuracy){
    //cout << M << endl;
    assert(error/trace < accuracy);
  }
  complex<double> tmp = 0.5 * (M[mu][e] + conj(M[e][mu]));
  M[mu][e] = tmp;
  M[e][mu] = conj(tmp);
  
  // real on-diagonals
  for(flavour f1=e; f1<=mu; f1++){
    error = abs(imag(M[f1][f1]));
    assert(error/trace < accuracy);
    M[f1][f1] = real(M[f1][f1]);
  }

  // density matrix probabilities
  double ad = abs(M[e][e ] * M[mu][mu]) / trace;
  double bc = abs(M[e][mu] * M[mu][e ]) / trace;
  double det = ad - bc;
  if(det < 0 && bc>0){
    assert( abs(det) < accuracy);
    M[e][mu] *= sqrt(ad/bc);
    M[mu][e] *= sqrt(ad/bc);
  }
}

void unitarize(MATRIX<complex<double>,2,2>& M, const double accuracy){
  // M = ( (a, b), (-e^Iphi b*, e^Iphi a*) )
  //   = ( (a, b), (c, d) )
  double a2 = real(M[e ][e ] * conj(M[e ][e ]) );
  double b2 = real(M[e ][mu] * conj(M[e ][mu]) );
  double c2 = real(M[mu][e ] * conj(M[mu][e ]) );
  double d2 = real(M[mu][mu] * conj(M[mu][mu]) );

  // aa* < 1
  if(a2 > 1.){
    assert( abs(a2-1.) < accuracy);
    M[e][e] /= sqrt(a2);
    a2 = real(M[e ][e ] * conj(M[e ][e ]) );
  }
  
  // aa* + bb* = 1
  assert( abs(a2 + b2 - 1.) < accuracy);
  M[e][mu] *= sqrt( max(0., 1.-a2) / b2 );
  b2 = real(M[e ][mu] * conj(M[e ][mu]) );
  
  // aa* = dd*
  assert( abs(a2-d2) < accuracy );
  M[mu][mu] *= sqrt(a2/d2);
  d2 = real(M[mu][mu] * conj(M[mu][mu]) );

  // bb* = cc*
  assert( abs(b2-c2) < accuracy );
  M[mu][e] *= sqrt(b2/c2);
  c2 = real(M[mu][e ] * conj(M[mu][e ]) );
  
  // det(M) = e^Iphi
  complex<double> newval;
  complex<double> eIphi = M[e][e]*M[mu][mu] - M[e][mu]*M[mu][e];
  double eIphiMag = sqrt(real( eIphi*conj(eIphi) ));
  assert( abs(eIphiMag - 1.) < accuracy);
  eIphi /= eIphiMag;
  newval = eIphi * conj(M[e][e ]);
  assert( abs(M[mu][mu]-newval) < accuracy);
  M[mu][mu] = newval;
  newval = -eIphi * conj(M[e][mu]);
  assert( abs(M[mu][e ]-newval) < accuracy);
  M[mu][e ] = newval;

  // crazy sanity check
  MATRIX<complex<double>,2,2> identity;
  identity = M*Adjoint(M);
  assert( abs( real(identity[e ][e ]*conj(identity[e ][e ])) - 1.) < accuracy);
  assert( abs( real(identity[mu][mu]*conj(identity[mu][mu])) - 1.) < accuracy);
  assert( abs( real(identity[e ][mu]*conj(identity[e ][mu]))     ) < accuracy);
  assert( abs( real(identity[mu][e ]*conj(identity[mu][e ]))     ) < accuracy);
}

double uniform(){
  return (float)rand() / (float)RAND_MAX;
}
double exponential_random(){
  return -log(uniform());
}

