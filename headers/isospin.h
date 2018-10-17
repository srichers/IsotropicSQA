//Yonglin Zhu
//zhuygln@gmail.com
//Calculate isospin vector for neutrino and antineutrino in two flavor model.
//Ref[Phys.Rev. D94 (2016) no.10, 105006,Phys.Rev. D86 (2012) 085015]
#ifndef _H_ISOSPIN
#define _H_ISOSPIN

#include <cmath>

vector<double> isospin(MATRIX<complex<double>,NF,NF> Sf, MATRIX<complex<double>,NF,NF> Sfbar){
  //int NE=E.size();
  vector<double> s(6);
  //vector<double> delta(NE);
  s[0] = real(Sf[e][e])*real(Sf[mu][e])+imag(Sf[e][e])*imag(Sf[mu][e]);
  s[1] = real(Sf[e][e])*imag(Sf[mu][e])-imag(Sf[e][e])*real(Sf[mu][e]);
  s[2] = (norm(Sf[e][e])-(pow(real(Sf[mu][e]),2)+pow(imag(Sf[e][mu]),2)))/2;
  s[3] =-(real(Sfbar[e][e])*real(Sfbar[mu][e])+imag(Sfbar[e][e])*imag(Sfbar[mu][e]));
  s[4] =  real(Sfbar[e][e])*imag(Sfbar[mu][e])-imag(Sfbar[e][e])*real(Sfbar[mu][e]);
  s[5] = -(norm(Sfbar[e][e])-(pow(real(Sfbar[mu][e]),2)+pow(imag(Sfbar[e][mu]),2)))/2;
  
  return s;
}

complex<double> pauli[4][2][2] = {
  { {0, 1}, {1, 0} }, // x
  { {0,-I}, {I, 0} }, // y
  { {1, 0}, {0,-1} }, // z
  { {1, 0}, {0, 1} }  // t
};

void pauli_decompose(const MATRIX<complex<double>,2,2>& M, double coefficients[4]){
  //assert(M[1][0] == conj(M[0][1]));
  //assert(imag(M[0][0]) == 0);
  //assert(imag(M[1][1]) == 0);

  coefficients[0] = 0.5 * (real(M[1][0]) + real(M[0][1]));
  coefficients[1] = 0.5 * (imag(M[1][0]) - imag(M[0][1]));
  coefficients[2] = 0.5 * (real(M[0][0]) - real(M[1][1]));
  coefficients[3] = 0.5 * (real(M[0][0]) + real(M[1][1]));
}

template<typename T>
void pauli_reconstruct(const T coefficients[4], MATRIX<complex<double>,2,2>& M){
  for(unsigned i=0; i<2; i++){
    for(unsigned j=0; j<2; j++){
      M[i][j] = 0;
      for(unsigned k=0; k<4; k++)
	M[i][j] += coefficients[k] * pauli[k][i][j];
    }
  }
}

void Hermitize(MATRIX<complex<double>,2,2>& M, double accuracy){
  double trace = real(M[e][e] + M[mu][mu]);
  assert(trace>0);

  // matching off-diagonals
  double error = abs(M[e][mu] - conj(M[mu][e]));
  if(error/trace >= accuracy){
    //cout << M << endl;
    assert(error/trace < accuracy);
  }
  complex<double> tmp = 0.5 * (M[e][mu] + conj(M[mu][e]));
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

void unitarize(MATRIX<complex<double>,2,2>& M, double accuracy){
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
#endif
