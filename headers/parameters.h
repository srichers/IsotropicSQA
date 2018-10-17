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

const int NM=2;
enum state { matter,antimatter};
state operator++(state &n,int){ state tmp=n; n=(state)( (int)n+1 ); return tmp;};

const int NF=2;
enum flavour { e, mu };
flavour operator++(flavour &n,int){ flavour tmp=n; n=(flavour)( (int)n+1 ); return tmp;};

// number of parametrs needed to describe neutrino S matrix
const int NY=6; 

// number of nuclei followed 
const int NN=1; 

// mass of mass state1, delta m^2 differences
double m1,dm21;

// number of energy bins, min and max energy
int NE;
double Emin, Emax;
vector<double> E, nu;
vector<double> dnu;

// the problem is broken up into seperate 'solutions'
const int NS=2; 
enum solution { msw, si};
solution operator++(solution &n,int){ solution tmp=n; n=(solution)( (int)n+1 ); return tmp;};

// vacuum eigenvalues
vector<vector<double> > kV;
int a1,a2;

// vacuum Hamiltonian and mixing matrices
vector<vector<MATRIX<complex<double>,NF,NF> > > HfV(NM);
vector<MATRIX<complex<double>,NF,NF> > UV(NM);

// vacuum values of the off-diagonal elements of the cofactor matrices
vector<vector<MATRIX<complex<double>,NF,NF> > > CV;

// mixing matrix element prefactors
vector<vector<vector<double> > > AV;

// initial mixing matrices, needs to be global
vector<vector<MATRIX<complex<double>,NF,NF> > > U0(NM); 

double theta12V;
vector<double> alphaV(NF), betaV(NF-1);
double c12V,s12V;

// neutron star radius
double Rnu;

// time snapshot
double t;

// vectors of mean energies, luminosities, temperatures, pinch paramaters
vector<vector<double> > meanE(NM,vector<double>(NF));
vector<vector<double> > L(NM,vector<double>(NF));
vector<vector<double> > pinch(NM,vector<double>(NF));

// pointers to initial spectra functions in found flux.h
vector<vector<double(*)(double)> > F0(NM,vector<double(*)(double)>(NF));
//



double M_2PI = 2.*M_PI;
complex<double> I = 1i;
// units, etc
namespace cgs{
  namespace units{
    double cm = 1.; //cm
    double eV = 1.60218e-12; //erg
  }
  namespace constants{
    double hbar = 1.05457266e-27; // erg s
    double c = 2.99792458e10; //cm/s
    double c2 = c*c;
    double c4 = c2*c2;
    double hbarc = hbar*c;
    double GF = 1.1663787e-5/*GeV^-2*//(1e9*1e9*units::eV*units::eV) * hbarc*hbarc*hbarc; //erg cm^3
    double Mp = 1.6726219e-24; // g
  }
}
