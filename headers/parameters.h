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

const int NE=25;

const int NM=2;
enum state { matter,antimatter};
state operator++(state &n,int){ state tmp=n; n=(state)( (int)n+1 ); return tmp;};

const int NF=2;
enum flavour { e, mu };
flavour operator++(flavour &n,int){ flavour tmp=n; n=(flavour)( (int)n+1 ); return tmp;};

// the problem is broken up into seperate 'solutions'
const int NS=2; 
enum solution { msw, si};
solution operator++(solution &n,int){ solution tmp=n; n=(solution)( (int)n+1 ); return tmp;};

// number of parametrs needed to describe neutrino S matrix
const int NY=6; 
const array<array<double,NY>,NS> YIdentity{array<double,NY>{M_PI/2., M_PI/2., 0., 1., 0., 0.},array<double,NY>{M_PI/2., M_PI/2., 0., 1., 0., 0.}};

const hsize_t max_dims[6] = {H5S_UNLIMITED, NM, NE, NF, NF, 2};
const hsize_t chunk_dims[6] = {1, NM, NE, NF, NF, 2};

const double M_2PI = 2.*M_PI;
const complex<double> I(0.0,1.0);
// units, etc
namespace cgs{
  namespace units{
    const double cm = 1.; //cm
    const double eV = 1.60218e-12; //erg
  }
  namespace constants{
    const double hbar = 1.05457266e-27; // erg s
    const double c = 2.99792458e10; //cm/s
    const double c2 = c*c;
    const double c4 = c2*c2;
    const double hbarc = hbar*c;
    const double GF = 1.1663787e-5/*GeV^-2*//(1e9*1e9*units::eV*units::eV) * hbarc*hbarc*hbarc; //erg cm^3
    const double Mp = 1.6726219e-24; // g
  }
}

// mass of mass state1, delta m^2 differences
const double m1 = 0./*eV*/ * cgs::units::eV/cgs::constants::c2; // g, mass of state 1
const double dm21 = 2.43e-3/*eV^2*/ *cgs::units::eV*cgs::units::eV/cgs::constants::c4; // g^2, delta m2^2-m1^2
const double a1 = dm21>0 ? -1. : +1.; // sign of eigenvalue 1 of vacuum Hamiltonian
const double a2 = dm21>0 ? +1. : -1.; // sign of eigenvalue 2 of vacuum Hamiltonian

// mixing angles
const double theta12V = 9. * M_PI/180.; // rad
const double c12V = cos(theta12V);
const double s12V = sin(theta12V);
const double alphaV[NF] = {0,0};
const double betaV[NF-1] = {0};
const double sin2thetaW = 0.23122;
