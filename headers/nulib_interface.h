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

#ifndef _NULIB_INTERFACE_H
#define _NULIB_INTERFACE_H

#include "H5Cpp.h"
#include <string>
#include <cstdlib>
//#include "mstl.h"

inline bool hdf5_dataset_exists(const char* filename, const char* datasetname){
  bool exists = true;

  // Temporarily turn off error printing
  H5E_auto2_t func;
  void* client_data;
  H5::Exception::getAutoPrint(func,&client_data);
  H5::Exception::dontPrint();

  // See if dataset exists
  H5::H5File file(filename, H5F_ACC_RDONLY);
  H5::DataSet dataset;
  try{
    dataset = file.openDataSet(datasetname);
  }
  catch(H5::FileIException& exception){
    exists = false;
  }

  // Turn error printing back on
  H5::Exception::setAutoPrint(func,client_data);
  file.close();

  return exists;
}

// module variables set in fortran NuLib code
extern int     __nulibtable_MOD_nulibtable_number_species;
extern int     __nulibtable_MOD_nulibtable_number_easvariables;
extern int     __nulibtable_MOD_nulibtable_number_groups;
extern int     __nulibtable_MOD_nulibtable_nrho;
extern int     __nulibtable_MOD_nulibtable_ntemp;
extern int     __nulibtable_MOD_nulibtable_nye;
extern int     __nulibtable_MOD_nulibtable_nitemp;
extern int     __nulibtable_MOD_nulibtable_nieta;
extern double* __nulibtable_MOD_nulibtable_energies;
extern double* __nulibtable_MOD_nulibtable_ewidths;
extern double* __nulibtable_MOD_nulibtable_ebottom;
extern double* __nulibtable_MOD_nulibtable_etop;
extern double* __nulibtable_MOD_nulibtable_logrho;
extern double* __nulibtable_MOD_nulibtable_logtemp;
extern double* __nulibtable_MOD_nulibtable_ye;
extern double* __nulibtable_MOD_nulibtable_logitemp;
extern double* __nulibtable_MOD_nulibtable_logieta;
extern double* __nulibtable_MOD_nulibtable_nu4scat;
extern double* __nulibtable_MOD_nulibtable_nu4pair;
extern double  __nulibtable_MOD_nulibtable_logtemp_min;
extern double  __nulibtable_MOD_nulibtable_logtemp_max;
extern double  __nulibtable_MOD_nulibtable_logrho_min;
extern double  __nulibtable_MOD_nulibtable_logrho_max;
extern double  __nulibtable_MOD_nulibtable_ye_min;
extern double  __nulibtable_MOD_nulibtable_ye_max;
extern double  __nulibtable_MOD_nulibtable_logitemp_min;
extern double  __nulibtable_MOD_nulibtable_logitemp_max;
extern double  __nulibtable_MOD_nulibtable_logieta_min;
extern double  __nulibtable_MOD_nulibtable_logieta_max;
extern int     __nulib_MOD_total_eos_variables;

// pointers to be freed after getting interaction rates
extern double* __nulibtable_MOD_nulibtable_emissivities;
extern double* __nulibtable_MOD_nulibtable_absopacity;
extern double* __nulibtable_MOD_nulibtable_scatopacity;
extern double* __nulibtable_MOD_nulibtable_delta;
extern double* __nulibtable_MOD_nulibtable_itable_phi0;
extern double* __nulibtable_MOD_nulibtable_itable_phi1;
extern double* __nulibtable_MOD_nulibtable_epannihiltable_phi0;
extern double* __nulibtable_MOD_nulibtable_epannihiltable_phi1;

// These are fortran functions and module variables in nulib.a
extern "C"{
  void nulibtable_range_species_range_energy_(
		  double*, //rho
		  double*, //temp
		  double*, //ye
		  double*, //eas_species_energy (3D array)
		  int*,    //number of species (3,5,6)
		  int*,    //number of groups
		  int*);   //number of easvariables (3)

  void nulibtable_single_species_range_energy_(
		  double*, //rho
		  double*, //temp
		  double*, //Ye
		  int*,    //species number
		  double*, //eas_energy (2D array)
		  int*,    //number of groups
		  int*);   //number of easvariables(3)

  void nulibtable_epannihil_single_species_range_energy_(
		  double* temp,  // MeV
		  double* eta,   // mu/kT
		  int* lns,      // species number
		  double* phi,   // phi[legendre-p/a index][this_group][anti-group]
		  int* ngroups1,
		  int* ngroups2,
		  int* n_phis);

  void nulibtable_inelastic_single_species_range_energy_(
		  double* temp,  // MeV
		  double* eta,   // mu/kT
		  int* lns,      // species number
		  double* phi,   // phi[legendre index][out group][in group]
		  int* ngroups1, // ng in
		  int* ngroups2, // ng out (should be same as eas_n1)
		  int* n_phis);   // number of legendre terms (=2)

  void nulibtable_reader_(char*,int*,int*,int*,int*,int*,int);
  void set_eos_variables_(double* eos_variables);
  void read_eos_table_(char* filename);
}

class EAS{
 public:
  int do_iscat, do_pair, do_delta, do_nu4scat, do_nu4pair;
  int ns, ng, nv;
  double munue = 0./0.; /*erg*/
  double temperature = 0./0.; /*MeV*/
  array<double,NE> E; // erg
  array<double,NE> nu, dnu; // Hz
  array<double, 4*NE*NM*NF> eas;
  array<double,NE*NE*NM*NF> escat_kernel0; // out-scattering
  array<double,NE*NE*NM*NF> escat_kernel1;
  array<double,NE*NE*NM*NF> pair_kernel0; // neutrino pair annihilation
  array<double,NE*NE*NM*NF> pair_kernel1;

  EAS(){}
  EAS(const string filename, const string eosfilename){
    // select what to read
    do_iscat = 0;
    do_pair = 0;
    do_delta = 0;
    do_nu4scat = 0;
    do_nu4pair = 0;
    if(hdf5_dataset_exists(filename.c_str(),"/scattering_delta")) do_delta = true;
    if(hdf5_dataset_exists(filename.c_str(),"/inelastic_phi0"))   do_iscat = true;
    if(hdf5_dataset_exists(filename.c_str(),"/epannihil_phi0"))   do_pair = true;
    if(hdf5_dataset_exists(filename.c_str(),"/nu4scat_kernel"))   do_nu4scat = true;
    if(hdf5_dataset_exists(filename.c_str(),"/nu4pair_kernel"))   do_nu4pair = true;
    cout << "TABLE ELEMENTS: " << do_delta << do_iscat << do_pair << do_nu4scat << do_nu4pair << endl;
    nulibtable_reader_((char*)filename.c_str(), &do_iscat, &do_pair, &do_delta, &do_nu4scat, &do_nu4pair, filename.length());
    read_eos_table_((char*)eosfilename.c_str());

    // resize arrays
    ns = __nulibtable_MOD_nulibtable_number_species;
    ng = __nulibtable_MOD_nulibtable_number_groups;
    nv = __nulibtable_MOD_nulibtable_number_easvariables;
    cout << ng << endl; cout << NE << endl;
    assert(NE==ng);
    assert(NM*NF == ns);

    // set energy grid
    cout << endl;
    cout<<"NE="<<ng << endl;
    for(int i=0;i<ng;i++){
      E[i]            = __nulibtable_MOD_nulibtable_energies[i]*1e6*cgs::units::eV; // erg
      nu[i]           = E[i] / (2.*M_PI*cgs::constants::hbar); // Hz
      double nubottom = __nulibtable_MOD_nulibtable_ebottom [i]*1e6*cgs::units::eV / (2.*M_PI*cgs::constants::hbar); // Hz
      double nutop    = __nulibtable_MOD_nulibtable_etop    [i]*1e6*cgs::units::eV / (2.*M_PI*cgs::constants::hbar); // Hz
      dnu[i]          = nutop - nubottom   ;
      cout << E[i]/(1.e6*cgs::units::eV) << " ";
    }
    cout.flush();

    cout << "#   logrho range: {" << __nulibtable_MOD_nulibtable_logrho_min << "," << __nulibtable_MOD_nulibtable_logrho_max << "} g/ccm" << endl;
    cout << "#   logT   range: {" << __nulibtable_MOD_nulibtable_logtemp_min << "," << __nulibtable_MOD_nulibtable_logtemp_max << "} MeV" << endl;
    cout << "#   Ye  range: {" << __nulibtable_MOD_nulibtable_ye_min << "," << __nulibtable_MOD_nulibtable_ye_max << "}" << endl;
    cout << "#   E   range: {" << __nulibtable_MOD_nulibtable_ebottom[0] << "," << __nulibtable_MOD_nulibtable_etop[ng-1] << "} MeV" << endl;
    cout << "#   n_species = " << ns << endl;
    cout << "#   n_rho   = " << __nulibtable_MOD_nulibtable_nrho << endl;
    cout << "#   n_T     = " << __nulibtable_MOD_nulibtable_ntemp << endl;
    cout << "#   n_Ye    = " << __nulibtable_MOD_nulibtable_nye << endl;
    cout << "#   n_E     = " << ng << endl;
  }

  void set(double rho_in, double T_in/*MeV*/, double Ye_in, const int do_interact){
    // condition inputs, get EOS info
    Ye_in = max(Ye_in,__nulibtable_MOD_nulibtable_ye_min);
    munue = get_munue(rho_in,T_in,Ye_in); // MeV
    temperature = T_in;
    double eta = get_eta(rho_in,T_in,Ye_in); // dimensionless
    int n_legendre = 2;

    if(do_interact){
    // BASIC EOS
    nulibtable_range_species_range_energy_(&rho_in, &T_in, &Ye_in, &eas.front(),
					   &__nulibtable_MOD_nulibtable_number_species,
					   &__nulibtable_MOD_nulibtable_number_groups,
					   &__nulibtable_MOD_nulibtable_number_easvariables);

    // INELASTIC SCATTERING KERNELS
    if(do_iscat){
      double phi[n_legendre][ng][ng]; //[a][j][i] = legendre index a, out group i, and in group j (ccm/s)
      for(int lns=1; lns<=__nulibtable_MOD_nulibtable_number_species; lns++){
	nulibtable_inelastic_single_species_range_energy_(&T_in, &eta, &lns, (double*)phi,
							  &ng,&ng,&n_legendre);
	for(int og=0; og<ng; og++)
	  for(int ig=0; ig<ng; ig++){
	    int ki = kernel_index(lns-1, ig, og);
	    assert(ki < escat_kernel0.size());
	    escat_kernel0[ki] = phi[0][og][ig];
	    escat_kernel1[ki] = phi[1][og][ig];
	  }
      }
    }

    // PAIR ANNIHILATION KERNELS
    if(do_pair){
      //[a][j][i] = out group i, and in group j (ccm/s)
      // a=0:Phi0p a=1:Phi0a a=2:Phi1p a=3:Phi1a
      int nvars = n_legendre*2;
      double phi[nvars][ng][ng];
      for(int lns=1; lns<=__nulibtable_MOD_nulibtable_number_species; lns++){
	nulibtable_epannihil_single_species_range_energy_(&T_in, &eta, &lns, (double*)phi,
							  &ng,&ng,&nvars);
	for(int og=0; og<ng; og++)
	  for(int ig=0; ig<ng; ig++){
	    int ki = kernel_index(lns-1, ig, og);
	    assert(ki < pair_kernel0.size());
	    pair_kernel0[ki] = phi[1][og][ig];
	    pair_kernel1[ki] = phi[3][og][ig];
	  }
      }
    }
    }

    // free memory
    free(__nulibtable_MOD_nulibtable_emissivities);
    free(__nulibtable_MOD_nulibtable_absopacity);
    free(__nulibtable_MOD_nulibtable_scatopacity);
    free(__nulibtable_MOD_nulibtable_delta);
    free(__nulibtable_MOD_nulibtable_itable_phi0);
    free(__nulibtable_MOD_nulibtable_itable_phi1);
    free(__nulibtable_MOD_nulibtable_epannihiltable_phi0);
    free(__nulibtable_MOD_nulibtable_epannihiltable_phi1);
  }

  double get_mue(const double rho /* g/ccm */, const double temp /*MeV*/, const double ye) const{ // MeV
    double eos_variables[__nulib_MOD_total_eos_variables];
    for(int i = 0; i<__nulib_MOD_total_eos_variables; i++)
    eos_variables[i] = 0;
    eos_variables[0] = rho;
    eos_variables[1] = temp;
    eos_variables[2] = ye;

    set_eos_variables_(eos_variables);
    double mue = eos_variables[10];
    return mue;
  }

  double get_munue(const double rho /* g/ccm */, const double temp /*MeV*/, const double ye) const{ // MeV
    double eos_variables[__nulib_MOD_total_eos_variables];
    for(int i=0; i<__nulib_MOD_total_eos_variables; i++) eos_variables[i] = 0;
    eos_variables[0] = rho;
    eos_variables[1] = temp;
    eos_variables[2] = ye;

    set_eos_variables_(eos_variables);
    double mue = eos_variables[10];
    double muhat = eos_variables[13];
    return (mue-muhat);
  }
  double get_eta(const double rho /* g/ccm */, const double temp /*MeV*/, const double ye) const{ // dimensionless
    double eos_variables[__nulib_MOD_total_eos_variables];
    for(int i=0; i<__nulib_MOD_total_eos_variables; i++) eos_variables[i] = 0;
    eos_variables[0] = rho;
    eos_variables[1] = temp;
    eos_variables[2] = ye;

    set_eos_variables_(eos_variables);
    double mue = eos_variables[10];
    return mue/eos_variables[1];
  }

  int index(const int is,const int ig,const int iv) const{
    return is + ig*ns + iv*ns*ng;
  }

  int kernel_index(const int is,const int igin, const int igout) const{
    return is + igin*ns + igout*ns*ng;
  }

  inline int nu4_kernel_index(const int ik, const int i1, const int i3) const{
    return i3 + ng*i1 + ng*ng*ik;
  }
  inline int nu4_bin2(const int ik, const int i1, const int i3) const{
    return i1+i3-ik;
  }

  double emis(const int is,const int ig) const{ // 1/cm/sr
    return abs(is,ig) * fermidirac(is,ig);
  }
  double abs(const int is,const int ig) const{ // 1/cm
    int ind = index(is,ig,1);
    assert(ind < eas.size());
    return eas[ind];
  }
  double scat(const int is,const int ig) const{ // 1/cm
    int ind = index(is,ig,2);
    assert(ind < eas.size());
    return eas[ind];
  }
  double delta(const int is,const int ig) const{ // 1/cm
    if(do_delta){
      int ind = index(is,ig,3);
      assert(ind < eas.size());
      return eas[ind];
    }
    else return 0;
  }
  double fermidirac(const int is, const int ig) const{
    double mu;
    if(is==0)      mu =  munue;
    else if(is==1) mu = -munue;
    else           mu = 0;
    mu *= 1e6*cgs::units::eV;
    double kTerg = temperature * 1e6*cgs::units::eV;
    double result = 1./(1. + exp((E[ig]-mu)/kTerg));
    return result;
  }
  double Phi0scat(const int is,const int igin, const int igout) const{ // cm^3/s/sr
    double result = 0;
    if(igin == igout)
      result += 2. * scat(is,igin)
	/(4.*M_PI*nu[igin]*nu[igin]*dnu[igin]/cgs::constants::c4);
    if(do_iscat)
      result += escat_kernel0[kernel_index(is,igin,igout)];
    return result;
  }
  double Phi1scat(const int is,const int igin, const int igout) const{ // cm^3/s/sr
    double result = 0;
    if(igin == igout)
      result += 2./3. * scat(is,igin)*delta(is,igin)
	/(4.*M_PI*nu[igin]*nu[igin]*dnu[igin]/cgs::constants::c4);
    if(do_iscat)
      result += escat_kernel1[kernel_index(is,igin,igout)];
    return result;
  }
  double Phi0pair(const int is,const int igin, const int igout) const{ // cm^3/s/sr
    double result = 0;
    if(do_pair)
      result += pair_kernel0[kernel_index(is,igin,igout)];
    return result;
  }
  double Phi1pair(const int is,const int igin, const int igout) const{ // cm^3/s/sr
    double result = 0;
    if(do_pair)
      result += pair_kernel1[kernel_index(is,igin,igout)];
    return result;
  }
};

#endif
