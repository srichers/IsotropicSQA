#include <cassert>

template<typename T, unsigned a, unsigned b>
  class MATRIX{
    T m[a][b];

 public:
    const T* operator[](const size_t i) const{
      return m[i];
    }
    T* operator[](const size_t i){
      return m[i];
    }

    MATRIX<T,a,b>(){
      for(unsigned i=0; i<a; i++)
	for(unsigned j=0; j<b; j++)
	  m[i][j] = 0;
    }
    
    template<typename T1, unsigned c>
    MATRIX<T,a,b> operator*(const MATRIX<T1,b,c>& right) const{
      MATRIX<T,a,c> result;
      for(unsigned i=0; i<a; i++)
	for(unsigned j=0; j<c; j++){
	  result[i][j] = 0;
	  for(unsigned k=0; k<b; k++){
	    result[i][j] += m[i][k]*right[k][j];
	  }
	}
      return result;
    }

    template<typename T1>
    MATRIX<T,a,b> operator+(const MATRIX<T1,a,b>& right) const{
      MATRIX<T,a,b> result;
      for(unsigned i=0; i<a; i++)
	for(unsigned j=0; j<b; j++)
	  result[i][j] = m[i][j] + right[i][j];
      return result;
    }

    template<typename T1>
    MATRIX<T,a,b> operator-(const MATRIX<T1,a,b>& right) const{
      MATRIX<T,a,b> result;
      for(unsigned i=0; i<a; i++)
	for(unsigned j=0; j<b; j++)
	  result[i][j] = m[i][j] - right[i][j];
      return result;
    }

    MATRIX<T,a,b> operator-() const{
      MATRIX<T,a,b> result;
      for(unsigned i=0; i<a; i++)
	for(unsigned j=0; j<b; j++)
	  result[i][j] = -m[i][j];
      return result;
    }

    template<typename T1>
    MATRIX<T,a,b> operator*=(const T1 input){
      for(unsigned i=0; i<a; i++)
	for(unsigned j=0; j<b; j++)
	  m[i][j] *= input;
      return *this;
    }

};

template<typename T, unsigned a, unsigned b>
ostream& operator<<(ostream& os, const MATRIX<T,a,b>& m){
  for(unsigned i=0; i<a; i++){
    os << "[";
    for(unsigned j=0; j<b; j++){
      os << m[i][j] << (j==b-1 ? "]" : ",");
    }
    os << endl;
  }
  return os;  
}

template<typename T, unsigned a, unsigned b>
MATRIX<T,a,b> Adjoint(MATRIX<T,a,b> input){
  assert(a==b);
  MATRIX<T,a,b> result;
  for(unsigned i=0; i<a; i++)
    for(unsigned j=0; j<b; j++)
      result[i][j] = conj(input[j][i]);
  return result;
}

template<typename T, unsigned a, unsigned b>
MATRIX<T,a,b> Conjugate(MATRIX<T,a,b> input){
  MATRIX<T,a,b> result;
  for(unsigned i=0; i<a; i++)
    for(unsigned j=0; j<b; j++)
      result[i][j] = conj(input[i][j]);
  return result;
}

template<typename T, unsigned a, unsigned b>
T Trace(MATRIX<T,a,b> input){
  assert(a==b);
  T result = 0;
  for(unsigned i=0; i<a; i++)
    result += input[i][i];
  return result;
}

double Sign(double input){
  return input>0 ? 1 : -1;
}

void RungeKuttaCashKarpParameters(int &NRK,int &NOrder,const double* &A,const double** &B,const double* &C,const double* &D)
     { NRK=6; NOrder=5;

       static const double a[]={ 0., 1./5., 3./10., 3./5., 1., 7./8. };
       static const double b0[]={};
       static const double b1[]={ 1./5. };
       static const double b2[]={ 3./40.,9./40. };
       static const double b3[]={ 3./10.,-9./10.,6./5. };
       static const double b4[]={ -11./54.,5./2.,-70./27.,35./27. };
       static const double b5[]={ 1631./55296.,175./512.,575./13824.,44275./110592.,253./4096. };
       static const double* b[]={ b0,b1,b2,b3,b4,b5 };
       static const double c[]={ 37./378.,0.,250./621.,125./594.,0.,512./1771. };
       static const double d[]={ 2825./27648.,0.,18575./48384.,13525./55296.,277./14336.,1./4. };

       A=a; B=b; C=c; D=d;
      } 
