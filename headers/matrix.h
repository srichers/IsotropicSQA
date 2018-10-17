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
