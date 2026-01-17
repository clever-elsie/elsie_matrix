#ifndef ELSIE_MATRIX_OPERATOR_MATRIX_MULPOW_HPP
#define ELSIE_MATRIX_OPERATOR_MATRIX_MULPOW_HPP

#include "../_image.hpp"

namespace elsie{

template<class T>
template<class U,class V>
matrix<T>& matrix<T>::fma(const matrix<U>&a,const matrix<V>&b){
  [[assume(dim().row==a.dim().row)]];
  [[assume(dim().col==b.dim().col)]];
  [[assume(a.dim().col==b.dim().row)]];
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator*=(const matrix<U>&rhs){
  [[assume(dim().col==rhs.dim().row)]];
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator*(const matrix<U>&lhs,const matrix<V>&rhs){
  [[assume(lhs.dim().col==rhs.dim().row)]];
}

template<class T>
matrix<T> matrix<T>::pow_impl(matrix<T>&&k, uint64_t b){
  [[assume(dim().row==dim().col)]];
  matrix<T> ret(k.dim(),T());
  for(size_t i=0;i<dim().row;++i)
    ret[i,i]=T(1);
  for(;b;b>>=1){
    if(b&1) ret=ret*k;
    k=k*k;
  }
}

template<class T>
void matrix<T>::pow(uint64_t b){
  *this=pow_impl(std::move(*this),b);
}

template<class T>
matrix<T> matrix<T>::pow(uint64_t b)const{
  return pow_impl(matrix<T>(this->copy()),b);
}

} // namespace elsie
#endif