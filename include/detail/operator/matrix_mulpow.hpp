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
  constexpr const size_t bs=64; // block size
  const auto[row,col]=dim();
  const size_t midCR=a.dim().col;
  size_t I,J,K,i,j,k,it,jt,kt;
  for(I=0,it=std::min(bs,row);I<row;I+=bs,it=std::min(it+bs,row))
  for(J=0,jt=std::min(bs,col);J<col;J+=bs,jt=std::min(jt+bs,col))
  for(K=0,kt=std::min(bs,midCR);K<midCR;K+=bs,kt=std::min(kt+bs,midCR))
    for(i=I;i<it;++i){
      const auto as=a[i];
      auto cs=(*this)[i];
      for(k=K;k<kt;++k){
        const auto&ak=as[k];
        const auto bs=b[k];
        for(j=J;j<jt;++j)
          cs[j]+=ak*bs[j];
      }
    }
  return *this;
}

template<class T>
template<class U>
matrix<T>& matrix<T>::operator*=(const matrix<U>&rhs){
  [[assume(dim().col==rhs.dim().row)]];
  return *this=*this*rhs;
}

template<class U,class V,class W=std::common_type_t<U,V>>
matrix<W> operator*(const matrix<U>&lhs,const matrix<V>&rhs){
  [[assume(lhs.dim().col==rhs.dim().row)]];
  matrix<W> ret(lhs.dim().row,rhs.dim().col,W());
  ret.fma(lhs,rhs);
  return ret;
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