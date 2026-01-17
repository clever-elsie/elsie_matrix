#ifndef ELSIE_MATRIX_IMAGE_HPP
#define ELSIE_MATRIX_IMAGE_HPP
#include <cstddef>
#include <cstdint>
#include <span>
#include <array>
#include <memory>
#include <string>

namespace elsie{

// コピーコンストラクタ，コピー代入はviewの作成
// moveは所有権を持っていればmove，viewの場合は単にコピーしてコピー元のviewを無効化
// deep copyは.copy()を使う
// 所有権の有無はプログラマに委ねる．assert用にowned()を用意する
template<class T>
class matrix{
  public:
  struct dimension{
    size_t row,col;
    dimension():row(0),col(0){}
    dimension(size_t row, size_t col):row(row),col(col){}
    dimension(const dimension& other):row(other.row),col(other.col){}
    dimension& operator=(const dimension& other){row=other.row; col=other.col; return *this;}
    size_t capacity()const{return row*col;}
    bool operator==(const dimension& other)const{return row==other.row&&col==other.col;}
    bool operator!=(const dimension& other)const{return !(*this==other);}
  };
  using pointer=T*;
  using const_pointer=const T*;
  using reference=T&;
  using const_reference=const T&;

  private:
  dimension dim_,leading_;
  // capacityは本当のオーバーフローをしないために，確保している配列の末尾までの距離を持つ
  size_t capacity_data; // { capacity:63, is_view:1 }
  pointer data_;

  private:
  static size_t calc_capacity_data(size_t capacity, bool is_view=false){
    return (capacity<<1) | (is_view?1:0);
  }
  size_t get_capacity()const{ return capacity_data>>1; }
  bool is_view()const{ return capacity_data&1; }

  // make_view用のコンストラクタ
  matrix(const matrix<T>& other, size_t row, size_t col, size_t row_offset, size_t col_offset);
  public:
  // constructor.destructor.operator=.hpp
  matrix();
  matrix(const dimension&dim);
  matrix(const dimension&dim, const T& init);
  matrix(const size_t row, const size_t col);
  matrix(const size_t row, const size_t col, const T& init);
  matrix(const matrix& other);
  matrix(matrix&& other);
  matrix& operator=(const matrix& other);
  matrix& operator=(matrix&& other);
  ~matrix();
  matrix<T> make_view(size_t row, size_t col, size_t row_offset, size_t col_offset)const;
  matrix<T> copy()const;
  matrix<T> copy(size_t row, size_t col, size_t row_offset, size_t col_offset)const;

  // access.hpp
  dimension dim()const;
  dimension leading_dim()const;
  reference operator[](size_t i, size_t j);
  const_reference operator[](size_t i, size_t j)const;
  T val(size_t i, size_t j)const;
  std::span<T> operator[](size_t i);
  std::span<const T> operator[](size_t i)const;
  bool owned()const{ return !is_view(); }
  private:
  pointer data(size_t i=0, size_t j=0)const{
    return data_+i*leading_.col+j;
  }
  public:
  
  // split.merge.hpp
  private:
  std::array<matrix,4> split()const;
  static void merge(matrix& dst,
    const size_t row_offset, const size_t col_offset,
    const matrix&m00,const matrix&m01,
    const matrix&m10,const matrix&m11);
  public:

  // operator/unary.hpp
  operator std::string()const;
  std::string string()const;
  void transpose();
  matrix<T> transpose()const;
  void negate();
  matrix<T> operator-()const;

  // operator/scalar.hpp
  template<class U> matrix<T>& operator*=(const U&);
  template<class U> matrix<T>& operator/=(const U&);
  template<class U,class V,class W>
  friend matrix<W> operator*(const U&,const matrix<V>&);
  template<class U,class V,class W>
  friend matrix<W> operator*(const matrix<U>&,const V&);
  template<class U,class V,class W>
  friend matrix<W> operator/(const matrix<U>&,const V&);

  // operator/matrix_addsub.hpp
  template<class U> matrix<T>& operator+=(const matrix<U>&);
  template<class U> matrix<T>& operator-=(const matrix<U>&);
  template<class U,class V,class W>
  friend matrix<W> operator+(const matrix<U>&,const matrix<V>&);
  template<class U,class V,class W>
  friend matrix<W> operator-(const matrix<U>&,const matrix<V>&);

  // operator/matrix_mulpow.hpp
  template<class U,class V> // *this+=A*B
  matrix<T>& fma(const matrix<U>&,const matrix<V>&);
  template<class U> matrix<T>& operator*=(const matrix<U>&);
  template<class U,class V,class W>
  friend matrix<W> operator*(const matrix<U>&,const matrix<V>&);
  void pow(uint64_t);
  matrix<T> pow(uint64_t)const;
  private:
  matrix<T> pow_impl(matrix<T>&&k,uint64_t b);
  public:

  // operation/reduction.hpp
};

} // namespace elsie
#endif