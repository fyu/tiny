#ifndef FURRY_COMMON_MEMORY_H
#define FURRY_COMMON_MEMORY_H

#include <cstdio>
#include <vector>
#include <unordered_map>

namespace furry {

template <typename T>
void Delete(std::vector<T*> *values) {
  for (T *p : *values) {
    delete p;
  }
  values->resize(0);
  values->shrink_to_fit();
}

template <typename K, typename T>
void Delete(std::unordered_map<K, T*> *values) {
  for (auto &v : *values) {
    delete v.second;
  }
  values->clear();
}

template<class T>
class ScopedPtr {
 public:
  ScopedPtr(T * p = 0); // never throws
  ScopedPtr(const ScopedPtr<T> &sp);
  ~ScopedPtr(); // never throws

  void Reset(T * p = 0); // never throws

  T& operator*() const; // never throws
  T* operator->() const; // never throws
  operator bool () const;
  operator T* () const;
  ScopedPtr<T>& operator = (const ScopedPtr<T>&) = delete;
  ScopedPtr<T>& operator = (ScopedPtr<T> &&p);
  ScopedPtr<T>& operator = (T *p);
  T* Get() const; // never throws
  T* get();

 private:
  T *p_ = nullptr;
};

template<class T>
class ScopedArray {
 public:
  explicit ScopedArray(T * p = nullptr); // never throws
  ~ScopedArray(); // never throws

  void Reset(T * p = nullptr); // never throws

  operator bool () const;
  operator T* () const;
  const T& operator[](std::ptrdiff_t i) const; // never throws
  T& operator[](std::ptrdiff_t i); // never throws
  // T* operator + (std::ptrdiff_t i);
  T* Get() const; // never throws

 private:
  T *p_ = nullptr;
};

////////////////////////////////////////////////////////////
// ScopedPtr inline
////////////////////////////////////////////////////////////

template <typename T> inline
ScopedPtr<T>::ScopedPtr(T *p) : p_(p) {
}

template <typename T> inline
ScopedPtr<T>::ScopedPtr(const ScopedPtr<T> &sp) {
  p_ = sp.p_;
}

template <typename T> inline
ScopedPtr<T>::~ScopedPtr() {
  Reset();
}

template <> inline
void ScopedPtr<std::FILE>::Reset(std::FILE *p) {
  if (p_ != nullptr)
    fclose(p_);
  p_ = p;
}

template <typename T> inline
void ScopedPtr<T>::Reset(T *p) {
  if (p_ != nullptr)
    delete p_;
  p_ = p;
}

template <typename T> inline
T& ScopedPtr<T>::operator * () const {
  return *p_;
}

template <typename T> inline
T* ScopedPtr<T>::operator -> () const {
  return p_;
}
template <typename T> inline
ScopedPtr<T>::operator bool () const {
  return p_ != nullptr;
}

template <typename T> inline
ScopedPtr<T>::operator T* () const {
  return p_;
}

template <typename T> inline
ScopedPtr<T>& ScopedPtr<T>::operator = (ScopedPtr<T> &&p) {
  if (p_ != p.p_) {
    Reset(p.p_);
    p.p_ = nullptr;
  }
  return *this;
}

template <typename T> inline
ScopedPtr<T>& ScopedPtr<T>::operator = (T *p) {
  if (p_ != p) {
    Reset(p);
  }
  return *this;
}

template <typename T> inline
T* ScopedPtr<T>::Get() const {
  return p_;
}

template <typename T> inline
T* ScopedPtr<T>::get() {
  return p_;
}

////////////////////////////////////////////////////////////
// ScopedArray inline
////////////////////////////////////////////////////////////

template <typename T> inline
ScopedArray<T>::ScopedArray(T *p) : p_(p) {
}

template <typename T> inline
ScopedArray<T>::~ScopedArray() {
  Reset();
}

template <typename T> inline
void ScopedArray<T>::Reset(T *p) {
  if (p_ != nullptr)
    delete [] p_;
  p_ = p;
}

template <typename T> inline
ScopedArray<T>::operator bool () const {
  return p_ != nullptr;
}

template <typename T> inline
ScopedArray<T>::operator T* () const {
  return p_;
}

template <typename T> inline
T& ScopedArray<T>::operator [] (std::ptrdiff_t i) {
  return p_[i];
}

template <typename T> inline
const T& ScopedArray<T>::operator [] (std::ptrdiff_t i) const {
  return p_[i];
}

// template <typename T> inline
// T* ScopedArray<T>::operator + (std::ptrdiff_t i) {
//   return p_ + i;
// }

template <typename T> inline
T* ScopedArray<T>::Get() const {
  return p_;
}

} // furry

#endif
