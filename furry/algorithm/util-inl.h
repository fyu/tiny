namespace furry {

template <typename T, typename CompFunc> inline
T Max(const std::vector<T> &v, CompFunc comp_func) {
  return *(std::max_element(v.begin(), v.end(), comp_func));
}

template <typename T>
T Mean(const std::vector<T> &v) {
  if (v.size() == 0)
    return T();
  else {
    return accumulate(v.begin() + 1, v.end(), *v.begin()) / v.size();
  }
}

template <typename T>
T Clamp(T v, const T &min, const T &max) {
  if (v < min)
    v = min;
  if (v > max)
    v = max;
  return v;
}

} // furry
