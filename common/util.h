#ifndef FURRY_COMMON_UTIL_H
#define FURRY_COMMON_UTIL_H

namespace furry {

class Indexable {
 public:
  Indexable() {}
  Indexable(const Indexable &indexable) : index_(indexable.index_) {}
  int index() const;
  void set_index(int i);
  void set_invalid();
  bool has_valid_index() const;
 private:
  int index_ = -1;
};

} // furry

#endif // FURRY_COMMON_UTIL_H
