#ifndef FURRY_COMMON_RANGE_H_
#define FURRY_COMMON_RANGE_H_

#include <iterator>

namespace furry {

class Range {
 private:
  class Iterator : public std::iterator<std::forward_iterator_tag, int> {
   public:
    Iterator(int value, const Range *range) : value_(value), range_(range) {}

    int operator * () const {
      return value_;
    }

    bool operator != (const Iterator &it) {
      return value_ != it.value_;
    }

    const Iterator& operator ++ () {
      value_ += range_->step_;
      return *this;
    }

   private:
    int value_;
    const Range *range_;
  };

 public:
  Range(int first, int last, int step)
      : first_(first), last_(last), step_(step) {
    // if (last_ < first_) step_ = -step;
    last_ = ((last_ - first_) / step_ + 1 * ((last_ - first_) % step > 0)) * step_ + first_;
  }

  Range(int first, int last) : first_(first), last_(last), step_(1) {
    if (last_ < first_) step_ = -1;
  }

  Range(int last) : first_(0), last_(last), step_(1) {
    if (last_ < first_) step_ = -1;
  }

  Range() : first_(0), last_(0), step_(1) {}

  Range(const Range &r)
      : first_(r.first_), last_(r.last_), step_(r.step_) {}

  Range::Iterator begin() const {
    return Range::Iterator(first_, this);
  }

  Range::Iterator end() const {
    return Range::Iterator(last_, this);
  }

  Range& set_first(int first) {
    first_ = first;
    return *this;
  }

  Range& set_last(int last) {
    last_ = last;
    return *this;
  }

  Range& set_step(int step) {
    step_ = step;
    return *this;
  }

  int first() const {
    return first_;
  }

  int last() const {
    return last_;
  }

  int step() const {
    return step_;
  }

  int num_steps() const {
    return (last_ - first_) / step_;
  }

 private:
  int first_;
  int last_;
  int step_;
};

} // furry

#endif
