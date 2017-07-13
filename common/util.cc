#include "furry/common/util.h"

namespace furry {

int Indexable::index() const {
  return index_;
}

void Indexable::set_index(int i) {
  index_ = i;
}

void Indexable::set_invalid() {
  index_ = -1;
}

bool Indexable::has_valid_index() const {
  return index_ >= 0;
}

} // furry
