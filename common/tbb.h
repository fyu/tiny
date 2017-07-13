#ifndef FURRY_COMMON_TBB_H
#define FURRY_COMMON_TBB_H

#include <functional>

namespace furry {

void ParFor(int begin, int end, std::function<void (int)> func);
void parfor(int begin, int end, std::function<void (int)> func);

}

#endif // FURRY_COMMON_TBB_H
