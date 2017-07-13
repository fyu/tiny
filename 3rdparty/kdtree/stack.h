#ifndef THIRD_PARTY_KDTREE_STACK_H__
#define THIRD_PARTY_KDTREE_STACK_H__


//   A kdtree implementation.

//   Contact: Ravi Kolluri
//            rkolluri@gmail.com

//   Copyright 2006 - Ravi Kolluri.


//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

//   Commercial use is absolutely prohibited.

//   See the accompanying LICENSE file for details.


//  A simple stack implementation used in the kdtree.
//  For efficiency, the size of the stack is fixed.
//  Use stl::stack if you need a more general stack class.

namespace kdtree {

template <typename Type, int Size>
class Stack {
 private:
  Type _data[Size];
  int _topindex;

 public:
  Stack() :
      _topindex(0)
  {}

  Type top() {
    return _data[_topindex - 1];
  }

  void pop() {
    if (_topindex > 0)
      _topindex--;
  }

  int size() {
    return _topindex;
  }

  void push(const Type& T) {
    _data[_topindex++] = T;
  }

  bool empty() {
    return (_topindex == 0);
  }

  void clear() {
    _topindex = 0;
  }
};

}  // namespace kdtree

#endif  // THIRD_PARTY_KDTREE_STACK_H__
