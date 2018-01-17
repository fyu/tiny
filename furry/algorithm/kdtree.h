#ifndef FURRY_ALG_KDTREE_H_
#define FURRY_ALG_KDTREE_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>

#include <Eigen/Dense>

#include "furry/algorithm/euclidean.h"
#include "furry/common/numeric.h"
#include "furry/common/eigen.h"
#include "furry/common/dist.h"

namespace furry {

template <int K, typename DataType = float, typename PointType = Eigen::Matrix<DataType, K, 1>>
class KdTree {
  // typedef Eigen::Matrix<DataType, K, 1> PointType;
 public:
  KdTree() {}
  ~KdTree() {
    Destroy();
  }

  KdTree(const std::vector<PointType>& points) {
    Build(points);
  }

  void Build(const std::vector<PointType> &points) {
    std::vector<PointType> points_;
    points_.reserve(points.size());
    for (auto it = points.begin(); it != points.end(); ++it) {
      points_.push_back(*it);
    }
    root_ = BuildTree_(points_.begin(), points_.end(), 0, NULL);
  }

  PointType NN(const PointType &p, double *nn_dist2 = NULL) {
    if (root_ == NULL) {
      std::cout << "Error: nothing in the kd tree\n";
      return p;
    }
    // find the partition
    Node *next_node = root_, *current_node;
    int axis = 0;
    int depth = -1;
    do {
      axis = ++depth % K;
      current_node = next_node;
      if (p[axis] < current_node->p[axis])
        next_node = current_node->left;
      else
        next_node = current_node->right;
    } while (next_node);

    // traverse the tree back
    double dist2, best_dist2 = Distance2(p, current_node->p);
    Node *best_node = current_node;
    next_node = current_node->parent;
    Node *sibling;
    while (next_node) {
      --depth;
      axis = depth % K;
      if (Distance2(p, next_node->p) < best_dist2) {
        best_dist2 = Distance2(p, next_node->p);
        best_node = next_node;
      }
      if (square(p[axis] - next_node->p[axis]) < best_dist2) {
        if (p[axis] < next_node->p[axis])
          sibling = next_node->right;
        else
          sibling = next_node->left;
        if (sibling != NULL) {
          auto node = NearestNeighbor_(p, sibling, &dist2);
          if (dist2 < best_dist2) {
            best_node = node;
            best_dist2 = dist2;
          }
        }
      }
      next_node = next_node->parent;
    }
    if (nn_dist2 != NULL)
      *nn_dist2 = best_dist2;
    return best_node->p;
  }

  template <typename OutputIterator>
  OutputIterator RangeSearch(const PointType &lower,
                   const PointType &upper,
                   OutputIterator result) {
    assert(LessThanAll(lower, upper));
    Cube bound = {lower, upper};
    Cube root_bound =
    // {-PointType::Constant(std::numeric_limits<DataType>::max()),
    //  PointType::Constant(std::numeric_limits<DataType>::max())};
        {std::numeric_limits<PointType>::lowest(),
         std::numeric_limits<PointType>::max()};
    return RangeSearch_(bound, root_, 0, root_bound, result);
  }

  void Destroy() {
    DeleteNode_(root_);
  }

  size_t size() {
    return Count_(root_);
  }
 private:
  struct Node {
    Node(const PointType &p_)
        : p(p_), left(nullptr), right(nullptr), parent(nullptr) {}
    PointType p;
    Node *left;
    Node *right;
    Node *parent;
  };

  struct Cube {
     PointType l; // lower bound of the cube
    PointType u; // upper bound of the cube
  };

  Node *root_ = NULL;

  template <typename InputIterator>
  Node* BuildTree_(InputIterator first,
                   InputIterator last,
                   int depth,
                   Node *parent) {
    if (first == last)
      return NULL;
    typedef typename InputIterator::value_type PointType_;
    int axis = depth % K;
    sort(first, last,
         [axis](const PointType_ &p0, const PointType_ &p1) {
           return p0[axis] < p1[axis];
         });
    int length = last - first;
    int median = length / 2;
    Node *root = new Node(*(first + median));
    // root->p = *(first + median);
    ++depth;
    root->left = BuildTree_(first, first + median, depth, root);
    root->right = BuildTree_(first + median + 1, last, depth, root);
    root->parent = parent;
    return root;
  }

  void DeleteNode_(Node *root) {
    if (root == NULL)
      return;
    DeleteNode_(root->left);
    DeleteNode_(root->right);
  }

  size_t Count_(Node *root) {
    if (root == NULL)
      return 0;
    return 1 + Count_(root->left) + Count_(root->right);
  }

  Node*  NearestNeighbor_(const PointType &p,
                          Node *root,
                          double *best_dist2) {
    Node *best_node = root;
    *best_dist2 = Distance2(p, root->p);
    Node *node;
    double dist2;
    if (root->left) {
      node = NearestNeighbor_(p, root->left, &dist2);
      if (dist2 < *best_dist2) {
        best_node = node;
        *best_dist2 = dist2;
      }
    }
    if (root->right) {
      node = NearestNeighbor_(p, root->right, &dist2);
      if (dist2 < *best_dist2) {
        best_node = node;
        *best_dist2 = dist2;
      }
    }
    return best_node;
  }

  template <typename OutputIterator> OutputIterator
  RangeSearch_(const Cube &bound,
               Node *root,
               int depth,
               Cube &root_bound,
               OutputIterator result) {
    if (root == NULL)
      return result;

    int axis = depth % K;
    ++depth;
    DataType lower = root_bound.l[axis];
    DataType upper = root_bound.u[axis];

    if (Contains_(bound, root->p))
      *result++ = root->p;

    // left side
    root_bound.u[axis] = root->p[axis];
    if (Contains_(bound, root_bound))
      result = Collect_(root->left, result);
    else if (Overlap_(bound, root_bound))
      result = RangeSearch_(bound, root->left, depth, root_bound, result);

    // right side
    root_bound.l[axis] = root->p[axis];
    root_bound.u[axis] = upper;
    if (Contains_(bound, root_bound))
      result = Collect_(root->right, result);
    else if (Overlap_(bound, root_bound))
      result = RangeSearch_(bound, root->right, depth, root_bound, result);

    root_bound.l[axis] = lower;
    //    root_bound.u[axis] = upper;

    return result;
  }

  bool LessEqualAll(const PointType &p0, const PointType &p1) {
    for (int i = 0; i < K; ++i) {
      if (p0[i] > p1[i]) return false;
    }
    return true;
  }

  bool LessThanAll(const PointType &p0, const PointType &p1) {
    for (int i = 0; i < K; ++i) {
      if (p0[i] >= p1[i]) return false;
    }
    return true;
  }


  bool Contains_(const Cube &c0, const Cube &c1) {
    // return (c0.l.array() <= c1.l.array()).all() &&
    //     (c1.u.array() <= c0.u.array()).all();
    return LessEqualAll(c0.l, c1.l) && LessEqualAll(c1.u, c0.u);
  }

  bool Contains_(const Cube &c,
                 const PointType &p) {
    // return (c.l.array() <= p.array()).all() &&
    //     (p.array() <= c.u.array()).all();
    return LessEqualAll(c.l, p) && LessEqualAll(p, c.u);
  }

  bool Overlap_(const Cube &c0, const Cube &c1) {
    // return (c0.l.array() < c1.u.array()).all() &&
    //     (c1.l.array() < c0.u.array()).all();
    return LessThanAll(c0.l, c1.u) && LessThanAll(c1.l, c0.u);
  }

  // in-order traversal
  template <typename OutputIterator>
  OutputIterator Collect_(Node *root, OutputIterator result) {
    if (root == NULL)
      return result;
    result = Collect_(root->left, result);
    *result++ = root->p;
    result = Collect_(root->right, result);
    return result;
  }
};

} // furry

#endif // FURRY_ALG_KDTREE_H_
