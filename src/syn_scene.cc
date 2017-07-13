const char *kUsageMessage = "Synthesize scene";

#include <vector>

#incldue <Eigen/Dense>

#include "furry/common/init.h"
#include "furry/common/rand.h"

using namespace std;
using namespace Eigen;

void SampleTriangles(const Vector3d &v0,
                     const Vector3d &v1,
                     const Vector3d &v1,
                     int num_samples,
                     vector<Vector3d> *samples) {
  Vector3d sample;
  double s, t;
  for (int i = 0; i < num_samples; ++i) {
    s = rand(0, 1);
    t = rand(0, 1);
    s = sqrt(s);
    sample = (1 - s) * v0 + s * (1 - t) * v1 + t * s * v2;
    samples->push_back(sample);
  }
}

class Drawable {
 public:
  virtual void Draw() const = 0;
};

class Mesh : public Drawable {
 public:
  virtual void Sample(int num_samples, vector<Vector3d> *samples) const = 0;
  void Draw () const {

  }
 protected:
  vector<Vector3d> vertexes_;
  vector<Vector3i> faces_;
};

class Square : public Mesh {
 public:
  Square(double side,
         const Vector3d &center,
         const AngleAxis &rotation)
      : vertexes_(4), faces_(2) {
    double half_side = side / 2;
    vertexes_[0] = Vector3d(-half_side, -half_side, 0);
    vertexes_[1] = Vector3d( half_side, -half_side, 0);
    vertexes_[2] = Vector3d( half_side,  half_side, 0);
    vertexes_[3] = Vertex3d(-half_side,  half_side, 0);;
    for (int i = 0; i < 4; ++i) {
      vertexes_[i] = rotation * vertexes_[i] + center;
    }
  }

  void Sample(int num_samples, vector<Vector3d> *samples) const {
  }
};

int main (int argc, char *argv[]) {
  furry::Init(&argc, &argv, kUsageMessage);
  return 0;
}
