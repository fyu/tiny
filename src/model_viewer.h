#ifndef FURRY_TINY_MODEL_VIEWER_H
#define FURRY_TINY_MODEL_VIEWER_H

#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Qt>
#include <QObject>
#include <QApplication>
#include <QGLWidget>
#include <QtOpenGL>

#if __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include "wrap/gui/trackball.h"

#include "model.h"

namespace furry {
namespace tiny {

class ModelViewer : public QGLWidget {
  Q_OBJECT

 public:
  ModelViewer(const Model *model, QWidget *parent = nullptr);

  void SetModel(const Model *model);

 protected:
  virtual QSize sizeHint() const;

  virtual void initializeGL();
  virtual void paintGL();
  virtual void mouseMoveEvent(QMouseEvent *event);
  virtual void mousePressEvent(QMouseEvent *event);
  virtual void mouseReleaseEvent(QMouseEvent *event);
  virtual void resizeGL(int w, int h);

  virtual void keyReleaseEvent(QKeyEvent * e);
  virtual void keyPressEvent(QKeyEvent * e);
  virtual void wheelEvent(QWheelEvent*e);

 protected:
  void DrawCamera(const Camera *camera);
  double FlipY(double y);
  bool SetCenter();

 protected:
  vcg::Trackball track_;
  const Model *model_;
  QSize window_size_;
  int mouse_[2] = { 0, 0 };
  Eigen::Vector3d max_coord_;
  Eigen::Vector3d min_coord_;
  Eigen::Vector3d center_;
  int point_size_ = 3;
  int line_width_ = 3;
};

} // tiny
} // furry

#endif // FURRY_TINY_MODEL_VIEWER_H
