#include "model_viewer.h"

#include <glog/logging.h>

#include "wrap/qt/trackball.h"

using namespace Eigen;

namespace furry {
namespace tiny {

ModelViewer::ModelViewer(const Model *model, QWidget *parent)
    : window_size_(1920, 1080), QGLWidget(parent) {
  SetModel(model);
  // Get rid of everything of the window

  if (SetCenter()) {
    setWindowFlags(Qt::Window |
                   Qt::FramelessWindowHint |
                   Qt::CustomizeWindowHint);
  }
}

void ModelViewer::SetModel(const Model *model) {
  model_ = model;
  model->GetTrackBoundingBox(&min_coord_, &max_coord_);
  min_coord_.z() = 0;
  center_ = model->GetMeanPoint();
}

QSize ModelViewer::sizeHint() const {
  return window_size_;
}

void ModelViewer::paintGL() {
  Vector3d eye(0, 0, -3);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glLoadIdentity();
  gluLookAt(eye.x(), eye.y(), eye.z(),
            center_.x(), center_.y(), center_.z(),
            0, -1, 0);

  track_.center = vcg::Point3f(center_.x(), center_.y(), center_.z());
  track_.radius = (center_ - eye).norm() / 2;
  track_.GetView();
  track_.Apply();

  // gluLookAt(0, 0, 0, 0, 0, 1, 0, -1, 0);
  glColor3f(1.0, 1.0, 1.0);
  glPointSize(point_size_);
  glBegin(GL_POINTS);
  for (int i = 0; i < model_->NumTracks(); ++i) {
    auto track = model_->GetTrack(i);
    glColor3ubv(track->GetColor().data());
    // LOG_IF(INFO, i % 40 == 0) << "Drawing " << track->GetPoint().transpose();
    glVertex3dv(track->GetPoint().data());
  }
  glEnd();

  DrawCamera(model_->GetCamera(0));

  track_.DrawPostApply();

  return;
}

void ModelViewer::initializeGL() {
  // Initialize background color
  glClearColor(0, 0, 0, 0);

  // Initialize lights
  // glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  // static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  // glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  // glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  // static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  // glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  // glEnable(GL_LIGHT0);
  // static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  // glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  // glEnable(GL_LIGHT1);
  // glEnable(GL_NORMALIZE);

  // Initialize graphics modes  
  glEnable(GL_DEPTH_TEST);
}

void ModelViewer::mouseMoveEvent(QMouseEvent *event) {
  int x = event->x();
  int y = event->y();
  auto buttons = event->buttons();
  // Invert y coordinate
  y = FlipY(y);

  // Compute mouse movement
  int dx = x - mouse_[0];
  int dy = y - mouse_[1];

  // Remember mouse position 
  mouse_[0] = x;
  mouse_[1] = y;

  auto e = event;
  if (e->buttons()) {
    track_.MouseMove(e->x(), FlipY(e->y()));
    updateGL();
  }
}

void ModelViewer::mousePressEvent(QMouseEvent *event) {
  mouse_[0] = event->x();
  mouse_[1] = FlipY(event->y());

  auto e = event;
  e->accept ();
  setFocus ();
  track_.MouseDown (e->x (), FlipY(e->y()), QT2VCG (e->button (), e->modifiers ()));
  updateGL ();
}

void ModelViewer::mouseReleaseEvent(QMouseEvent *event) {
  auto e = event;
  track_.MouseUp (e->x (), FlipY(e->y()), QT2VCG (e->button (), e->modifiers ()));
  updateGL ();
}

void ModelViewer::keyReleaseEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)
    track_.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)
    track_.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)
    track_.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}

void ModelViewer::keyPressEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)
    track_.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)
    track_.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)
    track_.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));
  if (e->key() == Qt::Key_Equal) ++point_size_;
  if (e->key() == Qt::Key_Minus) --point_size_;
  updateGL ();
}

void ModelViewer::wheelEvent (QWheelEvent * e)
{
  const int WHEEL_STEP = 120;
  track_.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
  updateGL ();
}

void ModelViewer::resizeGL(int w, int h) {
  window_size_.setWidth(w);
  window_size_.setHeight(h);

  // Resize window
  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90, (double)w / (double)h, 1e-3, 1e3);
  glMatrixMode(GL_MODELVIEW);
}

double ModelViewer::FlipY(double y) {
  return window_size_.height() - y - 1;
}

bool ModelViewer::SetCenter() {
  int w = window_size_.width();
  int h = window_size_.height();
  QDesktopWidget* desktop = QApplication::desktop();
  int sw = desktop->width();
  int sh = desktop->height();

  int x, y;
  if (sw > w && sh > h) {
    x = (sw - w) / 2;
    y = (sh - h) / 2;
    setGeometry(x, y, w, h);
    return true;
  } else {
    return false;
  }
}

void ModelViewer::DrawCamera(const Camera *camera) {
  Vector3d c = camera->GetCenter();
  Vector3d t = camera->GetTowards();
  Vector3d u = camera->GetUp();
  Vector3d s = t.cross(u);
  double fovx = camera->GetFovX();
  double fovy = camera->GetFovY();

  // LOG(INFO) << fovx << ' ' << fovy;
  // LOG(INFO) << '\n' << camera->R();

  double camera_size = 0.1 * (c - (center_) / 2).norm();
  double w = camera_size * tan(fovx / 2);
  double h = camera_size * tan(fovy / 2);

  Vector3d v[4];
  v[0] = c + t * camera_size + u * h + s * w;
  v[1] = c + t * camera_size + u * h - s * w;
  v[2] = c + t * camera_size - u * h - s * w;
  v[3] = c + t * camera_size - u * h + s * w;

  glLineWidth(line_width_);
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINES);
  glVertex3dv(c.data());
  glVertex3dv(v[0].data());
  glVertex3dv(c.data());
  glVertex3dv(v[1].data());
  glVertex3dv(c.data());
  glVertex3dv(v[2].data());
  glVertex3dv(c.data());
  glVertex3dv(v[3].data());

  glVertex3dv(v[0].data());
  glVertex3dv(v[1].data());
  glVertex3dv(v[1].data());
  glVertex3dv(v[2].data());
  glVertex3dv(v[2].data());
  glVertex3dv(v[3].data());
  glVertex3dv(v[3].data());
  glVertex3dv(v[0].data());
  glEnd();
}

} // tiny
} // furry
