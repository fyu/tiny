#ifndef FURRY_TINY_MVS_UI_H
#define FURRY_TINY_MVS_UI_H

#include <Qt>
#include <QObject>
#include <QApplication>
#include <QGLWidget>
#include <QtOpenGL>
#include <QGroupBox>

#include <Eigen/Dense>

#include "furry/gui/window_QT.h"
#include "mvs.h"

namespace furry {
namespace tiny {

class Plot2d : public CvWindow {
  Q_OBJECT
 public:
  Plot2d(QString title);
  void SetRange(const Eigen::Vector2f &x_range,
                const Eigen::Vector2f &y_range);

 public slots:
  void UpdatePoints(const std::vector<cv::Point2f> &points);
  void UpdatePoints(const std::vector<std::vector<cv::Point2f>> &points);

 signals:
  void PointClicked(const cv::Point2f &p);

 protected:
  void resizeEvent(QResizeEvent *event);

 protected:
  void DrawPoints();

 protected:
  std::vector<std::vector<cv::Point2f>> points_;
  cv::Mat canvas_;
  Eigen::Vector2f x_range_;
  Eigen::Vector2f y_range_;
  bool fixed_range_ = false;
};

class AdjustableSlider : public QWidget {
  Q_OBJECT
  public:
  AdjustableSlider(const std::string &title, bool int_type = false,
                   QWidget *parent = 0);

 public slots:
  void SetValueText(const QString & text);
  void SetMaxText(const QString & text);
  void SetMinText(const QString & text);
  void SetSliderValue(int value);

  void SetMaxValue(double max);
  void SetMinValue(double min);
  void SetValue(double value);

 signals:
  void ValueChanged(double value);
  void ValueEdited(double value);

 protected:
  void UpdateValue();

 private:
  QSlider *slider_;
  QLabel *value_label_;
  QLineEdit *value_text_;
  QLabel *max_label_;
  QLineEdit *max_text_;
  QLabel *min_label_;
  QLineEdit *min_text_;
  QLabel *title_;

  QDoubleValidator double_validator_;
  QIntValidator int_validator_;
  int num_steps_ = 100;
  double value_ = 0;
  double min_value_ = -10;
  double max_value_ = 10;
  bool int_type_ = false;
};

class MvsController : public QWidget {
  Q_OBJECT

  enum CombineMode {
    COMBINE_SINGLE,
    COMBINE_SUM,
  };

 public:
  MvsController(MvsModel* model, QWidget* parent = nullptr);

 public slots:
  void SetRefPoint(const cv::Point2f& point);
  void SetRefDepth(double d);
  void SetClickedPlotPoint(const cv::Point2f &point);
  void SetCombineMode(int m);

 private:
  void ShowPhotoConsistencyScores() const;

 private:
  MvsModel* model_;
  CvWindow* ref_window_;
  CvWindow* image_window_;
  CvWindow* patch_window_;
  Plot2d* plot_window_;

  int combine_mode_ = COMBINE_SINGLE;

  AdjustableSlider* depth_slider_;
};

} // tiny
} // furry

#endif // FURRY_TINY_MVS_UI_H
