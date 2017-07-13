#include "mvs_ui.h"

#include <iostream>
#include <algorithm>
#include <limits>
#include <functional>

#include <opencv2/core/core.hpp>
#include <Eigen/Dense>
#include <glog/logging.h>

#include "furry/common/str.h"

#include "util_ui.h"

using namespace cv;
using namespace std;
using namespace Eigen;

namespace furry {
namespace tiny {

Plot2d::Plot2d(QString title) : CvWindow(title, CV_WINDOW_NORMAL | CV_WINDOW_FREERATIO) {
  fixed_range_ = false;
  QSize size = myView->size();
  x_range_ = Vector2f(0, size.width() - 1);
  y_range_ = Vector2f(0, size.height() - 1);
  SetMouseCallback(
      // [&myView, &x_range_, &y_range_](int event, int x, int y, int, void*) {
      [&](int event, int x, int y, int, void*) {
      if (event == EVENT_LBUTTONDOWN) {
        QSize size = myView->size();
        Point2f p;
        p.x = x / (size.width() - 1.0f) * (x_range_[1] - x_range_[0]) +
            x_range_[0];
        p.y = (1.0f - y / (size.height() - 1.0f)) * (y_range_[1] - y_range_[0]) +
            y_range_[0];
        cout << "x_range: " << x_range_.transpose() << " Point: " << p << '\n';
        emit PointClicked(p);
      }
    });
  UpdatePoints(points_);
}

void Plot2d::SetRange(const Eigen::Vector2f &x_range,
                      const Eigen::Vector2f &y_range) {
  x_range_ = x_range;
  y_range_ = y_range;
  fixed_range_ = true;
}

void Plot2d::UpdatePoints(const vector<Point2f>& points) {
  points_.resize(1);
  points_[0] = points;
  DrawPoints();
}

void Plot2d::UpdatePoints(const vector<vector<Point2f>>& points) {
  points_ = points;
  DrawPoints();
}

void Plot2d::DrawPoints() {
  static const cv::Scalar colors[] = {kCvRed, kCvGreen, kCvBlue, kCvMagenta,
                                      kCvCyan, kCvYellow};
  static const int num_colors = 6;
  QSize size = myView->size();
  canvas_.create(size.height(), size.width(), CV_8UC3);
  canvas_ = CV_RGB(255, 255, 255);
  function<Point2f (const Point2f&)> change_point;
  Vector2f x_range, y_range;
  for (int group = 0; group < points_.size(); ++group) {
    // LOG(INFO) << "Draw points " << points_[0].size();
    if (points_[group].size() > 0) {
      if (fixed_range_) {
        x_range = x_range_;
        y_range = y_range_;
      } else {
        Vector2f lower(numeric_limits<float>::max(),
                       numeric_limits<float>::max());
        Vector2f upper(numeric_limits<float>::min(),
                       numeric_limits<float>::min());
        for (auto point : points_[group]) {
          if (point.x < lower[0]) {
            lower[0] = point.x;
          }
          if (point.x > upper[0]) {
            upper[0] = point.x;
          }
          if (point.y < lower[1]) {
            lower[1] = point.y;
          }
          if (point.y > upper[1]) {
            upper[1] = point.y;
          }
        }
        x_range[0] = lower[0];
        x_range[1] = upper[0];
        y_range[0] = lower[1];
        y_range[1] = upper[1];
      }
      LOG(INFO) << "Group: " << group << " Range: "
                << x_range.transpose() << ' '
                << y_range.transpose();
      x_range_ = x_range;
      y_range_ = y_range;
      change_point = [x_range, y_range, size](const Point2f &p) {
        return Point2f((p.x - x_range[0]) / (x_range[1] - x_range[0]) *
                       size.width(),
                       (1 - (p.y - y_range[0]) / (y_range[1] - y_range[0])) *
                       // (p.y - y_range[0]) / (y_range[1] - y_range[0]) *
                       size.height());
      };
      vector<Point2f> points_to_draw(points_[group].size());
      for (size_t i = 0; i < points_[group].size(); ++i) {
        points_to_draw[i] = change_point(points_[group][i]);
      }
      // LOG(INFO) << "Drawing "
      //     // << Join(points_to_draw.begin(), points_to_draw.end(), " ");
      //           << Join(points_[group].begin(), points_[group].end(), " ");
      for (size_t i = 1; i < points_[group].size(); ++i) {
        line(canvas_, points_to_draw[i-1], points_to_draw[i],
             colors[group % num_colors], 2);
      }
    }
  }

  UpdateImage(canvas_);
}

void Plot2d::resizeEvent(QResizeEvent *) {
  DrawPoints();
}

AdjustableSlider::AdjustableSlider(const std::string &title,
                                   bool int_type,
                                   QWidget *parent)
    : QWidget(parent), double_validator_(parent),
      int_validator_(parent), int_type_(int_type) {
  QFontMetrics metrics(QApplication::font());

  title_ = new QLabel(tr(title.c_str()));
  value_label_ = new QLabel(tr("Value"));
  value_text_ = new QLineEdit(QString::number(value_));
  max_label_ = new QLabel(tr("Max"));
  max_text_ = new QLineEdit(QString::number(max_value_));
  min_label_ = new QLabel(tr("Min"));
  min_text_ = new QLineEdit(QString::number(min_value_));

  slider_ = new QSlider(Qt::Horizontal);
  slider_->setFocusPolicy(Qt::StrongFocus);
  slider_->setTickPosition(QSlider::NoTicks);
  slider_->setTickInterval(10);
  slider_->setSingleStep(1);
  if (int_type_) {
    slider_->setMaximum(max_value_);
    slider_->setMinimum(min_value_);
  } else {
    slider_->setMaximum(num_steps_);
    slider_->setMinimum(0);
  }

  UpdateValue();

  QHBoxLayout *layout = new QHBoxLayout;
  layout->addWidget(title_);
  layout->addWidget(value_label_);
  layout->addWidget(value_text_);
  layout->addWidget(max_label_);
  layout->addWidget(max_text_);
  layout->addWidget(min_label_);
  layout->addWidget(min_text_);
  layout->addWidget(slider_);

  setLayout(layout);

  title_->setFixedWidth(metrics.width("88888888"));
  QLineEdit *texts[] = {value_text_, max_text_, min_text_};
  for (int i = 0; i < 3; ++i) {
    if (int_type_) {
      texts[i]->setValidator(&int_validator_);
    } else {
      texts[i]->setValidator(&double_validator_);
    }
    texts[i]->setAlignment(Qt::AlignCenter);
    texts[i]->setFixedWidth(metrics.width("88888888"));
  }

  connect(max_text_, SIGNAL(textEdited(const QString&)),
          this, SLOT(SetMaxText(const QString&)));
  connect(min_text_, SIGNAL(textEdited(const QString&)),
          this, SLOT(SetMinText(const QString&)));
  connect(value_text_, SIGNAL(textEdited(const QString&)),
          this, SLOT(SetValueText(const QString&)));
  connect(slider_, SIGNAL(sliderMoved(int)),
          this, SLOT(SetSliderValue(int)));
}

void AdjustableSlider::SetValueText(const QString & text) {
  SetValue(text.toDouble());
  emit ValueEdited(value_);
}

void AdjustableSlider::SetMaxText(const QString & text) {
  SetMaxValue(text.toDouble());
}

void AdjustableSlider::SetMinText(const QString & text) {
  SetMinValue(text.toDouble());
}

void AdjustableSlider::SetSliderValue(int value) {
  if (!int_type_) {
    value_ = (double)value / num_steps_ * (max_value_ - min_value_) + min_value_;
  } else {
    value_ = value;
  }
  value_text_->setText(QString::number(value_));
  emit ValueEdited(value_);
}

void AdjustableSlider::UpdateValue() {
  value_text_->setText(QString::number(value_));
  if (int_type_) {
    slider_->setValue(value_);
  } else {
    slider_->setValue((value_ - min_value_) / (max_value_ - min_value_) *
                      num_steps_);
  }
  emit ValueChanged(value_);
}

void AdjustableSlider::SetMaxValue(double max) {
  max_value_ = max;
  if (value_ > max_value_) {
    value_ = max_value_;
    // emit ValueEdited(value_);
  }
  max_text_->setText(QString::number(max_value_));
  if (int_type_) {
    slider_->setMaximum(max_value_);
  }
  UpdateValue();
}

void AdjustableSlider::SetMinValue(double min) {
  min_value_ = min;
  if (value_ < min_value_) {
    value_ = min_value_;
    // emit ValueEdited(value_);
  }
  min_text_->setText(QString::number(min_value_));
  if (int_type_) {
    slider_->setMinimum(min_value_);
  }
  UpdateValue();
}

void AdjustableSlider::SetValue(double value) {
  value_ = value;
  if (value_ < min_value_)
    value_ = min_value_;
  if (value_ > max_value_)
    value_ = max_value_;
  UpdateValue();
}

MvsController::MvsController(MvsModel* model, QWidget* parent)
    : QWidget(parent), model_(model) {

  setStyleSheet("background-color: white;");

  // patch_window_ = new CvWindow("patch");
  ref_window_ = new CvWindow("Reference");
  image_window_ = new CvWindow("Images");
  patch_window_ = new CvWindow("Patches");
  plot_window_ = new Plot2d("ncc");

  // plot_window_->SetRange(model->depth_range(), Vector2f(-1, 1));
  vector<Point2f> points;
  for (int i = 0; i < 20; ++i) {
    points.push_back(Point2f(i, i / 20.0f));
  }
  plot_window_->UpdatePoints(points);
  // patch_window_->resize(200, min(200 * model_->num_cameras(), 1200));

  ref_window_->UpdateImage(model_->GetRefImage());
  image_window_->UpdateImage(GridImages(model_->images()));
  patch_window_->UpdateImage(HListImages(model_->ExtractPatches()));

  ref_window_->SetMouseCallback(
      [&](int event, int x, int y, int, void*) {
        if( event != EVENT_LBUTTONDOWN )
          return;
        SetRefPoint(Point2f(x, y));
      });

  depth_slider_ = new AdjustableSlider("Depth");

  depth_slider_->SetMaxValue(20);
  depth_slider_->SetMinValue(0);
  depth_slider_->SetValue(model->ref_depth());
  QRadioButton *show_single_score = new QRadioButton(tr("&single"));
  QRadioButton *show_sum_score = new QRadioButton(tr("s&um"));
  show_single_score->setChecked(true);
  QVBoxLayout *score_layout = new QVBoxLayout;
  score_layout->addWidget(show_single_score);
  score_layout->addWidget(show_sum_score);
  QGroupBox *score_group = new QGroupBox(tr("Combine"));
  score_group->setLayout(score_layout);

  auto score_mapper = new QSignalMapper;
  score_mapper->setMapping(show_single_score, COMBINE_SINGLE);
  score_mapper->setMapping(show_sum_score, COMBINE_SUM);
  connect(show_single_score, SIGNAL(clicked()), score_mapper, SLOT(map()));
  connect(show_sum_score, SIGNAL(clicked()), score_mapper, SLOT(map()));
  connect(score_mapper, SIGNAL(mapped(int)),
          this, SLOT(SetCombineMode(int)));

  auto control_layout = new QVBoxLayout;
  control_layout->addWidget(depth_slider_);
  control_layout->addWidget(score_group);
  QWidget *control_container = new QWidget;
  control_container->setLayout(control_layout);

  // auto grid = new QGridLayout;
  // grid->addWidget(ref_window_, 0, 0);
  // grid->addWidget(image_window_, 0, 1);
  // grid->addWidget(patch_window_, 1, 0);
  // grid->addLayout(slider_layout, 1, 1);
  // grid->setRowMinimumHeight(1, 200);

  QSplitter *image_splitter = new QSplitter;
  image_splitter->setOrientation(Qt::Horizontal);
  image_splitter->addWidget(ref_window_);
  image_splitter->addWidget(image_window_);
  // image_splitter->addWidget(plot_window_);
  QSplitter *tool_splitter = new QSplitter;
  tool_splitter->setOrientation(Qt::Horizontal);
  tool_splitter->addWidget(patch_window_);
  tool_splitter->addWidget(control_container);
  QSplitter *section_splitter = new QSplitter;
  section_splitter->setOrientation(Qt::Vertical);
  section_splitter->addWidget(image_splitter);
  section_splitter->addWidget(tool_splitter);
  section_splitter->addWidget(plot_window_);

  QVBoxLayout *window_layout = new QVBoxLayout;
  window_layout->addWidget(section_splitter);

  // setLayout(slider_layout);
  setLayout(window_layout);

  connect(depth_slider_, SIGNAL(ValueEdited(double)),
          this, SLOT(SetRefDepth(double)));
  connect(plot_window_, SIGNAL(PointClicked(cv::Point2f)),
          this, SLOT(SetClickedPlotPoint(cv::Point2f)));
}

void MvsController::SetRefPoint(const cv::Point2f& point) {
  model_->SetRefPoint(point);
  auto marked_images = model_->MarkRefPointOnImages();
  ref_window_->UpdateImage(marked_images[model_->ref_camera_index()]);
  image_window_->UpdateImage(GridImages(marked_images));
  patch_window_->UpdateImage(HListImages(model_->ExtractPatches()));
  auto points = model_->GetPhotoConsistency();
  // LOG(INFO) << "Update points";
  switch (combine_mode_) {
    case COMBINE_SINGLE:
      plot_window_->UpdatePoints(points);
      break;
    case COMBINE_SUM:
      {
        vector<Point2f> sum_points(points[0].size());
          for (size_t i = 0; i < points[0].size(); ++i) {
            Point2f p(points[0][i].x, 0);
            for (size_t j = 0; j < points.size(); ++j) {
              p.y += points[j][i].y;
            }
            p.y /= points.size();
            sum_points[i] = p;
          }
          plot_window_->UpdatePoints(sum_points);
      }
      break;
  }

  // ShowPhotoConsistencyScores();
}

void MvsController::SetRefDepth(double d) {
  cout << "depth changed " << d << "\n";
  model_->set_ref_depth(d);
  image_window_->UpdateImage(GridImages(model_->MarkRefPointOnImages()));
  patch_window_->UpdateImage(HListImages(model_->ExtractPatches()));
  depth_slider_->SetValue(d);
  // ShowPhotoConsistencyScores();
}

void MvsController::ShowPhotoConsistencyScores() const {
  auto scores = model_->CalcPatchNcc();
  LOG(INFO) << "Ncc: " << Join(scores.begin(), scores.end(), " ");
}

void MvsController::SetClickedPlotPoint(const cv::Point2f &point) {
  SetRefDepth(point.x);
}

void MvsController::SetCombineMode(int m) {
  combine_mode_ = m;
  SetRefPoint(model_->GetRefPoint());
}

} // tiny
} // furry
