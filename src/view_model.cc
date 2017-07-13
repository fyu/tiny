const char *kUsageMessage = "View model";

#include "furry/common/init.h"

#include "model_viewer.h"

using namespace furry;
using namespace furry::tiny;


DEFINE_string(in_model, "", "");

int main(int argc, char *argv[]) {
  furry::Init(&argc, &argv, kUsageMessage);
  CHECK(!FLAGS_in_model.empty());

  Model model;
  model.ReadFile(FLAGS_in_model);

  QApplication app(argc, argv);
  ModelViewer model_viewer(&model);
  model_viewer.show();
  app.exec();

  return 0;
}
