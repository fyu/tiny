qt4_add_resources(GUI_WINDOW_QT_RCCS window_QT.qrc)
qt4_wrap_cpp(GUI_WINDOW_QT_MOC window_QT.h)

add_library(furry_gui_cvwindow ${GUI_WINDOW_QT_MOC} ${GUI_WINDOW_QT_RCCS} window_QT.cpp)
