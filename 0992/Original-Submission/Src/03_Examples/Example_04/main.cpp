#include <QtWidgets/QApplication>
#include "GUI/MainWindow.h"

using namespace cagd;

int main(int argc, char **argv)
{
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling, true);
    QApplication app(argc, argv);
    app.setAttribute(Qt::AA_UseDesktopOpenGL, true);
    app.setStyle("fusion");

    MainWindow mwnd;
    mwnd.showMaximized();

    return app.exec();
}
