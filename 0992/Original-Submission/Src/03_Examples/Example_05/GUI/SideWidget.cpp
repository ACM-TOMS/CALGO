#include "SideWidget.h"

using namespace cagd;

SideWidget::SideWidget(QWidget *parent):
        QWidget(parent)
{
    setupUi(this);

    QPalette p = rotate_x_slider->palette();

    p.setColor(QPalette::Highlight, QColor(255,50,10).lighter());

    rotate_x_slider->setPalette(p);

    p = rotate_y_slider->palette();

    p.setColor(QPalette::Highlight, QColor(50,255,10).lighter());

    rotate_y_slider->setPalette(p);

}
