#ifndef SIDEDWIDGET_H
#define SIDEDWIDGET_H

#include <QtWidgets/QWidget>
#include "ui_SideWidget.h"

namespace cagd
{
    class SideWidget: public QWidget, public Ui::SideWidget
    {
    public:
        // special and default constructor
        SideWidget(QWidget *parent = 0);
    };
}

#endif // SIDEDWIDGET_H
