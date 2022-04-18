#include "MainWindow.h"

using namespace cagd;

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent)
{
    setupUi(this);

/*

  the structure of the main window's central widget

 *---------------------------------------------------*
 |                 central widget                    |
 |                                                   |
 |  *---------------------------*-----------------*  |
 |  |     rendering context     |   scroll area   |  |
 |  |       OpenGL widget       | *-------------* |  |
 |  |                           | | side widget | |  |
 |  |                           | |             | |  |
 |  |                           | |             | |  |
 |  |                           | *-------------* |  |
 |  *---------------------------*-----------------*  |
 |                                                   |
 *---------------------------------------------------*

*/
    _side_widget = new SideWidget(this);

    _scroll_area = new QScrollArea(this);
    _scroll_area->setWidget(_side_widget);
    _scroll_area->setSizePolicy(_side_widget->sizePolicy());
    _scroll_area->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

    _gl_widget = new GLWidget(this);

    centralWidget()->setLayout(new QHBoxLayout());
    centralWidget()->layout()->addWidget(_gl_widget);
    centralWidget()->layout()->addWidget(_scroll_area);

    // creating a signal slot mechanism between the rendering context and the side widget
    connect(_side_widget->rotate_x_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(setAngleX(int)));
    connect(_side_widget->rotate_y_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(setAngleY(int)));
    connect(_side_widget->rotate_z_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(setAngleZ(int)));

    connect(_side_widget->zoom_factor_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(setZoomFactor(double)));

    connect(_side_widget->trans_x_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(setTransX(double)));
    connect(_side_widget->trans_y_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(setTransY(double)));
    connect(_side_widget->trans_z_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(setTransZ(double)));

    connect(_side_widget->apply_reflection_lines_check_box, SIGNAL(clicked(bool)), _gl_widget, SLOT(toggleReflectionLines(bool)));
    connect(_side_widget->compare_control_nets_check_box, SIGNAL(clicked(bool)), _gl_widget, SLOT(setVisibilityOfInitialControlNet(bool)));

    connect(_side_widget->show_initial_B_surface_radio_button, SIGNAL(toggled(bool)), _gl_widget, SLOT(setVisibilityOfInitialBSurface(bool)));
    connect(_side_widget->show_order_elevated_B_surface_radio_button, SIGNAL(toggled(bool)), _gl_widget, SLOT(setVisibilityOfOrderElevatedBSurface(bool)));
    connect(_side_widget->show_subdivisions_in_direction_u_radio_button, SIGNAL(toggled(bool)), _gl_widget, SLOT(setVisibilityOfSubdivisionsInDirectionU(bool)));
    connect(_side_widget->show_subdivisions_in_direction_v_radio_button, SIGNAL(toggled(bool)), _gl_widget, SLOT(setVisibilityOfSubdivisionsInDirectionV(bool)));
}

//--------------------------------
// implementation of private slots
//--------------------------------
void MainWindow::on_action_Quit_triggered()
{
    qApp->exit(0);
}
