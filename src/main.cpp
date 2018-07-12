
#include <iostream>
#include <QApplication>
#include "Simulation.h"
#include <QTextStream>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    simul *mafenetre = new simul;

    mafenetre->show();

    return app.exec();
}

