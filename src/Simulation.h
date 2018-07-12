#ifndef DEF_SIMULATION
#define DEF_SIMULATION
#include <QApplication>
#include <QPushButton>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QGroupBox>
#include <QTextEdit>
#include <QCheckBox>
#include <QDialog>
#include <QVector>
#include <QRadioButton>
#include "qcustomplot.h"
#include <QTextEdit>
#include <QTabWidget>
#include "muParser.h"

class simul : public QWidget
{
	Q_OBJECT

public:
		simul();

		public slots:
		
		void eval(); 
		void evalinfo();

private:
		QWidget *pageplot;
		QWidget *pagemesh; 
		QVBoxLayout *layoutV;
		QVBoxLayout *layoutplotmesh;
		QVBoxLayout *layoutplot;
		QHBoxLayout *layoutH3;
		QHBoxLayout *layoutH;
		QHBoxLayout *layoutH2;
		QHBoxLayout *layoutH4;
		QHBoxLayout *layoutBoundery1;
		QHBoxLayout *layoutBoundery2;
		QHBoxLayout *layoutBoundery3;
		QHBoxLayout *layoutBoundery4;
		QCustomPlot *customPlot;
		QCPColorScale *colorScale;
		QCustomPlot *customMesh;
		QFormLayout *formEDO;
		QFormLayout *formMesh;
		QFormLayout *formBp;
		QGroupBox *gEDO;
		QLineEdit *form1;
		QLineEdit *form1y;
		QLineEdit *form2;
		QLineEdit *form3; 
		QLineEdit *form4;
		QGroupBox *gCF;
		QVBoxLayout formCF;
		QLineEdit *form5;
		QLineEdit *form6;
		QLineEdit *form7;
		QLineEdit *form8;
		QLineEdit *form9;
		QLineEdit *form10;
		QLineEdit *form11;
		QLineEdit *form12;
		QLineEdit *form13;
		QLineEdit *form14;
		QGroupBox *gD;
		QGroupBox *glog;
		QGroupBox *gMesh;
		QGroupBox *gBp;
		QVBoxLayout *layoutlog; 
		QVBoxLayout *lD;
		QRadioButton *rb1;
		QRadioButton *rb2;
		QPushButton *but;
		QPushButton *butinfo;
		QTextEdit *log; 
		QDialog *dialog;
		mu::Parser *fmu;
		mu::Parser *axmu;
		mu::Parser *aymu;
		mu::Parser *bmu;
		mu::Parser *cmu;
		mu::Parser *B1mu;
		mu::Parser *B2mu;
		mu::Parser *B3mu;
		mu::Parser *B4mu;
		QPushButton *savemesh;
		QPushButton *saveplot;
		QPushButton *savedata;
		QPushButton *genmesh;
		QTabWidget *tabW;
		QComboBox *list1;
		QComboBox *list2;
		QComboBox *list3;
		QComboBox *list4;
		QProgressBar *pbar; 

		int iseval;
};


#endif 