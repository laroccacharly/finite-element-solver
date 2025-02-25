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
#include "arma_config.hpp"  // Custom Armadillo config to enable SuperLU
#include <QVector>
#include <algorithm>
#include "Simulation.h"
#include "qcustomplot.h"
#include "muParser.h"
#include <ctime>
#include <cstdio>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


simul::simul() : QWidget()

{

	// Create all the objects
	fmu = new mu::Parser;
	axmu = new mu::Parser;
	aymu = new mu::Parser;
	bmu = new mu::Parser;
	cmu = new mu::Parser;
	B1mu = new mu::Parser;
	B2mu = new mu::Parser;
	B3mu = new mu::Parser;
	B4mu = new mu::Parser;
	QTabWidget *tabW = new QTabWidget;
	QWidget *pageplot = new QWidget;
	QWidget *pagemesh = new QWidget;
	QVBoxLayout  *layoutV = new QVBoxLayout;
	QVBoxLayout  *layoutplot = new QVBoxLayout;
	QVBoxLayout  *layoutplotmesh = new QVBoxLayout;
	QHBoxLayout  *layoutH = new QHBoxLayout;
	QHBoxLayout  *layoutH2 = new QHBoxLayout;
	QHBoxLayout  *layoutH3 = new QHBoxLayout;
	QHBoxLayout  *layoutH4 = new QHBoxLayout;
	QHBoxLayout *layoutBoundery1 = new QHBoxLayout;
	QHBoxLayout *layoutBoundery2 = new QHBoxLayout;
	QHBoxLayout *layoutBoundery3 = new QHBoxLayout;
	QHBoxLayout *layoutBoundery4 = new QHBoxLayout;
	formEDO = new QFormLayout;
	formMesh = new QFormLayout;
	QVBoxLayout  *formCF = new QVBoxLayout;
	QVBoxLayout  *lD = new QVBoxLayout;
	layoutlog = new QVBoxLayout;
	formBp = new QFormLayout;

	but = new QPushButton("Evaluate");
	butinfo = new QPushButton("?");
	savemesh = new QPushButton("Export Mesh info");
	saveplot = new QPushButton("Save plot as .png");
	genmesh = new QPushButton("Generate a mesh");
	savedata = new QPushButton("Export plot data");

	butinfo->setMaximumWidth(20);
	customPlot = new QCustomPlot();
	customPlot->setMinimumSize(700, 500);
	colorScale = new QCPColorScale(customPlot);


	customMesh = new QCustomPlot();
	customMesh->setMinimumSize(700, 500);


	gEDO = new QGroupBox("PDE parameters");
	gCF = new QGroupBox("Boundary values ");
	gD = new QGroupBox("Degree of precision");
	glog = new QGroupBox("Log");
	gMesh = new QGroupBox("Meshing parameters");
	gBp = new QGroupBox("Boundary positions");

	form1 = new QLineEdit("1");
	form1y = new QLineEdit("1");
	form2 = new QLineEdit("3");
	form3 = new QLineEdit("0");
	form4 = new QLineEdit("x^2+y^2");

	form5 = new QLineEdit("1");
	form6 = new QLineEdit("1");
	form7 = new QLineEdit("1");
	form8 = new QLineEdit("1");
	form9 = new QLineEdit("0.5");
	form10 = new QLineEdit("5");

	form11 = new QLineEdit("0");
	form12 = new QLineEdit("1");
	form13 = new QLineEdit("0");
	form14 = new QLineEdit("1");

	log = new QTextEdit();
	log->setMinimumWidth(300);
	rb1 = new QRadioButton("Degree 1", this);
	rb2 = new QRadioButton("Degree 2", this);

	list1 = new QComboBox;
	list2 = new QComboBox;
	list3 = new QComboBox;
	list4 = new QComboBox;

	pbar = new QProgressBar; 

	int iseval;

	list1->addItem("g1(x,y)"); list2->addItem("g2(x,y)"); list3->addItem("g3(x,y)"); list4->addItem("g4(x,y)");
	list1->addItem("h1(x,y)"); list2->addItem("h2(x,y)"); list3->addItem("h3(x,y)"); list4->addItem("h4(x,y)");

	list1->setToolTip("y=ymin, x=[xmin, xmax]");
	list2->setToolTip("y=[ymin, ymax], x=xmax");
	list3->setToolTip("y=ymax, x=[xmin, xmax]");
	list4->setToolTip("y=[ymin, ymax], x=xmin");

	// UI construction

	formEDO->addRow("ax : ", form1);
	formEDO->addRow("ay : ", form1y);
	formEDO->addRow("b : ", form2);
	formEDO->addRow("c : ", form3);
	formEDO->addRow("f : ", form4);

	
	layoutBoundery1->addWidget(list1);
	layoutBoundery1->addWidget(form5);
	layoutBoundery2->addWidget(list2);
	layoutBoundery2->addWidget(form6);
	layoutBoundery3->addWidget(list3);
	layoutBoundery3->addWidget(form7);
	layoutBoundery4->addWidget(list4);
	layoutBoundery4->addWidget(form8);
	formCF->addLayout(layoutBoundery1);
	formCF->addLayout(layoutBoundery2);
	formCF->addLayout(layoutBoundery3);
	formCF->addLayout(layoutBoundery4);

	formBp->addRow("xmin = ", form11);
	formBp->addRow("xmax = ", form12);
	formBp->addRow("ymin = ", form13);
	formBp->addRow("ymax = ", form14);

	formMesh->addRow("Size of elements: ", form9);
	formMesh->addRow("Density of points : ", form10);

	lD->addWidget(rb1);
	lD->addWidget(rb2);

	layoutlog->addWidget(log);

	gEDO->setLayout(formEDO);
	gCF->setLayout(formCF);
	gD->setLayout(lD);
	glog->setLayout(layoutlog);
	gMesh->setLayout(formMesh);
	gBp->setLayout(formBp);

	layoutV->addWidget(gMesh);
	layoutV->addWidget(gEDO);
	layoutV->addWidget(gBp);
	layoutV->addWidget(gCF);
	layoutV->addWidget(gD);
	layoutH4->addWidget(but);
	layoutH4->addWidget(butinfo);
	layoutV->addLayout(layoutH4);
	layoutV->addWidget(pbar);
	layoutH->addWidget(glog);
	layoutH->addLayout(layoutV);

	layoutplotmesh->addWidget(customMesh);
	layoutplotmesh->addLayout(layoutH3);

	layoutplot->addWidget(customPlot);
	layoutplot->addLayout(layoutH2);
	pageplot->setLayout(layoutplot);
	pagemesh->setLayout(layoutplotmesh);
	tabW->addTab(pagemesh, "Mesh");
	tabW->addTab(pageplot, "Plot");

	layoutH->addWidget(tabW);
	rb1->setChecked(1);

	this->setLayout(layoutH);
		
	// Signal - Slot connections
	QWidget::connect(but, SIGNAL(clicked()), this, SLOT(eval()));
	QWidget::connect(butinfo, SIGNAL(clicked()), this, SLOT(evalinfo()));


}


void simul::evalinfo()
{
	QDesktopServices::openUrl(QUrl("doc.pdf")); // Open the documentation
}

void simul::eval() // All the simulation is in this function
{
	
	try 
	{
		// Initialisation 
		log->clear(); // Clear the log
		int deg = 1;
		std::cout << rb1->isChecked() << std::endl;
		if (rb1->isChecked()) // Check the degree of precision
		{
			deg = 1;
		}
		else
		{
			deg = 2;
		}

		if (iseval == 1) // Check if it is the first time the code is executed
		{
			// nothing
		}
		else
		{
			iseval = 0;
		}

		// Log string 
		QString logstr;
		logstr += "-- This section shows when each step ended. Values are in seconds -- \n";
		std::stringstream buffer;

		// Create Mesh with CGAL
		CDT cdt;

		double xmin, xmax, ymin, ymax;
		xmin = form11->text().toDouble(); ymin = form13->text().toDouble();
		ymax = form14->text().toDouble(); xmax = form12->text().toDouble();
		Vertex_handle va = cdt.insert(Point(xmin, ymin));
		Vertex_handle vb = cdt.insert(Point(xmin, ymax));
		Vertex_handle vc = cdt.insert(Point(xmax, ymax));
		Vertex_handle vd = cdt.insert(Point(xmax, ymin));

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, vd);
		cdt.insert_constraint(vd, va);

		clock_t mesh_start = clock(); // Start clock

		Mesher mesher(cdt);
		double crit1 = form9->text().toDouble();// size of elements
		double crit2 = 0.125; // Fixed value found in CGAL documentation. 
		double eps = crit1 / 1000;  // eps is used to avoid floating-point arithmetic error

		mesher.set_criteria(Criteria(crit2, crit1));
		mesher.refine_mesh();

		clock_t mesh_stop = clock(); // Stop clock
		buffer << "Refine mesh : " << double(mesh_stop - mesh_start) / CLOCKS_PER_SEC << std::endl; // Evaluate the time to mesh

		// Extract the points generated by the mesh

		CDT::Point_iterator itp;
		int nbpt = cdt.number_of_vertices(); // Number of points in the mesh
		arma::mat coor(nbpt, 2); // The position of each point is saved in the "coor" matrix
		QVector<double> cx(nbpt);
		QVector<double> cy(nbpt);

		int i = 0;

		for (itp = cdt.points_begin(); itp != cdt.points_end(); itp++)
		{
			coor(i, 0) = itp->x();
			coor(i, 1) = itp->y();
			cx[i] = coor(i, 0);
			cy[i] = coor(i, 1);
			i++;

		}

		clock_t point_stop = clock();
		buffer << "Extract the points from the mesh : " << double(point_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;

		// Display mesh with QCustomPlot
		customMesh->clearItems();
		customMesh->clearFocus();
		customMesh->clearMask();
		customMesh->clearPlottables();
		customMesh->clearGraphs();
		customMesh->addGraph(0);
		customMesh->graph(0)->setData(cx, cy);
		customMesh->graph(0)->setLineStyle(QCPGraph::lsNone);
		QCPScatterStyle::ScatterShape shape = QCPScatterStyle::ssCircle;
		customMesh->graph(0)->setScatterStyle(QCPScatterStyle(shape, 7));
		customMesh->rescaleAxes(true);
		customMesh->replot();

		// Identify the points of each face/element
		int nbel = cdt.number_of_faces(); // Number of elements in the mesh
		double f1, f2, f3;
		arma::mat fx(nbel, deg * 3);
		arma::mat fy(nbel, deg * 3);
		arma::colvec fsize(nbel); // The size of each element is saved in the matrix "fsize". The size of an element is the length of his longest edge. 
		CDT::Face_iterator itf;
		i = 0;

		for (itf = cdt.faces_begin(); itf != cdt.faces_end(); itf++)
		{

			fx(i, 0) = itf->vertex(0)->point().x();
			fy(i, 0) = itf->vertex(0)->point().y();

			fx(i, 1) = itf->vertex(CDT::ccw(0))->point().x();
			fy(i, 1) = itf->vertex(CDT::ccw(0))->point().y();

			fx(i, 2) = itf->vertex(CDT::ccw(CDT::ccw(0)))->point().x();
			fy(i, 2) = itf->vertex(CDT::ccw(CDT::ccw(0)))->point().y();


			// Size of each face
			f1 = std::sqrt(std::pow(fx(i, 0) - fx(i, 1), 2) + std::pow(fy(i, 0) - fy(i, 1), 2));
			f2 = std::sqrt(std::pow(fx(i, 1) - fx(i, 2), 2) + std::pow(fy(i, 1) - fy(i, 2), 2));
			f1 = std::sqrt(std::pow(fx(i, 0) - fx(i, 2), 2) + std::pow(fy(i, 0) - fy(i, 2), 2));
			fsize(i) = std::max(f1, std::max(f2, f3));

			i++;
		}

		double sizemax = arma::max(fsize);
		double sizemin = arma::min(fsize);
		
		clock_t face_stop = clock();
		buffer << "Identify the points of each face : " << double(face_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;


		// Create the connnectivity matrix 
		arma::Mat<int> connec(nbel, 3);
		i = 0;

		for (int k = 0; k < nbel; k++)
		{
			for (int l = 0; l < 3; l++)
			{

				for (int m = 0; m < nbpt; m++)
				{
					if (coor(i, 0) == fx(k, l) && coor(i, 1) == fy(k, l))
					{
						connec(k, l) = i;

						break;
					}

					i++;
				}
				i = 0;
			}

		}


		// Add additionnal nodes for 2nd degree of precision 
		int tok = 0;
		int tok2 = 0;
		arma::colvec edx, edy;
		double vx1, vx2, vy1, vy2;
		int posed;
		arma::mat coeq(6, 6);

		if (deg == 2)
		{
			i = 0;
			CDT::Edge_iterator itedge;
			for (itedge = cdt.edges_begin(); itedge != cdt.edges_end(); itedge++)
			{

				posed = itedge->second;
				vx1 = itedge->first->vertex(CDT::ccw(posed))->point().x();
				vx2 = itedge->first->vertex(CDT::cw(posed))->point().x();

				vy1 = itedge->first->vertex(CDT::ccw(posed))->point().y();
				vy2 = itedge->first->vertex(CDT::cw(posed))->point().y();
				edx.resize(i + 1);
				edy.resize(i + 1);
				edx(i) = (vx1 + vx2) / 2;
				edy(i) = (vy1 + vy2) / 2;

				i++;

			}

			arma::colvec npointx, npointy;
			npointx = edx;
			npointy = edy;

			// npointx and npointy are the matrix that contain the position of each secondary nodes. 

			// Now, these nodes are added to the matrix coor. Therefore, the matrix connec has to be modified. 
			arma::mat npoint = arma::join_horiz(npointx, npointy);
			coor = arma::join_vert(coor, npoint);
			arma::Mat<int> z3(nbel, 3); z3.fill(0);
			connec = arma::join_horiz(connec, z3);
			for (int k = 0; k < nbel; k++)
			{
				for (int i = 0; i < npointx.size(); i++)
				{
					if (((coor(connec(k, 0), 0) + coor(connec(k, 1), 0)) / 2) == npointx(i) && ((coor(connec(k, 0), 1) + coor(connec(k, 1), 1)) / 2) == npointy(i))
					{
						connec(k, 3) = i + nbpt;
					}
					else if (((coor(connec(k, 1), 0) + coor(connec(k, 2), 0)) / 2) == npointx(i) && ((coor(connec(k, 1), 1) + coor(connec(k, 2), 1)) / 2) == npointy(i))
					{
						connec(k, 4) = i + nbpt;
					}
					else if (((coor(connec(k, 0), 0) + coor(connec(k, 2), 0)) / 2) == npointx(i) && ((coor(connec(k, 0), 1) + coor(connec(k, 2), 1)) / 2) == npointy(i))
					{
						connec(k, 5) = i + nbpt;
					}
				}
			}


			// Parameters for interpolation functions

			coeq = {
				{1, -3, -3, 4, 2, 2},
				{0, -1, 0, 0, 2, 0},
				{0, 0, -1, 0, 0, 2},
				{0, 4, 0, -4, -4, 0},
				{0, 0, 0, 4, 0, 0},
				{0, 0, 4, -4, 0, -4}
			};

			// Replot the mesh with the new nodes
			QVector<double> px(npointx.size()), py(npointx.size());
			for (int i = 0; i < npointx.size(); i++)
			{
				px[i] = npointx(i);
				py[i] = npointy(i);

			}
			// Display mesh

			customMesh->addGraph();
			customMesh->graph(1)->setData(px, py);
			customMesh->graph(1)->setLineStyle(QCPGraph::lsNone);
			QCPScatterStyle::ScatterShape shape = QCPScatterStyle::ssCross;
			customMesh->graph(1)->setScatterStyle(QCPScatterStyle(shape, 7));
			customMesh->rescaleAxes(true);
			customMesh->replot();


		}

		clock_t connec_stop = clock();
		buffer << "Create connec and coor: " << double(connec_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;
		int nddl = coor.n_rows; // Number of nodes
	

		// muParser
		double fVal1, fVal2;
		double BVal1, BVal2;

		fmu->DefineVar("x", &fVal1);
		fmu->DefineVar("y", &fVal2);
		axmu->DefineVar("x", &fVal1);
		axmu->DefineVar("y", &fVal2);
		aymu->DefineVar("x", &fVal1);
		aymu->DefineVar("y", &fVal2);
		bmu->DefineVar("x", &fVal1);
		bmu->DefineVar("y", &fVal2);
		cmu->DefineVar("x", &fVal1);
		cmu->DefineVar("y", &fVal2);
		B1mu->DefineVar("x", &BVal1);
		B1mu->DefineVar("y", &BVal2);
		B2mu->DefineVar("x", &BVal1);
		B2mu->DefineVar("y", &BVal2);
		B3mu->DefineVar("x", &BVal1);
		B3mu->DefineVar("y", &BVal2);
		B4mu->DefineVar("x", &BVal1);
		B4mu->DefineVar("y", &BVal2);

		fmu->SetExpr(form4->text().toLocal8Bit().constData());
		axmu->SetExpr(form1->text().toLocal8Bit().constData());
		aymu->SetExpr(form1y->text().toLocal8Bit().constData());
		bmu->SetExpr(form2->text().toLocal8Bit().constData());
		cmu->SetExpr(form3->text().toLocal8Bit().constData());
		B1mu->SetExpr(form5->text().toLocal8Bit().constData());
		B2mu->SetExpr(form6->text().toLocal8Bit().constData());
		B3mu->SetExpr(form7->text().toLocal8Bit().constData());
		B4mu->SetExpr(form8->text().toLocal8Bit().constData());
		double feval, axeval, ayeval, beval, ceval, B1eval, B2eval, B3eval, B4eval, heval;

		// Dirichlet boundary condition

		tok = 0;
		tok2 = 0;
		QVector<double> test(2);
		arma::colvec ndsDDLConnus;
		arma::colvec ndsDDLInc;
		for (int i = 0; i < nddl; i++)
		{
			if ((coor(i, 0) == xmin && list4->currentIndex() == 0) || (coor(i, 0) == xmax && list2->currentIndex() == 0) || (coor(i, 1) == ymin && list1->currentIndex() == 0) || (coor(i, 1) == ymax && list3->currentIndex() == 0))
			{
				ndsDDLConnus.resize(tok + 1);
				ndsDDLConnus(tok) = i;

				tok++;
			}
			else
			{
				ndsDDLInc.resize(tok2 + 1);
				ndsDDLInc(tok2) = i;
				tok2++;
			}

		}


		// Create numer and adres

		arma::colvec numer(nddl);

		for (i = 0; i < ndsDDLInc.size(); i++)
		{
			numer(ndsDDLInc(i)) = i;
		}
		tok = 0;
		tok = numer.max() + 1;

		for (i = 0; i < ndsDDLConnus.size(); i++)
		{
			numer(ndsDDLConnus(i)) = tok;
			tok++;
		}


		arma::mat adres(nbel, deg * 3);

		for (int k = 0; k < nbel; k++)
		{
			for (int l = 0; l < deg * 3; l++)
			{
				adres(k, l) = numer(connec(k, l));

			}

		}

		clock_t table_stop = clock();
		buffer << "Create adres and numer: " << double(table_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;


		// Identify each point that belongs to the Dirichlet boundary condition
		int nbinc = ndsDDLInc.size();

		arma::colvec ddlConnus(ndsDDLConnus.size());
		
		for (int i = 0; i < ndsDDLConnus.size(); i++)
		{
			BVal1 = coor(ndsDDLConnus(i), 0);
			BVal2 = coor(ndsDDLConnus(i), 1);
			if (coor(ndsDDLConnus(i), 0) == xmin && list4->currentIndex() == 0)
			{
				ddlConnus(i) = B4mu->Eval();
			}
			else if (coor(ndsDDLConnus(i), 0) == xmax && list2->currentIndex() == 0)
			{
				ddlConnus(i) = B2mu->Eval();
			}
			else if (coor(ndsDDLConnus(i), 1) == ymin &&  list1->currentIndex() == 0)
			{
				ddlConnus(i) = B1mu->Eval();
			}
			else if (coor(ndsDDLConnus(i), 1) == ymax &&  list3->currentIndex() == 0)
			{
				ddlConnus(i) = B3mu->Eval();

			}


		}

		// A few parameters
		arma::rowvec zd, nd, e, we, xi, yi, wi;
		zd = {-1, 1, 0};
		nd = {-1, 0, 1};
		e = {-0.774596669241483, 0, 0.774596669241483};
		we = {0.555555555555556, 0.888888888888889, 0.555555555555556};
		xi = {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.816847572980459, 0.091576213509771, 0.091576213509771};
		yi = {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.091576213509771, 0.816847572980459, 0.091576213509771};
		wi = {0.1116907948390055, 0.1116907948390055, 0.1116907948390055, 0.0549758718276610, 0.0549758718276610, 0.0549758718276610};


		// Interpolation functions and elementary matrix
		int nC = deg * 3;
		arma::sp_mat A(nddl, nddl);
		A.zeros();
		arma::colvec F(nddl);
		F.zeros();
		arma::colvec S(nddl);
		S.zeros();
		arma::colvec xk(nC);
		arma::colvec yk(nC);
		double Jk, chi, eta, transx, transy, ax, ay;
		arma::mat G(2, 2);
		double s1n, s2n, s3n;

		arma::colvec IntS(nC), IntS1(nC), IntS2(nC), IntS3(nC);
		arma::colvec S1(3), S2(3), S3(3);
		arma::colvec psi(nC), IntF(nC);
		arma::mat IntA(nC, nC), Fk(nC, nbel);
		arma::cube Intz(nC, nC, nbel);
		arma::colvec v1, v2, v3;
		v1 = {0, -1};
		v2 = {-1, 0};
		v3 = {1, 1};
		int tokbar = 0;
		arma::colvec Fi(wi.size()), Aij(wi.size());


		for (int k = 0; k < nbel; k++)
		{
			pbar->setValue(std::floor(100*k/nbel));
			
			for (int i = 0; i < 3 * deg; i++)
			{
				xk(i) = coor(connec(k, i), 0);
				yk(i) = coor(connec(k, i), 1);

			}

			Jk = (xk(1) - xk(0))*(yk(2) - yk(0)) - (yk(0) - yk(1))*(xk(0) - xk(2));
			G(0, 0) = yk(2) - yk(0);
			G(0, 1) = yk(0) - yk(1);
			G(1, 0) = xk(0) - xk(2);
			G(1, 1) = xk(1) - xk(0);
			s1n = arma::norm(G*v1);
			s3n = arma::norm(G*v2);
			s2n = arma::norm(G*v3*(1 / sqrt(2)));


			for (int i = 0; i < nC; i++)
			{
				for (int j = 0; j < nC; j++)
				{
				
					for (int l = 0; l < we.size(); l++)

					{
						// Neumann boundary condition
						if (nC == 3)
						{
							for (int nat = 0; nat < 3; nat++)
							{


								if (nat == 0)
								{
									chi = (e(l) + 1) / 2; eta = 0;
								}
								else if (nat == 1)
								{
									chi = 1 - (e(l) + 1) / 2; eta = (e(l) + 1) / 2;
								}
								else if (nat == 2)
								{
									chi = 0; eta = 1 - (e(l) + 1) / 2;
								}
								psi(0) = 1 - chi - eta; psi(1) = chi; psi(2) = eta;
								transx = xk(0)*psi(0) + xk(1)*psi(1) + xk(2)*psi(2);
								transy = yk(0)*psi(0) + yk(1)*psi(1) + yk(2)*psi(2);

								BVal1 = transx;
								BVal2 = transy;

								if (std::abs(BVal1 - xmin) < eps  && list4->currentIndex() == 1)
								{
									heval = B4mu->Eval();
								}
								else if (std::abs(BVal1 - xmax) < eps && list2->currentIndex() == 1)
								{
									heval = B2mu->Eval();
								}
								else if (std::abs(BVal2 - ymin) < eps &&  list1->currentIndex() == 1)
								{
									heval = B1mu->Eval();
								}
								else if (std::abs(BVal2 - ymax) < eps   &&  list3->currentIndex() == 1)
								{
									heval = B3mu->Eval();

								}
								else
								{
									heval = 0;
								}


								if (nat == 0)
								{
									S1(l) = psi(i)*heval*s1n * 1 / 2;
								}
								else if (nat == 1)
								{
									S2(l) = psi(i)*heval*s2n * 1 / sqrt(2);
								}
								else if (nat == 2)
								{
									S3(l) = psi(i)*heval*s3n * 1 / 2;
								}

							}
						}
						else if (nC == 6)

						{

							for (int nat = 0; nat < 3; nat++)
							{


								if (nat == 0)
								{
									chi = (e(l) + 1) / 2; eta = 0;
								}
								else if (nat == 1)
								{
									chi = 1 - (e(l) + 1) / 2; eta = (e(l) + 1) / 2;
								}
								else if (nat == 2)
								{
									chi = 0; eta = 1 - (e(l) + 1) / 2;
								}
								arma::rowvec tov1(6);
								tov1 = {1, chi, eta, chi*eta, chi*chi, eta*eta};
								for (int m = 0; m < nC; m++)
								{
									psi(m) = arma::sum(coeq.row(m) % tov1);
								}
								transx = arma::sum(xk % psi);
								transy = arma::sum(yk % psi);
				


								BVal1 = transx;
								BVal2 = transy;

								if (std::abs(BVal1 - xmin) < eps  && list4->currentIndex() == 1)
								{
									heval = B4mu->Eval();
								}
								else if (std::abs(BVal1 - xmax) < eps && list2->currentIndex() == 1)
								{
									heval = B2mu->Eval();
								}
								else if (std::abs(BVal2 - ymin) < eps &&  list1->currentIndex() == 1)
								{
									heval = B1mu->Eval();
								}
								else if (std::abs(BVal2 - ymax) < eps   &&  list3->currentIndex() == 1)
								{
									heval = B3mu->Eval();

								}
								else
								{
									heval = 0;
								}

								if (nat == 0)
								{
									S1(l) = psi(i)*heval*s1n * 1 / 2;


								}
								else if (nat == 1)
								{
									S2(l) = psi(i)*heval*s2n * 1 / sqrt(2);
								}
								else if (nat == 2)
								{
									S3(l) = psi(i)*heval*s3n * 1 / 2;
								}





							}



						}
					}
					IntS1(i) = arma::sum(we % S1.t());
					IntS2(i) = arma::sum(we % S2.t());
					IntS3(i) = arma::sum(we % S3.t());

					IntS(i) = IntS1(i) + IntS2(i) + IntS3(i);

					for (int z = 0; z < wi.size(); z++)
					{

						chi = xi(z);
						eta = yi(z);

						if (nC == 3)
						{
							psi(0) = 1 - chi - eta;
							psi(1) = chi;
							psi(2) = eta;
							transx = xk(0)*psi(0) + xk(1)*psi(1) + xk(2)*psi(2);
							transy = yk(0)*psi(0) + yk(1)*psi(1) + yk(2)*psi(2);
						}
						else if (nC == 6)
						{
							nd.resize(6);
							zd.resize(6);
							arma::rowvec tov1(6), tov2(6), tov3(6);
							tov1 = {1, chi, eta, chi*eta, chi*chi, eta*eta};
							tov2 = {0, 1, 0, eta, 2 * chi, 0};
							tov3 = {0, 0, 1, chi, 0, 2 * eta};
							for (int m = 0; m < nC; m++)
							{
								psi(m) = arma::sum(coeq.row(m) % tov1);
								zd(m) = arma::sum(coeq.row(m) % tov2);
								nd(m) = arma::sum(coeq.row(m) % tov3);

							}

							transx = arma::sum(xk.t() % psi.t());
							transy = arma::sum(yk.t() % psi.t());

						}


						fVal1 = transx; 
						fVal2 = transy;
						axeval = axmu->Eval();
						ayeval = aymu->Eval();
						ax = axeval;
						ay = ayeval;

						feval = fmu->Eval();
						Fi(z) = feval*psi(i)*Jk;

						beval = bmu->Eval();
						ceval = cmu->Eval();
						Aij(z) = (beval*
							(1 / (Jk *Jk)*(((yk(2) - yk(0))*zd(j) + (yk(0) - yk(1))*nd(j))*((yk(2) - yk(0))*zd(i) + (yk(0) - yk(1))*nd(i))))) +
							ceval*psi(i)*psi(j) +
							1 / Jk*(psi(i))*(ax*((yk(2) - yk(0))*zd(j) + (yk(0) - yk(1))*nd(j)) + ay*((xk(0) - xk(2))*zd(j) + (xk(1) - xk(0))*nd(j)))*Jk;
			

					}

					IntF(i) = arma::sum(wi.t() % Fi);
					IntA(i, j) = arma::sum(wi.t() % Aij);

					A(adres(k, i), adres(k, j)) += IntA(i, j);


				}

				F(adres(k, i)) += IntF(i);
				S(adres(k, i)) += IntS(i);
			}


		}

		clock_t cal_stop = clock();
		buffer << "Create elementary matrix: " << double(cal_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;


		// Solve the system with spsolve
		arma::colvec Fc = F(arma::span(0, nbinc - 1));
		arma::colvec Sc = S(arma::span(0, nbinc - 1));
		arma::colvec UI = arma::spsolve(A(arma::span(0, nbinc - 1), arma::span(0, nbinc - 1)), Fc + Sc - A(arma::span(0, nbinc - 1), arma::span(nbinc, nddl - 1))*ddlConnus);
		arma::colvec Ut = arma::join_vert(UI, ddlConnus);

		clock_t sol_stop = clock();
		buffer << "Solve the system : " << double(sol_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;


		// Now we display the solution in the original plan
		// For each element, the solution is calculated n(n+1)/2 times where n= disc. The dots are uniformly distributed within the triangular element
		int disc = form10->text().toDouble() + 1; // resolution 

		int varx = 0;
		for (i = 1; i < disc; i++)
		{
			varx += i;
		}
		arma::mat coort(varx, 2);
		coort.fill(0);

		arma::colvec tok6(disc);
		tok6 = arma::linspace<arma::colvec>(0, 1 - (0 + 1 - 1.0) / (disc - 1), disc - (0 + 1 - 1));
		tok6.resize(disc - 1);
		coort = arma::join_horiz(tok6, tok6(0)*arma::ones(tok6.size(), 1));
		arma::colvec tok7;
		for (i = 1; i < disc - 1; i++)
		{
			tok7 = arma::linspace<arma::colvec>(0, 1 - (i + 1 - 1.0) / (disc - 1), disc - (i + 1 - 1));
			tok7.resize(tok7.size() - 1);
			coort = arma::join_vert(coort, arma::join_horiz(tok7, tok6(i)*arma::ones(tok7.size(), 1)));
		}

		coort = coort + 1.0 / ((disc - 1) * 2);

		double u1, u2, u3, u4, u5, u6, psi1, psi2, psi3, psi4, psi5, psi6, dx, dy;

		int tok9;
		dx = 1.5*sizemax / disc;
		dy = 1.5*sizemax / disc;
		int grx = (xmax - xmin) / dx;
		int gry = (ymax - ymin) / dy;

		arma::field<arma::colvec> rawdata(grx, gry);

		arma::colvec ue(varx* nbel), uposx(varx* nbel), uposy(varx* nbel);
		int tok8 = 0;
		for (int k = 0; k < nbel; k++)
		{
			for (int i = 0; i < 3 * deg; i++)
			{
				xk(i) = coor(connec(k, i), 0);
				yk(i) = coor(connec(k, i), 1);
			}

			Jk = (xk(1) - xk(0))*(yk(2) - yk(0)) - (yk(0) - yk(1))*(xk(0) - xk(2));
			G(0, 0) = yk(2) - yk(0);
			G(0, 1) = yk(0) - yk(1);
			G(1, 0) = xk(0) - xk(2);
			G(1, 1) = xk(1) - xk(0);



			u1 = Ut(numer(connec(k, 0)));
			u2 = Ut(numer(connec(k, 1)));
			u3 = Ut(numer(connec(k, 2)));

			if (nC == 6)
			{
				u4 = Ut(numer(connec(k, 3)));
				u5 = Ut(numer(connec(k, 4)));
				u6 = Ut(numer(connec(k, 5)));
			}

			for (int l = 0; l < varx; l++)
			{

				chi = coort(l, 0);
				eta = coort(l, 1);

				if (nC == 3)
				{
					psi1 = 1 - chi - eta;
					psi2 = chi;
					psi3 = eta;


					ue(tok8) = u1*psi1 + u2*psi2 + u3*psi3;
					uposx(tok8) = psi1*xk(0) + psi2*xk(1) + psi3*xk(2);
					uposy(tok8) = psi1*yk(0) + psi2*yk(1) + psi3*yk(2);
				}
				else if (nC == 6)

				{
					arma::rowvec tov1(6);
					tov1 = {1, chi, eta, chi*eta, chi*chi, eta*eta};
					psi1 = arma::sum(coeq.row(0) % tov1);
					psi2 = arma::sum(coeq.row(1) % tov1);
					psi3 = arma::sum(coeq.row(2) % tov1);
					psi4 = arma::sum(coeq.row(3) % tov1);
					psi5 = arma::sum(coeq.row(4) % tov1);
					psi6 = arma::sum(coeq.row(5) % tov1);


					ue(tok8) = u1*psi1 + u2*psi2 + u3*psi3 + u4*psi4 + u5*psi5 + u6*psi6;
					uposx(tok8) = psi1*xk(0) + psi2*xk(1) + psi3*xk(2) + psi4*xk(3) + psi5*xk(4) + psi6*xk(5);
					uposy(tok8) = psi1*yk(0) + psi2*yk(1) + psi3*yk(2) + psi4*yk(3) + psi5*yk(4) + psi6*yk(5);


				}


				for (int xp = 0; xp < grx; xp++)
				{
					for (int yp = 0; yp < gry; yp++)
					{

						if (uposx(tok8) >= (xmin + dx*xp) && uposx(tok8) <= (xmin + dx*(xp + 1)) &&
							uposy(tok8) >= (ymin + dy*yp) && uposy(tok8) <= (ymin + dy*(yp + 1)))
						{
							tok9 = rawdata(xp, yp).size();
							rawdata(xp, yp).resize(tok9 + 1);
							rawdata(xp, yp).at(tok9) = ue(tok8);

						}

					}
				}


				tok8++;

			}


		}
		clock_t return_stop = clock();
		buffer << "Return to the original plan : " << double(return_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;

		// configure axis rect:
		customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom); // this will also allow rescaling the color scale by dragging/zooming
		customPlot->axisRect()->setupFullAxesBox(true);
		customPlot->xAxis->setLabel("x");
		customPlot->yAxis->setLabel("y");

		// set up the QCPColorMap:

		QCPColorMap *colorMap = new QCPColorMap(customPlot->xAxis, customPlot->yAxis);
		customPlot->addPlottable(colorMap);
		int nx = grx;
		int ny = gry;
		colorMap->data()->setSize(nx, ny); // we want the color map to have nx * ny data points
		colorMap->data()->setRange(QCPRange(xmin, xmax), QCPRange(ymin, ymax)); 
		// now we assign some data, by accessing the QCPColorMapData instance of the color map:
		double x, y, z;
		for (int xIndex = 0; xIndex < nx; ++xIndex)
		{
			for (int yIndex = 0; yIndex < ny; ++yIndex)
			{

				z = arma::mean(rawdata(xIndex, yIndex));
				colorMap->data()->setCell(xIndex, yIndex, z);
			}
		}

		if (iseval==0)
		{

		customPlot->plotLayout()->addElement(0, 1, colorScale); 
		}
		colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)
		colorMap->setColorScale(colorScale); // associate the color map with the color scale
		colorScale->axis()->setLabel("Legend");
		//}
		// set the color gradient of the color map to one of the presets:
		colorMap->setGradient(QCPColorGradient::gpPolar);
		// we could have also created a QCPColorGradient instance and added own colors to
	
		// rescale the data dimension (color) such that all data points lie in the span visualized by the color gradient:
		QCPRange inter = QCPRange(arma::min(ue), arma::max(ue));
		colorMap->setDataRange(inter);
		colorMap->setInterpolate(1);

		// make sure the axis rect and color scale synchronize their bottom and top margins (so they line up):
		QCPMarginGroup *marginGroup = new QCPMarginGroup(customPlot);
		customPlot->axisRect()->setMarginGroup(QCP::msBottom | QCP::msTop, marginGroup);
		colorScale->setMarginGroup(QCP::msBottom | QCP::msTop, marginGroup);

		// rescale the key (x) and value (y) axes so the whole color map is visible:
		customPlot->rescaleAxes();

		customPlot->replot();

		clock_t disp_stop = clock();
		buffer << "Display the solution : " << double(disp_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;


		customPlot->savePng("plot.png"); // Save the plot


		// Save the data in .csv files
		std::ofstream myFlux("data.csv");
		std::ofstream myFluxconnec("Mesh_connec.csv");
		std::ofstream myFluxcoor("Mesh_coor.csv");
		if (myFlux && myFluxcoor && myFluxconnec)
		{
			 myFlux << "x" << "," <<"y" << ","<< "u" << std::endl; 
			 for (int i = 0; i < uposx.size(); i++)
			 {
				 myFlux << uposx(i) << "," << uposy(i) <<"," << ue(i) << std::endl;
			 }

			 for (int i = 0; i < nbel; i++)
			 {
				 if (nC == 3)
				 {
					 myFluxconnec << connec(i, 0) << "," << connec(i, 1) << "," << connec(i, 2) << std::endl;
				 }

				 if (nC == 6)
				 {
					 myFluxconnec << connec(i, 0) << "," << connec(i, 1) << "," << connec(i, 2) << "," << connec(i, 3) << "," << connec(i, 4) << "," << connec(i, 5) << std::endl;
				 }

			 }

			 for (int i = 0; i < nddl; i++)
			 {
				
					 myFluxcoor << coor(i, 0) << "," << coor(i, 1)  << std::endl;
				 
			 }

		}
		else
		{
			throw QString("Error with the .csv file");
		}

		
		clock_t data_stop = clock();
		buffer << "Export the data and the plot : " << double(data_stop - mesh_start) / CLOCKS_PER_SEC << std::endl;


		// Add info in the log
		buffer << "\n-- Meshing information -- \n" << std::endl;
		buffer << "Maximum element size : " << sizemax << std::endl << "Minimum element size : " << sizemin << std::endl;
		buffer << "More information about the Mesh can be found in Mesh_coor.csv and Mesh_connec.csv " << std::endl;
		buffer << "-- " << std::endl;

		QString bufferq = QString::fromStdString(buffer.str());
		logstr += bufferq;
		log->setText(logstr); // Display the log 
		iseval = 1;  // The code was executed 


	}
	catch (mu::Parser::exception_type &e)
	{
		
		QString errormu = QString::fromStdString(e.GetMsg());
		QMessageBox::information(this, "Error", errormu);
	}
	catch (QString & errorfile)
	{
		QMessageBox::information(this, "Error",errorfile);
	}
	catch (std::logic_error &e)
	{
		QMessageBox::information(this, "Error",e.what());
	}
	catch (std::runtime_error &e)
	{
		QMessageBox::information(this, "Error", "runtime_error");
	}
	catch (...)
	{
		QMessageBox::information(this, "Error", "Error unknown");
	}
	pbar->setValue(100);
	QMessageBox::information(this, "Information", "Done!");


}

