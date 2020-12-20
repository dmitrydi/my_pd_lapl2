//============================================================================
// Name        : my_pd_lapl2.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <limits>
#include "test_runner.h"
#include "profile.h"
#include "auxillary.h"
#include "tests.h"
#include "bessel.h"
#include "quadrature.h"
#include <algorithm>
#include <complex>
#include "matrix_func.h"
#include <Eigen/Dense>
#include "lapl_well.h"
#include "real_well.h"
using namespace std;

using Eigen::MatrixXd;

int main() {
	TestRunner tr;
	//RUN_TEST(tr, TestIbessk0::Precision_f);
	//RUN_TEST(tr, TestIbessk0::Precision_ab);
	//RUN_TEST(tr, TestLevin::Stability); //passed
	//RUN_TEST(tr, TestSexp::ExhaustiveCorrect); //passed
	//RUN_TEST(tr, TestSexp::ExhaustiveSpeed); //passed
	//RUN_TEST(tr, TestSexp::CheckSpurious);//passed
	//RUN_TEST(tr, TestSexp::Correct); //passed
	//RUN_TEST(tr, TestWijn::Correct); //passed
	//RUN_TEST(tr, TestiF2E::Correct); //passed
	//RUN_TEST(tr, TestIbessk0::PrintResults); //passed
	//RUN_TEST(tr, TestIbessk0::Speed);
	//RUN_TEST(tr, TestIbessk0::PrintQuad);
	//RUN_TEST(tr, TestIbessk0::Precision_ab);
	//RUN_TEST(tr, TestBessel::Speed);
	//RUN_TEST(tr, TestBessel::Precision);
	//RUN_TEST(tr, TestIntegrals::QrombOverIbessk0); //x <= 10 correct;
	//RUN_TEST(tr, TestIntegrals::Speed);

	//double u = 1.;
	int nseg = 20;
	double xed = 10;
	double xwd = 0.5*xed;
	double yed = 10;
	double ywd = 0.5*yed;
	double Fcd = 10.;
	//double td = 1.;
	RealWell well(xwd, xed, ywd, yed, Fcd);
	vector<double> tds = LogSpaced(0.1, 10., 10);
	for (auto td: tds) {
		cout <<"td = " << td << " pd = " << well.pd(td) << endl;
	}

}
