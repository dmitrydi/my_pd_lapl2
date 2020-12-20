/*
 * lapl_well.cpp
 *
 *  Created on: 20 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "lapl_well.h"

using namespace Eigen;

LaplWell::LaplWell(const double xwd_, const double xed_, const double ywd_, const double yed_, const double Fcd_):
	xwd(xwd_), xed(xed_), ywd(ywd_), yed(yed_), Fcd(Fcd_) {};
double LaplWell::pd(const double u) const {
	MatrixXd m = MakeMatrix(u, nseg, eps,
			xwd, xed, ywd, yed);
	VectorXd r = MakeRhs(u, nseg, Fcd);
	return m.colPivHouseholderQr().solve(r)(0);
};

