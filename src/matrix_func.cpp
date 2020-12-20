/*
 * matrix_func.cpp
 *
 *  Created on: 19 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "matrix_func.h"

MatrixXd MakeMatrix(const double u, const int nseg, const double eps,
		const double xwd, const double xed, const double ywd, const double yed) {
	MatrixXd ans(2*nseg + 1, 2*nseg + 1);
	double dx = 1./nseg;
	double yd = ywd;
	LaplFunc::iF1 _if1(u, ywd, yed, eps);
	LaplFunc::iF2E _if2E(u, xwd, xed, xed, 0., ywd, yed, eps);
	LaplFunc::i1F2H _i1f2H(u, xwd, xed, xed, 0., ywd, yed, eps);
	LaplFunc::i2F2H _i2f2H(u, ywd, 0.);
	double mult = -1.*LaplFunc::PI/xed;
	//+ _if2E(x1, x2, xd, yd) + _i1f2H(x1, x2, xd, yd) + _i2f2H(x1, x2, yd))
	for (int i = 0; i < 2*nseg; ++i) {
		ans(i,0) = 1.;
		ans(2*nseg, i+1) = 1.;
	}
	//{LOG_DURATION("iF1");
		for (int i = 0; i < 2*nseg; ++i) {
			double xd = xwd-1. + (i + 0.5)*dx;
			for (int j = 0; j < 2*nseg; ++j) {
				double x1 = -1. + j*dx;
				double x2 = x1 + dx;
				ans(i, j+1) = mult*_if1(x1, x2, yd);
			}
		}
	//}
	//{LOG_DURATION("iF2E");
		for (int i = 0; i < 2*nseg; ++i) {
			double xd = xwd-1. + (i + 0.5)*dx;
			for (int j = 0; j < 2*nseg; ++j) {
				double x1 = -1. + j*dx;
				double x2 = x1 + dx;
				ans(i, j+1) += mult*_if2E(x1, x2, xd, yd);
			}
		}
	////}
	//{LOG_DURATION("i1F2H");
		for (int i = 0; i < 2*nseg; ++i) {
				double xd = xwd-1. + (i + 0.5)*dx;
				for (int j = 0; j < 2*nseg; ++j) {
					double x1 = -1. + j*dx;
					double x2 = x1 + dx;
					ans(i, j+1) += mult*_i1f2H(x1, x2, xd, yd);
				}
			}
	//}
	//{LOG_DURATION("i2F2H");
		for (int i = 0; i < 2*nseg; ++i) {
				double xd = xwd-1. + (i + 0.5)*dx;
				for (int j = 0; j < 2*nseg; ++j) {
					double x1 = -1. + j*dx;
					double x2 = x1 + dx;
					ans(i, j+1) += mult*_i2f2H(x1, x2, yd);
				}
			}
	//}
	return ans;
}

VectorXd MakeRhs(const double u, const int nseg, const double Fcd) {
	VectorXd ans(2*nseg+1);
	double coef = LaplFunc::PI/Fcd/nseg/u;
	for (int i = 0; i < nseg; ++i) {
		ans(nseg+i) = coef*(i+0.5);
		ans(nseg-i-1) = ans[nseg+i];
	}
	ans(2*nseg) = 2.*nseg/u;
	return ans;
}




