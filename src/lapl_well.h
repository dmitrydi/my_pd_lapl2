/*
 * lapl_well.h
 *
 *  Created on: 20 дек. 2020 г.
 *      Author: Dmitry_Di
 */
#pragma once
#include "matrix_func.h"

class LaplWell {
public:
	LaplWell(const double xwd_, const double xed_, const double ywd_, const double yed_, const double Fcd_);
	double pd(const double u) const;
private:
	const double xwd, xed, ywd, yed, Fcd;
	const int nseg = 20;
	const double eps = 1.e-12;
};
