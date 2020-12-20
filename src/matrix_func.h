/*
 * matrix_func.h
 *
 *  Created on: 19 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once
#include "lapl_functions.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include "profile.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

MatrixXd MakeMatrix(const double u, const int nseg, const double eps,
		const double xwd, const double xed, const double ywd, const double yed);
VectorXd MakeRhs(const double u, const int nseg, const double Fcd);
