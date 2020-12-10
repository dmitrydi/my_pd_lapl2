/*
 * lapl_tests.h
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */
#pragma once

#include "lapl_functions.h"
#include "test_runner.h"
#include "auxillary.h"
#include "profile.h"
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <random>
#include <ctime>
#include <utility>

namespace TestLevin {
void Stability();
}

namespace TestWijn {
void Correct();
void Speed();
}

namespace TestSexp {
double Guaranteed(const double yed, const double ek);
void Correct();
void Speed();
void CheckSpurious();
void ExhaustiveCorrect();
void ExhaustiveSpeed();
}

namespace TestiF2E {
struct iF2EG_ret {
	double res;
	double eps;
	int64_t k;
};

iF2EG_ret Guaranteed(
		const double u,
		const double ksiwd,
		const double ksied,
		const double ksiede,
		const double alpha,
		const double ywd,
		const double yed,
		const double eps,
		double x1,
		double x2,
		const double ksid,
		const double yd);

void Correct();
}
