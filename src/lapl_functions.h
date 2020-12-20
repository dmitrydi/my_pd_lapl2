/*
 * lapl_functions.h
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include "levin.h"
#include "bessel.h"
#include "quadrature.h"
#include <cmath>
#include <cstdint>
#include <limits>


namespace LaplFunc {

static const double PI = 3.141592653589793;

static const double TINY = std::numeric_limits<double>::min();

class ek{
public:
	ek(const double u_, const double ksi_, const double alpha_);
	double operator()(const int64_t k) const;
private:
	const double u, ksi, alpha;
};

class sexp {
public:
	sexp(const double yed_, const double eps_, const LevinSpace::SumAlgo algo_ = LevinSpace::SumAlgo::d);
	double operator()(const double e) const;
	double get_k(const double e, const int64_t k) const;
	bool Uniform() const;
private:
	const double yed;
	const double eps;
	LevinSpace::SumAlgo algo;
};

class iF1{
public:
	iF1(
			const double u_,
			const double ywd_,
			const double yed_,
			const double eps_
	);
	double operator()(const double x1, const double x2, const double yd) const;
private:
	const double u, ywd, yed, squ, eps, _sexp;
};

class iF2E{
public:
	iF2E(
			const double u_,
			const double ksiwd_,
			const double ksied_,
			const double ksiede_,
			const double alpha_,
			const double ywd_,
			const double yed_,
			const double eps_
		);
	double get_k(const double x1, const double x2, const double ksid, const double yd, const int64_t k) const;
	double operator()(const double x1, const double x2, const double xd, const double yd) const;
	bool Uniform() const;
private:
	const double u;
	const double ksiwd;
	const double ksied;
	const double ksiede;
	const double alpha;
	const double ywd;
	const double yed;
	const double eps;
	ek _ek;
	sexp _sexpp;
};

class i1F2H {
public:
	i1F2H(
			const double u_,
			const double ksiwd_,
			const double ksied_,
			const double ksiede_,
			const double alpha_,
			const double ywd_,
			const double yed_,
			const double eps_
		);
	double get_k(const double x1, const double x2, const double ksid, const double yd, const double beta, const int64_t k) const;
	double operator()(const double x1, const double x2, const double xd, const double yd) const;
private:
	const double u;
	const double ksiwd;
	const double ksied;
	const double ksiede;
	const double alpha;
	const double ywd;
	const double yed;
	const double eps;
	const double squ;
	const double _beta_mult;
};

class i2F2H {
public:
	i2F2H(
			const double u_,
			const double ywd_,
			const double alpha_
			);
	double operator()(const double x1, const double x2, const double yd) const;
private:
	const double u;
	const double ywd;
	const double alpha;
	const double squ;
};

}
