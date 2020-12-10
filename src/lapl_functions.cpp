/*
 * lapl_functions.cpp
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "lapl_functions.h"

using namespace std;

static const double PI = 3.141592653589793;

namespace LaplFunc {

ek::ek(const double u_, const double ksi_, const double alpha_): u(u_), ksi(ksi_), alpha(alpha_){};

double ek::operator()(const int64_t k) const {
	double dum = static_cast<double>(k)*PI/ksi;
	return sqrt(u + dum*dum + alpha*alpha);
};

sexp::sexp(const double yed_, const double eps_, const LevinSpace::SumAlgo algo_): yed(yed_), eps(eps_), algo(algo_) {};

double sexp::operator()(const double e) const {
	auto fun = [e, this](int64_t m) { return this->get_k(e, m); };
	return LevinSpace::Sum(algo, fun, false, eps);
}

double sexp::get_k(const double e, const int64_t k) const {
	return exp(-2.*static_cast<double>(k)*yed*e);
}

bool sexp::Uniform() const {return true;};

iF2E::iF2E(
		const double u_,
		const double ksiwd_,
		const double ksied_,
		const double ksiede_,
		const double alpha_,
		const double ywd_,
		const double yed_,
		const double eps_
	): u(u_), ksiwd(ksiwd_), ksied(ksied_), ksiede(ksiede_), alpha(alpha_), ywd(ywd_), yed(yed_), eps(eps_) {};
double iF2E::get_k(const double x1, const double x2, const double ksid, const double yd, const int64_t k) const {
	double ek_ = ek(u, ksiede, alpha)(k);
	double sexp_ = sexp(yed, eps/100.)(ek_);
	double aydywd = std::abs(yd-ywd);
	double kpiOxed = k*PI/ksied;
	double ydPywd = yd + ywd;
	double res;
	res = 2./kpiOxed/ek_;
	res *= std::cos(kpiOxed*ksid);
	res *= std::sin(0.5*kpiOxed*(x2 - x1));
	res *= std::cos(0.5*kpiOxed*(2.*ksiwd + x1 + x2));
	res *= (std::exp(-ek_*(2.*yed - ydPywd)) + std::exp(-ek_*ydPywd) + std::exp(-ek_*(2.*yed-aydywd)))*(1 + sexp_)
			+ std::exp(-ek_*aydywd)*sexp_;
	return res;
}
double iF2E::operator()(double x1, double x2, const double ksid, const double yd) const {
	auto fun = [x1, x2, ksid, yd, this](int64_t m) { return this->get_k(x1, x2, ksid, yd, m); };
	return LevinSpace::LevinDSum(fun, false, eps);
}

bool iF2E::Uniform() const {return false;};

}



