/*
 * lapl_functions.cpp
 *
 *  Created on: 3 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "lapl_functions.h"

using namespace std;

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
	): u(u_), ksiwd(ksiwd_), ksied(ksied_), ksiede(ksiede_), alpha(alpha_), ywd(ywd_), yed(yed_), eps(eps_), _ek(u, ksiede, alpha),
			_sexpp(yed, eps/100.){};
double iF2E::get_k(const double x1, const double x2, const double ksid, const double yd, const int64_t k) const {
	double ek_ = _ek(k);//ek(u, ksiede, alpha)(k);
	double sexp_ = _sexpp(ek_);//sexp(yed, eps/100.)(ek_);
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
double iF2E::operator()(double x1, double x2, const double xd, const double yd) const {
	auto fun = [x1, x2, xd, yd, this](int64_t m) { return this->get_k(x1, x2, xd, yd, m); };
	double ans =0.;
	double term =0.;
	int64_t k = 0;
	do {
		term = fun(++k);
		ans += term;
	} while (10.*std::abs(term) > std::abs(ans)*eps);
	return ans;
}

bool iF2E::Uniform() const {return false;};

i1F2H::i1F2H(
		const double u_,
		const double ksiwd_,
		const double ksied_,
		const double ksiede_,
		const double alpha_,
		const double ywd_,
		const double yed_,
		const double eps_
		): u(u_), ksiwd(ksiwd_), ksied(ksied_), ksiede(ksiede_), alpha(alpha_), ywd(ywd_), yed(yed_), eps(eps_),
				squ(std::sqrt(u+alpha*alpha)), _beta_mult(ksied/ksiede/squ){
};
double i1F2H::get_k(const double x1, const double x2, const double ksid, const double yd, const double beta, const int64_t k) const{
	double t = ksiede*(ksid/ksied+beta*ksiwd/ksied-2*k);
	double t1 = t - ksiede*x2/ksied;
	double t2 = t - ksiede*x1/ksied;
	double adyd = std::abs(yd-ywd);
	if (adyd <= TINY) {
		double _mult = _beta_mult;
		if (t1 > 0.0) {
			return _mult*static_cast<double>(ibessk0q(squ*t2) - ibessk0q(squ*t1));
		} else if (t2 < 0) {
			return _mult*(ibessk0q(squ*std::abs(t1)) - ibessk0q(squ*std::abs(t2)));
		} else if (t1 <=0. && t2 >=0.) {
			return _mult*(ibessk0q(squ*std::abs(t1)) + ibessk0q(squ*std::abs(t2)));
		}
	} else {
		Bessik bess;
		auto func = [&bess, adyd, this](double x) {return bess.k0((this->squ)*std::sqrt(x*x+adyd*adyd));};
		return ksid/ksiede*qromb(func, t1, t2, eps);
	}
};
double i1F2H::operator()(const double x1, const double x2, const double xd, const double yd) const {
	auto fun1 = [x1, x2, xd, yd, this](int64_t k) { return this->get_k(x1, x2, xd, yd, 1, k); };
	auto fun2 = [x1, x2, xd, yd, this](int64_t k) { return this->get_k(x1, x2, xd, yd, 1, -k); };
	auto fun3 = [x1, x2, xd, yd, this](int64_t k) { return this->get_k(x1, x2, xd, yd, -1, k); };
	auto fun4 = [x1, x2, xd, yd, this](int64_t k) { return this->get_k(x1, x2, xd, yd, -1, -k); };
	double t1 = LevinSpace::LevinDSum(fun1, true, eps);
	double t2 = LevinSpace::LevinDSum(fun2, false, eps);
	double t3 = LevinSpace::LevinDSum(fun3, true, eps);
	double t4 = LevinSpace::LevinDSum(fun4, false, eps);
	return 0.5*ksiede/PI*(t1 + t2 + t3 + t4);
};

iF1::iF1(const double u_, const double ywd_, const double yed_, const double eps_): u(u_), ywd(ywd_), yed(yed_), squ(std::sqrt(u)),
		eps(eps_), _sexp(sexp(yed, eps/10.)(squ)) {};
double iF1::operator()(const double x1, const double x2, const double yd) const {
	using std::exp;
	double ans = 0.5*(x2-x1)/squ;
	double dy = std::abs(yd-ywd);
	double sumy = yd+ywd;
	ans *= exp(-squ*(2.*yed-sumy))+exp(-squ*sumy)+exp(-squ*(2.*yed-dy))+exp(-squ*dy);
	ans *= (1+_sexp);
	return ans;
};

i2F2H::i2F2H(const double u_, const double ywd_, const double alpha_): u(u_), ywd(ywd_), alpha(alpha_), squ(std::sqrt(u+alpha*alpha)) {};
double i2F2H::operator()(const double x1, const double x2, const double yd) const {
	return -0.5*std::exp(-squ*std::abs(yd-ywd))/squ*(x2-x1);
}
}



