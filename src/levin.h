/*
 * levin.h
 *
 *  Created on: 9 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once
#include "wijn.h"
#include <vector>
#include <cmath>
#include <limits>
#include <cstdint>
#include <sstream>
#include <iostream>
#include <sstream>

namespace LevinSpace {
enum SumAlgo {
	t,
	d,
	u,
	last_
};

static const double SMALL_CONST = std::numeric_limits<double>::min()*10.0;
static const double BIG_CONST = std::numeric_limits<double>::max();
static const double EPS = std::numeric_limits<double>::epsilon();
static const double TINY = std::numeric_limits<double>::min();
static const int LEVIN_MAXITER = 1000;

struct Levin {
	//Convergence acceleration of a sequence by the Levin transformation. Initialize by calling the
	//constructor with arguments nmax, an upper bound on the number of terms to be summed, and
	//epss, the desired accuracy. Then make successive calls to the function next, which returns
	//the current estimate of the limit of the sequence. The flag cnvgd is set when convergence is
	//detected.
	std::vector<double> numer,denom; //Numerator and denominator computed via (5.3.16).
	int n,ncv;
	bool cnvgd;
	double small,big; //Numbers near machine underflow and overflow limits.
	double eps,lastval,lasteps;
	Levin(int nmax, double epss) : numer(nmax), denom(nmax), n(0), ncv(0),
			cnvgd(0), eps(epss), lastval(0.) {
		small=SMALL_CONST;
		big=BIG_CONST;
	}
	double next(double sum, double omega, double beta=1.) {
	//Arguments: sum, the nth partial sum of the sequence; omega, the nth remainder estimate
	//!n, usually from (5.3.19); and the parameter beta, which should usually be set to 1, but
	//sometimes 0.5 works better. The current estimate of the limit of the sequence is returned.
		double val;
		if (std::abs(omega) < TINY) {
			val = sum;
		} else {
			int j;
			double fact,ratio,term;
			term=1.0/(beta+n);
			denom[n]=term/omega;
			numer[n]=sum*denom[n];

			if (n > 0) {
				ratio=(beta+n-1)*term;
				for (j=1;j<=n;j++) {
					fact=(n-j+beta)*term;
					numer[n-j]=numer[n-j+1]-fact*numer[n-j];
					denom[n-j]=denom[n-j+1]-fact*denom[n-j];
					term=term*ratio;
				}
			}
			n++;
			if (std::abs(denom[0]) < small) {
				val = lastval;
			} else {
				val = numer[0]/denom[0];
			}
		}
		lasteps = std::abs(val-lastval);
		if (lasteps <= eps*std::abs(val) || lasteps <= TINY) ncv++;
		if (ncv >= 2) cnvgd = 1;
		lastval = val;
		return lastval;
	}
};

template <typename F>
double LevinTSum(const F& func,const bool from_zero_i = false, const double eps = EPS) {
	Levin lev(LEVIN_MAXITER, eps);
	double sum = from_zero_i ? func(0): 0.;
	double omega;
	double ans;
	const double beta = 1.;
	int sign = 1;
	for (int64_t r = 1; r < LEVIN_MAXITER; ++r) {
		omega = static_cast<double>(sign) * Wijn::SumNeum(func, r);
		sum += omega;
		ans = lev.next(sum, omega, beta);
		if (lev.cnvgd) {
			return ans;
		}

		sign *= -1;
	}
	std::ostringstream os;
	os << "Levin::LevinTSum did not converge in " << LEVIN_MAXITER << " steps, " << " eps = " << std::scientific << std::abs(lev.lasteps/lev.lastval);
	os << " lasteps = " << std::scientific << std::abs(lev.lasteps) << " lastval = " << std::scientific << std::abs(lev.lastval);
	throw std::runtime_error(os.str());
}

template <typename F>
double LevinUSum(const F& func, const bool from_zero_i = false, const double eps = EPS) {
	Levin lev(LEVIN_MAXITER, eps);
	double sum = from_zero_i ? func(0): 0.;
	double omega;
	double ans;
	const double beta = 1.;
	int sign = 1;
	for (int64_t r = 1; r < LEVIN_MAXITER; ++r) {
		double mem = sign * Wijn::SumNeum(func, r);
		sum += mem;
		omega = (beta + r)*mem; //(beta + n)*an
		ans = lev.next(sum, omega, beta);
		if (lev.cnvgd) {
			return ans;
		}

		sign *= -1;
	}
	std::ostringstream os;
	os << "Levin::LevinUSum did not converge in " << LEVIN_MAXITER << " steps, " << " eps = " << std::scientific << lev.eps;
	throw std::runtime_error(os.str());
}

template <typename F>
double LevinDSum(const F& func, const bool from_zero_i = false, const double eps = EPS) {
	Levin lev(LEVIN_MAXITER, eps);
	double sum = from_zero_i ? func(0): 0.;
	double ans;
	const double beta = 1.;
	int sign = 1;
	double omega = sign*Wijn::SumNeum(func, 1);
	double d;
	for (int64_t r = 1; r < LEVIN_MAXITER; ++r) {
		d = omega;
		sum += d;
		sign *= -1;
		omega = sign * Wijn::SumNeum(func, r + 1);
		ans = lev.next(sum, omega, beta);
		if (lev.cnvgd) {
			return ans;
		}

	}
	std::ostringstream os;
	os << "Levin::LevinDSum did not converge in " << LEVIN_MAXITER << " steps, " << " eps = " << std::scientific << std::abs(lev.lasteps/lev.lastval);
	throw std::runtime_error(os.str());
}

template <typename F>
double Sum(const SumAlgo& algo, const F& func,const bool from_zero_i = false, const double eps = EPS) {
	switch (algo) {
	case SumAlgo::t:
		return LevinTSum(func, from_zero_i, eps);
	case SumAlgo::d:
		return LevinDSum(func, from_zero_i, eps);
	case SumAlgo::u:
		return LevinUSum(func, from_zero_i, eps);
	default:
		throw;
	}
}

}
