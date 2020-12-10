/*
 * wijn.h
 *
 *  Created on: 8 дек. 2020 г.
 *      Author: Dmitry_Di
 */
#pragma once

#include "naive_algos.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <iostream>

namespace Wijn {
static const int64_t MAXN = std::numeric_limits<int64_t>::max();
static const double SMALL = std::numeric_limits<double>::min();
static const double EPS = std::numeric_limits<double>::epsilon();

template <typename Func>
double SumVStack(const Func& f, int64_t r) {
	long double x = static_cast<long double>(MAXN)/static_cast<long double>(r);
	int nmax = static_cast<int>(floor(log2(x)));
	std::vector<double> vals;
	vals.reserve(nmax);
	int64_t mult = 1;
	for (int i = 0; i < nmax; ++i) {
		double v = mult*f(r*mult);
		if (std::isinf(v) || std::isnan(v)) break;
		vals.push_back(v);
		mult *= 2;
	}
	return Pairwise::SumStack(vals);
}

template <typename Func>
double SumVKahan(const Func& f, int64_t r) {
	long double x = static_cast<long double>(MAXN)/static_cast<long double>(r);
	int nmax = static_cast<int>(floor(log2(x)));
	std::vector<double> vals;
	vals.reserve(nmax);
	int64_t mult = 1;
	for (int i = 0; i < nmax; ++i) {
		double v = mult*f(r*mult);
		if (std::isinf(v) || std::isnan(v)) break;
		vals.push_back(v);
		mult *= 2;
	}
	return Compensated::KahanSum(vals);
}

template <typename Func>
double SumVNeum(const Func& f, int64_t r) {
	long double x = static_cast<long double>(MAXN)/static_cast<long double>(r);
	int nmax = static_cast<int>(floor(log2(x)));
	std::vector<double> vals;
	vals.reserve(nmax);
	int64_t mult = 1;
	for (int i = 0; i < nmax; ++i) {
		double v = mult*f(r*mult);
		if (std::isinf(v) || std::isnan(v)) break;
		vals.push_back(v);
		mult *= 2;
	}
	return Compensated::NeumaierSum(vals);
}

template <typename Func>
double SumKah(const Func& f, int64_t r) {
	long double x = static_cast<long double>(MAXN)/static_cast<long double>(r);
	int nmax = static_cast<int>(floor(log2(x)));
	int64_t mult = 1;
	Compensated::Kahan<double> Summator;
	for (int i = 0; i < nmax; ++i) {
		double v = mult*f(r*mult);
		if (std::isinf(v) || std::isnan(v)) break;
		Summator.Add(v);
		mult *= 2;
	}
	return Summator.Get();
}

template <typename Func>
double SumNeum(const Func& f, int64_t r) {
	long double x = static_cast<long double>(MAXN)/static_cast<long double>(r);
	int nmax = static_cast<int>(floor(log2(x)));
	int64_t mult = 1;
	Compensated::Neumaier<double> Summator;
	for (int i = 0; i < nmax; ++i) {
		double v = mult*f(r*mult);
		if (std::isinf(v) || std::isnan(v)) break;
		Summator.Add(v);
		mult *= 2;
	}
	return Summator.Get();
}

}
