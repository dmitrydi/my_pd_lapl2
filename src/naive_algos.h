/*
 * naive_algos.h
 *
 *  Created on: 8 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include <vector>
#include <stack>
#include <numeric>
#include <iterator>
#include <iostream>
#include <quadmath.h>


namespace Compensated {

template <typename T>
T KahanSum(const std::vector<T>& input) {
	T sum = 0.;
	T c = 0.;
	for (auto v: input) {
		T y = v - c;
		T t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	return sum;
}

template <typename T>
T NeumaierSum(const std::vector<T>& input) {
	T sum = 0.;
	T c = 0.;

	for (auto v: input) {
		T t = sum + v;
		if (std::abs(sum) >= std::abs(v)) {
			c += (sum - t) + v;
		} else {
			c += (v - t) + sum;
		}
		sum = t;
	}
	return sum + c;
}

template <typename T>
class Kahan {
public:
	Kahan(): sum(0.), c(0.) {};
	void Add(const T input) {
		T y = input - c;
		T t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	T Get() const { //to be called once at the very end of summation process
		return sum;
	}
private:
	T sum, c;
};

template <typename T>
class Neumaier {
public:
	Neumaier(): sum(0.), c(0.) {};
	void Add(const T input) {
		T t = sum + input;
		if (std::abs(sum) >= std::abs(input)) {
			c += (sum - t) + input;
		} else {
			c += (input - t) + sum;
		}
		sum = t;
	}
	T Get() const { //to be called once at the very end of summation process
		return sum + c;
	}
private:
	T sum, c;
};


class NeumaierQuad {
public:
	NeumaierQuad();
	void Add(const __float128 input);
	__float128 Get() const;
private:
	__float128 sum, c;
};

}

namespace Pairwise {

template <typename T>
T SumStack(const std::vector<T>& sequence) {
	std::stack<T> st;
	int j;
	double p, q;
	for (size_t i = 0; i < sequence.size(); ++i) {
		st.push(sequence[i]);
		j = i + 1;
		while (j % 2 == 0) {
			j /= 2;
			p = st.top();
			st.pop();
			q = st.top();
			st.pop();
			st.push(p+q);
		}
	}
	T total = 0.;
	while (st.size() > 0) {
		total += st.top();
		st.pop();
	}
	return total;
}

}
