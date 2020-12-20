/*
 * naive_algos.cpp
 *
 *  Created on: 15 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "naive_algos.h"

namespace Compensated {

NeumaierQuad::NeumaierQuad(): sum(0.0q), c(0.0q) {};
void NeumaierQuad::Add(const __float128 input) {
	__float128 t = sum + input;
	if (fabsq(sum) >= fabsq(input)) {
		c += (sum - t) + input;
		} else {
			c += (input - t) + sum;
		}
		sum = t;
};
__float128 NeumaierQuad::Get() const {
	return sum + c;
};

}
