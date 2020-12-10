//============================================================================
// Name        : my_pd_lapl2.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include "test_runner.h"
#include "profile.h"
#include "auxillary.h"
#include "tests.h"

using namespace std;

int main() {

	TestRunner tr;
	//RUN_TEST(tr, TestLevin::Stability); //passed
	//RUN_TEST(tr, TestSexp::ExhaustiveCorrect); //passed
	//RUN_TEST(tr, TestSexp::ExhaustiveSpeed); //passed
	//RUN_TEST(tr, TestSexp::CheckSpurious);//passed
	//RUN_TEST(tr, TestSexp::Correct); //passed
	//RUN_TEST(tr, TestWijn::Correct); //passed
	//RUN_TEST(tr, TestiF2E::Correct); //passed

	return 0;
}
