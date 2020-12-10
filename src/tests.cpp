/*
 * lapl_tests.cpp
 *
 *  Created on: 9 дек. 2020 г.
 *      Author: Dmitry_Di
 */

#include "tests.h"

namespace TestWijn {
using namespace std;
void Correct() {
	//double ek = 24.654;
	//double yed = 6.95193;
	struct Func {
		Func(double yed_, double ek_): uniform_conv(true), yed(yed_), ek(ek_) {};
		double operator()(int64_t m) const {
			if (m < 0) throw runtime_error("In TestWijn::Correct: negative index encountered");
			return exp(-2*m*yed*ek);
		}
		bool Uniform() const {return uniform_conv;};
		bool uniform_conv;
		double yed, ek;
	};

	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen_yed(rd()); //Standard mersenne_twister_engine seeded with rd()
	mt19937 gen_ek(rd());
	uniform_real_distribution<> dis_ek(0.1, 1.0);
	uniform_real_distribution<> dis_yed(1.0, 10.0);
	{ LOG_DURATION("TestWijn::Correct With uniform optimization")
		for (int i = 0; i < 1000; ++i) {
			double yed = dis_yed(gen_yed);
			double ek = dis_ek(gen_ek);
			Func fun(yed, ek);
			for (int64_t r = 1; r <= 1000; ++r) {
				double res = Wijn::SumNeum(fun, r);
				ASSERT(!isnan(res));
				ASSERT(!isinf(res));
			}
		}
	}
	{ LOG_DURATION("TestWijn::Correct Without uniform optimization")
		for (int i = 0; i < 1000; ++i) {
			double yed = dis_yed(gen_yed);
			double ek = dis_ek(gen_ek);
			Func fun(yed, ek);
			for (int64_t r = 1; r <= 1000; ++r) {
				double res = Wijn::SumNeum(fun, r);
				ASSERT(!isnan(res));
				ASSERT(!isinf(res));
			}
		}
	}
}
}

namespace TestSexp {
using namespace std;

double Guaranteed(const double yed, const double ek) {
// sums exp(-2*m*yed*ek) to almost numeric limit precision
	const double tiny = numeric_limits<double>::min();
	const double ln_tiny = log(tiny);
	const double dum = -0.5*ln_tiny/yed/ek;
	const int64_t max_member = static_cast<int64_t>(dum);
	Compensated::Neumaier<double> summator;
	for (int64_t i = max_member; i >= 1; --i) {
		double val = exp(-2*i*yed*ek);
		if (!isnan(val) && !isinf(val)) summator.Add(val);
	}
	return summator.Get();
}

void Correct() {
	const double TEST_SEXP_EPS = 1e-14;
	const double NTESTS = 1000;
	const int N = 20;
	const int KMAX = 10;
	const double xed_min = 1;
	const double xed_max  = 10;
	const double yed_min = 1;
	const double yed_max = 10;
	const double u_min = 1e-9;
	const double u_max = 1000;

	const auto xeds = LogSpaced(xed_min, xed_max, N);
	const auto yeds = LogSpaced(yed_min, yed_max, N);
	const auto us = LogSpaced(u_min, u_max, N);

	mt19937 rng;
	rng.seed(time(nullptr));
	uniform_int_distribution<std::mt19937::result_type> dist(0, N-1);
	uniform_int_distribution<std::mt19937::result_type> k_dist(1, KMAX);

	for (int i = 0; i < NTESTS; ++i) {
		auto x_ind = dist(rng);
		auto y_ind = dist(rng);
		auto u_ind = dist(rng);
		double u = us[u_ind];
		double xed = xeds[x_ind];
		double yed = yeds[y_ind];
		auto ek = LaplFunc::ek(u, xed, 0);
		auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS);
		int64_t k = k_dist(rng);
		double ek_ = ek(k);
		auto s = sexp(ek_);
		cerr << "ek = " << ek_ << " yed = " << yed << endl;
		ASSERT_CLOSE(s, Guaranteed(yed, ek_), TEST_SEXP_EPS);
	}
}

void ExhaustiveCorrect() {
	const double TEST_SEXP_EPS = 1e-14; // Actual minimum eps reachable
	const int N = 100;
	const double ekmin = 3.e-5;
	const double ekmax = 1000.;
	const double yedmin = 1.;
	const double yedmax = 1000.;
	const auto eks = LogSpaced(ekmin, ekmax, N);
	const auto yeds = LogSpaced(yedmin, yedmax, N);
	for (int i = LevinSpace::SumAlgo::t; i != LevinSpace::SumAlgo::last_; ++i) {
		LevinSpace::SumAlgo algo = static_cast<LevinSpace::SumAlgo>(i);
		switch (algo) {
		case LevinSpace::SumAlgo::t:
			cerr << "Test t-sum " << endl;
			break;
		case LevinSpace::SumAlgo::d:
			cerr << "Test d-sum " << endl;
			break;
		case LevinSpace::SumAlgo::u:
			cerr << "Test u-sum " << endl;
			break;
		}
		for (const auto yed: yeds) {
			for (const auto ek: eks) {
				cerr << "yed = " << yed << " ek = " << ek;
				auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS/2., algo);
				auto s = sexp(ek);
				cerr << " s = " << s;
				auto expected = Guaranteed(yed, ek);
				cerr << " eps = " << abs(s - expected)/expected << endl;
				ASSERT_CLOSE(s, expected, TEST_SEXP_EPS);
			}
		}
	}

}

void ExhaustiveSpeed() {
	const double TEST_SEXP_EPS = 1e-14; // Actual minimum eps reachable
	const int N = 300;
	const double ekmin = 3.e-5;
	const double ekmax = 1000.;
	const double yedmin = 1.;
	const double yedmax = 1000.;
	const auto eks = LogSpaced(ekmin, ekmax, N);
	const auto yeds = LogSpaced(yedmin, yedmax, N);
	double s;
	for (int i = LevinSpace::SumAlgo::t; i != LevinSpace::SumAlgo::last_; ++i) {
		LevinSpace::SumAlgo algo = static_cast<LevinSpace::SumAlgo>(i);
		ostringstream os;
		os << "Test sexp::sexp(), num evals " << N*N << ", Algo: ";
		switch (algo) {
		case LevinSpace::SumAlgo::t:
			os << "Test t-sum ";
			break;
		case LevinSpace::SumAlgo::d:
			os << "Test d-sum ";
			break;
		case LevinSpace::SumAlgo::u:
			os << "Test u-sum ";
			break;
		}
		{LOG_DURATION(os.str());
			for (const auto yed: yeds) {
				auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS/2.);
				for (const auto ek: eks) {
					s = sexp(ek);
				}
			}
			cerr << s;
		}
	}
}

void CheckSpurious() {
	const double TEST_SEXP_EPS = 1e-14;
	double yed = 6.15848;
	double ek = 31.9799;
	auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS);
	auto s = sexp(ek);
	cerr << "ek = " << ek << " yed = " << yed << endl;
	ASSERT_CLOSE(s, Guaranteed(yed, ek), TEST_SEXP_EPS);
}

void Speed() {
	const double TEST_SEXP_EPS = 1e-14;
	const int NTESTS = 1000'000;
	const int N = 20;
	const int KMAX = 1;
	const double xed_min = 1;
	const double xed_max  = 1000;
	const double yed_min = 1;
	const double yed_max = 1000;
	const double u_min = 1e-9;
	const double u_max = 1000;

	const auto xeds = LogSpaced(xed_min, xed_max, N);
	const auto yeds = LogSpaced(yed_min, yed_max, N);
	const auto us = LogSpaced(u_min, u_max, N);

	mt19937 rng;
	rng.seed(time(nullptr));
	uniform_int_distribution<std::mt19937::result_type> dist(0, N-1);
	uniform_int_distribution<std::mt19937::result_type> k_dist(1, KMAX);

	double s;

	ostringstream os;

	os << "Test_sexp_speed, Num tests: ";
	os << fixed << to_string(NTESTS);
	os << ", duration: ";

	{LOG_DURATION(os.str());
		for (int i = 0; i < NTESTS; ++i) {
			auto x_ind = dist(rng);
			auto y_ind = dist(rng);
			auto u_ind = dist(rng);
			double u = us[u_ind];
			double xed = xeds[x_ind];
			double yed = yeds[y_ind];
			auto ek = LaplFunc::ek(u, xed, 0);
			auto sexp = LaplFunc::sexp(yed, TEST_SEXP_EPS/100.);
			int64_t k = k_dist(rng);
			s = sexp(ek(k));
		}
	}
	cerr << s << endl;
}

}

namespace TestiF2E {
using namespace std;

iF2EG_ret iF2EGuaranteed(
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
		const double yd) {
	auto if2e = LaplFunc::iF2E(u, ksiwd, ksied, ksiede, alpha, ywd, yed, eps);
	double psum;
	double sum = 0;
	int64_t i = 0;
	const int64_t blk_size = 100000;
	do {
		psum = 0;
		++i;
		Compensated::Neumaier<double> summator;
		for (int64_t j = i*blk_size; j > (i-1)*blk_size; --j) {
			summator.Add(if2e.get_k(x1, x2, ksid, yd, j));
		}
		psum = summator.Get();
		sum += psum;
	} while (std::abs(psum) > numeric_limits<double>::min()*10);
	return {sum, psum, i*blk_size};
}

void Correct() {
	const double u = 1.e-9;
	const double ksiwd = 500;
	const double ksied = 1000;
	const double ksiede = 1000;
	const double alpha = 0;
	const double ywd = 500;
	const double yed = 1000;
	const double eps = 1e-14;
	double x1 = -1.;
	double x2 = x1 + 1./20.;
	const double ksid = ksiwd;
	const double yd = ywd;
	double ans = LaplFunc::iF2E(u, ksiwd, ksied, ksiede, alpha, ywd, yed, eps/100.)(x1, x2, ksid, yd);
	auto expected = iF2EGuaranteed(
			u,
			ksiwd,
			ksied,
			ksiede,
			alpha,
			ywd,
			yed,
			eps/100.,
			x1,
			x2,
			ksid,
			yd);
	cerr <<" k = " << expected.k << " ans = " << expected.res << endl;
	ASSERT_CLOSE(ans, expected.res, eps);
}
}







