#include "SASPSO.hpp"
#include <memory>

// test function
// f(x) = x1^2 + x2^2
double f(const RealVector<2>& x) {
	return x[0] * x[0] + x[1] * x[1];
}

int main() {
	Problem<2> problem(f, RealVector<2>({ -1.0, -1.0 }), RealVector<2>({ 1.0, 1.0 }));
	std::unique_ptr<Optimizer<2>> opt = std::make_unique<SASPSO<2>>(problem, 100, 100, 1e-6, 1.0, 1.0, 1.0);
	opt->initialize();
	opt->optimize();
}