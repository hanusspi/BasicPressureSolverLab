#include "kernel.h"
#include <iostream>

double pressureSolver::kernel::a = 0.0;
double pressureSolver::kernel::h = 0.0;

double pressureSolver::kernel::a_coh = 0.0;
double pressureSolver::kernel::c = 0.0;
double pressureSolver::kernel::sub_coh = 0.0;
double pressureSolver::kernel::a_adh = 0.0;

void pressureSolver::kernel::setAlpha(const double h_in)
{
	h = h_in;
	a = 3.0 / (2.0 * PI * h * h * h);
}

constexpr double pressureSolver::kernel::cubic_spline(const double q)
{
	assert(q >= 0.0);

	if (q < 1.0) {
		return ((2.0 / 3.0) - (q * q) + 0.5 * (q * q * q));
	}
	else if (q < 2.0) {
		return ((1.0 / 6.0) * (2.0 - q) * (2.0 - q) * (2.0 - q));
	}
	else {
		return 0.0;
	}
}

constexpr double pressureSolver::kernel::cubic_grad_spline(const double q)
{
	assert(q >= 0.0);

	if (q < 1.0) {
		return -2.0 * q + 1.5 * q * q;
	}
	else if (q < 2.0) {
		const double tmp = 2.0 - q;
		return -0.5 * tmp * tmp;
	}
	else {
		return 0.0;
	}
}

double pressureSolver::kernel::W(const Eigen::Vector3d& x)
{
	const double q = (x).norm() / h;
	if (std::isnan(q))
		return 0;
	return a * cubic_spline(q);
}

Eigen::Vector3d pressureSolver::kernel::W_grad(const Eigen::Vector3d& x)
{
	double r = (x).norm() + epsilon;
	const double q = r / h;
	const double gradientMag = a * cubic_grad_spline(q);
	auto tmo = gradientMag * ((x) / (r * h));
	if (tmo.hasNaN())
		return Eigen::Vector3d(0,0,0);
	return gradientMag * ((x) / (r * h));
}

double pressureSolver::kernel::W_cohesion(const Eigen::Vector3d& x)
{
	double r = (x).norm();
	if (r < c / 2.0 && 0 <= r)
		return a_coh * (2 * (c - r) * (c - r) * (c - r) * r * r * r - sub_coh);
	else if (r < c && c / 2.0 <= r)
		return a_coh * (c - r) * (c - r) * (c - r) * r * r * r;
	else
		return 0.0;
}

double pressureSolver::kernel::W_adhesion(const Eigen::Vector3d& x)
{
	double r = (x).norm();
	if (c / 2.0 <= r && r <= c)
		return a_adh * pow((-4 * r * r / c) + 6 * r - 2 * c, 0.25);
	else
		return 0.0;
}