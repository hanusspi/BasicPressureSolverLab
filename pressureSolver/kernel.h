#pragma once
#include <cassert>
#include <Eigen/Dense>

#define epsilon 1e-5

namespace pressureSolver
{
	namespace kernel
	{
		constexpr double PI = 3.14159265358979323846;

		extern double a;
		extern double h;

		extern double c;
		extern double a_coh;
		extern double a_adh;
		extern double sub_coh; // for subtraction term in cohesion kernel

		void setAlpha(const double h);

		constexpr double cubic_spline(const double q);
		constexpr double cubic_grad_spline(const double q);

		double W(const Eigen::Vector3d& x);
		Eigen::Vector3d W_grad(const Eigen::Vector3d& x);

		double W_cohesion(const Eigen::Vector3d& x);
		double W_adhesion(const Eigen::Vector3d& x);

		//Eigen::Vector3d W_grad_num(const Eigen::Vector3d& x);
	};
};