#pragma once

#include <vector>
#include "../TrivariateFunction.h"
#include "../../util/MathUtils.h"

/**
 * 3D-spline function.
 *
 */
class Tricubic_Function : Trivariate_Function
{
private:
	/** Number of points. */
	static constexpr short N{ 4 };
	/** Coefficients */
	static std::vector<std::vector<std::vector<double>>> my_a;

public:
	/**
	 * @param aV List of spline coefficients.
	 */
	Tricubic_Function(std::vector<double>& aV)
	{
		my_a = std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>(N, std::vector<double>(N)));
		for (int i{}; i < N; i++)
		{
			for (int j{}; j < N; j++)
			{
				for (int k{}; k < N; k++)
				{
					my_a[i][j][k] = aV[i + N * (j + N * k)];
				}
			}
		}
	}

	/**
	 * @param x x-coordinate of the interpolation point.
	 * @param y y-coordinate of the interpolation point.
	 * @param z z-coordinate of the interpolation point.
	 * @return the interpolated value.
	 * @ if {@code x}, {@code y} or
	 * {@code z} are not in the interval {@code [0, 1]}.
	 */
	 //override
	double value(const double& x, const double& y, const double& z)
	{
		//Math_Utils::check_range_inclusive(x, 0, 1);
		//Math_Utils::check_range_inclusive(y, 0, 1);
		//Math_Utils::check_range_inclusive(z, 0, 1);

		const double x2 = x * x;
		const double x3 = x2 * x;
		const std::vector<double> pX = { 1, x, x2, x3 };

		const double y2 = y * y;
		const double y3 = y2 * y;
		const std::vector<double> pY = { 1, y, y2, y3 };

		const double z2 = z * z;
		const double z3 = z2 * z;
		const std::vector<double> pZ = { 1, z, z2, z3 };

		double result{};
		for (int i{}; i < N; i++)
		{
			for (int j{}; j < N; j++)
			{
				for (int k{}; k < N; k++)
				{
					result += my_a[i][j][k] * pX[i] * pY[j] * pZ[k];
				}
			}
		}

		return result;
	}
};