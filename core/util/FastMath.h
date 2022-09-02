#pragma once
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

 /*
  * This is not the original file distributed by the Apache Software Foundation
  * It has been modified by the Hipparchus project
  */
  //package org.hipparchus.util;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Runtime_Exception;
#include <type_traits>
#include <vector>
#include "../CalculusFieldElement.hpp"

/**
 * Faster, more accurate, portable alternative to {@link Math} and
 * {@link StrictMath} for large scale computation.
 * <p>
 * FastMath is a drop-in replacement for both Math and StrictMath. This
 * means that for any method in Math (say {@code Math.sin(x)} or
 * {@code Math.cbrt(y)}), user can directly change the class and use the
 * methods as is (using {@code std::sin(x)} or {@code std::cbrt(y)}
 * in the previous example).
 * <p>
 * FastMath speed is achieved by relying heavily on optimizing compilers
 * to native code present in many JVMs today and use of large tables.
 * The larger tables are lazily initialized on first use, so that the setup
 * time does not penalize methods that don't need them.
 * <p>
 * Note that FastMath is
 * extensively used inside Hipparchus, so by calling some algorithms, * the overhead when the the tables need to be initialized will occur
 * regardless of the end-user calling FastMath methods directly or not.
 * Performance figures for a specific JVM and hardware can be evaluated by
 * running the FastMathTestPerformance tests in the test directory of the source
 * distribution.
 * <p>
 * FastMath accuracy should be mostly independent of the JVM as it relies only
 * on IEEE-754 basic operations and on embedded tables. Almost all operations
 * are accurate to about 0.5 ulp throughout the domain range. This statement, * of course is only a rough global observed behavior, it is <em>not</em> a
 * guarantee for <em>every</em> double numbers input (see William Kahan's <a
 * href="http://en.wikipedia.org/wiki/Rounding#The_table-maker.27s_dilemma">Table
 * Maker's Dilemma</a>).
 * <p>
 * FastMath additionally : the following methods not found in Math/StrictMath:
 * <ul>
 * <li>{@link #asinhstatic_cast<double>(}</li>
 * <li>{@link #acoshstatic_cast<double>(}</li>
 * <li>{@link #atanhstatic_cast<double>(}</li>
 * </ul>
 * The following methods are found in Math/StrictMath since 1.6 only, they are provided
 * by FastMath even in 1.5 Java virtual machines
 * <ul>
 * <li>{@link #copy_sign(double, double)}</li>
 * <li>{@link #get_exponentstatic_cast<double>(}</li>
 * <li>{@link #next_after(double,double)}</li>
 * <li>{@link #next_upstatic_cast<double>(}</li>
 * <li>{@link #scalb(double, int)}</li>
 * <li>{@link #copy_sign(float, float)}</li>
 * <li>{@link #get_exponent(float)}</li>
 * <li>{@link #next_after(float,double)}</li>
 * <li>{@link #next_up(float)}</li>
 * <li>{@link #scalb(float, int)}</li>
 * </ul>
 */
class FastMath
{
public:
	/** Archimede's constant PI, ratio of circle circumference to diameter. */
	static const double PI = 105414357.0 / 33554432.0 + 1.984187159361080883e-9;

	/** Napier's constant e, base of the natural logarithm. */
	static const double E = 2850325.0 / 1048576.0 + 8.254840070411028747e-8;

	/** Index of exp(0) in the array of integer exponentials. */
	static const int EXP_INT_TABLE_MAX_INDEX = 750;
	/** Length of the array of integer exponentials. */
	static const int EXP_INT_TABLE_LEN = EXP_INT_TABLE_MAX_INDEX * 2;
	/** Logarithm table length. */
	static const int LN_MANT_LEN = 1024;
	/** Exponential fractions table length. */
	static const int EXP_FRAC_TABLE_LEN = 1025; // 0, 1/1024, ... 1024/1024

private:
	/** StrictMath.log(Double.MAX_VALUE): {@value} */
	static const double LOG_MAX_VALUE = StrictMath.log(Double.MAX_VALUE);

	/** Indicator for tables initialization.
	 * <p>
	 * This compile-time constant should be set to true only if one explicitly
	 * wants to compute the tables at class loading time instead of using the
	 * already computed ones provided as literal arrays below.
	 * </p>
	 */
	static constexpr bool RECOMPUTE_TABLES_AT_RUNTIME{};

	/** log(2) (high bits). */
	static constexpr double LN_2_A{ 0.693147063255310059 };

	/** log(2) (low bits). */
	static constexpr double LN_2_B{ 1.17304635250823482e-7 };

	/** Coefficients for log, when input 0.99 < x < 1.01. */
	static const std::vector<std::vector<double>> LN_QUICK_COEF =
	{
		{1.0, 5.669184079525E-24},
		{-0.25, -0.25},
		{0.3333333134651184, 1.986821492305628E-8},
		{-0.25, -6.663542893624021E-14},
		{0.19999998807907104, 1.1921056801463227E-8},
		{-0.1666666567325592, -7.800414592973399E-9},
		{0.1428571343421936, 5.650007086920087E-9},
		{-0.12502530217170715, -7.44321345601866E-11},
		{0.11113807559013367, 9.219544613762692E-9}
	};

	/** Coefficients for log in the range of 1.0 < x < 1.0 + 2^-10. */
	static const std::vector<std::vector<double>> LN_HI_PREC_COEF =
	{
		{1.0, -6.032174644509064E-23}, {-0.25, -0.25}, {0.3333333134651184, 1.9868161777724352E-8}, {-0.2499999701976776, -2.957007209750105E-8}, {0.19999954104423523, 1.5830993332061267E-10}, {-0.16624879837036133, -2.6033824355191673E-8}
	};

	/** Sine, Cosine, Tangent tables are for 0, 1/8, 2/8, ... 13/8 = PI/2 approx. */

	/** Sine table (high bits). */
	static const std::vector<double> SINE_TABLE_A =
	{
		+0.0,
		+0.1246747374534607d,
		+0.24740394949913025d,
		+0.366272509098053d,
		+0.4794255495071411d,
		+0.5850973129272461d,
		+0.6816387176513672d,
		+0.7675435543060303d,
		+0.8414709568023682d,
		+0.902267575263977d,
		+0.9489846229553223d,
		+0.9808930158615112d,
		+0.9974949359893799d,
		+0.9985313415527344d
	};

	/** Sine table (low bits). */
	static const std::vector<double> SINE_TABLE_B =
	{
		+0.0,
		-4.068233003401932E-9d,
		+9.755392680573412E-9d,
		+1.9987994582857286E-8d,
		-1.0902938113007961E-8d,
		-3.9986783938944604E-8d,
		+4.23719669792332E-8d,
		-5.207000323380292E-8d,
		+2.800552834259E-8d,
		+1.883511811213715E-8d,
		-3.5997360512765566E-9d,
		+4.116164446561962E-8d,
		+5.0614674548127384E-8d,
		-1.0129027912496858E-9d
	};

	/** Cosine table (high bits). */
	static const std::vector<double> COSINE_TABLE_A =
	{
		+1.0,
		+0.9921976327896118d,
		+0.9689123630523682d,
		+0.9305076599121094d,
		+0.8775825500488281d,
		+0.8109631538391113d,
		+0.7316888570785522d,
		+0.6409968137741089d,
		+0.5403022766113281d,
		+0.4311765432357788d,
		+0.3153223395347595d,
		+0.19454771280288696d,
		+0.07073719799518585d,
		-0.05417713522911072d
	};

	/** Cosine table (low bits). */
	static const std::vector<double> COSINE_TABLE_B =
	{
		+0.0,
		+3.4439717236742845E-8d,
		+5.865827662008209E-8d,
		-3.7999795083850525E-8d,
		+1.184154459111628E-8d,
		-3.43338934259355E-8d,
		+1.1795268640216787E-8d,
		+4.438921624363781E-8d,
		+2.925681159240093E-8d,
		-2.6437112632041807E-8d,
		+2.2860509143963117E-8d,
		-4.813899778443457E-9d,
		+3.6725170580355583E-9d,
		+2.0217439756338078E-10d
	};

	/** Tangent table, used by atan() (high bits). */
	static const std::vector<double> TANGENT_TABLE_A =
	{
		+0.0,
		+0.1256551444530487d,
		+0.25534194707870483d,
		+0.3936265707015991d,
		+0.5463024377822876d,
		+0.7214844226837158d,
		+0.9315965175628662d,
		+1.1974215507507324d,
		+1.5574076175689697d,
		+2.092571258544922d,
		+3.0095696449279785d,
		+5.041914939880371d,
		+14.101419448852539d,
		-18.430862426757812d
	};

	/** Tangent table, used by atan() (low bits). */
	static const std::vector<double> TANGENT_TABLE_B =
	{
		+0.0,
		-7.877917738262007E-9d,
		-2.5857668567479893E-8d,
		+5.2240336371356666E-9d,
		+5.206150291559893E-8d,
		+1.8307188599677033E-8d,
		-5.7618793749770706E-8d,
		+7.848361555046424E-8d,
		+1.0708593250394448E-7d,
		+1.7827257129423813E-8d,
		+2.893485277253286E-8d,
		+3.1660099222737955E-7d,
		+4.983191803254889E-7d,
		-3.356118100840571E-7d
	};

	/** Bits of 1/(2*pi), need for reducePayneHanek(). */
	static const std::vector<long> RECIP_2PI =
	{
		(0x28be60dbL << 32) | 0x9391054aL, (0x7f09d5f4L << 32) | 0x7d4d3770L, (0x36d8a566L << 32) | 0x4f10e410L, (0x7f9458eaL << 32) | 0xf7aef158L, (0x6dc91b8eL << 32) | 0x909374b8L, (0x01924bbaL << 32) | 0x82746487L, (0x3f877ac7L << 32) | 0x2c4a69cfL, (0xba208d7dL << 32) | 0x4baed121L, (0x3a671c09L << 32) | 0xad17df90L, (0x4e64758eL << 32) | 0x60d4ce7dL, (0x272117e2L << 32) | 0xef7e4a0eL, (0xc7fe25ffL << 32) | 0xf7816603L, (0xfbcbc462L << 32) | 0xd6829b47L, (0xdb4d9fb3L << 32) | 0xc9f2c26dL, (0xd3d18fd9L << 32) | 0xa797fa8bL, (0x5d49eeb1L << 32) | 0xfaf97c5eL, (0xcf41ce7dL << 32) | 0xe294a4baL, 0x9afed7ecL << 32 };

	/** Bits of pi/4, need for reducePayneHanek(). */
	static const std::vector<long> PI_O_4_BITS =
	{
		(0xc90fdaa2L << 32) | 0x2168c234L, (0xc4c6628bL << 32) | 0x80dc1cd1L };

	/** Eighths.
	 * This is used by sinQ, because its faster to do a table lookup than
	 * a multiply in this time-critical routine
	 */
	static const std::vector<double> EIGHTHS =
	{
		0,
		0.125,
		0.25,
		0.375,
		0.5,
		0.625,
		0.75,
		0.875,
		1.0,
		1.125,
		1.25,
		1.375,
		1.5,
		1.625
	};

	/** Table of 2^((n+2)/3) */
	static const std::vector<double> CBRTTWO =
	{
		0.6299605249474366,
		0.7937005259840998,
		1.0,
		1.2599210498948732,
		1.5874010519681994
	};

	/*
	 *  There are 52 bits in the mantissa of a double.
	 *  For additional precision, the code splits double numbers into two parts, *  by clearing the low order 30 bits if possible, and then performs the arithmetic
	 *  on each half separately.
	 */

	 /**
	  * 0x40000000 - used to split a double into two parts, both with the low order bits cleared.
	  * Equivalent to 2^30.
	  */
	static constexpr long HEX_40000000{ 0x40000000L }; // 1073741824L

	/** Mask used to clear low order 30 bits */
	static constexpr long MASK_30BITS{ -1L - (HEX_40000000 - 1) }; // 0xFFFFFFFFC0000000L;

	/** Mask used to clear the non-sign part of an int. */
	static constexpr int MASK_NON_SIGN_INT{ 0x7fffffff };

	/** Mask used to clear the non-sign part of a long. */
	static constexpr long MASK_NON_SIGN_LONG{ 0x7fffffffffffffffl };

	/** Mask used to extract exponent from double bits. */
	static constexpr long MASK_DOUBLE_EXPONENT{ 0x7ff0000000000000L };

	/** Mask used to extract mantissa from double bits. */
	static constexpr long MASK_DOUBLE_MANTISSA{ 0x000fffffffffffffL };

	/** Mask used to add implicit high order bit for normalized double. */
	static constexpr long IMPLICIT_HIGH_BIT{ 0x0010000000000000L };

	/** 2^52 - double numbers this large must be integral (no fraction) or NaN or Infinite */
	static constexpr double TWO_POWER_52{ 4503599627370496.0 };

	/** Constant: {@value}. */
	static constexpr double F_1_3{ 1.0 / 3 };
	/** Constant: {@value}. */
	static constexpr double F_1_5{ 1.0 / 5.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_7{ 1.0 / 7.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_9{ 1.0 / 9.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_11{ 1.0 / 11.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_13{ 1.0 / 13.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_15{ 1.0 / 15. };
	/** Constant: {@value}. */
	static constexpr double F_1_17{ 1.0 / 17. };
	/** Constant: {@value}. */
	static constexpr double F_3_4{ 3.0 / 4.0 };
	/** Constant: {@value}. */
	static constexpr double F_15_16{ 15.0 / 16.0 };
	/** Constant: {@value}. */
	static constexpr double F_13_14{ 13.0 / 14.0 };
	/** Constant: {@value}. */
	static constexpr double F_11_12{ 11.0 / 12.0 };
	/** Constant: {@value}. */
	static constexpr double F_9_10{ 9.0 / 10.0 };
	/** Constant: {@value}. */
	static constexpr double F_7_8{ 7.0 / 8.0 };
	/** Constant: {@value}. */
	static constexpr double F_5_6{ 5.0 / 6.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_2{ 1.0 / 2.0 };
	/** Constant: {@value}. */
	static constexpr double F_1_4{ 1.0 / 4.0 };

	/**
	 * Private Constructor
	 */
	FastMath() = default

		// Generic helper methods

		/**
		 * Get the high order bits from the mantissa.
		 * Equivalent to adding and subtracting HEX_40000 but also works for very large numbers
		 *
		 * @param d the value to split
		 * @return the high order part of the mantissa
		 */
		static double doubleHigh_part(const double& d)
	{
		if (d > -Precision.SAFE_MIN && d < Precision.SAFE_MIN)
		{
			return d; // These are un-normalised - don't try to convert
		}
		long xl = Double.double_to_raw_long_bits(d); // can take raw bits because just gonna convert it back
		xl &= MASK_30BITS; // Drop low order bits
		return Double.long_bits_to_double(xl);
	}

public:
	/** Compute the square root of a number.
	 * <p><b>Note:</b> this implementation currently delegates to {@link Math#sqrt}
	 * @param a number on which evaluation is done
	 * @return square root of a
	 */
	static double sqrt(const double& a) const
	{
		return std::sqrt(a);
	}

	/** Compute the hyperbolic cosine of a number.
	 * @param x number on which evaluation is done
	 * @return hyperbolic cosine of x
	 */
	static double cosh(const double& x)
	{
		if (std::isnan(x))
		{
			return x;
		}

		// cosh[z] = (exp(z) + exp(-z))/2

		// for numbers with magnitude 20 or so, // exp(-z) can be ignored in comparison with exp(z)

		if (x > 20)
		{
			if (x >= LOG_MAX_VALUE)
			{
				// Avoid overflow (MATH-905).
				const double t = exp(0.5 * x);
				return (0.5 * t) * t;
			}
			return 0.5 * exp(x);
		}
		if (x < -20)
		{
			if (x <= -LOG_MAX_VALUE)
			{
				// Avoid overflow (MATH-905).
				const double t = exp(-0.5 * x);
				return (0.5 * t) * t;
			}
			return 0.5 * exp(-x);
		}

		auto hiPrec = std::vector<double>(2);
		if (x < 0.0)
		{
			x = -x;
		}
		exp(x, 0.0, hiPrec);

		double ya = hiPrec[0] + hiPrec[1];
		double yb = -(ya - hiPrec[0] - hiPrec[1]);

		double temp = ya * HEX_40000000;
		double yaa = ya + temp - temp;
		double yab = ya - yaa;

		// recip = 1/y
		double recip = 1.0 / ya;
		temp = recip * HEX_40000000;
		double recipa = recip + temp - temp;
		double recipb = recip - recipa;

		// Correct for rounding in division
		recipb += (1.0 - yaa * recipa - yaa * recipb - yab * recipa - yab * recipb) * recip;
		// Account for yb
		recipb += -yb * recip * recip;

		// y = y + 1/y
		temp = ya + recipa;
		yb += -(temp - ya - recipa);
		ya = temp;
		temp = ya + recipb;
		yb += -(temp - ya - recipb);
		ya = temp;

		double result = ya + yb;
		result *= 0.5;
		return result;
	}

	/** Compute the hyperbolic sine of a number.
	 * @param x number on which evaluation is done
	 * @return hyperbolic sine of x
	 */
	static double sinh(const double& x)
	{
		if (std::isnan(x))
		{
			return x;
		}
		bool negate{};
		// sinh[z] = (exp(z) - exp(-z) / 2

		// for values of z larger than about 20, // exp(-z) can be ignored in comparison with exp(z)

		if (x > 20)
		{
			if (x >= LOG_MAX_VALUE)
			{
				// Avoid overflow (MATH-905).
				const double t = exp(0.5 * x);
				return (0.5 * t) * t;
			}
			return 0.5 * exp(x);
		}

		if (x < -20)
		{
			if (x <= -LOG_MAX_VALUE)
			{
				// Avoid overflow (MATH-905).
				const double t = exp(-0.5 * x);
				return (-0.5 * t) * t;
			}
			return -0.5 * exp(-x);
		}

		if (x == 0)
		{
			return x;
		}

		if (x < 0.0)
		{
			x = -x;
			negate = true;
		}

		auto hiPrec = std::vector<double>(2);
		double result;

		if (x > 0.25)
		{
			exp(x, 0.0, hiPrec);

			double ya = hiPrec[0] + hiPrec[1];
			double yb = -(ya - hiPrec[0] - hiPrec[1]);

			double temp = ya * HEX_40000000;
			double yaa = ya + temp - temp;
			double yab = ya - yaa;

			// recip = 1/y
			double recip = 1.0 / ya;
			temp = recip * HEX_40000000;
			double recipa = recip + temp - temp;
			double recipb = recip - recipa;

			// Correct for rounding in division
			recipb += (1.0 - yaa * recipa - yaa * recipb - yab * recipa - yab * recipb) * recip;
			// Account for yb
			recipb += -yb * recip * recip;

			recipa = -recipa;
			recipb = -recipb;

			// y = y - 1/y
			temp = ya + recipa;
			yb += -(temp - ya - recipa);
			ya = temp;
			temp = ya + recipb;
			yb += -(temp - ya - recipb);
			ya = temp;

			result = ya + yb;
			result *= 0.5;
		}
		else
		{
			expm1(x, hiPrec);

			double ya = hiPrec[0] + hiPrec[1];
			double yb = -(ya - hiPrec[0] - hiPrec[1]);

			/* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
			double denom = 1.0 + ya;
			double denomr = 1.0 / denom;
			double denomb = -(denom - 1.0 - ya) + yb;
			double ratio = ya * denomr;
			double temp = ratio * HEX_40000000;
			double ra = ratio + temp - temp;
			double rb = ratio - ra;

			temp = denom * HEX_40000000;
			double za = denom + temp - temp;
			double zb = denom - za;

			rb += (ya - za * ra - za * rb - zb * ra - zb * rb) * denomr;

			// Adjust for yb
			rb += yb * denomr;                        // numerator
			rb += -ya * denomb * denomr * denomr;   // denominator

			// y = y - 1/y
			temp = ya + ra;
			yb += -(temp - ya - ra);
			ya = temp;
			temp = ya + rb;
			yb += -(temp - ya - rb);
			ya = temp;

			result = ya + yb;
			result *= 0.5;
		}

		return negate
			? -result
			: result;
	}

	/**
	 * Combined hyperbolic sine and hyperbolic cosine function.
	 *
	 * @param x Argument.
	 * @return [sinh(x), cosh(x)]
	 */
	static Sinh_Cosh sinh_cosh(const double& x)
	{
		bool negate{};
		if (std::isnan(x))
		{
			return Sinh_Cosh(x, x);
		}

		// sinh[z] = (exp(z) - exp(-z) / 2
		// cosh[z] = (exp(z) + exp(-z))/2

		// for values of z larger than about 20, // exp(-z) can be ignored in comparison with exp(z)

		if (x > 20)
		{
			const double e;
			if (x >= LOG_MAX_VALUE)
			{
				// Avoid overflow (MATH-905).
				const double t = exp(0.5 * x);
				e = (0.5 * t) * t;
			}
			else
			{
				e = 0.5 * exp(x);
			}
			return Sinh_Cosh(e, e);
		}
		if (x < -20)
		{
			const double e;
			if (x <= -LOG_MAX_VALUE)
			{
				// Avoid overflow (MATH-905).
				const double t = exp(-0.5 * x);
				e = (-0.5 * t) * t;
			}
			else
			{
				e = -0.5 * exp(-x);
			}
			return Sinh_Cosh(e, -e);
		}

		if (x == 0)
		{
			return Sinh_Cosh(x, 1.0);
		}

		if (x < 0.0)
		{
			x = -x;
			negate = true;
		}

		auto hiPrec = std::vector<double>(2);
		double resultM;
		double resultP;

		if (x > 0.25)
		{
			exp(x, 0.0, hiPrec);

			const double ya = hiPrec[0] + hiPrec[1];
			const double yb = -(ya - hiPrec[0] - hiPrec[1]);

			double temp = ya * HEX_40000000;
			double yaa = ya + temp - temp;
			double yab = ya - yaa;

			// recip = 1/y
			double recip = 1.0 / ya;
			temp = recip * HEX_40000000;
			double recipa = recip + temp - temp;
			double recipb = recip - recipa;

			// Correct for rounding in division
			recipb += (1.0 - yaa * recipa - yaa * recipb - yab * recipa - yab * recipb) * recip;
			// Account for yb
			recipb += -yb * recip * recip;

			// y = y - 1/y
			temp = ya - recipa;
			double ybM = yb - (temp - ya + recipa);
			double ya_m = temp;
			temp = ya_m - recipb;
			ybM += -(temp - ya_m + recipb);
			ya_m = temp;
			resultM = ya_m + ybM;
			resultM *= 0.5;

			// y = y + 1/y
			temp = ya + recipa;
			double ybP = yb - (temp - ya - recipa);
			double yaP = temp;
			temp = yaP + recipb;
			ybP += -(temp - yaP - recipb);
			yaP = temp;
			resultP = yaP + ybP;
			resultP *= 0.5;
		}
		else
		{
			expm1(x, hiPrec);

			const double ya = hiPrec[0] + hiPrec[1];
			const double yb = -(ya - hiPrec[0] - hiPrec[1]);

			/* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
			double denom = 1.0 + ya;
			double denomr = 1.0 / denom;
			double denomb = -(denom - 1.0 - ya) + yb;
			double ratio = ya * denomr;
			double temp = ratio * HEX_40000000;
			double ra = ratio + temp - temp;
			double rb = ratio - ra;

			temp = denom * HEX_40000000;
			double za = denom + temp - temp;
			double zb = denom - za;

			rb += (ya - za * ra - za * rb - zb * ra - zb * rb) * denomr;

			// Adjust for yb
			rb += yb * denomr;                        // numerator
			rb += -ya * denomb * denomr * denomr;   // denominator

			// y = y - 1/y
			temp = ya + ra;
			double ybM = yb - (temp - ya - ra);
			double ya_m = temp;
			temp = ya_m + rb;
			ybM += -(temp - ya_m - rb);
			ya_m = temp;
			resultM = ya_m + ybM;
			resultM *= 0.5;

			// y = y + 1/y + 2
			temp = ya - ra;
			double ybP = yb - (temp - ya + ra);
			double yaP = temp;
			temp = yaP - rb;
			ybP += -(temp - yaP + rb);
			yaP = temp;
			resultP = yaP + ybP + 2;
			resultP *= 0.5;
		}

		if (negate)
		{
			resultM = -resultM;
		}

		return Sinh_Cosh(resultM, resultP);
	}

	/**
	 * Combined hyperbolic sine and hyperbolic cosine function.
	 *
	 * @param x Argument.
	 * @param <T> the type of the field element
	 * @return [sinh(x), cosh(x)]
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Sinh_Cosh<T> sinh_cosh(T x)
	{
		return x.sinh_cosh();
	}

	/** Compute the hyperbolic tangent of a number.
	 * @param x number on which evaluation is done
	 * @return hyperbolic tangent of x
	 */
	static double tanh(const double& x)
	{
		bool negate{};

		if (std::isnan(x))
		{
			return x;
		}

		// tanh[z] = sinh[z] / cosh[z]
		// = (exp(z) - exp(-z)) / (exp(z) + exp(-z))
		// = (exp(2x) - 1) / (exp(2x) + 1)

		// for magnitude > 20, sinh[z] == cosh[z] in double precision

		if (x > 20.0)
		{
			return 1.0;
		}

		if (x < -20)
		{
			return -1.0;
		}

		if (x == 0)
		{
			return x;
		}

		if (x < 0.0)
		{
			x = -x;
			negate = true;
		}

		double result;
		if (x >= 0.5)
		{
			auto hiPrec = std::vector<double>(2);
			// tanh(x) = (exp(2x) - 1) / (exp(2x) + 1)
			exp(x * 2.0, 0.0, hiPrec);

			double ya = hiPrec[0] + hiPrec[1];
			double yb = -(ya - hiPrec[0] - hiPrec[1]);

			/* Numerator */
			double na = -1.0 + ya;
			double nb = -(na + 1.0 - ya);
			double temp = na + yb;
			nb += -(temp - na - yb);
			na = temp;

			/* Denominator */
			double da = 1.0 + ya;
			double db = -(da - 1.0 - ya);
			temp = da + yb;
			db += -(temp - da - yb);
			da = temp;

			temp = da * HEX_40000000;
			double daa = da + temp - temp;
			double dab = da - daa;

			// ratio = na/da
			double ratio = na / da;
			temp = ratio * HEX_40000000;
			double ratioa = ratio + temp - temp;
			double ratiob = ratio - ratioa;

			// Correct for rounding in division
			ratiob += (na - daa * ratioa - daa * ratiob - dab * ratioa - dab * ratiob) / da;

			// Account for nb
			ratiob += nb / da;
			// Account for db
			ratiob += -db * na / da / da;

			result = ratioa + ratiob;
		}
		else
		{
			auto hiPrec = std::vector<double>(2);
			// tanh(x) = expm1(2x) / (expm1(2x) + 2)
			expm1(x * 2.0, hiPrec);

			double ya = hiPrec[0] + hiPrec[1];
			double yb = -(ya - hiPrec[0] - hiPrec[1]);

			/* Numerator */
			double na = ya;
			double nb = yb;

			/* Denominator */
			double da = 2.0 + ya;
			double db = -(da - 2.0 - ya);
			double temp = da + yb;
			db += -(temp - da - yb);
			da = temp;

			temp = da * HEX_40000000;
			double daa = da + temp - temp;
			double dab = da - daa;

			// ratio = na/da
			double ratio = na / da;
			temp = ratio * HEX_40000000;
			double ratioa = ratio + temp - temp;
			double ratiob = ratio - ratioa;

			// Correct for rounding in division
			ratiob += (na - daa * ratioa - daa * ratiob - dab * ratioa - dab * ratiob) / da;

			// Account for nb
			ratiob += nb / da;
			// Account for db
			ratiob += -db * na / da / da;

			result = ratioa + ratiob;
		}

		if (negate)
		{
			result = -result;
		}

		return result;
	}

	/** Compute the inverse hyperbolic cosine of a number.
	 * @param a number on which evaluation is done
	 * @return inverse hyperbolic cosine of a
	 */
	static double acosh(const double& a)
	{
		return std::log(a + std::sqrt(a * a - 1));
	}

	/** Compute the inverse hyperbolic sine of a number.
	 * @param a number on which evaluation is done
	 * @return inverse hyperbolic sine of a
	 */
	static double asinh(const double& a)
	{
		bool negative{};
		if (a < 0)
		{
			negative = true;
			a = -a;
		}

		double absAsinh;
		if (a > 0.167)
		{
			absAsinh = std::log(std::sqrt(a * a + 1) + a);
		}
		else
		{
			const double a2 = a * a;
			if (a > 0.097)
			{
				absAsinh = a * (1 - a2 * (F_1_3 - a2 * (F_1_5 - a2 * (F_1_7 - a2 * (F_1_9 - a2 * (F_1_11 - a2 * (F_1_13 - a2 * (F_1_15 - a2 * F_1_17 * F_15_16) * F_13_14) * F_11_12) * F_9_10) * F_7_8) * F_5_6) * F_3_4) * F_1_2);
			}
			else if (a > 0.036)
			{
				absAsinh = a * (1 - a2 * (F_1_3 - a2 * (F_1_5 - a2 * (F_1_7 - a2 * (F_1_9 - a2 * (F_1_11 - a2 * F_1_13 * F_11_12) * F_9_10) * F_7_8) * F_5_6) * F_3_4) * F_1_2);
			}
			else if (a > 0.0036)
			{
				absAsinh = a * (1 - a2 * (F_1_3 - a2 * (F_1_5 - a2 * (F_1_7 - a2 * F_1_9 * F_7_8) * F_5_6) * F_3_4) * F_1_2);
			}
			else
			{
				absAsinh = a * (1 - a2 * (F_1_3 - a2 * F_1_5 * F_3_4) * F_1_2);
			}
		}

		return negative ? -absAsinh : absAsinh;
	}

	/** Compute the inverse hyperbolic tangent of a number.
	 * @param a number on which evaluation is done
	 * @return inverse hyperbolic tangent of a
	 */
	static double atanh(const double& a)
	{
		bool negative{};
		if (a < 0)
		{
			negative = true;
			a = -a;
		}

		double absAtanh;
		if (a > 0.15)
		{
			absAtanh = 0.5 * std::log((1 + a) / (1 - a));
		}
		else
		{
			const double a2 = a * a;
			if (a > 0.087)
			{
				absAtanh = a * (1 + a2 * (F_1_3 + a2 * (F_1_5 + a2 * (F_1_7 + a2 * (F_1_9 + a2 * (F_1_11 + a2 * (F_1_13 + a2 * (F_1_15 + a2 * F_1_17))))))));
			}
			else if (a > 0.031)
			{
				absAtanh = a * (1 + a2 * (F_1_3 + a2 * (F_1_5 + a2 * (F_1_7 + a2 * (F_1_9 + a2 * (F_1_11 + a2 * F_1_13))))));
			}
			else if (a > 0.003)
			{
				absAtanh = a * (1 + a2 * (F_1_3 + a2 * (F_1_5 + a2 * (F_1_7 + a2 * F_1_9))));
			}
			else
			{
				absAtanh = a * (1 + a2 * (F_1_3 + a2 * F_1_5));
			}
		}

		return negative
			? -absAtanh
			: absAtanh;
	}

	/** Compute the signum of a number.
	 * The signum is -1 for negative numbers, +1 for positive numbers and 0 otherwise
	 * @param a number on which evaluation is done
	 * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
	 */
	static double signum(const double& a)
	{
		return (a < 0.0) ? -1.0 : ((a > 0.0) ? 1.0 : a); // return +0.0/-0.0/NaN depending on a
	}

	/** Compute the signum of a number.
	 * The signum is -1 for negative numbers, +1 for positive numbers and 0 otherwise
	 * @param a number on which evaluation is done
	 * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
	 */
	static float signum(const float& a)
	{
		// return +0.0/-0.0/NaN depending on a
		return (a < 0.0f)
			? -1.0f
			: ((a > 0.0f) ? 1.0f : a);
	}

	/** Compute next number towards positive infinity.
	 * @param a number to which neighbor should be computed
	 * @return neighbor of a towards positive infinity
	 */
	static double next_up(const double& a)
	{
		return next_after(a, INFINITY);
	}

	/** Compute next number towards positive infinity.
	 * @param a number to which neighbor should be computed
	 * @return neighbor of a towards positive infinity
	 */
	static float next_up(const float& a)
	{
		return next_after(a, Float.POSITIVE_INFINITY);
	}

	/** Compute next number towards negative infinity.
	 * @param a number to which neighbor should be computed
	 * @return neighbor of a towards negative infinity
	 */
	static double next_down(const double& a)
	{
		return next_after(a, -INFINITY);
	}

	/** Compute next number towards negative infinity.
	 * @param a number to which neighbor should be computed
	 * @return neighbor of a towards negative infinity
	 */
	static float next_down(const float& a)
	{
		return next_after(a, Float.NEGATIVE_INFINITY);
	}

	/** Returns a pseudo-random number between 0.0 and 1.0.
	 * <p><b>Note:</b> this implementation currently delegates to {@link Math#random}
	 * @return a random number between 0.0 and 1.0
	 */
	static double random()
	{
		return Math.random();
	}

	/**
	 * Exponential function.
	 *
	 * Computes exp(x), function result is nearly rounded.   It will be correctly
	 * rounded to the theoretical value for 99.9% of input values, otherwise it will
	 * have a 1 ULP error.
	 *
	 * Method:
	 *    Lookup intVal = exp(int(x))
	 *    Lookup fracVal = exp(int(x-int(x) / 1024.0) * 1024.0 );
	 *    Compute z as the exponential of the remaining bits by a polynomial minus one
	 *    exp(x) = intVal * fracVal * (1 + z)
	 *
	 * Accuracy:
	 *    Calculation is done with 63 bits of precision, so result should be correctly
	 *    rounded for 99.9% of input values, with less than 1 ULP error otherwise.
	 *
	 * @param x   a double
	 * @return double e<sup>x</sup>
	 */
	static double exp(const double& x)
	{
		return exp(x, 0.0, NULL);
	}

private:
	/**
	 * Internal helper method for exponential function.
	 * @param x original argument of the exponential function
	 * @param extra extra bits of precision on input (To Be Confirmed)
	 * @param hiPrec extra bits of precision on output (To Be Confirmed)
	 * @return exp(x)
	 */
	static double exp(const double& x, const double& extra, const std::vector<double>& hiPrec)
	{
		double int_part_a;
		double int_part_b;
		int intVal = static_cast<int>(x;

		/* Lookup exp(floor(x)).
		 * int_part_a will have the upper 22 bits, int_part_b will have the lower
		 * 52 bits.
		 */
		if (x < 0.0)
		{
			// We don't check against intVal here as conversion of large negative double values
			// may be affected by a JIT bug. Subsequent comparisons can safely use intVal
			if (x < -746d)
			{
				if (hiPrec != NULL)
				{
					hiPrec[0] = 0.0;
					hiPrec[1] = 0.0;
				}
				return 0.0;
			}

			if (intVal < -709)
			{
				/* This will produce a subnormal output */
				const double result = exp(x + 40.19140625, extra, hiPrec) / 285040095144011776.0;
				if (hiPrec != NULL)
				{
					hiPrec[0] /= 285040095144011776.0;
					hiPrec[1] /= 285040095144011776.0;
				}
				return result;
			}

			if (intVal == -709)
			{
				/* exp(1.494140625) is nearly a machine number... */
				const double result = exp(x + 1.494140625, extra, hiPrec) / 4.455505956692756620;
				if (hiPrec != NULL)
				{
					hiPrec[0] /= 4.455505956692756620;
					hiPrec[1] /= 4.455505956692756620;
				}
				return result;
			}

			intVal--;
		}
		else
		{
			if (intVal > 709)
			{
				if (hiPrec != NULL)
				{
					hiPrec[0] = INFINITY;
					hiPrec[1] = 0.0;
				}
				return INFINITY;
			}
		}

		int_part_a = ExpIntTable.EXP_INT_TABLE_A[EXP_INT_TABLE_MAX_INDEX + intVal];
		int_part_b = ExpIntTable.EXP_INT_TABLE_B[EXP_INT_TABLE_MAX_INDEX + intVal];

		/* Get the fractional part of x, find the greatest multiple of 2^-10 less than
		 * x and look up the exp function of it.
		 * frac_part_a will have the upper 22 bits, frac_part_b the lower 52 bits.
		 */
		const int intFrac = static_cast<int>(((x - intVal) * 1024.0);
		const double frac_part_a = Exp_fracTable.EXP_FRAC_TABLE_A[intFrac];
		const double frac_part_b = Exp_fracTable.EXP_FRAC_TABLE_B[intFrac];

		/* epsilon is the difference in x from the nearest multiple of 2^-10.  It
		 * has a value in the range 0 <= epsilon < 2^-10.
		 * Do the subtraction from x as the last step to avoid possible loss of precision.
		 */
		const double epsilon = x - (intVal + intFrac / 1024.0);

		/* Compute z = exp(epsilon) - 1.0 via a minimax polynomial.  z has
	   full double precision (52 bits).  sin_ce z < 2^-10, we will have
	   62 bits of precision when combined with the constant 1.  This will be
	   used in the last addition below to get proper rounding. */

	   /* Remez generated polynomial.  Converges on the interval [0, 2^-10], error
	  is less than 0.5 ULP */
		double z = 0.04168701738764507;
		z = z * epsilon + 0.1666666505023083;
		z = z * epsilon + 0.5000000000042687;
		z = z * epsilon + 1.0;
		z = z * epsilon + -3.940510424527919E-20;

		/* Compute (int_part_a+int_part_b) * (frac_part_a+frac_part_b) by binomial
	   expansion.
	   tempA is exact since int_part_a and int_part_b only have 22 bits each.
	   tempB will have 52 bits of precision.
		 */
		double tempA = int_part_a * frac_part_a;
		double tempB = int_part_a * frac_part_b + int_part_b * frac_part_a + int_part_b * frac_part_b;

		/* Compute the result.  (1+z)(tempA+tempB).  Order of operations is
	   important.  For accuracy add by increasing size.  tempA is exact and
	   much larger than the others.  If there are extra bits specified from the
	   pow() function, use them. */
		const double tempC = tempB + tempA;

		// If tempC is positive infinite, the evaluation below could result in NaN, // because z could be negative at the same time.
		if (tempC == INFINITY)
		{
			return INFINITY;
		}

		const double result;
		if (extra != 0.0)
		{
			result = tempC * extra * z + tempC * extra + tempC * z + tempB + tempA;
		}
		else
		{
			result = tempC * z + tempB + tempA;
		}

		if (hiPrec != NULL)
		{
			// If requesting high precision
			hiPrec[0] = tempA;
			hiPrec[1] = tempC * extra * z + tempC * extra + tempC * z + tempB;
		}

		return result;
	}

	/** Compute exp(x) - 1
	 * @param x number to compute shifted exponential
	 * @return exp(x) - 1
	 */
	public static double expm1(const double& x)
	{
		return expm1(x, NULL);
	}

	/** Internal helper method for expm1
	 * @param x number to compute shifted exponential
	 * @param hi_prec_out receive high precision result for -1.0 < x < 1.0
	 * @return exp(x) - 1
	 */
	private static double expm1(const double& x, const std::vector<double>& hi_prec_out)
	{
		if (std::isnan(x) || x == 0.0) { // NaN or zero
			return x;
		}

		if (x <= -1.0 || x >= 1.0)
		{
			// If not between +/- 1.0
			//return exp(x) - 1.0;
			auto hiPrec = std::vector<double>(2);
			exp(x, 0.0, hiPrec);
			if (x > 0.0)
			{
				return -1.0 + hiPrec[0] + hiPrec[1];
			}
			else
			{
				const double ra = -1.0 + hiPrec[0];
				double rb = -(ra + 1.0 - hiPrec[0]);
				rb += hiPrec[1];
				return ra + rb;
			}
		}

		double baseA;
		double baseB;
		double epsilon;
		bool negative = false;

		if (x < 0.0)
		{
			x = -x;
			negative = true;
		}

		{
			int intFrac = static_cast<int>((x * 1024.0);
			double tempA = Exp_fracTable.EXP_FRAC_TABLE_A[intFrac] - 1.0;
			double tempB = Exp_fracTable.EXP_FRAC_TABLE_B[intFrac];

			double temp = tempA + tempB;
			tempB = -(temp - tempA - tempB);
			tempA = temp;

			temp = tempA * HEX_40000000;
			baseA = tempA + temp - temp;
			baseB = tempB + (tempA - baseA);

			epsilon = x - intFrac / 1024.0;
		}

		/* Compute expm1(epsilon) */
		double zb = 0.008336750013465571;
		zb = zb * epsilon + 0.041666663879186654;
		zb = zb * epsilon + 0.16666666666745392;
		zb = zb * epsilon + 0.49999999999999994;
		zb *= epsilon;
		zb *= epsilon;

		double za = epsilon;
		double temp = za + zb;
		zb = -(temp - za - zb);
		za = temp;

		temp = za * HEX_40000000;
		temp = za + temp - temp;
		zb += za - temp;
		za = temp;

		/* Combine the parts.   expm1(a+b) = expm1(a) + expm1(b) + expm1(a)*expm1(b) */
		double ya = za * baseA;
		//double yb = za*baseB + zb*baseA + zb*baseB;
		temp = ya + za * baseB;
		double yb = -(temp - ya - za * baseB);
		ya = temp;

		temp = ya + zb * baseA;
		yb += -(temp - ya - zb * baseA);
		ya = temp;

		temp = ya + zb * baseB;
		yb += -(temp - ya - zb * baseB);
		ya = temp;

		//ya = ya + za + baseA;
		//yb = yb + zb + baseB;
		temp = ya + baseA;
		yb += -(temp - baseA - ya);
		ya = temp;

		temp = ya + za;
		//yb += (ya > za) ? -(temp - ya - za) : -(temp - za - ya);
		yb += -(temp - ya - za);
		ya = temp;

		temp = ya + baseB;
		//yb += (ya > baseB) ? -(temp - ya - baseB) : -(temp - baseB - ya);
		yb += -(temp - ya - baseB);
		ya = temp;

		temp = ya + zb;
		//yb += (ya > zb) ? -(temp - ya - zb) : -(temp - zb - ya);
		yb += -(temp - ya - zb);
		ya = temp;

		if (negative)
		{
			/* Compute expm1(-x) = -expm1(x) / (expm1(x) + 1) */
			double denom = 1.0 + ya;
			double denomr = 1.0 / denom;
			double denomb = -(denom - 1.0 - ya) + yb;
			double ratio = ya * denomr;
			temp = ratio * HEX_40000000;
			const double ra = ratio + temp - temp;
			double rb = ratio - ra;

			temp = denom * HEX_40000000;
			za = denom + temp - temp;
			zb = denom - za;

			rb += (ya - za * ra - za * rb - zb * ra - zb * rb) * denomr;

			// f(x) = x/1+x
			// Compute f'(x)
			// Product rule:  d(uv) = du*v + u*dv
			// Chain rule:  d(f(g(x)) = f'(g(x))*f(g'(x))
			// d(1/x) = -1/(x*x)
			// d(1/1+x) = -1/( (1+x)^2) *  1 =  -1/((1+x)*(1+x))
			// d(x/1+x) = -x/((1+x)(1+x)) + 1/1+x = 1 / ((1+x)(1+x))

			// Adjust for yb
			rb += yb * denomr;                      // numerator
			rb += -ya * denomb * denomr * denomr;   // denominator

			// negate
			ya = -ra;
			yb = -rb;
		}

		if (hi_prec_out != NULL)
		{
			hi_prec_out[0] = ya;
			hi_prec_out[1] = yb;
		}

		return ya + yb;
	}

	/**
	 * Natural logarithm.
	 *
	 * @param x   a double
	 * @return log(x)
	 */
	public static double log(const double& x)
	{
		return log(x, NULL);
	}

	/**
	 * Internal helper method for natural logarithm function.
	 * @param x original argument of the natural logarithm function
	 * @param hiPrec extra bits of precision on output (To Be Confirmed)
	 * @return log(x)
	 */
	private static double log(const double& x, const std::vector<double>& hiPrec)
	{
		if (x == 0)   // Handle special case of +0/-0
		{
			return -INFINITY;
		}
		long bits = Double.double_to_raw_long_bits(x);

		/* Handle special cases of negative input, and NaN */
		if (((bits & 0x8000000000000000L) != 0 || std::isnan(x)) && x != 0.0)
		{
			if (hiPrec != NULL)
			{
				hiPrec[0] = std::numeric_limits<double>::quiet_NaN();
			}

			return std::numeric_limits<double>::quiet_NaN();
		}

		/* Handle special cases of Positive infinity. */
		if (x == INFINITY)
		{
			if (hiPrec != NULL)
			{
				hiPrec[0] = INFINITY;
			}

			return INFINITY;
		}

		/* Extract the exponent */
		int exp = static_cast<int>((bits >> 52) - 1023;

		if ((bits & 0x7ff0000000000000L) == 0)
		{
			// Subnormal!
			if (x == 0)
			{
				// Zero
				if (hiPrec != NULL)
				{
					hiPrec[0] = -INFINITY;
				}

				return -INFINITY;
			}

			/* Normalize the subnormal number. */
			bits <<= 1;
			while ((bits & 0x0010000000000000L) == 0)
			{
				--exp;
				bits <<= 1;
			}
		}

		if ((exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && hiPrec == NULL)
		{
			/* The normal method doesn't work well in the range [0.99, 1.01], so call do a straight
		   polynomial expansion in higer precision. */

		   /* Compute x - 1.0 and split it */
			double xa = x - 1.0;
			double tmp = xa * HEX_40000000;
			double aa = xa + tmp - tmp;
			double ab = xa - aa;
			xa = aa;
			double xb = ab;

			const std::vector<double> ln_coef_last = LN_QUICK_COEF[LN_QUICK_COEF.size() - 1];
			double ya = ln_coef_last[0];
			double yb = ln_coef_last[1];

			for (int i = LN_QUICK_COEF.size() - 2; i >= 0; i--)
			{
				/* Multiply a = y * x */
				aa = ya * xa;
				ab = ya * xb + yb * xa + yb * xb;
				/* split, so now y = a */
				tmp = aa * HEX_40000000;
				ya = aa + tmp - tmp;
				yb = aa - ya + ab;

				/* Add  a = y + lnQuickCoef */
				const std::vector<double> ln_coef_i = LN_QUICK_COEF[i];
				aa = ya + ln_coef_i[0];
				ab = yb + ln_coef_i[1];
				/* Split y = a */
				tmp = aa * HEX_40000000;
				ya = aa + tmp - tmp;
				yb = aa - ya + ab;
			}

			/* Multiply a = y * x */
			aa = ya * xa;
			ab = ya * xb + yb * xa + yb * xb;
			/* split, so now y = a */
			tmp = aa * HEX_40000000;
			ya = aa + tmp - tmp;
			yb = aa - ya + ab;

			return ya + yb;
		}

		// lnm is a log of a number in the range of 1.0 - 2.0, so 0 <= lnm < ln(2)
		const std::vector<double> lnm = lnMant.LN_MANT[static_cast<int>(((bits & 0x000ffc0000000000L) >> 42)];

		/*
	double epsilon = x / Double.long_bits_to_double(bits & 0xfffffc0000000000L);

	epsilon -= 1.0;
		 */

		 // y is the most significant 10 bits of the mantissa
		 //double y = Double.long_bits_to_double(bits & 0xfffffc0000000000L);
		 //double epsilon = (x - y) / y;
		const double epsilon = (bits & 0x3ffffffffffL) / (TWO_POWER_52 + (bits & 0x000ffc0000000000L));

		double lnza;
		double lnzb = 0.0;

		if (hiPrec != NULL)
		{
			/* split epsilon -> x */
			double tmp = epsilon * HEX_40000000;
			double aa = epsilon + tmp - tmp;
			double ab = epsilon - aa;
			double xa = aa;
			double xb = ab;

			/* Need a more accurate epsilon, so adjust the division. */
			const double numer = bits & 0x3ffffffffffL;
			const double denom = TWO_POWER_52 + (bits & 0x000ffc0000000000L);
			aa = numer - xa * denom - xb * denom;
			xb += aa / denom;

			/* Remez polynomial evaluation */
			const std::vector<double> ln_coef_last = LN_HI_PREC_COEF[LN_HI_PREC_COEF.size() - 1];
			double ya = ln_coef_last[0];
			double yb = ln_coef_last[1];

			for (int i = LN_HI_PREC_COEF.size() - 2; i >= 0; i--)
			{
				/* Multiply a = y * x */
				aa = ya * xa;
				ab = ya * xb + yb * xa + yb * xb;
				/* split, so now y = a */
				tmp = aa * HEX_40000000;
				ya = aa + tmp - tmp;
				yb = aa - ya + ab;

				/* Add  a = y + lnHiPrecCoef */
				const std::vector<double> ln_coef_i = LN_HI_PREC_COEF[i];
				aa = ya + ln_coef_i[0];
				ab = yb + ln_coef_i[1];
				/* Split y = a */
				tmp = aa * HEX_40000000;
				ya = aa + tmp - tmp;
				yb = aa - ya + ab;
			}

			/* Multiply a = y * x */
			aa = ya * xa;
			ab = ya * xb + yb * xa + yb * xb;

			/* split, so now lnz = a */
			/*
	  tmp = aa * 1073741824.0;
	  lnza = aa + tmp - tmp;
	  lnzb = aa - lnza + ab;
			 */
			lnza = aa + ab;
			lnzb = -(lnza - aa - ab);
		}
		else
		{
			/* High precision not required.  Eval Remez polynomial
		 using standard double precision */
			lnza = -0.16624882440418567;
			lnza = lnza * epsilon + 0.19999954120254515;
			lnza = lnza * epsilon + -0.2499999997677497;
			lnza = lnza * epsilon + 0.3333333333332802;
			lnza = lnza * epsilon + -0.5;
			lnza = lnza * epsilon + 1.0;
			lnza *= epsilon;
		}

		/* Relative sizes:
		 * lnzb     [0, 2.33E-10]
		 * lnm[1]   [0, 1.17E-7]
		 * ln2B*exp [0, 1.12E-4]
		 * lnza      [0, 9.7E-4]
		 * lnm[0]   [0, 0.692]
		 * ln2A*exp [0, 709]
		 */

		 /* Compute the following sum:
		  * lnzb + lnm[1] + ln2B*exp + lnza + lnm[0] + ln2A*exp;
		  */

		  //return lnzb + lnm[1] + ln2B*exp + lnza + lnm[0] + ln2A*exp;
		double a = LN_2_A * exp;
		double b = 0.0;
		double c = a + lnm[0];
		double d = -(c - a - lnm[0]);
		a = c;
		b += d;

		c = a + lnza;
		d = -(c - a - lnza);
		a = c;
		b += d;

		c = a + LN_2_B * exp;
		d = -(c - a - LN_2_B * exp);
		a = c;
		b += d;

		c = a + lnm[1];
		d = -(c - a - lnm[1]);
		a = c;
		b += d;

		c = a + lnzb;
		d = -(c - a - lnzb);
		a = c;
		b += d;

		if (hiPrec != NULL)
		{
			hiPrec[0] = a;
			hiPrec[1] = b;
		}

		return a + b;
	}

	/**
	 * Computes log(1 + x).
	 *
	 * @param x Number.
	 * @return {@code log(1 + x)}.
	 */
	public static double log1p(const double& x)
	{
		if (x == -1)
		{
			return -INFINITY;
		}

		if (x == INFINITY)
		{
			return INFINITY;
		}

		if (x > 1e-6 ||
			x < -1e-6)
		{
			const double xpa = 1 + x;
			const double xpb = -(xpa - 1 - x);

			const std::vector<double> hiPrec = std::vector<double>(2);
			const double lores = log(xpa, hiPrec);
			if (std::isinf(lores)) { // Don't allow this to be converted to NaN
				return lores;
			}

			// Do a taylor series expansion around xpa:
			//   f(x+y) = f(x) + f'(x) y + f''(x)/2 y^2
			const double fx1 = xpb / xpa;
			const double epsilon = 0.5 * fx1 + 1;
			return epsilon * fx1 + hiPrec[1] + hiPrec[0];
		}
		else
		{
			// Value is small |x| < 1e6, do a Taylor series centered on 1.
			const double y = (x * F_1_3 - F_1_2) * x + 1;
			return y * x;
		}
	}

	/** Compute the base 10 logarithm.
	 * @param x a number
	 * @return log10(x)
	 */
	public static double log10(const double& x)
	{
		auto hiPrec = std::vector<double>(2);

		const double lores = log(x, hiPrec);
		if (std::isinf(lores)) { // don't allow this to be converted to NaN
			return lores;
		}

		const double tmp = hiPrec[0] * HEX_40000000;
		const double lna = hiPrec[0] + tmp - tmp;
		const double lnb = hiPrec[0] - lna + hiPrec[1];

		const double rln10a = 0.4342944622039795;
		const double rln10b = 1.9699272335463627E-8;

		return rln10b * lnb + rln10b * lna + rln10a * lnb + rln10a * lna;
	}

	/**
	 * Computes the <a href="http://mathworld.wolfram.com/Logarithm.html">
	 * logarithm</a> in a given base.
	 *
	 * Returns {@code NaN} if either argument is negative.
	 * If {@code base} is 0 and {@code x} is positive, 0 is returned.
	 * If {@code base} is positive and {@code x} is 0, * {@code -INFINITY} is returned.
	 * If both arguments are 0, the result is {@code NaN}.
	 *
	 * @param base Base of the logarithm, must be greater than 0.
	 * @param x Argument, must be greater than 0.
	 * @return the value of the logarithm, i.e. the number {@code y} such that
	 * <code>base<sup>y</sup> = x</code>.
	 */
	public static double log(const double& base, const double& x)
	{
		return log(x) / log(base);
	}

	/**
	 * Power function.  Compute x^y.
	 *
	 * @param x   a double
	 * @param y   a double
	 * @return double
	 */
	public static double pow(const double& x, const double y)
	{
		if (y == 0)
		{
			// y = -0 or y = +0
			return 1.0;
		}
		else
		{
			const long yBits = Double.double_to_raw_long_bits(y);
			const int  y_raw_exp = static_cast<int>(((yBits & MASK_DOUBLE_EXPONENT) >> 52);
			const long y_raw_mantissa = yBits & MASK_DOUBLE_MANTISSA;
			const long x_bits = Double.double_to_raw_long_bits(x);
			const int  x_raw_exp = static_cast<int>(((x_bits & MASK_DOUBLE_EXPONENT) >> 52);
			const long x_raw_mantissa = x_bits & MASK_DOUBLE_MANTISSA;

			if (y_raw_exp > 1085)
			{
				// y is either a very large integral value that does not fit in a long or it is a special number

				if ((y_raw_exp == 2047 && y_raw_mantissa != 0) ||
					(x_raw_exp == 2047 && x_raw_mantissa != 0))
				{
					// NaN
					return std::numeric_limits<double>::quiet_NaN();
				}
				else if (x_raw_exp == 1023 && x_raw_mantissa == 0)
				{
					// x = -1.0 or x = +1.0
					if (y_raw_exp == 2047)
					{
						// y is infinite
						return std::numeric_limits<double>::quiet_NaN();
					}
					else
					{
						// y is a large even integer
						return 1.0;
					}
				}
				else
				{
					// the absolute value of x is either greater or smaller than 1.0

					// if y_raw_exp == 2047 and mantissa is 0, y = -infinity or y = +infinity
					// if 1085 < y_raw_exp < 2047, y is simply a large number, however, due to limited
					// accuracy, at this magnitude it behaves just like infinity with regards to x
					if ((y > 0) ^ (x_raw_exp < 1023))
					{
						// either y = +infinity (or large engouh) and abs(x) > 1.0
						// or     y = -infinity (or large engouh) and abs(x) < 1.0
						return INFINITY;
					}
					else
					{
						// either y = +infinity (or large engouh) and abs(x) < 1.0
						// or     y = -infinity (or large engouh) and abs(x) > 1.0
						return +0.0;
					}
				}
			}
			else
			{
				// y is a regular non-zero number

				if (y_raw_exp >= 1023)
				{
					// y may be an integral value, which should be handled specifically
					const long y_full_mantissa = IMPLICIT_HIGH_BIT | y_raw_mantissa;
					if (y_raw_exp < 1075)
					{
						// normal number with negative shift that may have a fractional part
						const long integralMask = (-1L) << (1075 - y_raw_exp);
						if ((y_full_mantissa & integralMask) == y_full_mantissa)
						{
							// all fractional bits are 0, the number is really integral
							const long l = y_full_mantissa >> (1075 - y_raw_exp);
							return std::pow(x, (y < 0) ? -l : l);
						}
					}
					else
					{
						// normal number with positive shift, always an integral value
						// we know it fits in a primitive long because y_raw_exp > 1085 has been handled above
						const long l = y_full_mantissa << (y_raw_exp - 1075);
						return std::pow(x, (y < 0) ? -l : l);
					}
				}

				// y is a non-integral value

				if (x == 0)
				{
					// x = -0 or x = +0
					// the integer powers have already been handled above
					return y < 0 ? INFINITY : +0.0;
				}
				else if (x_raw_exp == 2047)
				{
					if (x_raw_mantissa == 0)
					{
						// x = -infinity or x = +infinity
						return (y < 0) ? +0.0 : INFINITY;
					}
					else
					{
						// NaN
						return std::numeric_limits<double>::quiet_NaN();
					}
				}
				else if (x < 0)
				{
					// the integer powers have already been handled above
					return std::numeric_limits<double>::quiet_NaN();
				}
				else
				{
					// this is the general case, for regular fractional numbers x and y

					// Split y into ya and yb such that y = ya+yb
					const double tmp = y * HEX_40000000;
					const double ya = (y + tmp) - tmp;
					const double yb = y - ya;

					/* Compute ln(x) */
					auto lns = std::vector<double>(2);
					const double lores = log(x, lns);
					if (std::isinf(lores)) // don't allow this to be converted to NaN
					{
						return lores;
					}

					double lna = lns[0];
					double lnb = lns[1];

					/* resplit lns */
					const double tmp1 = lna * HEX_40000000;
					const double tmp2 = (lna + tmp1) - tmp1;
					lnb += lna - tmp2;
					lna = tmp2;

					// y*ln(x) = (aa+ab)
					const double& aa = lna * ya;
					const double& ab = lna * yb + lnb * ya + lnb * yb;

					lna = aa + ab;
					lnb = -(lna - aa - ab);

					double z = 1.0 / 120.0;
					z = z * lnb + (1.0 / 24.0);
					z = z * lnb + (1.0 / 6.0);
					z = z * lnb + 0.5;
					z = z * lnb + 1.0;
					z *= lnb;

					return exp(lna, z, NULL);
				}
			}
		}
	}

	/**
	 * Raise a double to an int power.
	 *
	 * @param d Number to raise.
	 * @param e Exponent.
	 * @return d<sup>e</sup>
	 */
	public static double pow(const double& d, const int& e)
	{
		return pow(d, static_cast<long>(e);
	}

	/**
	 * Raise a double to a long power.
	 *
	 * @param d Number to raise.
	 * @param e Exponent.
	 * @return d<sup>e</sup>
	 */
	public static double pow(const double& d, const long& e)
	{
		if (e == 0)
		{
			return 1.0;
		}
		if (e > 0)
		{
			return Split(d).pow(e).full;
		}
		return Split(d).reciprocal().pow(-e).full;
	}

	/** Class operator on double numbers split into one 26 bits number and one 27 bits number. */
	private static class Split
	{
	public:
		/** Split version of NaN. */
		static const Split NAN = Split(Double.NaN, 0);

		/** Split version of positive infinity. */
		static const Split POSITIVE_INFINITY = Split(INFINITY, 0);

		/** Split version of negative infinity. */
		static const Split NEGATIVE_INFINITY = Split(-INFINITY, 0);

	private:
		/** Full number. */
		const double my_full;

		/** High order bits. */
		const double my_high;

		/** Low order bits. */
		private const double my_low;

		/** Simple constructor.
		 * @param x number to split
		 */
		Split(const double& x)
			:
			my_full{ x },
			my_high{ Double.long_bits_to_double(Double.double_to_raw_long_bits(x) & ((-1L) << 27)) },
			my_low{ x - high }
		{
		}

		/** Simple constructor.
		 * @param high high order bits
		 * @param low low order bits
		 */
		Split(const double& high, const double& low)
		{
			Split(high == 0.0
				? (
					low == 0.0 && Double.double_to_raw_long_bits(high) == long.MIN_VALUE /* negative zero */
					? -0.0
					: low
					)
				: high + low, high, low);
		}

		/** Simple constructor.
		 * @param full full number
		 * @param high high order bits
		 * @param low low order bits
		 */
		Split(const double& full, const double& high, const double& low)
			:
			my_full{ full },
			my_high{ high },
			my_low{ low }
		{
		}

		/** Multiply the instance by another one.
		 * @param b other instance to multiply by
		 * @return product
		 */
		public Split multiply(const Split& b)
		{
			// beware the following expressions must NOT be simplified, they rely on floating point arithmetic properties
			const Split  mulBasic = Split(full * b.full);
			const double mulError = low * b.low - (((mulBasic.full - high * b.high) - low * b.high) - high * b.low);
			return Split(mulBasic.high, mulBasic.low + mulError);
		}

		/** Compute the reciprocal of the instance.
		 * @return reciprocal of the instance
		 */
		public Split reciprocal()
		{
			const double& approximateInv = 1.0 / full;
			const Split  splitInv = Split(approximateInv);

			// if 1.0/d were computed perfectly, remultiplying it by d should give 1.0
			// we want to estimate the error so we can fix the low order bits of approximateInvLow
			// beware the following expressions must NOT be simplified, they rely on floating point arithmetic properties
			const Split product = multiply(splitInv);
			const double error = (product.high - 1) + product.low;

			// better accuracy estimate of reciprocal
			return std::isnan(error) ? splitInv : Split(splitInv.high, splitInv.low - error / full);
		}

		/** Computes this^e.
		 * @param e exponent (beware, here it MUST be > 0; the only exclusion is long.MIN_VALUE)
		 * @return d^e, split in high and low bits
		 */
		private Split pow(const long& e)
		{
			// prepare result
			Split result = Split(1);

			// d^(2p)
			Split d2p = Split(full, high, low);

			for (long p{ e }; p != 0; p >>>= 1)
			{
				if ((p & 0x1) != 0)
				{
					// accurate multiplication result = result * d^(2p) using Veltkamp TwoProduct algorithm
					result = result.multiply(d2p);
				}

				// accurate squaring d^(2(p+1)) = d^(2p) * d^(2p) using Veltkamp TwoProduct algorithm
				d2p = d2p.multiply(d2p);
			}

			if (std::isnan(result.full))
			{
				if (std::isnan(full))
				{
					return Split.NAN;
				}
				else
				{
					// some intermediate numbers exceeded capacity, // and the low order bits became NaN (because infinity - infinity = NaN)
					if (std::abs(full) < 1)
					{
						return Split(std::copysign(0.0, full), 0.0);
					}
					else if (full < 0 && (e & 0x1) == 1)
					{
						return Split.NEGATIVE_INFINITY;
					}
					else
					{
						return Split.POSITIVE_INFINITY;
					}
				}
			}
			else
			{
				return result;
			}
		}
	}

	/**
	 *  Computes sin(x) - x, where |x| < 1/16.
	 *  Use a Remez polynomial approximation.
	 *  @param x a number smaller than 1/16
	 *  @return sin(x) - x
	 */
	private static double polySine(const double& x)

	{
		double x2 = x * x;

		double p = 2.7553817452272217E-6;
		p = p * x2 + -1.9841269659586505E-4;
		p = p * x2 + 0.008333333333329196;
		p = p * x2 + -0.16666666666666666;
		//p *= x2;
		//p *= x;
		p = p * x2 * x;

		return p;
	}

	/**
	 *  Computes cos(x) - 1, where |x| < 1/16.
	 *  Use a Remez polynomial approximation.
	 *  @param x a number smaller than 1/16
	 *  @return cos(x) - 1
	 */
	private static double polyCosine(const double& x)
	{
		double x2 = x * x;

		double p = 2.479773539153719E-5;
		p = p * x2 + -0.0013888888689039883;
		p = p * x2 + 0.041666666666621166;
		p = p * x2 + -0.49999999999999994;
		p *= x2;

		return p;
	}

	/**
	 *  Compute sine over the first quadrant (0 < x < pi/2).
	 *  Use combination of table lookup and rational polynomial expansion.
	 *  @param xa number from which sine is requested
	 *  @param xb extra bits for x (may be 0.0)
	 *  @return sin(xa + xb)
	 */
	private static double sinQ(const double& xa, const double& xb)
	{
		int idx = static_cast<int>(((xa * 8.0) + 0.5);
		const double epsilon = xa - EIGHTHS[idx]; //idx*0.125;

		// Table lookups
		const double sintA = SINE_TABLE_A[idx];
		const double sintB = SINE_TABLE_B[idx];
		const double costA = COSINE_TABLE_A[idx];
		const double costB = COSINE_TABLE_B[idx];

		// Polynomial eval of sin(epsilon), cos(epsilon)
		double sinEpsA = epsilon;
		double sin_eps_b = polySine(epsilon);
		const double cos_eps_a = 1.0;
		const double cos_eps_b = polyCosine(epsilon);

		// Split epsilon   xa + xb = x
		const double temp = sinEpsA * HEX_40000000;
		double temp2 = (sinEpsA + temp) - temp;
		sin_eps_b += sinEpsA - temp2;
		sinEpsA = temp2;

		/* Compute sin(x) by angle addition formula */
		double result;

		/* Compute the following sum:
		 *
		 * result = sintA + costA*sinEpsA + sintA*cos_eps_b + costA*sin_eps_b +
		 *          sintB + costB*sinEpsA + sintB*cos_eps_b + costB*sin_eps_b;
		 *
		 * Ranges of elements
		 *
		 * xxxtA   0            PI/2
		 * xxxtB   -1.5e-9      1.5e-9
		 * sinEpsA -0.0625      0.0625
		 * sin_eps_b -6e-11       6e-11
		 * cos_eps_a  1.0
		 * cos_eps_b  0           -0.0625
		 *
		 */

		 //result = sintA + costA*sinEpsA + sintA*cos_eps_b + costA*sin_eps_b +
		 //          sintB + costB*sinEpsA + sintB*cos_eps_b + costB*sin_eps_b;

		 //result = sintA + sintA*cos_eps_b + sintB + sintB * cos_eps_b;
		 //result += costA*sinEpsA + costA*sin_eps_b + costB*sinEpsA + costB * sin_eps_b;
		double a{};
		double b{};

		double t = sintA;
		double c = a + t;
		double d = -(c - a - t);
		a = c;
		b += d;

		t = costA * sinEpsA;
		c = a + t;
		d = -(c - a - t);
		a = c;
		b += d;

		b = b + sintA * cos_eps_b + costA * sin_eps_b;
		/*
	t = sintA*cos_eps_b;
	c = a + t;
	d = -(c - a - t);
	a = c;
	b = b + d;

	t = costA*sin_eps_b;
	c = a + t;
	d = -(c - a - t);
	a = c;
	b = b + d;
		 */

		b = b + sintB + costB * sinEpsA + sintB * cos_eps_b + costB * sin_eps_b;
		/*
	t = sintB;
	c = a + t;
	d = -(c - a - t);
	a = c;
	b = b + d;

	t = costB*sinEpsA;
	c = a + t;
	d = -(c - a - t);
	a = c;
	b = b + d;

	t = sintB*cos_eps_b;
	c = a + t;
	d = -(c - a - t);
	a = c;
	b = b + d;

	t = costB*sin_eps_b;
	c = a + t;
	d = -(c - a - t);
	a = c;
	b = b + d;
		 */

		if (xb != 0.0)
		{
			t = ((costA + costB) * (cos_eps_a + cos_eps_b) -
				(sintA + sintB) * (sinEpsA + sin_eps_b)) * xb;  // approximate cosine*xb
			c = a + t;
			d = -(c - a - t);
			a = c;
			b += d;
		}

		result = a + b;

		return result;
	}

	/**
	 * Compute cosine in the first quadrant by subtracting input from PI/2 and
	 * then calling sinQ.  This is more accurate as the input approaches PI/2.
	 *  @param xa number from which cosine is requested
	 *  @param xb extra bits for x (may be 0.0)
	 *  @return cos(xa + xb)
	 */
	private static double cosQ(const double& xa, const double& xb)
	{
		constexpr double pi2a{ 1.5707963267948966 };
		constexpr double pi2b{ 6.123233995736766E-17 };

		const double a = pi2a - xa;
		double b = -(a - pi2a + xa);
		b += pi2b - xb;

		return sinQ(a, b);
	}

	/**
	 *  Compute tangent (or cotangent) over the first quadrant.   0 < x < pi/2
	 *  Use combination of table lookup and rational polynomial expansion.
	 *  @param xa number from which sine is requested
	 *  @param xb extra bits for x (may be 0.0)
	 *  @param cotanFlag if true, compute the cotangent instead of the tangent
	 *  @return tan(xa+xb) (or cotangent, depending on cotanFlag)
	 */
	private static double tanQ(const double& xa, const double& xb, bool cotanFlag)
	{
		int idx = static_cast<int>(((xa * 8.0) + 0.5);
		const double epsilon = xa - EIGHTHS[idx]; //idx*0.125;

		// Table lookups
		const double sintA = SINE_TABLE_A[idx];
		const double sintB = SINE_TABLE_B[idx];
		const double costA = COSINE_TABLE_A[idx];
		const double costB = COSINE_TABLE_B[idx];

		// Polynomial eval of sin(epsilon), cos(epsilon)
		double sinEpsA = epsilon;
		double sin_eps_b = polySine(epsilon);
		const double cos_eps_a = 1.0;
		const double cos_eps_b = polyCosine(epsilon);

		// Split epsilon   xa + xb = x
		double temp = sinEpsA * HEX_40000000;
		double temp2 = (sinEpsA + temp) - temp;
		sin_eps_b += sinEpsA - temp2;
		sinEpsA = temp2;

		/* Compute sin(x) by angle addition formula */

		/* Compute the following sum:
		 *
		 * result = sintA + costA*sinEpsA + sintA*cos_eps_b + costA*sin_eps_b +
		 *          sintB + costB*sinEpsA + sintB*cos_eps_b + costB*sin_eps_b;
		 *
		 * Ranges of elements
		 *
		 * xxxtA   0            PI/2
		 * xxxtB   -1.5e-9      1.5e-9
		 * sinEpsA -0.0625      0.0625
		 * sin_eps_b -6e-11       6e-11
		 * cos_eps_a  1.0
		 * cos_eps_b  0           -0.0625
		 *
		 */

		 //result = sintA + costA*sinEpsA + sintA*cos_eps_b + costA*sin_eps_b +
		 //          sintB + costB*sinEpsA + sintB*cos_eps_b + costB*sin_eps_b;

		 //result = sintA + sintA*cos_eps_b + sintB + sintB * cos_eps_b;
		 //result += costA*sinEpsA + costA*sin_eps_b + costB*sinEpsA + costB * sin_eps_b;
		double a = 0;
		double b = 0;

		// Compute sine
		double t = sintA;
		double c = a + t;
		double d = -(c - a - t);
		a = c;
		b += d;

		t = costA * sinEpsA;
		c = a + t;
		d = -(c - a - t);
		a = c;
		b += d;

		b += sintA * cos_eps_b + costA * sin_eps_b;
		b += sintB + costB * sinEpsA + sintB * cos_eps_b + costB * sin_eps_b;

		double sina = a + b;
		double sinb = -(sina - a - b);

		// Compute cosine

		a = 0.0;
		b = 0.0;

		t = costA * cos_eps_a;
		c = a + t;
		d = -(c - a - t);
		a = c;
		b += d;

		t = -sintA * sinEpsA;
		c = a + t;
		d = -(c - a - t);
		a = c;
		b += d;

		b += costB * cos_eps_a + costA * cos_eps_b + costB * cos_eps_b;
		b -= sintB * sinEpsA + sintA * sin_eps_b + sintB * sin_eps_b;

		double cosa = a + b;
		double cosb = -(cosa - a - b);

		if (cotanFlag)
		{
			double tmp;
			tmp = cosa; cosa = sina; sina = tmp;
			tmp = cosb; cosb = sinb; sinb = tmp;
		}

		/* estimate and correct, compute 1.0/(cosa+cosb) */
		/*
	double est = (sina+sinb)/(cosa+cosb);
	double err = (sina - cosa*est) + (sinb - cosb*est);
	est += err/(cosa+cosb);
	err = (sina - cosa*est) + (sinb - cosb*est);
		 */

		 // f(x) = 1/x,   f'(x) = -1/x^2

		double est = sina / cosa;

		/* Split the estimate to get more accurate read on division rounding */
		temp = est * HEX_40000000;
		double esta = (est + temp) - temp;
		double estb = est - esta;

		temp = cosa * HEX_40000000;
		double cosaa = (cosa + temp) - temp;
		double cosab = cosa - cosaa;

		//double err = (sina - est*cosa)/cosa;  // Correction for division rounding
		double err = (sina - esta * cosaa - esta * cosab - estb * cosaa - estb * cosab) / cosa;  // Correction for division rounding
		err += sinb / cosa;                     // Change in est due to sinb
		err += -sina * cosb / cosa / cosa;    // Change in est due to cosb

		if (xb != 0.0)
		{
			// tan' = 1 + tan^2      cot' = -(1 + cot^2)
			// Approximate impact of xb
			double xbadj = xb + est * est * xb;
			if (cotanFlag)
			{
				xbadj = -xbadj;
			}

			err += xbadj;
		}

		return est + err;
	}

	/** Reduce the input argument using the Payne and Hanek method.
	 *  This is good for all inputs 0.0 < x < inf
	 *  Output is remainder after dividing by PI/2
	 *  The result array should contain 3 numbers.
	 *  result[0] is the integer portion, so mod 4 this gives the quadrant.
	 *  result[1] is the upper bits of the remainder
	 *  result[2] is the lower bits of the remainder
	 *
	 * @param x number to reduce
	 * @param result placeholder where to put the result
	 */
	private static void reducePayneHanek(const double& x, const std::vector<double>& result)
	{
		/* Convert input double to bits */
		long inbits = Double.double_to_raw_long_bits(x);
		int exponent = static_cast<int>(((inbits >> 52) & 0x7ff) - 1023;

		/* Convert to fixed point representation */
		inbits &= 0x000fffffffffffffL;
		inbits |= 0x0010000000000000L;

		/* Normalize input to be between 0.5 and 1.0 */
		exponent++;
		inbits <<= 11;

		/* Based on the exponent, get a shifted copy of recip2pi */
		long shpi0;
		long shpiA;
		long shpiB;
		int idx = exponent >> 6;
		int shift = exponent - (idx << 6);

		if (shift != 0)
		{
			shpi0 = (idx == 0) ? 0 : (RECIP_2PI[idx - 1] << shift);
			shpi0 |= RECIP_2PI[idx] >> > (64 - shift);
			shpiA = (RECIP_2PI[idx] << shift) | (RECIP_2PI[idx + 1] >> > (64 - shift));
			shpiB = (RECIP_2PI[idx + 1] << shift) | (RECIP_2PI[idx + 2] >> > (64 - shift));
		}
		else
		{
			shpi0 = (idx == 0) ? 0 : RECIP_2PI[idx - 1];
			shpiA = RECIP_2PI[idx];
			shpiB = RECIP_2PI[idx + 1];
		}

		/* Multiply input by shpiA */
		long a = inbits >> > 32;
		long b = inbits & 0xffffffffL;

		long c = shpiA >> > 32;
		long d = shpiA & 0xffffffffL;

		long ac = a * c;
		long bd = b * d;
		long bc = b * c;
		long ad = a * d;

		long prodB = bd + (ad << 32);
		long prodA = ac + (ad >> > 32);

		bool bita = (bd & 0x8000000000000000L) != 0;
		bool bitb = (ad & 0x80000000L) != 0;
		bool bitsum = (prodB & 0x8000000000000000L) != 0;

		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prodA++;
		}

		bita = (prodB & 0x8000000000000000L) != 0;
		bitb = (bc & 0x80000000L) != 0;

		prodB += bc << 32;
		prodA += bc >> > 32;

		bitsum = (prodB & 0x8000000000000000L) != 0;

		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prodA++;
		}

		/* Multiply input by shpiB */
		c = shpiB >> > 32;
		d = shpiB & 0xffffffffL;
		ac = a * c;
		bc = b * c;
		ad = a * d;

		/* Collect terms */
		ac += (bc + ad) >> > 32;

		bita = (prodB & 0x8000000000000000L) != 0;
		bitb = (ac & 0x8000000000000000L) != 0;
		prodB += ac;
		bitsum = (prodB & 0x8000000000000000L) != 0;
		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prodA++;
		}

		/* Multiply by shpi0 */
		c = shpi0 >> > 32;
		d = shpi0 & 0xffffffffL;

		bd = b * d;
		bc = b * c;
		ad = a * d;

		prodA += bd + ((bc + ad) << 32);

		/*
		 * prodA, prodB now contain the remainder as a fraction of PI.  We want this as a fraction of
		 * PI/2, so use the following steps:
		 * 1.) multiply by 4.
		 * 2.) do a fixed point muliply by PI/4.
		 * 3.) Convert to floating point.
		 * 4.) Multiply by 2
		 */

		 /* This identifies the quadrant */
		int intPart = static_cast<int>((prodA >> > 62);

		/* Multiply by 4 */
		prodA <<= 2;
		prodA |= prodB >> > 62;
		prodB <<= 2;

		/* Multiply by PI/4 */
		a = prodA >> > 32;
		b = prodA & 0xffffffffL;

		c = PI_O_4_BITS[0] >> > 32;
		d = PI_O_4_BITS[0] & 0xffffffffL;

		ac = a * c;
		bd = b * d;
		bc = b * c;
		ad = a * d;

		long prod2B = bd + (ad << 32);
		long prod2A = ac + (ad >> > 32);

		bita = (bd & 0x8000000000000000L) != 0;
		bitb = (ad & 0x80000000L) != 0;
		bitsum = (prod2B & 0x8000000000000000L) != 0;

		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prod2A++;
		}

		bita = (prod2B & 0x8000000000000000L) != 0;
		bitb = (bc & 0x80000000L) != 0;

		prod2B += bc << 32;
		prod2A += bc >> > 32;

		bitsum = (prod2B & 0x8000000000000000L) != 0;

		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prod2A++;
		}

		/* Multiply input by pio4bits[1] */
		c = PI_O_4_BITS[1] >> > 32;
		d = PI_O_4_BITS[1] & 0xffffffffL;
		ac = a * c;
		bc = b * c;
		ad = a * d;

		/* Collect terms */
		ac += (bc + ad) >> > 32;

		bita = (prod2B & 0x8000000000000000L) != 0;
		bitb = (ac & 0x8000000000000000L) != 0;
		prod2B += ac;
		bitsum = (prod2B & 0x8000000000000000L) != 0;
		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prod2A++;
		}

		/* Multiply inputB by pio4bits[0] */
		a = prodB >> > 32;
		b = prodB & 0xffffffffL;
		c = PI_O_4_BITS[0] >> > 32;
		d = PI_O_4_BITS[0] & 0xffffffffL;
		ac = a * c;
		bc = b * c;
		ad = a * d;

		/* Collect terms */
		ac += (bc + ad) >> > 32;

		bita = (prod2B & 0x8000000000000000L) != 0;
		bitb = (ac & 0x8000000000000000L) != 0;
		prod2B += ac;
		bitsum = (prod2B & 0x8000000000000000L) != 0;
		/* Carry */
		if ((bita && bitb) ||
			((bita || bitb) && !bitsum))
		{
			prod2A++;
		}

		/* Convert to double */
		double tmpA = (prod2A >> > 12) / TWO_POWER_52;  // High order 52 bits
		double tmpB = (((prod2A & 0xfffL) << 40) + (prod2B >> > 24)) / TWO_POWER_52 / TWO_POWER_52; // Low bits

		double sumA = tmpA + tmpB;
		double sum_b = -(sumA - tmpA - tmpB);

		/* Multiply by PI/2 and return */
		result[0] = intPart;
		result[1] = sumA * 2.0;
		result[2] = sum_b * 2.0;
	}

	/**
	 * Sine function.
	 *
	 * @param x Argument.
	 * @return sin(x)
	 */
	public static double sin(const double& x)
	{
		bool negative = false;
		int quadrant = 0;
		double xa;
		double xb = 0.0;

		/* Take absolute value of the input */
		xa = x;
		if (x < 0)
		{
			negative = true;
			xa = -xa;
		}

		/* Check for zero and negative zero */
		if (xa == 0.0)
		{
			long bits = Double.double_to_raw_long_bits(x);
			if (bits < 0)
			{
				return -0.0;
			}
			return 0.0;
		}

		if (xa != xa || xa == INFINITY)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		/* Perform any argument reduction */
		if (xa > 3294198.0)
		{
			// PI * (2**20)
			// Argument too big for CodyWaite reduction.  Must use
			// PayneHanek.
			auto reduceResults = std::vector<double>(3];
			reducePayneHanek(xa, reduceResults);
			quadrant = (static_cast<int>(reduceResults[0]) & 3;
			xa = reduceResults[1];
			xb = reduceResults[2];
		}
		else if (xa > 1.5707963267948966)
		{
			const CodyWaite cw = CodyWaite(xa);
			quadrant = cw.get_k() & 3;
			xa = cw.get_remA();
			xb = cw.get_remB();
		}

		if (negative)
		{
			quadrant ^= 2;  // Flip bit 1
		}

		switch (quadrant)
		{
		case 0:
			return sinQ(xa, xb);
		case 1:
			return cosQ(xa, xb);
		case 2:
			return -sinQ(xa, xb);
		case 3:
			return -cosQ(xa, xb);
		default:
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	/**
	 * Cosine function.
	 *
	 * @param x Argument.
	 * @return cos(x)
	 */
	public static double cos(const double& x)
	{
		int quadrant = 0;

		/* Take absolute value of the input */
		double xa = x;
		if (x < 0)
		{
			xa = -xa;
		}

		if (xa != xa || xa == INFINITY)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		/* Perform any argument reduction */
		double xb = 0;
		if (xa > 3294198.0)
		{
			// PI * (2**20)
			// Argument too big for CodyWaite reduction.  Must use
			// PayneHanek.
			auto reduceResults = std::vector<double>(3];
			reducePayneHanek(xa, reduceResults);
			quadrant = (static_cast<int>(reduceResults[0]) & 3;
			xa = reduceResults[1];
			xb = reduceResults[2];
		}
		else if (xa > 1.5707963267948966)
		{
			const CodyWaite cw = CodyWaite(xa);
			quadrant = cw.get_k() & 3;
			xa = cw.get_remA();
			xb = cw.get_remB();
		}

		//if (negative)
		//  quadrant = (quadrant + 2) % 4;

		switch (quadrant)
		{
		case 0:
			return cosQ(xa, xb);
		case 1:
			return -sinQ(xa, xb);
		case 2:
			return -cosQ(xa, xb);
		case 3:
			return sinQ(xa, xb);
		default:
			return std::numeric_limits<double>::quiet_NaN();
		}
	}

	/**
	 * Combined Sine and Cosine function.
	 *
	 * @param x Argument.
	 * @return [sin(x), cos(x)]
	 */
	public static Sin_Cos sin_cos(const double& x)
	{
		bool negative = false;
		int quadrant = 0;
		double xa;
		double xb = 0.0;

		/* Take absolute value of the input */
		xa = x;
		if (x < 0)
		{
			negative = true;
			xa = -xa;
		}

		/* Check for zero and negative zero */
		if (xa == 0.0)
		{
			long bits = Double.double_to_raw_long_bits(x);
			if (bits < 0)
			{
				return Sin_Cos(-0.0, 1.0);
			}
			return Sin_Cos(0.0, 1.0);
		}

		if (xa != xa || xa == INFINITY)
		{
			return Sin_Cos(Double.NaN, NAN);
		}

		/* Perform any argument reduction */
		if (xa > 3294198.0)
		{
			// PI * (2**20)
			// Argument too big for CodyWaite reduction.  Must use
			// PayneHanek.
			auto reduceResults = std::vector<double>(3);
			reducePayneHanek(xa, reduceResults);
			quadrant = (static_cast<int>(reduceResults[0]) & 3;
			xa = reduceResults[1];
			xb = reduceResults[2];
		}
		else if (xa > 1.5707963267948966)
		{
			const CodyWaite cw = CodyWaite(xa);
			quadrant = cw.get_k() & 3;
			xa = cw.get_remA();
			xb = cw.get_remB();
		}

		switch (quadrant)
		{
		case 0:
			return Sin_Cos(negative ? -sinQ(xa, xb) : sinQ(xa, xb), cosQ(xa, xb));
		case 1:
			return Sin_Cos(negative ? -cosQ(xa, xb) : cosQ(xa, xb), -sinQ(xa, xb));
		case 2:
			return Sin_Cos(negative ? sinQ(xa, xb) : -sinQ(xa, xb), -cosQ(xa, xb));
		case 3:
			return Sin_Cos(negative ? cosQ(xa, xb) : -cosQ(xa, xb), sinQ(xa, xb));
		default:
			return Sin_Cos(Double.NaN, NAN);
		}
	}

	/**
	 * Combined Sine and Cosine function.
	 *
	 * @param x Argument.
	 * @param <T> the type of the field element
	 * @return [sin(x), cos(x)]
	 * @since 1.4
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static Field_Sin_Cos<T> sin_cos(T x)
	{
		return x.sin_cos();
	}

	/**
	 * Tangent function.
	 *
	 * @param x Argument.
	 * @return tan(x)
	 */
	public static double tan(const double& x)
	{
		bool negative = false;
		int quadrant = 0;

		/* Take absolute value of the input */
		double xa = x;
		if (x < 0)
		{
			negative = true;
			xa = -xa;
		}

		/* Check for zero and negative zero */
		if (xa == 0.0)
		{
			long bits = Double.double_to_raw_long_bits(x);
			if (bits < 0)
			{
				return -0.0;
			}
			return 0.0;
		}

		if (xa != xa || xa == INFINITY)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		/* Perform any argument reduction */
		double xb = 0;
		if (xa > 3294198.0)
		{
			// PI * (2**20)
			// Argument too big for CodyWaite reduction.  Must use
			// PayneHanek.
			auto reduceResults = std::vector<double>(3);
			reducePayneHanek(xa, reduceResults);
			quadrant = (static_cast<int>(reduceResults[0]) & 3;
			xa = reduceResults[1];
			xb = reduceResults[2];
		}
		else if (xa > 1.5707963267948966)
		{
			const CodyWaite cw = CodyWaite(xa);
			quadrant = cw.get_k() & 3;
			xa = cw.get_remA();
			xb = cw.get_remB();
		}

		if (xa > 1.5)
		{
			// Accuracy suffers between 1.5 and PI/2
			const double pi2a = 1.5707963267948966;
			const double pi2b = 6.123233995736766E-17;

			const double& a = pi2a - xa;
			double b = -(a - pi2a + xa);
			b += pi2b - xb;

			xa = a + b;
			xb = -(xa - a - b);
			quadrant ^= 1;
			negative ^= true;
		}

		double result;
		if ((quadrant & 1) == 0)
		{
			result = tanQ(xa, xb, false);
		}
		else
		{
			result = -tanQ(xa, xb, true);
		}

		if (negative)
		{
			result = -result;
		}

		return result;
	}

	/**
	 * Arctangent function
	 *  @param x a number
	 *  @return atan(x)
	 */
	public static double atan(const double& x)
	{
		return atan(x, 0.0, false);
	}

	/** Internal helper function to compute arctangent.
	 * @param xa number from which arctangent is requested
	 * @param xb extra bits for x (may be 0.0)
	 * @param left_plane if true, result angle must be put in the left half plane
	 * @return atan(xa + xb) (or angle shifted by {@code PI} if left_plane is true)
	 */
	private static double atan(const double& xa, const double& xb, bool left_plane)
	{
		if (xa == 0.0) { // Matches +/- 0.0; return correct sign
			return left_plane ? copy_sign(std::numbers::pi, xa) : xa;
		}

		const bool negate;
		if (xa < 0)
		{
			// negative
			xa = -xa;
			xb = -xb;
			negate = true;
		}
		else
		{
			negate = false;
		}

		if (xa > 1.633123935319537E16) { // Very large input
			return (negate ^ left_plane) ? (-std::numbers::pi * F_1_2) : (std::numbers::pi * F_1_2);
		}

		/* Estimate the closest tabulated arctan value, compute eps = xa-tangentTable */
		const int idx;
		if (xa < 1)
		{
			idx = static_cast<int>((((-1.7168146928204136 * xa * xa + 8.0) * xa) + 0.5);
		}
		else
		{
			const double one_over_xa = 1 / xa;
			idx = static_cast<int>((-((-1.7168146928204136 * one_over_xa * one_over_xa + 8.0) * one_over_xa) + 13.07);
		}

		const double ttA = TANGENT_TABLE_A[idx];
		const double ttB = TANGENT_TABLE_B[idx];

		double epsA = xa - ttA;
		double epsB = -(epsA - xa + ttA);
		epsB += xb - ttB;

		double temp = epsA + epsB;
		epsB = -(temp - epsA - epsB);
		epsA = temp;

		/* Compute eps = eps / (1.0 + xa*tangent) */
		temp = xa * HEX_40000000;
		double ya = xa + temp - temp;
		double yb = xb + xa - ya;
		xa = ya;
		xb += yb;

		//if (idx > 8 || idx == 0)
		if (idx == 0)
		{
			/* If the slope of the arctan is gentle enough (< 0.45), this approximation will suffice */
			//double denom = 1.0 / (1.0 + xa*tangent_Table_A[idx] + xb*tangent_Table_A[idx] + xa*tangentTableB[idx] + xb*tangentTableB[idx]);
			const double denom = 1.0 / (1.0 + (xa + xb) * (ttA + ttB));
			//double denom = 1.0 / (1.0 + xa*tangent_Table_A[idx]);
			ya = epsA * denom;
			yb = epsB * denom;
		}
		else
		{
			double temp2 = xa * ttA;
			double za = 1.0 + temp2;
			double zb = -(za - 1.0 - temp2);
			temp2 = xb * ttA + xa * ttB;
			temp = za + temp2;
			zb += -(temp - za - temp2);
			za = temp;

			zb += xb * ttB;
			ya = epsA / za;

			temp = ya * HEX_40000000;
			const double yaa = (ya + temp) - temp;
			const double yab = ya - yaa;

			temp = za * HEX_40000000;
			const double zaa = (za + temp) - temp;
			const double zab = za - zaa;

			/* Correct for rounding in division */
			yb = (epsA - yaa * zaa - yaa * zab - yab * zaa - yab * zab) / za;

			yb += -epsA * zb / za / za;
			yb += epsB / za;
		}

		epsA = ya;
		epsB = yb;

		/* Evaluate polynomial */
		const double epsA2 = epsA * epsA;

		/*
	yb = -0.09001346640161823;
	yb = yb * epsA2 + 0.11110718400605211;
	yb = yb * epsA2 + -0.1428571349122913;
	yb = yb * epsA2 + 0.19999999999273194;
	yb = yb * epsA2 + -0.33333333333333093;
	yb = yb * epsA2 * epsA;
		 */

		yb = 0.07490822288864472;
		yb = yb * epsA2 - 0.09088450866185192;
		yb = yb * epsA2 + 0.11111095942313305;
		yb = yb * epsA2 - 0.1428571423679182;
		yb = yb * epsA2 + 0.19999999999923582;
		yb = yb * epsA2 - 0.33333333333333287;
		yb = yb * epsA2 * epsA;

		ya = epsA;

		temp = ya + yb;
		yb = -(temp - ya - yb);
		ya = temp;

		/* Add in effect of epsB.   atan'(x) = 1/(1+x^2) */
		yb += epsB / (1.0 + epsA * epsA);

		const double eighths = EIGHTHS[idx];

		//result = yb + eighths[idx] + ya;
		double za = eighths + ya;
		double zb = -(za - eighths - ya);
		temp = za + yb;
		zb += -(temp - za - yb);
		za = temp;

		double result = za + zb;

		if (left_plane)
		{
			// Result is in the left plane
			const double resultb = -(result - za - zb);
			const double pia = 1.5707963267948966 * 2;
			const double pib = 6.123233995736766E-17 * 2;

			za = pia - result;
			zb = -(za - pia + result);
			zb += pib - resultb;

			result = za + zb;
		}

		if (negate ^ left_plane)
		{
			result = -result;
		}

		return result;
	}

	/**
	 * Two arguments arctangent function
	 * @param y ordinate
	 * @param x abscissa
	 * @return phase angle of point (x,y) between {@code -PI} and {@code PI}
	 */
	public static double atan2(const double& y, const double& x)
	{
		if (std::isnan(x) || std::isnan(y))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (y == 0)
		{
			const double result = x * y;
			const double invx = 1.0 / x;

			if (invx == 0) { // X is infinite
				if (x > 0)
				{
					return y; // return +/- 0.0
				}
				else
				{
					return copy_sign(std::numbers::pi, y);
				}
			}

			if (x < 0 || invx < 0)
			{
				return copy_sign(std::numbers::pi, y);
			}
			else
			{
				return result;
			}
		}

		// y cannot now be zero

		if (y == INFINITY)
		{
			if (x == INFINITY)
			{
				return std::numbers::pi * F_1_4;
			}

			if (x == -INFINITY)
			{
				return std::numbers::pi * F_3_4;
			}

			return std::numbers::pi * F_1_2;
		}

		if (y == -INFINITY)
		{
			if (x == INFINITY)
			{
				return -std::numbers::pi * F_1_4;
			}

			if (x == -INFINITY)
			{
				return -std::numbers::pi * F_3_4;
			}

			return -std::numbers::pi * F_1_2;
		}

		if (x == INFINITY)
		{
			return copy_sign(0d, y);
		}

		if (x == -INFINITY)

		{
			return copy_sign(std::numbers::pi, y);
		}

		// Neither y nor x can be infinite or NAN here

		if (x == 0)
		{
			return copy_sign(std::numbers::pi * F_1_2, y);
		}

		// Compute ratio r = y/x
		const double r = y / x;
		if (std::isinf(r)) { // bypass calculations that can create NaN
			return atan(r, 0, x < 0);
		}

		double ra = doubleHigh_part(r);
		double rb = r - ra;

		// Split x
		const double xa = doubleHigh_part(x);
		const double xb = x - xa;

		rb += (y - ra * xa - ra * xb - rb * xa - rb * xb) / x;

		const double temp = ra + rb;
		rb = -(temp - ra - rb);
		ra = temp;

		if (ra == 0) { // Fix up the sign so atan works correctly
			ra = copy_sign(0d, y);
		}

		// Call atan
		return atan(ra, rb, x < 0);
	}

	/** Compute the arc sine of a number.
	 * @param x number on which evaluation is done
	 * @return arc sine of x
	 */
	public static double asin(const double& x)
	{
		if (std::isnan(x))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (x > 1.0 || x < -1.0)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (x == 1.0)
		{
			return std::numbers::pi / 2.0;
		}

		if (x == -1.0)
		{
			return -std::numbers::pi / 2.0;
		}

		if (x == 0.0) { // Matches +/- 0.0; return correct sign
			return x;
		}

		/* Compute asin(x) = atan(x/sqrt(1-x*x)) */

		/* Split x */
		double temp = x * HEX_40000000;
		const double xa = x + temp - temp;
		const double xb = x - xa;

		/* Square it */
		double ya = xa * xa;
		double yb = xa * xb * 2.0 + xb * xb;

		/* Subtract from 1 */
		ya = -ya;
		yb = -yb;

		double za = 1.0 + ya;
		double zb = -(za - 1.0 - ya);

		temp = za + yb;
		zb += -(temp - za - yb);
		za = temp;

		/* Square root */
		double y;
		y = sqrt(za);
		temp = y * HEX_40000000;
		ya = y + temp - temp;
		yb = y - ya;

		/* Extend precision of sqrt */
		yb += (za - ya * ya - 2 * ya * yb - yb * yb) / (2.0 * y);

		/* Contribution of zb to sqrt */
		double dx = zb / (2.0 * y);

		// Compute ratio r = x/y
		double r = x / y;
		temp = r * HEX_40000000;
		double ra = r + temp - temp;
		double rb = r - ra;

		rb += (x - ra * ya - ra * yb - rb * ya - rb * yb) / y;  // Correct for rounding in division
		rb += -x * dx / y / y;  // Add in effect additional bits of sqrt.

		temp = ra + rb;
		rb = -(temp - ra - rb);
		ra = temp;

		return atan(ra, rb, false);
	}

	/** Compute the arc cosine of a number.
	 * @param x number on which evaluation is done
	 * @return arc cosine of x
	 */
	public static double acos(const double& x)
	{
		if (std::isnan(x))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (x > 1.0 || x < -1.0)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (x == -1.0)
		{
			return std::numbers::pi;
		}

		if (x == 1.0)
		{
			return 0.0;
		}

		if (x == 0)
		{
			return std::numbers::pi / 2.0;
		}

		/* Compute acos(x) = atan(sqrt(1-x*x)/x) */

		/* Split x */
		double temp = x * HEX_40000000;
		const double xa = x + temp - temp;
		const double xb = x - xa;

		/* Square it */
		double ya = xa * xa;
		double yb = xa * xb * 2.0 + xb * xb;

		/* Subtract from 1 */
		ya = -ya;
		yb = -yb;

		double za = 1.0 + ya;
		double zb = -(za - 1.0 - ya);

		temp = za + yb;
		zb += -(temp - za - yb);
		za = temp;

		/* Square root */
		double y = sqrt(za);
		temp = y * HEX_40000000;
		ya = y + temp - temp;
		yb = y - ya;

		/* Extend precision of sqrt */
		yb += (za - ya * ya - 2 * ya * yb - yb * yb) / (2.0 * y);

		/* Contribution of zb to sqrt */
		yb += zb / (2.0 * y);
		y = ya + yb;
		yb = -(y - ya - yb);

		// Compute ratio r = y/x
		double r = y / x;

		// Did r overflow?
		if (std::isinf(r)) { // x is effectively zero
			return std::numbers::pi / 2; // so return the appropriate value
		}

		double ra = doubleHigh_part(r);
		double rb = r - ra;

		rb += (y - ra * xa - ra * xb - rb * xa - rb * xb) / x;  // Correct for rounding in division
		rb += yb / x;  // Add in effect additional bits of sqrt.

		temp = ra + rb;
		rb = -(temp - ra - rb);
		ra = temp;

		return atan(ra, rb, x < 0);
	}

	/** Compute the cubic root of a number.
	 * @param x number on which evaluation is done
	 * @return cubic root of x
	 */
	public static double cbrt(const double& x)
	{
		/* Convert input double to bits */
		long inbits = Double.double_to_raw_long_bits(x);
		int exponent = static_cast<int>(((inbits >> 52) & 0x7ff) - 1023;
		bool subnormal = false;

		if (exponent == -1023)
		{
			if (x == 0)
			{
				return x;
			}

			/* Subnormal, so normalize */
			subnormal = true;
			x *= 1.8014398509481984E16;  // 2^54
			inbits = Double.double_to_raw_long_bits(x);
			exponent = static_cast<int>(((inbits >> 52) & 0x7ff) - 1023;
		}

		if (exponent == 1024)
		{
			// Nan or infinity.  Don't care which.
			return x;
		}

		/* Divide the exponent by 3 */
		int exp3 = exponent / 3;

		/* p2 will be the nearest power of 2 to x with its exponent divided by 3 */
		double p2 = Double.long_bits_to_double((inbits & 0x8000000000000000L) |
			static_cast<long>((((exp3 + 1023) & 0x7ff)) << 52);

		/* This will be a number between 1 and 2 */
		const double mant = Double.long_bits_to_double((inbits & 0x000fffffffffffffL) | 0x3ff0000000000000L);

		/* Estimate the cube root of mant by polynomial */
		double est = -0.010714690733195933;
		est = est * mant + 0.0875862700108075;
		est = est * mant + -0.3058015757857271;
		est = est * mant + 0.7249995199969751;
		est = est * mant + 0.5039018405998233;

		est *= CBRTTWO[exponent % 3 + 2];

		// est should now be good to about 15 bits of precision.   Do 2 rounds of
		// Newton's method to get closer,  this should get us full double precision
		// Scale down x for the purpose of doing newtons method.  This avoids over/under flows.
		const double xs = x / (p2 * p2 * p2);
		est += (xs - est * est * est) / (3 * est * est);
		est += (xs - est * est * est) / (3 * est * est);

		// Do one round of Newton's method in extended precision to get the last bit right.
		double temp = est * HEX_40000000;
		double ya = est + temp - temp;
		double yb = est - ya;

		double za = ya * ya;
		double zb = ya * yb * 2.0 + yb * yb;
		temp = za * HEX_40000000;
		double temp2 = za + temp - temp;
		zb += za - temp2;
		za = temp2;

		zb = za * yb + ya * zb + zb * yb;
		za *= ya;

		double na = xs - za;
		double nb = -(na - xs + za);
		nb -= zb;

		est += (na + nb) / (3 * est * est);

		/* Scale by a power of two, so this is exact. */
		est *= p2;

		if (subnormal)
		{
			est *= 3.814697265625E-6;  // 2^-18
		}

		return est;
	}

	/**
	 *  Convert degrees to radians, with error of less than 0.5 ULP
	 *  @param x angle in degrees
	 *  @return x converted into radians
	 */
	public static double to_radians(const double& x)

	{
		if (std::isinf(x) || x == 0.0) { // Matches +/- 0.0; return correct sign
			return x;
		}

		// These are PI/180 split into high and low order bits
		const double facta = 0.01745329052209854;
		const double factb = 1.997844754509471E-9;

		double xa = doubleHigh_part(x);
		double xb = x - xa;

		double result = xb * factb + xb * facta + xa * factb + xa * facta;
		if (result == 0)
		{
			result *= x; // ensure correct sign if calculation underflows
		}
		return result;
	}

	/**
	 *  Convert radians to degrees, with error of less than 0.5 ULP
	 *  @param x angle in radians
	 *  @return x converted into degrees
	 */
	public static double to_degrees(const double& x)

	{
		if (std::isinf(x) || x == 0.0) { // Matches +/- 0.0; return correct sign
			return x;
		}

		// These are 180/PI split into high and low order bits
		const double facta = 57.2957763671875;
		const double factb = 3.145894820876798E-6;

		double xa = doubleHigh_part(x);
		double xb = x - xa;

		return xb * factb + xb * facta + xa * factb + xa * facta;
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @return abs(x)
	 */
	public static int abs(const int& x)
	{
		const int i = x >> > 31;
		return (x ^ (~i + 1)) + i;
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @return abs(x)
	 */
	public static long abs(const long& x)
	{
		const long l = x >> > 63;
		// l is one if x negative zero else
		// ~l+1 is zero if x is positive, -1 if x is negative
		// x^(~l+1) is x is x is positive, ~x if x is negative
		// add around
		return (x ^ (~l + 1)) + l;
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @return abs(x), or an exception for {@code std::numeric_limits<int>::min()}
	 */
	public static int absExact(const int& x)
	{
		if (x == std::numeric_limits<int>::min())
		{
			throw Arithmetic_Exception();
		}
		return abs(x);
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @return abs(x), or an exception for {@code long.MIN_VALUE}
	 * @since 2.0
	 */
	public static long absExact(const long& x)
	{
		if (x == long.MIN_VALUE)
		{
			throw Arithmetic_Exception();
		}
		return abs(x);
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @return abs(x)
	 * @since 2.0
	 */
	public static float abs(const float& x)
	{
		return Float.int_bits_to_float(MASK_NON_SIGN_INT & Float.float_to_raw_int_bits(x));
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @return abs(x)
	 */
	public static double abs(const double& x)
	{
		return Double.long_bits_to_double(MASK_NON_SIGN_LONG & Double.double_to_raw_long_bits(x));
	}

	/**
	 * Negates the argument.
	 * @param x number from which opposite value is requested
	 * @return -x, or an exception for {@code std::numeric_limits<int>::min()}
	 * @since 2.0
	 */
	public static int negateExact(const int& x)
	{
		if (x == std::numeric_limits<int>::min())
		{
			throw Arithmetic_Exception();
		}
		return -x;
	}

	/**
	 * Negates the argument.
	 * @param x number from which opposite value is requested
	 * @return -x, or an exception for {@code long.MIN_VALUE}
	 * @since 2.0
	 */
	public static long negateExact(const long& x)
	{
		if (x == long.MIN_VALUE)
		{
			throw Arithmetic_Exception();
		}
		return -x;
	}

	/**
	 * Compute least significant bit (Unit in Last Position) for a number.
	 * @param x number from which ulp is requested
	 * @return ulp(x)
	 */
	public static double ulp(const double& x)
	{
		if (std::isinf(x))
		{
			return INFINITY;
		}
		return abs(x - Double.long_bits_to_double(Double.double_to_raw_long_bits(x) ^ 1));
	}

	/**
	 * Compute least significant bit (Unit in Last Position) for a number.
	 * @param x number from which ulp is requested
	 * @return ulp(x)
	 */
	public static float ulp(float x)
	{
		if (Float.std::isinfinite(x))
		{
			return Float.POSITIVE_INFINITY;
		}
		return abs(x - Float.int_bits_to_float(Float.float_to_raw_int_bits(x) ^ 1));
	}

	/**
	 * Multiply a double number by a power of 2.
	 * @param d number to multiply
	 * @param n power of 2
	 * @return d &times; 2<sup>n</sup>
	 */
	public static double scalb(const double d, const int& n)
	{
		// first simple and fast handling when 2^n can be represented using normal numbers
		if ((n > -1023) && (n < 1024))
		{
			return d * Double.long_bits_to_double((static_cast<long>((n + 1023)) << 52);
		}

		// handle special cases
		if (std::isnan(d) || std::isinf(d) || (d == 0))
		{
			return d;
		}
		if (n < -2098)
		{
			return (d > 0) ? 0.0 : -0.0;
		}
		if (n > 2097)
		{
			return (d > 0) ? INFINITY : -INFINITY;
		}

		// decompose d
		const long bits = Double.double_to_raw_long_bits(d);
		const long sign = bits & 0x8000000000000000L;
		int  exponent = (static_cast<int>((bits >> > 52)) & 0x7ff;
		long mantissa = bits & 0x000fffffffffffffL;

		// compute scaled exponent
		int scaledExponent = exponent + n;

		if (n < 0)
		{
			// we are really in the case n <= -1023
			if (scaledExponent > 0)
			{
				// both the input and the result are normal numbers, we only adjust the exponent
				return Double.long_bits_to_double(sign | ((static_cast<long>(scaledExponent) << 52) | mantissa);
			}
			else if (scaledExponent > -53)
			{
				// the input is a normal number and the result is a subnormal number

				// recover the hidden mantissa bit
				mantissa |= 1L << 52;

				// scales down complete mantissa, hence losing least significant bits
				const long most_significant_lost_bit = mantissa & (1L << (-scaledExponent));
				mantissa >>>= 1 - scaledExponent;
				if (most_significant_lost_bit != 0)
				{
					// we need to add 1 bit to round up the result
					mantissa++;
				}
				return Double.long_bits_to_double(sign | mantissa);
			}
			else
			{
				// no need to compute the mantissa, the number scales down to 0
				return (sign == 0L) ? 0.0 : -0.0;
			}
		}
		else
		{
			// we are really in the case n >= 1024
			if (exponent == 0)
			{
				// the input number is subnormal, normalize it
				while ((mantissa >> > 52) != 1)
				{
					mantissa <<= 1;
					--scaledExponent;
				}
				++scaledExponent;
				mantissa &= 0x000fffffffffffffL;

				if (scaledExponent < 2047)
				{
					return Double.long_bits_to_double(sign | ((static_cast<long>(scaledExponent) << 52) | mantissa);
				}
				else
				{
					return (sign == 0L) ? INFINITY : -INFINITY;
				}
			}
			else if (scaledExponent < 2047)
			{
				return Double.long_bits_to_double(sign | ((static_cast<long>(scaledExponent) << 52) | mantissa);
			}
			else
			{
				return (sign == 0L) ? INFINITY : -INFINITY;
			}
		}
	}

	/**
	 * Multiply a float number by a power of 2.
	 * @param f number to multiply
	 * @param n power of 2
	 * @return f &times; 2<sup>n</sup>
	 */
	public static float scalb(const float f, const int& n)
	{
		// first simple and fast handling when 2^n can be represented using normal numbers
		if ((n > -127) && (n < 128))
		{
			return f * Float.int_bits_to_float((n + 127) << 23);
		}

		// handle special cases
		if (Float.is_nan(f) || Float.std::isinfinite(f) || (f == 0f))
		{
			return f;
		}
		if (n < -277)
		{
			return (f > 0) ? 0.0f : -0.0f;
		}
		if (n > 276)
		{
			return (f > 0) ? Float.POSITIVE_INFINITY : Float.NEGATIVE_INFINITY;
		}

		// decompose f
		const int bits = Float.floatToIntBits(f);
		const int sign = bits & 0x80000000;
		int  exponent = (bits >> > 23) & 0xff;
		int mantissa = bits & 0x007fffff;

		// compute scaled exponent
		int scaledExponent = exponent + n;

		if (n < 0)
		{
			// we are really in the case n <= -127
			if (scaledExponent > 0)
			{
				// both the input and the result are normal numbers, we only adjust the exponent
				return Float.int_bits_to_float(sign | (scaledExponent << 23) | mantissa);
			}
			else if (scaledExponent > -24)
			{
				// the input is a normal number and the result is a subnormal number

				// recover the hidden mantissa bit
				mantissa |= 1 << 23;

				// scales down complete mantissa, hence losing least significant bits
				const int most_significant_lost_bit = mantissa & (1 << (-scaledExponent));
				mantissa >>>= 1 - scaledExponent;
				if (most_significant_lost_bit != 0)
				{
					// we need to add 1 bit to round up the result
					mantissa++;
				}
				return Float.int_bits_to_float(sign | mantissa);
			}
			else
			{
				// no need to compute the mantissa, the number scales down to 0
				return (sign == 0) ? 0.0f : -0.0f;
			}
		}
		else
		{
			// we are really in the case n >= 128
			if (exponent == 0)
			{
				// the input number is subnormal, normalize it
				while ((mantissa >> > 23) != 1)
				{
					mantissa <<= 1;
					--scaledExponent;
				}
				++scaledExponent;
				mantissa &= 0x007fffff;

				if (scaledExponent < 255)
				{
					return Float.int_bits_to_float(sign | (scaledExponent << 23) | mantissa);
				}
				else
				{
					return (sign == 0) ? Float.POSITIVE_INFINITY : Float.NEGATIVE_INFINITY;
				}
			}
			else if (scaledExponent < 255)
			{
				return Float.int_bits_to_float(sign | (scaledExponent << 23) | mantissa);
			}
			else
			{
				return (sign == 0) ? Float.POSITIVE_INFINITY : Float.NEGATIVE_INFINITY;
			}
		}
	}

	/**
	 * Get the next machine representable number after a number, moving
	 * in the direction of another number.
	 * <p>
	 * The ordering is as follows (increasing):
	 * <ul>
	 * <li>-INFINITY</li>
	 * <li>-MAX_VALUE</li>
	 * <li>-MIN_VALUE</li>
	 * <li>-0.0</li>
	 * <li>+0.0</li>
	 * <li>+MIN_VALUE</li>
	 * <li>+MAX_VALUE</li>
	 * <li>+INFINITY</li>
	 * <li></li>
	 * <p>
	 * If arguments compare equal, then the second argument is returned.
	 * <p>
	 * If {@code direction} is greater than {@code d}, * the smallest machine representable number strictly greater than
	 * {@code d} is returned; if less, then the largest representable number
	 * strictly less than {@code d} is returned.</p>
	 * <p>
	 * If {@code d} is infinite and direction does not
	 * bring it back to finite numbers, it is returned unchanged.</p>
	 *
	 * @param d base number
	 * @param direction (the only important thing is whether
	 * {@code direction} is greater or smaller than {@code d})
	 * @return the next machine representable number in the specified direction
	 */
	public static double next_after(const double& d, double direction)
	{
		// handling of some important special cases
		if (std::isnan(d) || std::isnan(direction))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		else if (d == direction)
		{
			return direction;
		}
		else if (std::isinf(d))
		{
			return (d < 0) ? -Double.MAX_VALUE : Double.MAX_VALUE;
		}
		else if (d == 0)
		{
			return (direction < 0) ? -Double.MIN_VALUE : Double.MIN_VALUE;
		}
		// special cases MAX_VALUE to infinity and  MIN_VALUE to 0
		// are handled just as normal numbers
		// can use raw bits since already dealt with infinity and NaN
		const long bits = Double.double_to_raw_long_bits(d);
		const long sign = bits & 0x8000000000000000L;
		if ((direction < d) ^ (sign == 0L))
		{
			return Double.long_bits_to_double(sign | ((bits & 0x7fffffffffffffffL) + 1));
		}
		else
		{
			return Double.long_bits_to_double(sign | ((bits & 0x7fffffffffffffffL) - 1));
		}
	}

	/**
	 * Get the next machine representable number after a number, moving
	 * in the direction of another number.
	 * <p>
	 * The ordering is as follows (increasing):
	 * <ul>
	 * <li>-INFINITY</li>
	 * <li>-MAX_VALUE</li>
	 * <li>-MIN_VALUE</li>
	 * <li>-0.0</li>
	 * <li>+0.0</li>
	 * <li>+MIN_VALUE</li>
	 * <li>+MAX_VALUE</li>
	 * <li>+INFINITY</li>
	 * <li></li>
	 * <p>
	 * If arguments compare equal, then the second argument is returned.
	 * <p>
	 * If {@code direction} is greater than {@code f}, * the smallest machine representable number strictly greater than
	 * {@code f} is returned; if less, then the largest representable number
	 * strictly less than {@code f} is returned.</p>
	 * <p>
	 * If {@code f} is infinite and direction does not
	 * bring it back to finite numbers, it is returned unchanged.</p>
	 *
	 * @param f base number
	 * @param direction (the only important thing is whether
	 * {@code direction} is greater or smaller than {@code f})
	 * @return the next machine representable number in the specified direction
	 */
	public static float next_after(const float f, const double direction)
	{
		// handling of some important special cases
		if (std::isnan(f) || std::isnan(direction))
		{
			return Float.NaN;
		}
		else if (f == direction)
		{
			return (float)direction;
		}
		else if (Float.std::isinfinite(f))
		{
			return (f < 0f) ? -Float.MAX_VALUE : Float.MAX_VALUE;
		}
		else if (f == 0f)
		{
			return (direction < 0) ? -Float.MIN_VALUE : Float.MIN_VALUE;
		}
		// special cases MAX_VALUE to infinity and  MIN_VALUE to 0
		// are handled just as normal numbers

		const int bits = Float.floatToIntBits(f);
		const int sign = bits & 0x80000000;
		if ((direction < f) ^ (sign == 0))
		{
			return Float.int_bits_to_float(sign | ((bits & 0x7fffffff) + 1));
		}
		else
		{
			return Float.int_bits_to_float(sign | ((bits & 0x7fffffff) - 1));
		}
	}

	/** Get the largest whole number smaller than x.
	 * @param x number from which floor is requested
	 * @return a double number f such that f is an integer f &lt;= x &lt; f + 1.0
	 */
	public static double floor(const double& x)
	{
		long y;

		if (std::isnan(x))
		{
			return x;
		}

		if (x >= TWO_POWER_52 || x <= -TWO_POWER_52)
		{
			return x;
		}

		y = static_cast<long>(x;
		if (x < 0 && y != x)
		{
			y--;
		}

		if (y == 0)
		{
			return x * y;
		}

		return y;
	}

	/** Get the smallest whole number larger than x.
	 * @param x number from which ceil is requested
	 * @return a double number c such that c is an integer c - 1.0 &lt; x &lt;= c
	 */
	public static double ceil(const double& x)
	{
		double y;

		if (std::isnan(x))
		{
			return x;
		}

		y = floor(x);
		if (y == x)
		{
			return y;
		}

		y += 1.0;

		if (y == 0)
		{
			return x * y;
		}

		return y;
	}

	/** Get the whole number that is the nearest to x, or the even one if x is exactly half way between two integers.
	 * @param x number from which nearest whole number is requested
	 * @return a double number r such that r is an integer r - 0.5 &lt;= x &lt;= r + 0.5
	 */
	public static double rint(const double& x)
	{
		double y = floor(x);
		double d = x - y;

		if (d > 0.5)
		{
			if (y == -1.0)
			{
				return -0.0; // Preserve sign of operand
			}
			return y + 1.0;
		}
		if (d < 0.5)
		{
			return y;
		}

		/* half way, round to even */
		long z = static_cast<long>(y;
		return (z & 1) == 0 ? y : y + 1.0;
	}

	/** Get the closest long to x.
	 * @param x number from which closest long is requested
	 * @return closest long to x
	 */
	public static long round(const double& x)
	{
		const long bits = Double.double_to_raw_long_bits(x);
		const int biasedExp = (static_cast<int>((bits >> 52)) & 0x7ff;
		// Shift to get rid of bits past comma except first one: will need to
		// 1-shift to the right to end up with correct magnitude.
		const int shift = (52 - 1 + Double.MAX_EXPONENT) - biasedExp;
		if ((shift & -64) == 0)
		{
			// shift in [0,63], so unbiased exp in [-12,51].
			long extendedMantissa = 0x0010000000000000L | (bits & 0x000fffffffffffffL);
			if (bits < 0)
			{
				extendedMantissa = -extendedMantissa;
			}
			// If value is positive and first bit past comma is 0, rounding
			// to lower integer, else to upper one, which is what "+1" and
			// then ">>1" do.
			return ((extendedMantissa >> shift) + 1L) >> 1;
		}
		// +-Infinity, NaN, or a mathematical integer.
		return static_cast<long>(x;
	}

	/** Get the closest int to x.
	 * @param x number from which closest int is requested
	 * @return closest int to x
	 */
	public static int round(const float& x)
	{
		const int bits = Float.float_to_raw_int_bits(x);
		const int biasedExp = (bits >> 23) & 0xff;
		// Shift to get rid of bits past comma except first one: will need to
		// 1-shift to the right to end up with correct magnitude.
		const int shift = (23 - 1 + Float.MAX_EXPONENT) - biasedExp;
		if ((shift & -32) == 0)
		{
			// shift in [0,31], so unbiased exp in [-9,22].
			int extendedMantissa = 0x00800000 | (bits & 0x007fffff);
			if (bits < 0)
			{
				extendedMantissa = -extendedMantissa;
			}
			// If value is positive and first bit past comma is 0, rounding
			// to lower integer, else to upper one, which is what "+1" and
			// then ">>1" do.
			return ((extendedMantissa >> shift) + 1) >> 1;
		}
		else
		{
			// +-Infinity, NaN, or a mathematical integer.
			return static_cast<int>(x;
		}
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @return a if a is lesser or equal to b, b otherwise
	 */
	public static int min(const int& a, const int& b)
	{
		return (a <= b) ? a : b;
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @return a if a is lesser or equal to b, b otherwise
	 */
	public static long min(const long& a, const long& b)
	{
		return (a <= b) ? a : b;
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @return a if a is lesser or equal to b, b otherwise
	 */
	public static float min(const float& a, const float& b)
	{
		if (a > b)
		{
			return b;
		}
		if (a < b)
		{
			return a;
		}
		/* if either arg is NaN, return NaN */
		if (a != b)
		{
			return Float.NaN;
		}
		/* min(+0.0,-0.0) == -0.0 */
		/* 0x80000000 == Float.float_to_raw_int_bits(-0.0) */
		int bits = Float.float_to_raw_int_bits(a);
		if (bits == 0x80000000)
		{
			return a;
		}
		return b;
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @return a if a is lesser or equal to b, b otherwise
	 */
	public static double min(const double& a, const double& b)
	{
		if (a > b)
		{
			return b;
		}
		if (a < b)
		{
			return a;
		}
		/* if either arg is NaN, return NaN */
		if (a != b)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		/* min(+0.0,-0.0) == -0.0 */
		/* 0x8000000000000000L == Double.double_to_raw_long_bits(-0.0) */
		long bits = Double.double_to_raw_long_bits(a);
		if (bits == 0x8000000000000000L)
		{
			return a;
		}
		return b;
	}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @return b if a is lesser or equal to b, a otherwise
	 */
	public static int max(const int& a, const int& b)
	{
		return (a <= b) ? b : a;
	}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @return b if a is lesser or equal to b, a otherwise
	 */
	public static long max(const long& a, const long& b)
	{
		return (a <= b) ? b : a;
	}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @return b if a is lesser or equal to b, a otherwise
	 */
	public static float max(const float& a, const float& b)
	{
		if (a > b)
		{
			return a;
		}
		if (a < b)
		{
			return b;
		}
		/* if either arg is NaN, return NaN */
		if (a != b)
		{
			return Float.NaN;
		}
		/* min(+0.0,-0.0) == -0.0 */
		/* 0x80000000 == Float.float_to_raw_int_bits(-0.0) */
		int bits = Float.float_to_raw_int_bits(a);
		if (bits == 0x80000000)
		{
			return b;
		}
		return a;
	}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @return b if a is lesser or equal to b, a otherwise
	 */
	public static double max(const double& a, const double& b)
	{
		if (a > b)
		{
			return a;
		}
		if (a < b)
		{
			return b;
		}
		/* if either arg is NaN, return NaN */
		if (a != b)
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		/* min(+0.0,-0.0) == -0.0 */
		/* 0x8000000000000000L == Double.double_to_raw_long_bits(-0.0) */
		long bits = Double.double_to_raw_long_bits(a);
		if (bits == 0x8000000000000000L)
		{
			return b;
		}
		return a;
	}

	/**
	 * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
	 * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)<br/>
	 * avoiding intermediate overflow or underflow.
	 *
	 * <ul>
	 * <li> If either argument is infinite, then the result is positive infinity.</li>
	 * <li> else, if either argument is NaN then the result is NaN.</li>
	 * </ul>
	 *
	 * @param x a value
	 * @param y a value
	 * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
	 */
	public static double hypot(const double& x, const double y)
	{
		if (std::isinf(x) || std::isinf(y))
		{
			return INFINITY;
		}
		else if (std::isnan(x) || std::isnan(y))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		else
		{
			const int exp_x = get_exponent(x);
			const int exp_y = get_exponent(y);
			if (exp_x > exp_y + 27)
			{
				// y is neglectible with respect to x
				return abs(x);
			}
			else if (exp_y > exp_x + 27)
			{
				// x is neglectible with respect to y
				return abs(y);
			}
			else
			{
				// find an intermediate scale to avoid both overflow and underflow
				const int middle_exp = (exp_x + exp_y) / 2;

				// scale parameters without losing precision
				const double scaled_x = scalb(x, -middle_exp);
				const double scaled_y = scalb(y, -middle_exp);

				// compute scaled hypotenuse
				const double scaled_h = sqrt(scaled_x * scaled_x + scaled_y * scaled_y);

				// remove scaling
				return scalb(scaled_h, middle_exp);
			}
		}
	}

	/**
	 * Computes the remainder as prescribed by the IEEE 754 standard.
	 * The remainder value is mathematically equal to {@code x - y*n}
	 * where {@code n} is the mathematical integer closest to the exact mathematical value
	 * of the quotient {@code x/y}.
	 * If two mathematical integers are equally close to {@code x/y} then
	 * {@code n} is the integer that is even.
	 * <p>
	 * <ul>
	 * <li>If either operand is NaN, the result is NaN.</li>
	 * <li>If the result is not NaN, the sign of the result equals the sign of the dividend.</li>
	 * <li>If the dividend is an infinity, or the divisor is a zero, or both, the result is NaN.</li>
	 * <li>If the dividend is finite and the divisor is an infinity, the result equals the dividend.</li>
	 * <li>If the dividend is a zero and the divisor is finite, the result equals the dividend.</li>
	 * </ul>
	 * @param dividend the number to be divided
	 * @param divisor the number by which to divide
	 * @return the remainder, rounded
	 */
	public static double IEE_Eremainder(const double dividend, const double divisor)
	{
		if (get_exponent(dividend) == 1024 || get_exponent(divisor) == 1024 || divisor == 0.0)
		{
			// we are in one of the special cases
			if (std::isinf(divisor) && !std::isinf(dividend))
			{
				return dividend;
			}
			else
			{
				return std::numeric_limits<double>::quiet_NaN();
			}
		}
		else
		{
			// we are in the general case
			const double n = std::rint(dividend / divisor);
			const double remainder = std::isinf(n) ? 0.0 : dividend - divisor * n;
			return (remainder == 0) ? std::copysign(remainder, dividend) : remainder;
		}
	}

	/** Convert a long to interger, detecting overflows
	 * @param n number to convert to int
	 * @return integer with same valie as n if no overflows occur
	 * @exception Math_Runtime_Exception if n cannot fit into an int
	 */
	public static int toIntExact(const long n) Math_Runtime_Exception
	{
		if (n < std::numeric_limits<int>::min() || n > std::numeric_limits<int>::max())
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW);
		}
		return static_cast<int>(n;
	}

	/** Increment a number, detecting overflows.
	 * @param n number to increment
	 * @return n+1 if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static int incrementExact(const int& n) Math_Runtime_Exception
	{
		if (n == std::numeric_limits<int>::max())
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_ADDITION, n, 1);
		}

		return n + 1;
	}

	/** Increment a number, detecting overflows.
	 * @param n number to increment
	 * @return n+1 if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static long incrementExact(const long n) Math_Runtime_Exception
	{
		if (n == long.MAX_VALUE)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_ADDITION, n, 1);
		}

		return n + 1;
	}

	/** Decrement a number, detecting overflows.
	 * @param n number to decrement
	 * @return n-1 if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static int decrementExact(const int& n) Math_Runtime_Exception
	{
		if (n == std::numeric_limits<int>::min())
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_SUBTRACTION, n, 1);
		}

		return n - 1;
	}

	/** Decrement a number, detecting overflows.
	 * @param n number to decrement
	 * @return n-1 if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static long decrementExact(const long n) Math_Runtime_Exception
	{
		if (n == long.MIN_VALUE)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_SUBTRACTION, n, 1);
		}

		return n - 1;
	}

	/** Add two numbers, detecting overflows.
	 * @param a first number to add
	 * @param b second number to add
	 * @return a+b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static int add_exact(const int& a, const int& b) Math_Runtime_Exception
	{
		// compute sum
		const int sum = a + b;

		// check for overflow
		if ((a ^ b) >= 0 && (sum ^ b) < 0)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_ADDITION, a, b);
		}

		return sum;
	}

	/** Add two numbers, detecting overflows.
	 * @param a first number to add
	 * @param b second number to add
	 * @return a+b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static long add_exact(const long& a, const long& b) Math_Runtime_Exception
	{
		// compute sum
		const long sum = a + b;

		// check for overflow
		if ((a ^ b) >= 0 && (sum ^ b) < 0)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_ADDITION, a, b);
		}

		return sum;
	}

	/** Subtract two numbers, detecting overflows.
	 * @param a first number
	 * @param b second number to subtract from a
	 * @return a-b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static int subtractExact(const int& a, const int& b)
	{
		// compute subtraction
		const int sub = a - b;

		// check for overflow
		if ((a ^ b) < 0 && (sub ^ b) >= 0)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_SUBTRACTION, a, b);
		}

		return sub;
	}

	/** Subtract two numbers, detecting overflows.
	 * @param a first number
	 * @param b second number to subtract from a
	 * @return a-b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static long subtractExact(const long& a, const long& b)
	{
		// compute subtraction
		const long sub = a - b;

		// check for overflow
		if ((a ^ b) < 0 && (sub ^ b) >= 0)
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_SUBTRACTION, a, b);
		}

		return sub;
	}

	/** Multiply two numbers, detecting overflows.
	 * @param a first number to multiply
	 * @param b second number to multiply
	 * @return a*b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static int multiply_exact(const int& a, const int& b)
	{
		if (((b > 0) && (a > std::numeric_limits<int>::max() / b || a < std::numeric_limits<int>::min() / b)) ||
			((b < -1) && (a > std::numeric_limits<int>::min() / b || a < std::numeric_limits<int>::max() / b)) ||
			((b == -1) && (a == std::numeric_limits<int>::min())))
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_MULTIPLICATION, a, b);
		}
		return a * b;
	}

	/** Multiply two numbers, detecting overflows.
	 * @param a first number to multiply
	 * @param b second number to multiply
	 * @return a*b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 * @since 1.3
	 */
	public static long multiply_exact(const long& a, const int& b)
	{
		return multiply_exact(a, static_cast<long>(b);
	}

	/** Multiply two numbers, detecting overflows.
	 * @param a first number to multiply
	 * @param b second number to multiply
	 * @return a*b if no overflows occur
	 * @exception Math_Runtime_Exception if an overflow occurs
	 */
	public static long multiply_exact(const long& a, const long& b)
	{
		if (((b > 0l) && (a > long.MAX_VALUE / b || a < long.MIN_VALUE / b)) ||
			((b < -1l) && (a > long.MIN_VALUE / b || a < long.MAX_VALUE / b)) ||
			((b == -1l) && (a == long.MIN_VALUE)))
		{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_MULTIPLICATION, a, b);
		}
		return a * b;
	}

	/** Multiply two integers and give an exact result without overflow.
	 * @param a first factor
	 * @param b second factor
	 * @return a * b exactly
	 * @since 1.3
	 */
	public static long multiplyFull(const int& a, const int& b)
	{
		return (static_cast<long>(a) * (static_cast<long>(b);
	}

	/** Multiply two long integers and give the 64 most significant bits of the result.
	 * <p>
	 * Beware that as Java primitive long are always considered to be signed, there are some
	 * intermediate values {@code a} and {@code b} for which {@code a * b} exceeds {@code long.MAX_VALUE}
	 * but this method will still return 0l. This happens for example for {@code a = 2\xc2\xb3\xc2\xb9} and
	 * {@code b = 2\xc2\xb3\xc2\xb2} as {@code a * b = 2\xe2\x81\xb6\xc2\xb3 = long.MAX_VALUE + 1}, so it exceeds the max value
	 * for a long, but still fits in 64 bits, so this method correctly returns 0l in this case, * but multiplication result would be considered negative (and in fact equal to {@code long.MIN_VALUE}
	 * </p>
	 * @param a first factor
	 * @param b second factor
	 * @return a * b / 2<sup>64</sup>
	 * @since 1.3
	 */
	public static long multiplyHigh(const long& a, const long& b)
	{
		// all computations below are performed on unsigned numbers because we start
		// by using logical shifts (and not arithmetic shifts). We will therefore
		// need to take care of sign before returning
		// a negative long n between -2\xe2\x81\xb6\xc2\xb3 and -1, interpreted as an unsigned long
		// corresponds to 2\xe2\x81\xb6\xe2\x81\xb4 + n (which is between 2\xe2\x81\xb6\xc2\xb3 and 2\xe2\x81\xb6\xe2\x81\xb4-1)
		// so if this number is multiplied by p, what we really compute
		// is (2\xe2\x81\xb6\xe2\x81\xb4 + n) * p = 2\xe2\x81\xb6\xe2\x81\xb4 * p + n * p, therefore the part above 64 bits
		// will have an extra term p that we will need to remove
		const long tobeRemoved = ((a < 0) ? b : 0) + ((b < 0) ? a : 0);

		// split numbers in two 32 bits parts
		const long a_high = a >> > 32;
		const long a_low = a & 0xFFFFFFFFl;
		const long b_high = b >> > 32;
		const long b_low = b & 0xFFFFFFFFl;

		// ab = a_high * b_high * 2\xe2\x81\xb6\xe2\x81\xb4 + (a_high * b_low + a_low * b_high) * 2\xc2\xb3\xc2\xb2 + a_low * b_low
		const long hh = a_high * b_high;
		const long hl1 = a_high * b_low;
		const long hl2 = a_low * b_high;
		const long ll = a_low * b_low;

		// adds up everything in the above 64 bit part, taking care to avoid overflow
		const long hlHigh = (hl1 >> > 32) + (hl2 >> > 32);
		const long hlLow = (hl1 & 0xFFFFFFFFl) + (hl2 & 0xFFFFFFFFl);
		const long carry = (hlLow + (ll >> > 32)) >> > 32;
		const long unsignedResult = hh + hlHigh + carry;

		return unsignedResult - tobeRemoved;
	}

	/** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code  b < 0}.
	 * <p>
	 * This methods returns the same value as integer division when
	 * a and b are same signs, but returns a different value when
	 * they are opposite (i.e. q is negative).
	 *
	 * @param a dividend
	 * @param b divisor
	 * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code  b < 0}
	 * @exception Math_Runtime_Exception if b == 0
	 * @see #floorMod(int, int)
	 */
	public static int floorDiv(const int& a, const int& b)
	{
		if (b == 0)
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
		}

		const int m = a % b;
		return ((a ^ b) >= 0 || m == 0)
			? a / b	// a an b have same sign, or division is exact
			: (a / b) - 1; // a and b have opposite signs and division is not exact
	}

	/** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
	 * <p>
	 * This methods returns the same value as integer division when
	 * a and b are same signs, but returns a different value when
	 * they are opposite (i.e. q is negative).
	 *
	 * @param a dividend
	 * @param b divisor
	 * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
	 * @exception Math_Runtime_Exception if b == 0
	 * @see #floorMod(long, int)
	 * @since 1.3
	 */
	public static long floorDiv(const long& a, const int& b)
	{
		return floorDiv(a, static_cast<long>(b);
	}

	/** Finds q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
	 * <p>
	 * This methods returns the same value as integer division when
	 * a and b are same signs, but returns a different value when
	 * they are opposite (i.e. q is negative).
	 *
	 * @param a dividend
	 * @param b divisor
	 * @return q such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
	 * @exception Math_Runtime_Exception if b == 0
	 * @see #floorMod(long, long)
	 */
	public static long floorDiv(const long& a, const long& b) Math_Runtime_Exception
	{
		if (b == 0l)
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
		}

		const long m = a % b;
		return ((a ^ b) >= 0l || m == 0l)
			? a / b // a an b have same sign, or division is exact
			: (a / b) - 1l; // a and b have opposite signs and division is not exact
	}

	/** Finds r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
	 * <p>
	 * This methods returns the same value as integer modulo when
	 * a and b are same signs, but returns a different value when
	 * they are opposite (i.e. q is negative).
	 * </p>
	 * @param a dividend
	 * @param b divisor
	 * @return r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
	 * @exception Math_Runtime_Exception if b == 0
	 * @see #floorDiv(int, int)
	 */
	public static int floorMod(const int& a, const int& b) Math_Runtime_Exception
	{
		if (b == 0)
		{
			throw std::exception("not implemented");
			//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
		}

		const int m = a % b;
		return ((a ^ b) >= 0 || m == 0)
			? m // a an b have same sign, or division is exact
			: b + m; // a and b have opposite signs and division is not exact
	}
}

/** Finds r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
 * <p>
 * This methods returns the same value as integer modulo when
 * a and b are same signs, but returns a different value when
 * they are opposite (i.e. q is negative).
 * </p>
 * @param a dividend
 * @param b divisor
 * @return r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
 * @exception Math_Runtime_Exception if b == 0
 * @see #floorDiv(long, int)
 * @since 1.3
 */
public static int floorMod(const long& a, const int& b)
{
	return static_cast<int>(floorMod(a, static_cast<long>(b);
}

/** Finds r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}.
 * <p>
 * This methods returns the same value as integer modulo when
 * a and b are same signs, but returns a different value when
 * they are opposite (i.e. q is negative).
 * </p>
 * @param a dividend
 * @param b divisor
 * @return r such that {@code a = q b + r} with {@code 0 <= r < b} if {@code b > 0} and {@code b < r <= 0} if {@code b < 0}
 * @exception Math_Runtime_Exception if b == 0
 * @see #floorDiv(long, long)
 */
public static long floorMod(const long& a, const long& b)
{
	if (b == 0l)
	{
		throw std::exception("not implemented");
		//throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR);
	}

	const long m = a % b;
	return ((a ^ b) >= 0l || m == 0l)
		? m // a an b have same sign, or division is exact
		: b + m; // a and b have opposite signs and division is not exact
}
	}

	/**
	 * Returns the first argument with the sign of the second argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param magnitude the value to return
	 * @param sign the sign for the returned value
	 * @return the magnitude with the same sign as the {@code sign} argument
	 */
	public static double copy_sign(const double& magnitude, const double& sign)
	{
		// The highest order bit is going to be zero if the
		// highest order bit of m and s is the same and one otherwise.
		// So (m^s) will be positive if both m and s have the same sign
		// and negative otherwise.
		const long m = Double.double_to_raw_long_bits(magnitude); // don't care about NaN
		const long s = Double.double_to_raw_long_bits(sign);
		if ((m ^ s) >= 0)
		{
			return magnitude;
		}
		return -magnitude; // flip sign
	}

	/**
	 * Returns the first argument with the sign of the second argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param magnitude the value to return
	 * @param sign the sign for the returned value
	 * @return the magnitude with the same sign as the {@code sign} argument
	 */
	public static float copy_sign(const float& magnitude, const float& sign)
	{
		// The highest order bit is going to be zero if the
		// highest order bit of m and s is the same and one otherwise.
		// So (m^s) will be positive if both m and s have the same sign
		// and negative otherwise.
		const int m = Float.float_to_raw_int_bits(magnitude);
		const int s = Float.float_to_raw_int_bits(sign);
		if ((m ^ s) >= 0)
		{
			return magnitude;
		}
		return -magnitude; // flip sign
	}

	/**
	 * Return the exponent of a double number, removing the bias.
	 * <p>
	 * For double numbers of the form 2<sup>x</sup>, the unbiased
	 * exponent is exactly x.
	 * </p>
	 * @param d number from which exponent is requested
	 * @return exponent for d in IEEE754 representation, without bias
	 */
	public static int get_exponent(const double& d)
	{
		// NaN and Infinite will return 1024 anyhow so can use raw bits
		return static_cast<int>(((Double.double_to_raw_long_bits(d) >> > 52) & 0x7ff) - 1023;
	}

	/**
	 * Return the exponent of a float number, removing the bias.
	 * <p>
	 * For float numbers of the form 2<sup>x</sup>, the unbiased
	 * exponent is exactly x.
	 * </p>
	 * @param f number from which exponent is requested
	 * @return exponent for d in IEEE754 representation, without bias
	 */
	public static int get_exponent(const float& f)
	{
		// NaN and Infinite will return the same exponent anyhow so can use raw bits
		return ((Float.float_to_raw_int_bits(f) >> > 23) & 0xff) - 127;
	}

	/** Compute Fused-multiply-add operation a * b + c.
	 * <p>
	 * This method was introduced in the regular {@code Math} and {@code StrictMath}
	 * methods with Java 9, and then added to Hipparchus for consistency. However, * a more general method was available in Hipparchus that also allow to repeat
	 * this computation across several terms: {@link Math_Arrays#linear_combination(std::vector<double>, std::vector<double>)}.
	 * The linear combination method should probably be preferred in most cases.
	 * </p>
	 * @param a first factor
	 * @param b second factor
	 * @param c additive term
	 * @return a * b + c, using extended precision in the multiplication
	 * @see Math_Arrays#linear_combination(std::vector<double>, std::vector<double>)
	 * @see Math_Arrays#linear_combination(double, double, double, double)
	 * @see Math_Arrays#linear_combination(double, double, double, double, double, double)
	 * @see Math_Arrays#linear_combination(double, double, double, double, double, double, double, double)
	 * @since 1.3
	 */
	public static double fma(const double& a, const double b, const double& c)
	{
		return Math_Arrays::linear_combination(a, b, 1.0, c);
	}

	/** Compute Fused-multiply-add operation a * b + c.
	 * <p>
	 * This method was introduced in the regular {@code Math} and {@code StrictMath}
	 * methods with Java 9, and then added to Hipparchus for consistency. However, * a more general method was available in Hipparchus that also allow to repeat
	 * this computation across several terms: {@link Math_Arrays#linear_combination(std::vector<double>, std::vector<double>)}.
	 * The linear combination method should probably be preferred in most cases.
	 * </p>
	 * @param a first factor
	 * @param b second factor
	 * @param c additive term
	 * @return a * b + c, using extended precision in the multiplication
	 * @see Math_Arrays#linear_combination(std::vector<double>, std::vector<double>)
	 * @see Math_Arrays#linear_combination(double, double, double, double)
	 * @see Math_Arrays#linear_combination(double, double, double, double, double, double)
	 * @see Math_Arrays#linear_combination(double, double, double, double, double, double, double, double)
	 */
	public static float fma(const float& a, const float b, const float c)
	{
		return (float)Math_Arrays::linear_combination(a, b, 1.0, c);
	}

	/** Compute the square root of a number.
	 * @param a number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return square root of a
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T sqrt(const T a)
	{
		return a.sqrt();
	}

	/** Compute the hyperbolic cosine of a number.
	 * @param x number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return hyperbolic cosine of x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T cosh(const T& x)
	{
		return x.cosh();
	}

	/** Compute the hyperbolic sine of a number.
	 * @param x number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return hyperbolic sine of x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T sinh(const T& x)
	{
		return x.sinh();
	}

	/** Compute the hyperbolic tangent of a number.
	 * @param x number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return hyperbolic tangent of x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T tanh(const T& x)
	{
		return x.tanh();
	}

	/** Compute the inverse hyperbolic cosine of a number.
	 * @param a number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return inverse hyperbolic cosine of a
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T acosh(const T a)
	{
		return a.acosh();
	}

	/** Compute the inverse hyperbolic sine of a number.
	 * @param a number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return inverse hyperbolic sine of a
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T asinh(const T a)
	{
		return a.asinh();
	}

	/** Compute the inverse hyperbolic tangent of a number.
	 * @param a number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return inverse hyperbolic tangent of a
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T atanh(const T a)
	{
		return a.atanh();
	}

	/** Compute the sign of a number.
	 * The sign is -1 for negative numbers, +1 for positive numbers and 0 otherwise, * for std::complex<double> number, it is extended on the unit circle (equivalent to z/|z|, * with special handling for 0 and NaN)
	 * @param a number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return -1.0, -0.0, +0.0, +1.0 or NaN depending on sign of a
	 * @since 2.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T sign(const T a)
	{
		return a.sign();
	}

	/**
	 * Exponential function.
	 *
	 * Computes exp(x), function result is nearly rounded.   It will be correctly
	 * rounded to the theoretical value for 99.9% of input values, otherwise it will
	 * have a 1 ULP error.
	 *
	 * Method:
	 *    Lookup intVal = exp(int(x))
	 *    Lookup fracVal = exp(int(x-int(x) / 1024.0) * 1024.0 );
	 *    Compute z as the exponential of the remaining bits by a polynomial minus one
	 *    exp(x) = intVal * fracVal * (1 + z)
	 *
	 * Accuracy:
	 *    Calculation is done with 63 bits of precision, so result should be correctly
	 *    rounded for 99.9% of input values, with less than 1 ULP error otherwise.
	 *
	 * @param x   a double
	 * @param <T> the type of the field element
	 * @return double e<sup>x</sup>
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T exp(const T& x)
	{
		return x.exp();
	}

	/** Compute exp(x) - 1
	 * @param x number to compute shifted exponential
	 * @param <T> the type of the field element
	 * @return exp(x) - 1
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T expm1(const T& x)
	{
		return x.expm1();
	}

	/**
	 * Natural logarithm.
	 *
	 * @param x   a double
	 * @param <T> the type of the field element
	 * @return log(x)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T log(const T& x)
	{
		return x.log();
	}

	/**
	 * Computes log(1 + x).
	 *
	 * @param x Number.
	 * @param <T> the type of the field element
	 * @return {@code log(1 + x)}.
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T log1p(const T& x)
	{
		return x.log1p();
	}

	/** Compute the base 10 logarithm.
	 * @param x a number
	 * @param <T> the type of the field element
	 * @return log10(x)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T log10(const T& x)
	{
		return x.log10();
	}

	/**
	 * Power function.  Compute x<sup>y</sup>.
	 *
	 * @param x   a double
	 * @param y   a double
	 * @param <T> the type of the field element
	 * @return x<sup>y</sup>
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T pow(const T& x, const T& y)
	{
		return x.pow(y);
	}

	/**
	 * Power function.  Compute x<sup>y</sup>.
	 *
	 * @param x   a double
	 * @param y   a double
	 * @param <T> the type of the field element
	 * @return x<sup>y</sup>
	 * @since 1.7
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T pow(const T& x, const double y)
	{
		return x.pow(y);
	}

	/**
	 * Raise a double to an int power.
	 *
	 * @param d Number to raise.
	 * @param e Exponent.
	 * @param <T> the type of the field element
	 * @return d<sup>e</sup>
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T pow(T d, const int& e)
	{
		return d.pow(e);
	}

	/**
	 * Sine function.
	 *
	 * @param x Argument.
	 * @param <T> the type of the field element
	 * @return sin(x)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T sin(const T& x)
	{
		return x.sin();
	}

	/**
	 * Cosine function.
	 *
	 * @param x Argument.
	 * @param <T> the type of the field element
	 * @return cos(x)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T cos(const T& x)
	{
		return x.cos();
	}

	/**
	 * Tangent function.
	 *
	 * @param x Argument.
	 * @param <T> the type of the field element
	 * @return tan(x)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T tan(const T& x)
	{
		return x.tan();
	}

	/**
	 * Arctangent function
	 *  @param x a number
	 * @param <T> the type of the field element
	 *  @return atan(x)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T atan(const T& x)
	{
		return x.atan();
	}

	/**
	 * Two arguments arctangent function
	 * @param y ordinate
	 * @param x abscissa
	 * @param <T> the type of the field element
	 * @return phase angle of point (x,y) between {@code -PI} and {@code PI}
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T atan2(const T y, const T x)
	{
		return y.atan2(x);
	}

	/** Compute the arc sine of a number.
	 * @param x number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return arc sine of x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T asin(const T& x)
	{
		return x.asin();
	}

	/** Compute the arc cosine of a number.
	 * @param x number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return arc cosine of x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T acos(const T& x)
	{
		return x.acos();
	}

	/** Compute the cubic root of a number.
	 * @param x number on which evaluation is done
	 * @param <T> the type of the field element
	 * @return cubic root of x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T cbrt(const T& x)
	{
		return x.cbrt();
	}

	/**
	 * Norm.
	 * @param x number from which norm is requested
	 * @param <T> the type of the field element
	 * @return norm(x)
	 * @since 2.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static double norm(const T& x)
	{
		return x.norm();
	}

	/**
	 * Absolute value.
	 * @param x number from which absolute value is requested
	 * @param <T> the type of the field element
	 * @return abs(x)
	 * @since 2.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T abs(const T& x)
	{
		return x.abs();
	}

	/**
	 *  Convert degrees to radians, with error of less than 0.5 ULP
	 *  @param x angle in degrees
	 *  @param <T> the type of the field element
	 *  @return x converted into radians
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T to_radians(T x)
	{
		return x.to_radians();
	}

	/**
	 *  Convert radians to degrees, with error of less than 0.5 ULP
	 *  @param x angle in radians
	 *  @param <T> the type of the field element
	 *  @return x converted into degrees
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T to_degrees(T x)
	{
		return x.to_degrees();
	}

	/**
	 * Multiply a double number by a power of 2.
	 * @param d number to multiply
	 * @param n power of 2
	 * @param <T> the type of the field element
	 * @return d &times; 2<sup>n</sup>
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T scalb(const T& d, const int& n)
	{
		return d.scalb(n);
	}

	/**
	 * Compute least significant bit (Unit in Last Position) for a number.
	 * @param x number from which ulp is requested
	 * @param <T> the type of the field element
	 * @return ulp(x)
	 * @since 2.0
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T ulp(const T& x)
	{
		return (std::isinf(x.get_real()))
			? x.new_instance(INFINITY)
			: x.ulp();
	}

	/** Get the largest whole number smaller than x.
	 * @param x number from which floor is requested
	 * @param <T> the type of the field element
	 * @return a double number f such that f is an integer f &lt;= x &lt; f + 1.0
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T floor(const T& x)
	{
		return x.floor();
	}

	/** Get the smallest whole number larger than x.
	 * @param x number from which ceil is requested
	 * @param <T> the type of the field element
	 * @return a double number c such that c is an integer c - 1.0 &lt; x &lt;= c
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T ceil(const T& x)
	{
		return x.ceil();
	}

	/** Get the whole number that is the nearest to x, or the even one if x is exactly half way between two integers.
	 * @param x number from which nearest whole number is requested
	 * @param <T> the type of the field element
	 * @return a double number r such that r is an integer r - 0.5 &lt;= x &lt;= r + 0.5
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T rint(const T& x)
	{
		return x.rint();
	}

	/** Get the closest long to x.
	 * @param x number from which closest long is requested
	 * @param <T> the type of the field element
	 * @return closest long to x
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static long round(const T& x)
	{
		return x.round();
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @param <T> the type of the field element
	 * @return a if a is lesser or equal to b, b otherwise
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T min(const T& a, const T& b)
	{
		const double aR = a.get_real();
		const double bR = b.get_real();
		if (aR < bR)
		{
			return a;
		}
		if (bR < aR)
		{
			return b;
		}
		// either the numbers are equal, or one of them is a NaN
		return std::isnan(aR) ? a : b;
	}

	/** Compute the minimum of two values
	 * @param a first value
	 * @param b second value
	 * @param <T> the type of the field element
	 * @return a if a is lesser or equal to b, b otherwise
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T min(const T& a, const double& b)
	{
		const double aR = a.get_real();
		if (aR < b)
		{
			return a;
		}
		if (b < aR)
		{
			return a.get_field().get_zero().add(b);
		}
		// either the numbers are equal, or one of them is a NaN
		return std::isnan(aR)
			? a
			: a.get_field().get_zero().add(b);
	}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @param <T> the type of the field element
	 * @return b if a is lesser or equal to b, a otherwise
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T max(const T& a, const T& b)
	{
		const double aR = a.get_real();
		const double bR = b.get_real();
		if (aR < bR)
		{
			return b;
		}
		if (bR < aR)
		{
			return a;
		}
		// either the numbers are equal, or one of them is a NaN
		return std::isnan(aR)
			? a
			: b;
	}

	/** Compute the maximum of two values
	 * @param a first value
	 * @param b second value
	 * @param <T> the type of the field element
	 * @return b if a is lesser or equal to b, a otherwise
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T max(const T& a, const double& b)
	{
		const double aR = a.get_real();
		if (aR < b)
		{
			return a.get_field().get_zero().add(b);
		}
		if (b < aR)
		{
			return a;
		}
		// either the numbers are equal, or one of them is a NaN
		return std::isnan(aR)
			? a
			: a.get_field().get_zero().add(b);
	}

	/**
	 * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
	 * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)<br/>
	 * avoiding intermediate overflow or underflow.
	 *
	 * <ul>
	 * <li> If either argument is infinite, then the result is positive infinity.</li>
	 * <li> else, if either argument is NaN then the result is NaN.</li>
	 * </ul>
	 *
	 * @param x a value
	 * @param y a value
	 * @param <T> the type of the field element
	 * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T hypot(const T& x, const T& y)
	{
		return x.hypot(y);
	}

	/**
	 * Computes the remainder as prescribed by the IEEE 754 standard.
	 * The remainder value is mathematically equal to {@code x - y*n}
	 * where {@code n} is the mathematical integer closest to the exact mathematical value
	 * of the quotient {@code x/y}.
	 * If two mathematical integers are equally close to {@code x/y} then
	 * {@code n} is the integer that is even.
	 * <p>
	 * <ul>
	 * <li>If either operand is NaN, the result is NaN.</li>
	 * <li>If the result is not NaN, the sign of the result equals the sign of the dividend.</li>
	 * <li>If the dividend is an infinity, or the divisor is a zero, or both, the result is NaN.</li>
	 * <li>If the dividend is finite and the divisor is an infinity, the result equals the dividend.</li>
	 * <li>If the dividend is a zero and the divisor is finite, the result equals the dividend.</li>
	 * </ul>
	 * @param dividend the number to be divided
	 * @param divisor the number by which to divide
	 * @param <T> the type of the field element
	 * @return the remainder, rounded
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T IEE_Eremainder(const T dividend, const double divisor)
	{
		return dividend.remainder(divisor);
	}

	/**
	 * Computes the remainder as prescribed by the IEEE 754 standard.
	 * The remainder value is mathematically equal to {@code x - y*n}
	 * where {@code n} is the mathematical integer closest to the exact mathematical value
	 * of the quotient {@code x/y}.
	 * If two mathematical integers are equally close to {@code x/y} then
	 * {@code n} is the integer that is even.
	 * <p>
	 * <ul>
	 * <li>If either operand is NaN, the result is NaN.</li>
	 * <li>If the result is not NaN, the sign of the result equals the sign of the dividend.</li>
	 * <li>If the dividend is an infinity, or the divisor is a zero, or both, the result is NaN.</li>
	 * <li>If the dividend is finite and the divisor is an infinity, the result equals the dividend.</li>
	 * <li>If the dividend is a zero and the divisor is finite, the result equals the dividend.</li>
	 * </ul>
	 * @param dividend the number to be divided
	 * @param divisor the number by which to divide
	 * @param <T> the type of the field element
	 * @return the remainder, rounded
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T IEE_Eremainder(const T dividend, const T divisor)
	{
		return dividend.remainder(divisor);
	}

	/**
	 * Returns the first argument with the sign of the second argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param magnitude the value to return
	 * @param sign the sign for the returned value
	 * @param <T> the type of the field element
	 * @return the magnitude with the same sign as the {@code sign} argument
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T copy_sign(T magnitude, T sign)
	{
		return magnitude.copy_sign(sign);
	}

	/**
	 * Returns the first argument with the sign of the second argument.
	 * A NaN {@code sign} argument is treated as positive.
	 *
	 * @param magnitude the value to return
	 * @param sign the sign for the returned value
	 * @param <T> the type of the field element
	 * @return the magnitude with the same sign as the {@code sign} argument
	 * @since 1.3
	 */
	template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	public static T copy_sign(T magnitude, double sign)
	{
		return magnitude.copy_sign(sign);
	}

	//    /**
	//     * Print out contents of arrays, and check the length.
	//     * <p>used to generate the preset arrays originally.</p>
	//     * @param a unused
	//     */
	//    public static void main(const std::vector<std::string>& a)
	{
		//        Fast_Math_Calc.printarray(System.out, "EXP_INT_TABLE_A", EXP_INT_TABLE_LEN, ExpIntTable.EXP_INT_TABLE_A);
		//        Fast_Math_Calc.printarray(System.out, "EXP_INT_TABLE_B", EXP_INT_TABLE_LEN, ExpIntTable.EXP_INT_TABLE_B);
		//        Fast_Math_Calc.printarray(System.out, "EXP_FRAC_TABLE_A", EXP_FRAC_TABLE_LEN, Exp_fracTable.EXP_FRAC_TABLE_A);
		//        Fast_Math_Calc.printarray(System.out, "EXP_FRAC_TABLE_B", EXP_FRAC_TABLE_LEN, Exp_fracTable.EXP_FRAC_TABLE_B);
		//        Fast_Math_Calc.printarray(System.out, "LN_MANT",LN_MANT_LEN, lnMant.LN_MANT);
		//        Fast_Math_Calc.printarray(System.out, "SINE_TABLE_A", SINE_TABLE_LEN, SINE_TABLE_A);
		//        Fast_Math_Calc.printarray(System.out, "SINE_TABLE_B", SINE_TABLE_LEN, SINE_TABLE_B);
		//        Fast_Math_Calc.printarray(System.out, "COSINE_TABLE_A", SINE_TABLE_LEN, COSINE_TABLE_A);
		//        Fast_Math_Calc.printarray(System.out, "COSINE_TABLE_B", SINE_TABLE_LEN, COSINE_TABLE_B);
		//        Fast_Math_Calc.printarray(System.out, "TANGENT_TABLE_A", SINE_TABLE_LEN, TANGENT_TABLE_A);
		//        Fast_Math_Calc.printarray(System.out, "TANGENT_TABLE_B", SINE_TABLE_LEN, TANGENT_TABLE_B);
		//    }

			/** Enclose large data table in nested static class so it's only loaded on first access. */
		private static class ExpIntTable
		{
		private:
			/** Exponential evaluated at integer values, * exp(x) =  exp_int_table_a[x + EXP_INT_TABLE_MAX_INDEX] + exp_int_table_b[x+EXP_INT_TABLE_MAX_INDEX].
			 */
			static const std::vector<double> EXP_INT_TABLE_A;
			/** Exponential evaluated at integer values, * exp(x) =  exp_int_table_a[x + EXP_INT_TABLE_MAX_INDEX] + exp_int_table_b[x+EXP_INT_TABLE_MAX_INDEX]
			 */
			static const std::vector<double> EXP_INT_TABLE_B;

		public:
			static
			{
				if (RECOMPUTE_TABLES_AT_RUNTIME)
				{
					EXP_INT_TABLE_A = std::vector<double>(FastMath.EXP_INT_TABLE_LEN];
					EXP_INT_TABLE_B = std::vector<double>(FastMath.EXP_INT_TABLE_LEN];

					auto tmp = std::vector<double>(2);
					auto recip = std::vector<double>(2);

					// Populate expIntTable
					for (int i{}; i < FastMath.EXP_INT_TABLE_MAX_INDEX; i++)
					{
						Fast_Math_Calc.expint(i, tmp);
						EXP_INT_TABLE_A[i + FastMath.EXP_INT_TABLE_MAX_INDEX] = tmp[0];
						EXP_INT_TABLE_B[i + FastMath.EXP_INT_TABLE_MAX_INDEX] = tmp[1];

						if (i != 0)
						{
							// Negative integer powers
							Fast_Math_Calc.split_reciprocal(tmp, recip);
							EXP_INT_TABLE_A[FastMath.EXP_INT_TABLE_MAX_INDEX - i] = recip[0];
							EXP_INT_TABLE_B[FastMath.EXP_INT_TABLE_MAX_INDEX - i] = recip[1];
						}
					}
				}
				else
				{
					EXP_INT_TABLE_A = Fast_Math_Literal_Arrays.load_exp_int_a();
					EXP_INT_TABLE_B = Fast_Math_Literal_Arrays.load_exp_int_b();
				}
			}
		}

		/** Enclose large data table in nested static class so it's only loaded on first access. */
		private static class Exp_fracTable
		{
		private:
			/** Exponential over the range of 0 - 1 in increments of 2^-10
			 * exp(x/1024) =  exp_frac_table_a[x] + exp_frac_table_b[x].
			 * 1024 = 2^10
			 */
			static const std::vector<double> EXP_FRAC_TABLE_A;
			/** Exponential over the range of 0 - 1 in increments of 2^-10
			 * exp(x/1024) =  exp_frac_table_a[x] + exp_frac_table_b[x].
			 */
			static const std::vector<double> EXP_FRAC_TABLE_B;

		public:
			static
			{
				if (RECOMPUTE_TABLES_AT_RUNTIME)
				{
					EXP_FRAC_TABLE_A = std::vector<double>(FastMath.EXP_FRAC_TABLE_LEN];
					EXP_FRAC_TABLE_B = std::vector<double>(FastMath.EXP_FRAC_TABLE_LEN];

					auto tmp = std::vector<double>(2);

					// Populate exp_fracTable
					const double factor = 1.0 / (EXP_FRAC_TABLE_LEN - 1);
					for (int i{}; i < EXP_FRAC_TABLE_A.size(); i++)
					{
						Fast_Math_Calc.slowexp(i * factor, tmp);
						EXP_FRAC_TABLE_A[i] = tmp[0];
						EXP_FRAC_TABLE_B[i] = tmp[1];
					}
				}
				else
				{
					EXP_FRAC_TABLE_A = Fast_Math_Literal_Arrays.load_exp_frac_a();
					EXP_FRAC_TABLE_B = Fast_Math_Literal_Arrays.load_exp_frac_b();
				}
			}
		}

		/** Enclose large data table in nested static class so it's only loaded on first access. */
		private static class lnMant
		{
		private:
			/** Extended precision logarithm table over the range 1 - 2 in increments of 2^-10. */
			static std::vector<std::vector<double>> LN_MANT;

			static
			{
				if (RECOMPUTE_TABLES_AT_RUNTIME)
				{
					LN_MANT = std::vector<std::vector<double>>(FastMath.LN_MANT_LEN);

					// Populate lnMant table
					for (int i{}; i < LN_MANT.size(); i++)
					{
						const double d = Double.long_bits_to_double(((static_cast<long>(i) << 42) | 0x3ff0000000000000L);
						LN_MANT[i] = Fast_Math_Calc.slow_log(d);
					}
				}
				else
				{
					LN_MANT = Fast_Math_Literal_Arrays.load_ln_mant();
				}
			}
		}

		/** Enclose the Cody/Waite reduction (used in "sin", "cos" and "tan"). */
		private static class CodyWaite
		{
		private:
			/** k */
			int my_constK;
			/** remA */
			double my_const_rem_a;
			/** remB */
			double my_const_rem_b;

		public:
			/**
			 * @param xa Argument.
			 */
			CodyWaite(const double& xa)
			{
				// Estimate k.
				//k = static_cast<int>((xa / 1.5707963267948966);
				int k = static_cast<int>((xa * 0.6366197723675814);

				// Compute remainder.
				double remA;
				double remB;
				while (true)
				{
					double a = -k * 1.570796251296997;
					remA = xa + a;
					remB = -(remA - xa - a);

					a = -k * 7.549789948768648E-8;
					double b = remA;
					remA = a + b;
					remB += -(remA - b - a);

					a = -k * 6.123233995736766E-17;
					b = remA;
					remA = a + b;
					remB += -(remA - b - a);

					if (remA > 0)
					{
						break;
					}

					// Remainder is negative, so decrement k and try again.
					// This should only happen if the input is very close
					// to an even multiple of pi/2.
					--k;
				}

				my_constK = k;
				my_const_rem_a = remA;
				my_const_rem_b = remB;
			}

			/**
			 * @return k
			 */
			int get_k() const
			{
				return my_constK;
			}
			/**
			 * @return remA
			 */
			double get_remA() const
			{
				return my_const_rem_a;
			}
			/**
			 * @return remB
			 */
			double get_remB() const
			{
				return my_const_rem_b;
			}
		}
	};