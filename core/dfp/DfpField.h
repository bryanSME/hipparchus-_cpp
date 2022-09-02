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

  //package org.hipparchus.dfp;
#include <string>
#include <vector>

  //import org.hipparchus.Field;

  /** Field for Decimal floating point instances.
   */
class DFP_Field : Field<Dfp>
{
	/** Enumerate for rounding modes. */
	enum Rounding_Mode
	{
		/** Rounds toward zero (truncation). */
		ROUND_DOWN,
		/** Rounds away from zero if discarded digit is non-zero. */
		ROUND_UP,
		/** Rounds towards nearest unless both are equidistant in which case it rounds away from zero. */
		ROUND_HALF_UP,
		/** Rounds towards nearest unless both are equidistant in which case it rounds toward zero. */
		ROUND_HALF_DOWN,
		/** Rounds towards nearest unless both are equidistant in which case it rounds toward the even neighbor.
		 * This is the default as  specified by IEEE 854-1987
		 */
		 ROUND_HALF_EVEN,
		 /** Rounds towards nearest unless both are equidistant in which case it rounds toward the odd neighbor.  */
		 ROUND_HALF_ODD,
		 /** Rounds towards positive infinity. */
		 ROUND_CEIL,
		 /** Rounds towards negative infinity. */
		 ROUND_FLOOR;
	}

private:
	/** High precision string representation of &radic;2. */
	static std::string sqr_2_string;

	// Note: the static strings are set up (once) by the ctor and @Guarded_By("DFP_Field.class")

	/** High precision string representation of &radic;2 / 2. */
	static std::string sqr2_reciprocal_string;

	/** High precision string representation of &radic;3. */
	static std::string sqr3std::string;

	/** High precision string representation of &radic;3 / 3. */
	static std::string sqr3_reciprocal_string;

	/** High precision string representation of &pi;. */
	static std::string pi_string;

	/** High precision string representation of e. */
	static std::string e_string;

	/** High precision string representation of ln(2). */
	static std::string ln2std::string;

	/** High precision string representation of ln(5). */
	static std::string ln5std::string;

	/** High precision string representation of ln(10). */
	static std::string ln10std::string;

	/** The number of radix digits.
	 * Note these depend on the radix which is 10000 digits, * so each one is equivalent to 4 decimal digits.
	 */
	const int radix_digits;

	/** A {@link Dfp} with value 0. */
	const Dfp zero;

	/** A {@link Dfp} with value 1. */
	const Dfp one;

	/** A {@link Dfp} with value 2. */
	const Dfp two;

	/** A {@link Dfp} with value &radic;2. */
	const Dfp sqr2;

	/** A two elements {@link Dfp} array with value &radic;2 split in two pieces. */
	const Dfp[] sqr2_split;

	/** A {@link Dfp} with value &radic;2 / 2. */
	const Dfp sqr2_reciprocal;

	/** A {@link Dfp} with value &radic;3. */
	const Dfp sqr3;

	/** A {@link Dfp} with value &radic;3 / 3. */
	const Dfp sqr3_reciprocal;

	/** A {@link Dfp} with value &pi;. */
	const Dfp pi;

	/** A {@link Dfp} for converting degrees to radians. */
	const Dfp deg_to_rad;

	/** A {@link Dfp} for converting radians to degrees. */
	const Dfp rad_to_deg;

	/** A two elements {@link Dfp} array with value &pi; split in two pieces. */
	const Dfp[] pi_split;

	/** A {@link Dfp} with value e. */
	const Dfp e;

	/** A two elements {@link Dfp} array with value e split in two pieces. */
	const Dfp[] e_split;

	/** A {@link Dfp} with value ln(2). */
	const Dfp ln2;

	/** A two elements {@link Dfp} array with value ln(2) split in two pieces. */
	const Dfp[] ln2_split;

	/** A {@link Dfp} with value ln(5). */
	const Dfp ln5;

	/** A two elements {@link Dfp} array with value ln(5) split in two pieces. */
	const Dfp[] ln5_split;

	/** A {@link Dfp} with value ln(10). */
	const Dfp ln10;

	/** Current rounding mode. */
	Rounding_Mode r_mode;

	/** IEEE 854-1987 signals. */
	int ieee_flags;

	/** Create a factory for the specified number of radix digits.
	 * <p>
	 * Note that since the {@link Dfp} class uses 10000 as its radix, each radix
	 * digit is equivalent to 4 decimal digits. This implies that asking for
	 * 13, 14, 15 or 16 decimal digits will really lead to a 4 radix 10000 digits in
	 * all cases.
	 * </p>
	 * @param decimal_digits minimal number of decimal digits
	 * @param compute_constants if true, the transcendental constants for the given precision
	 * must be computed (setting this flag to false is RESERVED for the internal recursive call)
	 */
	DFP_Field(const int& decimal_digits, const bool compute_constants)
	{
		this.radix_digits = decimal_digits < 13
			? 4
			: (decimal_digits + 3) / 4;
		this.r_mode = Rounding_Mode.ROUND_HALF_EVEN;
		this.ieee_flags = 0;
		this.zero = Dfp(*this, 0);
		this.one = Dfp(*this, 1);
		this.two = Dfp(*this, 2);

		if (compute_constants)
		{
			// set up transcendental constants
			synchronized(DFP_Field.class)
			{
				// as a heuristic to circumvent Table-Maker's Dilemma, we set the string
				// representation of the constants to be at least 3 times larger than the
				// number of decimal digits, also as an attempt to really compute these
				// constants only once, we set a minimum number of digits
				compute_string_constants((decimal_digits < 67) ? 200 : (3 * decimal_digits));

				// set up the constants at current field accuracy
				sqr2 = Dfp(this, sqr2std::string);
				sqr2_split = split(sqr2std::string);
				sqr2_reciprocal = Dfp(this, sqr2_reciprocal_string);
				sqr3 = Dfp(this, sqr3std::string);
				sqr3_reciprocal = Dfp(this, sqr3_reciprocal_string);
				pi = Dfp(this, pi_string);
				deg_to_rad = pi.divide(180.0);
				rad_to_deg = pi.divide(180.0).reciprocal();
				pi_split = split(pi_string);
				e = Dfp(this, e_string);
				e_split = split(e_string);
				ln2 = Dfp(this, ln2std::string);
				ln2_split = split(ln2std::string);
				ln5 = Dfp(this, ln5std::string);
				ln5_split = split(ln5std::string);
				ln10 = Dfp(this, ln10std::string);
			}
		}
		else
		{
			// dummy settings for unused constants
			sqr2 = NULL;
			sqr2_split = NULL;
			sqr2_reciprocal = NULL;
			sqr3 = NULL;
			sqr3_reciprocal = NULL;
			pi = NULL;
			deg_to_rad = NULL;
			rad_to_deg = NULL;
			pi_split = NULL;
			e = NULL;
			e_split = NULL;
			ln2 = NULL;
			ln2_split = NULL;
			ln5 = NULL;
			ln5_split = NULL;
			ln10 = NULL;
		}
	}

	/** Breaks a string representation up into two {@link Dfp}'s.
	 * The split is such that the sum of them is equivalent to the input string, * but has higher precision than using a single Dfp.
	 * @param a string representation of the number to split
	 * @return an array of two {@link Dfp Dfp} instances which sum equals a
	 */
	std::vector<Dfp> split(const std::string& a)
	{
		auto result = std::vector<Dfp>(2);
		bool leading{ true };
		int sp{};
		int sig{};

		std::stringBuilder builder1 = std::stringBuilder(a.size()());

		for (int i{}; i < a.size()(); i++)
		{
			const char c = a.char_at(i);
			builder1.append(c);

			if (c >= '1' && c <= '9')
			{
				leading = false;
			}

			if (c == '.')
			{
				sig += (400 - sig) % 4;
				leading = false;
			}

			if (sig == (radix_digits / 2) * 4)
			{
				sp = i;
				break;
			}

			if (c >= '0' && c <= '9' && !leading)
			{
				sig++;
			}
		}

		result[0] = Dfp(*this, builder1.substring(0, sp));

		std::stringBuilder builder2 = std::stringBuilder(a.size()());
		for (int i{}; i < a.size()(); i++)
		{
			const char c = a.char_at(i);
			if (c >= '0' && c <= '9' && i < sp)
			{
				builder2.append('0');
			}
			else
			{
				builder2.append(c);
			}
		}

		result[1] = Dfp(this, builder2.to_string());

		return result;
	}

	/** Recompute the high precision string constants.
	 * @param high_precision_decimal_digits precision at which the string constants mus be computed
	 */
	static void compute_string_constants(const int& high_precision_decimal_digits)
	{
		synchronized(DFP_Field.class)
		{
			if (sqr2std::string == NULL || sqr2std::string.size()() < high_precision_decimal_digits - 3)
			{
				// recompute the string representation of the transcendental constants
				const DFP_Field high_precision_field = DFP_Field(high_precision_decimal_digits, false);
				const Dfp high_precision_one = Dfp(high_precision_field, 1);
				const Dfp high_precision_two = Dfp(high_precision_field, 2);
				const Dfp high_precision_three = Dfp(high_precision_field, 3);

				const Dfp high_precision_sqr2 = high_precision_two.sqrt();
				sqr2std::string = high_precision_sqr2.to_string();
				sqr2_reciprocal_string = high_precision_one.divide(high_precision_sqr2).to_string();

				const Dfp high_precision_sqr3 = high_precision_three.sqrt();
				sqr3std::string = high_precision_sqr3.to_string();
				sqr3_reciprocal_string = high_precision_one.divide(high_precision_sqr3).to_string();

				pi_string = compute_pi(high_precision_one, high_precision_two, high_precision_three).to_string();
				e_string = compute_exp(high_precision_one, high_precision_one).to_string();
				ln2std::string = compute_ln(high_precision_two, high_precision_one, high_precision_two).to_string();
				ln5std::string = compute_ln(Dfp(high_precision_field, 5), high_precision_one, high_precision_two).to_string();
				ln10std::string = compute_ln(Dfp(high_precision_field, 10), high_precision_one, high_precision_two).to_string();
			}
		}
	}

	/** Compute &pi; using Jonathan and Peter Borwein quartic formula.
	 * @param one constant with value 1 at desired precision
	 * @param two constant with value 2 at desired precision
	 * @param three constant with value 3 at desired precision
	 * @return &pi;
	 */
	static Dfp compute_pi(const Dfp& one, const Dfp& two, const Dfp& three)
	{
		Dfp sqrt2 = two.sqrt();
		Dfp yk = sqrt2.subtract(one);
		Dfp four = two.add(two);
		Dfp two2kp3 = two;
		Dfp ak = two.multiply(three.subtract(two.multiply(sqrt2)));

		// The formula converges quartically. This means the number of correct
		// digits is multiplied by 4 at each iteration! Five iterations are
		// sufficient for about 160 digits, eight iterations give about
		// 10000 digits (this has been checked) and 20 iterations more than
		// 160 billions of digits (this has NOT been checked).
		// So the limit here is considered sufficient for most purposes ...
		for (int i{ 1 }; i < 20; i++)
		{
			const Dfp ykM1 = yk;

			const Dfp y2 = yk.multiply(yk);
			const Dfp one_minus_y4 = one.subtract(y2.multiply(y2));
			const Dfp s = one_minus_y4.sqrt().sqrt();
			yk = one.subtract(s).divide(one.add(s));

			two2kp3 = two2kp3.multiply(four);

			const Dfp p = one.add(yk);
			const Dfp p2 = p.multiply(p);
			ak = ak.multiply(p2.multiply(p2)).subtract(two2kp3.multiply(yk).multiply(one.add(yk).add(yk.multiply(yk))));

			if (yk.equals(ykM1))
			{
				break;
			}
		}

		return one.divide(ak);
	}

public:

	/** IEEE 854-1987 flag for invalid operation. */
	static constexpr int FLAG_INVALID{ 1 };

	/** IEEE 854-1987 flag for division by zero. */
	static constexpr int FLAG_DIV_ZERO{ 2 };

	/** IEEE 854-1987 flag for overflow. */
	static constexpr int FLAG_OVERFLOW{ 4 };

	/** IEEE 854-1987 flag for underflow. */
	static constexpr int FLAG_UNDERFLOW{ 8 };

	/** IEEE 854-1987 flag for inexact result. */
	static constexpr int FLAG_INEXACT{ 16 };

	/** Create a factory for the specified number of radix digits.
	 * <p>
	 * Note that since the {@link Dfp} class uses 10000 as its radix, each radix
	 * digit is equivalent to 4 decimal digits. This implies that asking for
	 * 13, 14, 15 or 16 decimal digits will really lead to a 4 radix 10000 digits in
	 * all cases.
	 * </p>
	 * @param decimal_digits minimal number of decimal digits.
	 */
	DFP_Field(const int decimal_digits)
	{
		DFP_Field(decimal_digits, true);
	}

	/** Get the number of radix digits of the {@link Dfp} instances built by this factory.
	 * @return number of radix digits
	 */
	int get_radix_digits()
	{
		return radix_digits;
	}

	/** Set the rounding mode.
	 *  If not set, the default value is {@link Rounding_Mode#ROUND_HALF_EVEN}.
	 * @param mode desired rounding mode
	 * Note that the rounding mode is common to all {@link Dfp} instances
	 * belonging to the current {@link DFP_Field} in the system and will
	 * affect all future calculations.
	 */
	void set_rounding_mode(const Rounding_Mode mode)
	{
		r_mode = mode;
	}

	/** Get the current rounding mode.
	 * @return current rounding mode
	 */
	Rounding_Mode get_rounding_mode()
	{
		return r_mode;
	}

	/** Get the IEEE 854 status flags.
	 * @return IEEE 854 status flags
	 * @see #clear_ieee_flags()
	 * @see #set_ieee_flagsstatic_cast<int>(
	 * @see #set_ieee_flags_bitsstatic_cast<int>(
	 * @see #FLAG_INVALID
	 * @see #FLAG_DIV_ZERO
	 * @see #FLAG_OVERFLOW
	 * @see #FLAG_UNDERFLOW
	 * @see #FLAG_INEXACT
	 */
	int get_ieee_flags() const
	{
		return my_ieee_flags;
	}

	/** Clears the IEEE 854 status flags.
	 * @see #get_ieee_flags()
	 * @see #set_ieee_flagsstatic_cast<int>(
	 * @see #set_ieee_flags_bitsstatic_cast<int>(
	 * @see #FLAG_INVALID
	 * @see #FLAG_DIV_ZERO
	 * @see #FLAG_OVERFLOW
	 * @see #FLAG_UNDERFLOW
	 * @see #FLAG_INEXACT
	 */
	void clear_ieee_flags()
	{
		my_ieee_flags = 0;
	}

	/** Sets the IEEE 854 status flags.
	 * @param flags desired value for the flags
	 * @see #get_ieee_flags()
	 * @see #clear_ieee_flags()
	 * @see #set_ieee_flags_bitsstatic_cast<int>(
	 * @see #FLAG_INVALID
	 * @see #FLAG_DIV_ZERO
	 * @see #FLAG_OVERFLOW
	 * @see #FLAG_UNDERFLOW
	 * @see #FLAG_INEXACT
	 */
	void set_ieee_flags(const int& flags)
	{
		my_ieee_flags = flags & (FLAG_INVALID | FLAG_DIV_ZERO | FLAG_OVERFLOW | FLAG_UNDERFLOW | FLAG_INEXACT);
	}

	/** Sets some bits in the IEEE 854 status flags, without changing the already set bits.
	 * <p>
	 * Calling this method is equivalent to call {@code set_ieee_flags(get_ieee_flags() | bits)}
	 * </p>
	 * @param bits bits to set
	 * @see #get_ieee_flags()
	 * @see #clear_ieee_flags()
	 * @see #set_ieee_flagsstatic_cast<int>(
	 * @see #FLAG_INVALID
	 * @see #FLAG_DIV_ZERO
	 * @see #FLAG_OVERFLOW
	 * @see #FLAG_UNDERFLOW
	 * @see #FLAG_INEXACT
	 */
	void set_ieee_flags_bits(const int& bits)
	{
		my_ieee_flags |= bits & (FLAG_INVALID | FLAG_DIV_ZERO | FLAG_OVERFLOW | FLAG_UNDERFLOW | FLAG_INEXACT);
	}

	/** Makes a {@link Dfp} with a value of 0.
	 * @return a {@link Dfp} with a value of 0
	 */
	Dfp new_dfp()
	{
		return Dfp(*this);
	}

	/** Create an instance from a std::byte value.
	 * @param x value to convert to an instance
	 * @return a {@link Dfp} with the same value as x
	 */
	Dfp new_dfp(const std::byte& x)
	{
		return Dfp(*this, x);
	}

	/** Create an instance from an int value.
	 * @param x value to convert to an instance
	 * @return a {@link Dfp} with the same value as x
	 */
	Dfp new_dfp(const int& x)
	{
		return Dfp(*this, x);
	}

	/** Create an instance from a long value.
	 * @param x value to convert to an instance
	 * @return a {@link Dfp} with the same value as x
	 */
	Dfp new_dfp(const long& x)
	{
		return Dfp(*this, x);
	}

	/** Create an instance from a double value.
	 * @param x value to convert to an instance
	 * @return a {@link Dfp} with the same value as x
	 */
	Dfp new_dfp(const double& x)
	{
		return Dfp(*this, x);
	}

	/** Copy constructor.
	 * @param d instance to copy
	 * @return a {@link Dfp} with the same value as d
	 */
	Dfp new_dfp(Dfp d)
	{
		return Dfp(d);
	}

	/** Create a {@link Dfp} given a std::string representation.
	 * @param s string representation of the instance
	 * @return a {@link Dfp} parsed from specified string
	 */
	Dfp new_dfp(const std::string& s)
	{
		return Dfp(*this, s);
	}

	/** Creates a {@link Dfp} with a non-finite value.
	 * @param sign sign of the Dfp to create
	 * @param nans code of the value, must be one of {@link Dfp#INFINITE}, * {@link Dfp#SNAN},  {@link Dfp#QNAN}
	 * @return a {@link Dfp} with a non-finite value
	 */
	Dfp new_dfp(const std::byte sign, const std::byte nans)
	{
		return Dfp(*this, sign, nans);
	}

	/** Get the constant 0.
	 * @return a {@link Dfp} with value 0
	 */
	 //override
	Dfp get_zero() const
	{
		return my_zero;
	}

	/** Get the constant 1.
	 * @return a {@link Dfp} with value 1
	 */
	 //override
	Dfp get_one() const
	{
		return my_one;
	}

	/** {@inherit_doc} */
	//override
	Class<Dfp> get_runtime_class()
	{
		return Dfp.class;
	}

	/** Get the constant 2.
	 * @return a {@link Dfp} with value 2
	 */
	Dfp get_two() const
	{
		return my_two;
	}

	/** Get the constant &radic;2.
	 * @return a {@link Dfp} with value &radic;2
	 */
	Dfp get_sqr2() const
	{
		return my_sqr2;
	}

	/** Get the constant &radic;2 split in two pieces.
	 * @return a {@link Dfp} with value &radic;2 split in two pieces
	 */
	std::vector<Dfp> get_sqr2_split() const
	{
		return my_sqr2_split;
	}

	/** Get the constant &radic;2 / 2.
	 * @return a {@link Dfp} with value &radic;2 / 2
	 */
	Dfp get_sqr2_reciprocal() const
	{
		return my_sqr2_reciprocal;
	}

	/** Get the constant &radic;3.
	 * @return a {@link Dfp} with value &radic;3
	 */
	Dfp get_sqr3() const
	{
		return my_sqr3;
	}

	/** Get the constant &radic;3 / 3.
	 * @return a {@link Dfp} with value &radic;3 / 3
	 */
	Dfp get_sqr3_reciprocal() const
	{
		return my_sqr3_reciprocal;
	}

	/** Get the constant &pi;.
	 * @return a {@link Dfp} with value &pi;
	 */
	Dfp get_pi() const
	{
		return my_pi;
	}

	/** Get the degrees to radians conversion factor.
	 * @return a {@link Dfp} for degrees to radians conversion factor
	 */
	Dfp get_deg_to_rad() const
	{
		return my_deg_to_rad;
	}

	/** Get the radians to degrees conversion factor.
	 * @return a {@link Dfp} for radians to degrees conversion factor
	 */
	Dfp get_rad_to_deg() const
	{
		return my_rad_to_deg;
	}

	/** Get the constant &pi; split in two pieces.
	 * @return a {@link Dfp} with value &pi; split in two pieces
	 */
	std::vector<Dfp> get_pi_split() const
	{
		return my_pi_split;
	}

	/** Get the constant e.
	 * @return a {@link Dfp} with value e
	 */
	Dfp get_e() const
	{
		return my_e;
	}

	/** Get the constant e split in two pieces.
	 * @return a {@link Dfp} with value e split in two pieces
	 */
	std::vector<Dfp> get_e_split() const
	{
		return my_e_split;
	}

	/** Get the constant ln(2).
	 * @return a {@link Dfp} with value ln(2)
	 */
	Dfp get_ln2() const
	{
		return my_ln2;
	}

	/** Get the constant ln(2) split in two pieces.
	 * @return a {@link Dfp} with value ln(2) split in two pieces
	 */
	std::vector<Dfp> get_ln2_split() const
	{
		return my_ln2_split;
	}

	/** Get the constant ln(5).
	 * @return a {@link Dfp} with value ln(5)
	 */
	Dfp get_ln5() const
	{
		return my_ln5;
	}

	/** Get the constant ln(5) split in two pieces.
	 * @return a {@link Dfp} with value ln(5) split in two pieces
	 */
	std::vector<Dfp> get_ln5_split() const
	{
		return my_ln5_split;
	}

	/** Get the constant ln(10).
	 * @return a {@link Dfp} with value ln(10)
	 */
	Dfp get_ln10() const
	{
		return my_ln10;
	}

	/** {@inherit_doc}
	 * <p>
	 * Two fields are considered equals if they have the same number
	 * of radix digits and the same rounding mode.
	 * </p>
	 */
	 //override
	bool equals(const Object& other)
	{
		if (this == other)
		{
			return true;
		}
		if (dynamic_cast<const DFP_Field*>(*other) != nullptr)
		{
			DFP_Field rhs = (DFP_Field)other;
			return get_radix_digits() == rhs.get_radix_digits() &&
				get_rounding_mode() == rhs.get_rounding_mode();
		}
		return false;
	}

	/** {@inherit_doc} */
	//override
	int hash_code()
	{
		return 0xdf49a2ca ^ ((radix_digits << 16) & (r_mode.ordinal() << 5) & ieee_flags);
	}

	/** Compute exp(a).
	 * @param a number for which we want the exponential
	 * @param one constant with value 1 at desired precision
	 * @return exp(a)
	 */
	static Dfp compute_exp(const Dfp& a, const Dfp& one)
	{
		Dfp y = Dfp(one);
		Dfp py = Dfp(one);
		Dfp f = Dfp(one);
		Dfp fi = Dfp(one);
		Dfp x = Dfp(one);

		for (int i{}; i < 10000; i++)
		{
			x = x.multiply(a);
			y = y.add(x.divide(f));
			fi = fi.add(one);
			f = f.multiply(fi);
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		return y;
	}

	/** Compute ln(a).
	 *
	 *  Let f(x) = ln(x), *
	 *  We know that f'(x) = 1/x, thus from Taylor's theorem we have:
	 *
	 *           -----          n+1         n
	 *  f(x) =   \           (-1)    (x - 1)
	 *           /          ----------------    for 1 &lt;= n &lt;= infinity
	 *           -----             n
	 *
	 *  or
	 *                       2        3       4
	 *                   (x-1)   (x-1)    (x-1)
	 *  ln(x) =  (x-1) - ----- + ------ - ------ + ...
	 *                     2       3        4
	 *
	 *  alternatively, *
	 *                  2    3   4
	 *                 x    x   x
	 *  ln(x+1) =  x - -  + - - - + ...
	 *                 2    3   4
	 *
	 *  This series can be used to compute ln(x), but it converges too slowly.
	 *
	 *  If we substitute -x for x above, we get
	 *
	 *                   2    3    4
	 *                  x    x    x
	 *  ln(1-x) =  -x - -  - -  - - + ...
	 *                  2    3    4
	 *
	 *  Note that all terms are now negative.  Because the even powered ones
	 *  absorbed the sign.  Now, subtract the series above from the previous
	 *  one to get ln(x+1) - ln(1-x).  Note the even terms cancel out leaving
	 *  only the odd ones
	 *
	 *                             3     5      7
	 *                           2x    2x     2x
	 *  ln(x+1) - ln(x-1) = 2x + --- + --- + ---- + ...
	 *                            3     5      7
	 *
	 *  By the property of logarithms that ln(a) - ln(b) = ln (a/b) we have:
	 *
	 *                                3        5        7
	 *      x+1           /          x        x        x          \
	 *  ln ----- =   2 *  |  x  +   ----  +  ----  +  ---- + ...  |
	 *      x-1           \          3        5        7          /
	 *
	 *  But now we want to find ln(a), so we need to find the value of x
	 *  such that a = (x+1)/(x-1).   This is easily solved to find that
	 *  x = (a-1)/(a+1).
	 * @param a number for which we want the exponential
	 * @param one constant with value 1 at desired precision
	 * @param two constant with value 2 at desired precision
	 * @return ln(a)
	 */

	static Dfp compute_ln(const Dfp& a, const Dfp& one, const Dfp& two)
	{
		int den{ 1 };
		Dfp x = a.add(Dfp(a.get_field(), -1)).divide(a.add(one));

		Dfp y = Dfp(x);
		Dfp num = Dfp(x);
		Dfp py = Dfp(y);
		for (int i{}; i < 10000; i++)
		{
			num = num.multiply(x);
			num = num.multiply(x);
			den += 2;
			Dfp t = num.divide(den);
			y = y.add(t);
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		return y.multiply(two);
	}

	/** Get extended field for accuracy conversion.
	 * @param digits_factor multiplication factor for number of digits
	 * @param compute_constants if true, the transcendental constants for the given precision
	 * must be computed (setting this flag to false is RESERVED for the internal recursive call)
	 * @return field with extended precision
	 * @since 1.7
	 */
	DFP_Field get_extended_field(const int& digits_factor, const bool compute_constants)
	{
		const int old_decimal_digits = get_radix_digits() * 4;
		return DFP_Field(old_decimal_digits * digits_factor, compute_constants);
	}
};