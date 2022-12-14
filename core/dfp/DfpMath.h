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

  /** Mathematical routines for use with {@link Dfp}.
   * The constants are defined in {@link DFP_Field}
   */
class Dfp_Math
{
private:
	/** Name for traps triggered by pow. */
	static const std::string POW_TRAP{ "pow" };

	/**
	 * Private Constructor.
	 */
	Dfp_Math() = default;

protected:
	/** Breaks a string representation up into two dfp's.
	 * <p>The two dfp are such that the sum of them is equivalent
	 * to the input string, but has higher precision than using a
	 * single dfp. This is useful for improving accuracy of
	 * exponentiation and critical multiplies.
	 * @param field field to which the Dfp must belong
	 * @param a string representation to split
	 * @return an array of two {@link Dfp} which sum is a
	 */
	static std::vector<Dfp> split(const DFP_Field& field, const std::string& a)
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

			if (sig == (field.get_radix_digits() / 2) * 4)
			{
				sp = i;
				break;
			}

			if (c >= '0' && c <= '9' && !leading)
			{
				sig++;
			}
		}

		result[0] = field.new_dfp(builder1.substring(0, sp));

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

		result[1] = field.new_dfp(builder2.to_string());

		return result;
	}

	/** Splits a {@link Dfp} into 2 {@link Dfp}'s such that their sum is equal to the input {@link Dfp}.
	 * @param a number to split
	 * @return two elements array containing the split number
	 */
	static std::vector<Dfp> split(const Dfp& a)
	{
		auto result = std::vector<Dfp>(2);
		const Dfp shift = a.multiply(a.power10_k(a.get_radix_digits() / 2));
		result[0] = a.add(shift).subtract(shift);
		result[1] = a.subtract(result[0]);
		return result;
	}

	/** Multiply two numbers that are split in to two pieces that are
	 *  meant to be added together.
	 *  Use binomial multiplication so ab = a0 b0 + a0 b1 + a1 b0 + a1 b1
	 *  Store the first term in result0, the rest in result1
	 *  @param a first factor of the multiplication, in split form
	 *  @param b second factor of the multiplication, in split form
	 *  @return a &times; b, in split form
	 */
	static std::vector<Dfp> split_mult(const std::vector<Dfp>& a, const std::vector<Dfp>& b)
	{
		auto result = std::vector<Dfp>(2);

		result[1] = a[0].get_zero();
		result[0] = a[0].multiply(b[0]);

		/* If result[0] is infinite or zero, don't compute result[1].
		 * Attempting to do so may produce NaNs.
		 */

		if (result[0].classify() == Dfp.INFINITE || result[0].equals(result[1]))
		{
			return result;
		}

		result[1] = a[0].multiply(b[1]).add(a[1].multiply(b[0])).add(a[1].multiply(b[1]));

		return result;
	}

	/** Divide two numbers that are split in to two pieces that are meant to be added together.
	 * Inverse of split multiply above:
	 *  (a+b) / (c+d) = (a/c) + ( (bc-ad)/(c**2+cd) )
	 *  @param a dividend, in split form
	 *  @param b divisor, in split form
	 *  @return a / b, in split form
	 */
	static std::vector<Dfp> split_div(const std::vector<Dfp>& a, const std::vector<Dfp>& b)
	{
		std::vector<Dfp> result;

		result = std::vector<Dfp>(2);

		result[0] = a[0].divide(b[0]);
		result[1] = a[1].multiply(b[0]).subtract(a[0].multiply(b[1]));
		result[1] = result[1].divide(b[0].multiply(b[0]).add(b[0].multiply(b[1])));

		return result;
	}

	/** Raise a split base to the a power.
	 * @param base number to raise
	 * @param a power
	 * @return base<sup>a</sup>
	 */
	static Dfp split_pow(const std::vector<Dfp>& base, const int& a)
	{
		bool invert{};

		auto r = std::vector<Dfp>(2);

		auto result = std::vector<Dfp>(2);
		result[0] = base[0].get_one();
		result[1] = base[0].get_zero();

		if (a == 0)
		{
			// Special case a = 0
			return result[0].add(result[1]);
		}

		if (a < 0)
		{
			// If a is less than zero
			invert = true;
			a = -a;
		}

		// Exponentiate by successive squaring
		do
		{
			r[0] = Dfp(base[0]);
			r[1] = Dfp(base[1]);
			int trial = 1;

			int prevtrial;
			while (true)
			{
				prevtrial = trial;
				trial *= 2;
				if (trial > a)
				{
					break;
				}
				r = split_mult(r, r);
			}

			trial = prevtrial;

			a -= trial;
			result = split_mult(result, r);
		} while (a >= 1);

		result[0] = result[0].add(result[1]);

		if (invert)
		{
			result[0] = base[0].get_one().divide(result[0]);
		}

		return result[0];
	}

	/** Computes e to the given power.
	 * Where -1 &lt; a &lt; 1.  Use the classic Taylor series.  1 + x**2/2! + x**3/3! + x**4/4!  ...
	 * @param a power at which e should be raised
	 * @return e<sup>a</sup>
	 */
	static Dfp exp_internal(const Dfp& a)
	{
		Dfp y = a.get_one();
		Dfp x = a.get_one();
		Dfp fact = a.get_one();
		Dfp py = Dfp(y);

		for (int i{ 1 }; i < 90; i++)
		{
			x = x.multiply(a);
			fact = fact.divide(i);
			y = y.add(x.multiply(fact));
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		return y;
	}

	/** Computes the natural log of a number between 0 and 2.
	 *  Let f(x) = ln(x), *
	 *  We know that f'(x) = 1/x, thus from Taylor's theorum we have:
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
	 * @param a number from which logarithm is requested, in split form
	 * @return log(a)
	 */
	static std::vector<Dfp> log_internal(const std::vector<Dfp>& a)
	{
		/* Now we want to compute x = (a-1)/(a+1) but this is prone to
		 * loss of precision.  So instead, compute x = (a/4 - 1/4) / (a/4 + 1/4)
		 */
		Dfp t = a[0].divide(4).add(a[1].divide(4));
		Dfp x = t.add(a[0].new_instance("-0.25")).divide(t.add(a[0].new_instance("0.25")));

		Dfp y = Dfp(x);
		Dfp num = Dfp(x);
		Dfp py = Dfp(y);
		int den = 1;
		for (int i{}; i < 10000; i++)
		{
			num = num.multiply(x);
			num = num.multiply(x);
			den += 2;
			t = num.divide(den);
			y = y.add(t);
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		y = y.multiply(a[0].get_two());

		return split(y);
	}

	/** Computes sin(a)  Used when 0 &lt; a &lt; pi/4.
	 * Uses the classic Taylor series.  x - x**3/3! + x**5/5!  ...
	 * @param a number from which sine is desired, in split form
	 * @return sin(a)
	 */
	static Dfp sin_internal(const std::vector<Dfp>& a)
	{
		Dfp c = a[0].add(a[1]);
		Dfp y = c;
		c = c.multiply(c);
		Dfp x = y;
		Dfp fact = a[0].get_one();
		Dfp py = Dfp(y);

		for (int i = 3; i < 90; i += 2)
		{
			x = x.multiply(c);
			x = x.negate();

			fact = fact.divide((i - 1) * i);  // 1 over fact
			y = y.add(x.multiply(fact));
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		return y;
	}

	/** Computes cos(a)  Used when 0 &lt; a &lt; pi/4.
	 * Uses the classic Taylor series for cosine.  1 - x**2/2! + x**4/4!  ...
	 * @param a number from which cosine is desired, in split form
	 * @return cos(a)
	 */
	static Dfp cos_internal(const std::vector<Dfp>& a)
	{
		const Dfp one = a[0].get_one();

		Dfp x = one;
		Dfp y = one;
		Dfp c = a[0].add(a[1]);
		c = c.multiply(c);

		Dfp fact = one;
		Dfp py = Dfp(y);

		for (int i{ 2 }; i < 90; i += 2)
		{
			x = x.multiply(c);
			x = x.negate();

			fact = fact.divide((i - 1) * i);  // 1 over fact

			y = y.add(x.multiply(fact));
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		return y;
	}

	/** computes the arc-tangent of the argument.
	 * @param a number from which arc-tangent is desired
	 * @return atan(a)
	 */
	static Dfp atan_internal(const Dfp& a)
	{
		Dfp y = Dfp(a);
		Dfp x = Dfp(y);
		Dfp py = Dfp(y);

		for (int i{ 3 }; i < 90; i += 2)
		{
			x = x.multiply(a);
			x = x.multiply(a);
			x = x.negate();
			y = y.add(x.divide(i));
			if (y.equals(py))
			{
				break;
			}
			py = Dfp(y);
		}

		return y;
	}

public:

	/** Raises base to the power a by successive squaring.
	 * @param base number to raise
	 * @param a power
	 * @return base<sup>a</sup>
	 */
	static Dfp pow(const Dfp& base, const int& a)
	{
		bool invert{};

		Dfp result = base.get_one();

		if (a == 0)
		{
			// Special case
			return result;
		}

		if (a < 0)
		{
			invert = true;
			a = -a;
		}

		// Exponentiate by successive squaring
		do
		{
			Dfp r = Dfp(base);
			Dfp prevr;
			int trial{ 1 };
			int prevtrial;

			do
			{
				prevr = Dfp(r);
				prevtrial = trial;
				r = r.multiply(r);
				trial *= 2;
			} while (a > trial);

			r = prevr;
			trial = prevtrial;

			a -= trial;
			result = result.multiply(r);
		} while (a >= 1);

		if (invert)
		{
			result = base.get_one().divide(result);
		}

		return base.new_instance(result);
	}

	/** Computes e to the given power.
	 * a is broken into two parts, such that a = n+m  where n is an integer.
	 * We use pow() to compute e<sup>n</sup> and a Taylor series to compute
	 * e<sup>m</sup>.  We return e*<sup>n</sup> &times; e<sup>m</sup>
	 * @param a power at which e should be raised
	 * @return e<sup>a</sup>
	 */
	static Dfp exp(const Dfp& a)
	{
		const Dfp inta = a.rint();
		const Dfp fraca = a.subtract(inta);

		const int ia = inta.int_value();
		if (ia > 2147483646)
		{
			// return +Infinity
			return a.new_instance((byte)1, Dfp.INFINITE);
		}

		if (ia < -2147483646)
		{
			// return 0;
			return a.new_instance();
		}

		const Dfp einta = split_pow(a.get_field().get_e_split(), ia);
		const Dfp efraca = exp_internal(fraca);

		return einta.multiply(efraca);
	}

	/** Returns the natural logarithm of a.
	 * a is first split into three parts such that  a = (10000^h)(2^j)k.
	 * ln(a) is computed by ln(a) = ln(5)*h + ln(2)*(h+j) + ln(k)
	 * k is in the range 2/3 &lt; k &lt; 4/3 and is passed on to a series expansion.
	 * @param a number from which logarithm is requested
	 * @return log(a)
	 */
	static Dfp log(const Dfp& a)
	{
		int lr;
		Dfp x;
		int ix;
		int p2 = 0;

		// Check the arguments somewhat here
		if (a.equals(a.get_zero()) || a.less_than(a.get_zero()) || a.is_nan())
		{
			// negative, zero or NaN
			a.get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			return a.dotrap(DFP_Field.FLAG_INVALID, "ln", a, a.new_instance((byte)1, Dfp.QNAN));
		}

		if (a.classify() == Dfp.INFINITE)
		{
			return a;
		}

		x = Dfp(a);
		lr = x.log10_k();

		x = x.divide(pow(a.new_instance(10000), lr));  /* This puts x in the range 0-10000 */
		ix = x.floor().int_value();

		while (ix > 2)
		{
			ix >>= 1;
			p2++;
		}

		std::vector<Dfp> spx = split(x);
		auto spy = std::vector<Dfp>(2);
		spy[0] = pow(a.get_two(), p2);          // use spy[0] temporarily as a divisor
		spx[0] = spx[0].divide(spy[0]);
		spx[1] = spx[1].divide(spy[0]);

		spy[0] = a.new_instance("1.33333");    // Use spy[0] for comparison
		while (spx[0].add(spx[1]).greater_than(spy[0]))
		{
			spx[0] = spx[0].divide(2);
			spx[1] = spx[1].divide(2);
			p2++;
		}

		// X is now in the range of 2/3 < x < 4/3
		Dfp[] spz = log_internal(spx);

		spx[0] = a.new_instance(std::stringBuilder().append(p2 + 4 * lr).to_string());
		spx[1] = a.get_zero();
		spy = split_mult(a.get_field().get_ln2_split(), spx);

		spz[0] = spz[0].add(spy[0]);
		spz[1] = spz[1].add(spy[1]);

		spx[0] = a.new_instance(std::stringBuilder().append(4 * lr).to_string());
		spx[1] = a.get_zero();
		spy = split_mult(a.get_field().get_ln5_split(), spx);

		spz[0] = spz[0].add(spy[0]);
		spz[1] = spz[1].add(spy[1]);

		return a.new_instance(spz[0].add(spz[1]));
	}

	/** Computes x to the y power.<p>
	 *
	 *  Uses the following method:<p>
	 *
	 *  <ol>
	 *  <li> Set u = rint(y), v = y-u
	 *  <li> Compute a = v * ln(x)
	 *  <li> Compute b = rint( a/ln(2) )
	 *  <li> Compute c = a - b*ln(2)
	 *  <li> x<sup>y</sup> = x<sup>u</sup>  *   2<sup>b</sup> * e<sup>c</sup>
	 *  </ol>
	 *  if |y| &gt; 1e8, then we compute by exp(y*ln(x))   <p>
	 *
	 *  <b>Special Cases</b><p>
	 *  <ul>
	 *  <li>  if y is 0.0 or -0.0 then result is 1.0
	 *  <li>  if y is 1.0 then result is x
	 *  <li>  if y is NaN then result is NaN
	 *  <li>  if x is NaN and y is not zero then result is NaN
	 *  <li>  if |x| &gt; 1.0 and y is +Infinity then result is +Infinity
	 *  <li>  if |x| &lt; 1.0 and y is -Infinity then result is +Infinity
	 *  <li>  if |x| &gt; 1.0 and y is -Infinity then result is +0
	 *  <li>  if |x| &lt; 1.0 and y is +Infinity then result is +0
	 *  <li>  if |x| = 1.0 and y is +/-Infinity then result is NaN
	 *  <li>  if x = +0 and y &gt; 0 then result is +0
	 *  <li>  if x = +Inf and y &lt; 0 then result is +0
	 *  <li>  if x = +0 and y &lt; 0 then result is +Inf
	 *  <li>  if x = +Inf and y &gt; 0 then result is +Inf
	 *  <li>  if x = -0 and y &gt; 0, finite, not odd integer then result is +0
	 *  <li>  if x = -0 and y &lt; 0, finite, and odd integer then result is -Inf
	 *  <li>  if x = -Inf and y &gt; 0, finite, and odd integer then result is -Inf
	 *  <li>  if x = -0 and y &lt; 0, not finite odd integer then result is +Inf
	 *  <li>  if x = -Inf and y &gt; 0, not finite odd integer then result is +Inf
	 *  <li>  if x &lt; 0 and y &gt; 0, finite, and odd integer then result is -(|x|<sup>y</sup>)
	 *  <li>  if x &lt; 0 and y &gt; 0, finite, and not integer then result is NaN
	 *  </ul>
	 *  @param x base to be raised
	 *  @param y power to which base should be raised
	 *  @return x<sup>y</sup>
	 */
	static Dfp pow(const Dfp& x, const Dfp& y)
	{
		// make sure we don't mix number with different precision
		if (x.get_field().get_radix_digits() != y.get_field().get_radix_digits())
		{
			x.get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = x.new_instance(x.get_zero());
			result.nans = Dfp.QNAN;
			return x.dotrap(DFP_Field.FLAG_INVALID, POW_TRAP, x, result);
		}

		const Dfp zero = x.get_zero();
		const Dfp one = x.get_one();
		const Dfp two = x.get_two();
		bool invert = false;
		int ui;

		/* Check for special cases */
		if (y.equals(zero))
		{
			return x.new_instance(one);
		}

		if (y.equals(one))
		{
			if (x.is_nan())
			{
				// Test for NaNs
				x.get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
				return x.dotrap(DFP_Field.FLAG_INVALID, POW_TRAP, x, x);
			}
			return x;
		}

		if (x.is_nan() || y.is_nan())
		{
			// Test for NaNs
			x.get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			return x.dotrap(DFP_Field.FLAG_INVALID, POW_TRAP, x, x.new_instance((byte)1, Dfp.QNAN));
		}

		// X == 0
		if (x.equals(zero))
		{
			if (Dfp.copysign(one, x).greater_than(zero))
			{
				// X == +0
				if (y.greater_than(zero))
				{
					return x.new_instance(zero);
				}
				else
				{
					return x.new_instance(x.new_instance((byte)1, Dfp.INFINITE));
				}
			}
			else
			{
				// X == -0
				if (y.classify() == Dfp.FINITE && y.rint().equals(y) && !y.remainder(two).equals(zero))
				{
					// If y is odd integer
					if (y.greater_than(zero))
					{
						return x.new_instance(zero.negate());
					}
					else
					{
						return x.new_instance(x.new_instance((byte)-1, Dfp.INFINITE));
					}
				}
				else
				{
					// Y is not odd integer
					if (y.greater_than(zero))
					{
						return x.new_instance(zero);
					}
					else
					{
						return x.new_instance(x.new_instance((byte)1, Dfp.INFINITE));
					}
				}
			}
		}

		if (x.less_than(zero))
		{
			// Make x positive, but keep track of it
			x = x.negate();
			invert = true;
		}

		if (x.greater_than(one) && y.classify() == Dfp.INFINITE)
		{
			if (y.greater_than(zero))
			{
				return y;
			}
			else
			{
				return x.new_instance(zero);
			}
		}

		if (x.less_than(one) && y.classify() == Dfp.INFINITE)
		{
			if (y.greater_than(zero))
			{
				return x.new_instance(zero);
			}
			else
			{
				return x.new_instance(Dfp.copysign(y, one));
			}
		}

		if (x.equals(one) && y.classify() == Dfp.INFINITE)
		{
			x.get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			return x.dotrap(DFP_Field.FLAG_INVALID, POW_TRAP, x, x.new_instance((byte)1, Dfp.QNAN));
		}

		if (x.classify() == Dfp.INFINITE)
		{
			// x = +/- inf
			if (invert)
			{
				// negative infinity
				if (y.classify() == Dfp.FINITE && y.rint().equals(y) && !y.remainder(two).equals(zero))
				{
					// If y is odd integer
					return y.greater_than(zero)
						? x.new_instance(x.new_instance((byte)-1, Dfp.INFINITE))
						: x.new_instance(zero.negate());
				}
				else
				{
					// Y is not odd integer
					return y.greater_than(zero)
						? x.new_instance(x.new_instance((byte)1, Dfp.INFINITE))
						: x.new_instance(zero);
				}
			}
			else
			{
				// positive infinity
				return y.greater_than(zero)
					? x
					: x.new_instance(zero);
			}
		}

		if (invert && !y.rint().equals(y))
		{
			x.get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			return x.dotrap(DFP_Field.FLAG_INVALID, POW_TRAP, x, x.new_instance((byte)1, Dfp.QNAN));
		}

		// End special cases

		Dfp r;
		if (y.less_than(x.new_instance(100000000)) && y.greater_than(x.new_instance(-100000000)))
		{
			const Dfp u = y.rint();
			ui = u.int_value();

			const Dfp v = y.subtract(u);

			if (v.unequal(zero))
			{
				const Dfp a = v.multiply(log(x));
				const Dfp b = a.divide(x.get_field().get_ln2()).rint();

				const Dfp c = a.subtract(b.multiply(x.get_field().get_ln2()));
				r = split_pow(split(x), ui);
				r = r.multiply(pow(two, b.int_value()));
				r = r.multiply(exp(c));
			}
			else
			{
				r = split_pow(split(x), ui);
			}
		}
		else
		{
			// very large exponent.  |y| > 1e8
			r = exp(log(x).multiply(y));
		}

		if (invert && y.rint().equals(y) && !y.remainder(two).equals(zero))
		{
			// if y is odd integer
			r = r.negate();
		}

		return x.new_instance(r);
	}

	/** computes the sine of the argument.
	 * @param a number from which sine is desired
	 * @return sin(a)
	 */
	static Dfp sin(const Dfp& a)
	{
		const Dfp pi = a.get_field().get_pi();
		const Dfp zero = a.get_field().get_zero();
		bool neg = false;

		/* First reduce the argument to the range of +/- PI */
		Dfp x = a.remainder(pi.multiply(2));

		/* if x < 0 then apply identity sin(-x) = -sin(x) */
		/* This puts x in the range 0 < x < PI            */
		if (x.less_than(zero))
		{
			x = x.negate();
			neg = true;
		}

		/* sin_ce sine(x) = sine(pi - x) we can reduce the range to
		 * 0 < x < pi/2
		 */

		if (x.greater_than(pi.divide(2)))
		{
			x = pi.subtract(x);
		}

		Dfp y;
		if (x.less_than(pi.divide(4)))
		{
			y = sin_internal(split(x));
		}
		else
		{
			auto c = std::vector<Dfp>(2);
			const std::vector<Dfp>& pi_split = a.get_field().get_pi_split();
			c[0] = pi_split[0].divide(2).subtract(x);
			c[1] = pi_split[1].divide(2);
			y = cos_internal(c);
		}

		if (neg)
		{
			y = y.negate();
		}

		return a.new_instance(y);
	}

	/** computes the cosine of the argument.
	 * @param a number from which cosine is desired
	 * @return cos(a)
	 */
	static Dfp cos(const Dfp& a)
	{
		const Dfp pi = a.get_field().get_pi();
		const Dfp zero = a.get_field().get_zero();
		bool neg = false;

		/* First reduce the argument to the range of +/- PI */
		Dfp x = a.remainder(pi.multiply(2));

		/* if x < 0 then apply identity cos(-x) = cos(x) */
		/* This puts x in the range 0 < x < PI           */
		if (x.less_than(zero))
		{
			x = x.negate();
		}

		/* sin_ce cos(x) = -cos(pi - x) we can reduce the range to
		 * 0 < x < pi/2
		 */

		if (x.greater_than(pi.divide(2)))
		{
			x = pi.subtract(x);
			neg = true;
		}

		Dfp y;
		if (x.less_than(pi.divide(4)))
		{
			auto c = std::vector<Dfp>(2);
			c[0] = x;
			c[1] = zero;

			y = cos_internal(c);
		}
		else
		{
			auto c = std::vector<Dfp>(2);
			const std::vector<Dfp>& pi_split = a.get_field().get_pi_split();
			c[0] = pi_split[0].divide(2).subtract(x);
			c[1] = pi_split[1].divide(2);
			y = sin_internal(c);
		}

		if (neg)
		{
			y = y.negate();
		}

		return a.new_instance(y);
	}

	/** computes the tangent of the argument.
	 * @param a number from which tangent is desired
	 * @return tan(a)
	 */
	static Dfp tan(const Dfp& a)
	{
		return sin(a).divide(cos(a));
	}

	/** computes the arc tangent of the argument
	 *
	 *  Uses the typical taylor series
	 *
	 *  but may reduce arguments using the following identity
	 * tan(x+y) = (tan(x) + tan(y)) / (1 - tan(x)*tan(y))
	 *
	 * since tan(PI/8) = sqrt(2)-1, *
	 * atan(x) = atan( (x - sqrt(2) + 1) / (1+x*sqrt(2) - x) + PI/8.0
	 * @param a number from which arc-tangent is desired
	 * @return atan(a)
	 */
	static Dfp atan(const Dfp& a)
	{
		const Dfp zero = a.get_field().get_zero();
		const Dfp one = a.get_field().get_one();
		const std::vector<Dfp>& sqr2_split = a.get_field().get_sqr2_split();
		const std::vector<Dfp>& pi_split = a.get_field().get_pi_split();
		bool recp{};
		bool neg{};
		bool sub{};

		const Dfp ty = sqr2_split[0].subtract(one).add(sqr2_split[1]);

		Dfp x = Dfp(a);
		if (x.less_than(zero))
		{
			neg = true;
			x = x.negate();
		}

		if (x.greater_than(one))
		{
			recp = true;
			x = one.divide(x);
		}

		if (x.greater_than(ty))
		{
			auto sty = std::vector<Dfp>(2);
			sub = true;

			sty[0] = sqr2_split[0].subtract(one);
			sty[1] = sqr2_split[1];

			Dfp[] xs = split(x);

			Dfp[] ds = split_mult(xs, sty);
			ds[0] = ds[0].add(one);

			xs[0] = xs[0].subtract(sty[0]);
			xs[1] = xs[1].subtract(sty[1]);

			xs = split_div(xs, ds);
			x = xs[0].add(xs[1]);

			//x = x.subtract(ty).divide(dfp.one.add(x.multiply(ty)));
		}

		Dfp y = atan_internal(x);

		if (sub)
		{
			y = y.add(pi_split[0].divide(8)).add(pi_split[1].divide(8));
		}

		if (recp)
		{
			y = pi_split[0].divide(2).subtract(y).add(pi_split[1].divide(2));
		}

		if (neg)
		{
			y = y.negate();
		}

		return a.new_instance(y);
	}

	/** computes the arc-sine of the argument.
	 * @param a number from which arc-sine is desired
	 * @return asin(a)
	 */
	static Dfp asin(const Dfp& a)
	{
		return atan(a.divide(a.get_one().subtract(a.multiply(a)).sqrt()));
	}

	/** computes the arc-cosine of the argument.
	 * @param a number from which arc-cosine is desired
	 * @return acos(a)
	 */
	static Dfp acos(const Dfp& a)
	{
		Dfp result;
		bool negative{};

		if (a.less_than(a.get_zero()))
		{
			negative = true;
		}

		a = Dfp.copysign(a, a.get_one());  // absolute value

		result = atan(a.get_one().subtract(a.multiply(a)).sqrt().divide(a));

		if (negative)
		{
			result = a.get_field().get_pi().subtract(result);
		}

		return a.new_instance(result);
	}
};