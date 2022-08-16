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

  //import java.util.Arrays;

  //import org.hipparchus.Calculus_Field_Element;
  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Field_Sin_Cos;
  //import org.hipparchus.util.Field_Sinh_Cosh;
  //import org.hipparchus.util.Math_Utils;

  /**
   *  Decimal floating point library for Java
   *
   *  <p>Another floating point class.  This one is built using radix 10000
   *  which is 10<sup>4</sup>, so its almost decimal.</p>
   *
   *  <p>The design goals here are:
   *  <ol>
   *    <li>Decimal math, or close to it</li>
   *    <li>Settable precision (but no mix between numbers using different settings)</li>
   *    <li>Portability.  Code should be kept as portable as possible.</li>
   *    <li>Performance</li>
   *    <li>Accuracy  - Results should always be +/- 1 ULP for basic
   *         algebraic operation</li>
   *    <li>Comply with IEEE 854-1987 as much as possible.
   *         (See IEEE 854-1987 notes below)</li>
   *  </ol></p>
   *
   *  <p>Trade offs:
   *  <ol>
   *    <li>Memory foot print.  I'm using more memory than necessary to
   *         represent numbers to get better performance.</li>
   *    <li>Digits are bigger, so rounding is a greater loss.  So, if you
   *         really need 12 decimal digits, better use 4 base 10000 digits
   *         there can be one partially filled.</li>
   *  </ol></p>
   *
   *  <p>Numbers are represented  in the following form:
   *  <pre>
   *  n  =  sign &times; mant &times; (radix)<sup>exp</sup>;</p>
   *  </pre>
   *  where sign is &plusmn;1, mantissa represents a fractional number between
   *  zero and one.  mant[0] is the least significant digit.
   *  exp is in the range of -32767 to 32768</p>
   *
   *  <p>IEEE 854-1987  Notes and differences</p>
   *
   *  <p>IEEE 854 requires the radix to be either 2 or 10.  The radix here is
   *  10000, so that requirement is not met, but  it is possible that a
   *  subclassed can be made to make it behave as a radix 10
   *  number.  It is my opinion that if it looks and behaves as a radix
   *  10 number then it is one and that requirement would be met.</p>
   *
   *  <p>The radix of 10000 was chosen because it should be faster to operate
   *  on 4 decimal digits at once instead of one at a time.  Radix 10 behavior
   *  can be realized by adding an additional rounding step to ensure that
   *  the number of decimal digits represented is constant.</p>
   *
   *  <p>The IEEE standard specifically leaves out internal data encoding, *  so it is reasonable to conclude that such a subclass of this radix
   *  10000 system is merely an encoding of a radix 10 system.</p>
   *
   *  <p>IEEE 854 also specifies the existence of "sub-normal" numbers.  This
   *  class does not contain any such entities.  The most significant radix
   *  10000 digit is always non-zero.  Instead, we support "gradual underflow"
   *  by raising the underflow flag for numbers less with exponent less than
   *  exp_min, but don't flush to zero until the exponent reaches MIN_EXP-digits.
   *  Thus the smallest number we can represent would be:
   *  1E(-(MIN_EXP-digits-1)*4),  eg, for digits=5, MIN_EXP=-32767, that would
   *  be 1e-131092.</p>
   *
   *  <p>IEEE 854 defines that the implied radix point lies just to the right
   *  of the most significant digit and to the left of the remaining digits.
   *  This implementation puts the implied radix point to the left of all
   *  digits including the most significant one.  The most significant digit
   *  here is the one just to the right of the radix point.  This is a fine
   *  detail and is really only a matter of definition.  Any side effects of
   *  this can be rendered invisible by a subclass.</p>
   * @see DFP_Field
   */
class Dfp : Calculus_Field_Element<Dfp>
{
	/** The radix, or base of this system.  Set to 10000 */
	public static const int RADIX = 10000;

	/** The minimum exponent before underflow is signaled.  Flush to zero
	 *  occurs at minExp-DIGITS */
	public static const int MIN_EXP = -32767;

	/** The maximum exponent before overflow is signaled and results flushed
	 *  to infinity */
	public static const int MAX_EXP = 32768;

	/** The amount under/overflows are scaled by before going to trap handler */
	public static const int ERR_SCALE = 32760;

	/** Indicator value for normal finite numbers. */
	public static const std::byte FINITE = 0;

	/** Indicator value for Infinity. */
	public static const std::byte INFINITE = 1;

	/** Indicator value for signaling NaN. */
	public static const std::byte SNAN = 2;

	/** Indicator value for quiet NaN. */
	public static const std::byte QNAN = 3;

	/** std::string for NaN representation. */
	private static const std::string NAN_STRING = "NaN";

	/** std::string for positive infinity representation. */
	private static const std::string POS_INFINITY_STRING = "Infinity";

	/** std::string for negative infinity representation. */
	private static const std::string NEG_INFINITY_STRING = "-Infinity";

	/** Name for traps triggered by addition. */
	private static const std::string ADD_TRAP = "add";

	/** Name for traps triggered by multiplication. */
	private static const std::string MULTIPLY_TRAP = "multiply";

	/** Name for traps triggered by division. */
	private static const std::string DIVIDE_TRAP = "divide";

	/** Name for traps triggered by square root. */
	private static const std::string SQRT_TRAP = "sqrt";

	/** Name for traps triggered by alignment. */
	private static const std::string ALIGN_TRAP = "align";

	/** Name for traps triggered by truncation. */
	private static const std::string TRUNC_TRAP = "trunc";

	/** Name for traps triggered by next_after. */
	private static const std::string NEXT_AFTER_TRAP = "next_after";

	/** Name for traps triggered by less_than. */
	private static const std::string LESS_THAN_TRAP = "less_than";

	/** Name for traps triggered by greater_than. */
	private static const std::string GREATER_THAN_TRAP = "greater_than";

	/** Name for traps triggered by new_instance. */
	private static const std::string NEW_INSTANCE_TRAP = "new_instance";

	/** Multiplication factor for number of digits used to compute linear combinations. */
	private static const int LINEAR_COMBINATION_DIGITS_FACTOR = 2;

	/** Mantissa. */
	protected std::vector<int> mant;

	/** Sign bit: 1 for positive, -1 for negative. */
	protected std::byte sign;

	/** Exponent. */
	protected int exp;

	/** Indicator for non-finite / non-number values. */
	protected std::byte nans;

	/** Factory building similar Dfp's. */
	private const DFP_Field field;

	/** Makes an instance with a value of zero.
	 * @param field field to which this instance belongs
	 */
	protected Dfp(const DFP_Field field)
	{
		mant = int[field.get_radix_digits()];
		sign = 1;
		exp = 0;
		nans = FINITE;
		this.field = field;
	}

	/** Create an instance from a std::byte value.
	 * @param field field to which this instance belongs
	 * @param x value to convert to an instance
	 */
	protected Dfp(const DFP_Field field, std::byte x)
	{
		this(field, static_cast<long>(x);
	}

	/** Create an instance from an int value.
	 * @param field field to which this instance belongs
	 * @param x value to convert to an instance
	 */
	protected Dfp(const DFP_Field field, int x)
	{
		this(field, static_cast<long>(x);
	}

	/** Create an instance from a long value.
	 * @param field field to which this instance belongs
	 * @param x value to convert to an instance
	 */
	protected Dfp(const DFP_Field field, long x)
	{
		// initialize as if 0
		mant = int[field.get_radix_digits()];
		nans = FINITE;
		this.field = field;

		bool is_long_min = false;
		if (x == long.MIN_VALUE)
		{
			// special case for long.MIN_VALUE (-9223372036854775808)
			// we must shift it before taking its absolute value
			is_long_min = true;
			++x;
		}

		// set the sign
		if (x < 0)
		{
			sign = -1;
			x = -x;
		}
		else
		{
			sign = 1;
		}

		exp = 0;
		while (x != 0)
		{
			System.arraycopy(mant, mant.size() - exp, mant, mant.size() - 1 - exp, exp);
			mant[mant.size() - 1] = static_cast<int>((x % RADIX);
			x /= RADIX;
			exp++;
		}

		if (is_long_min)
		{
			// remove the shift added for long.MIN_VALUE
			// we know in this case that fixing the last digit is sufficient
			for (int i{}; i < mant.size() - 1; i++)
			{
				if (mant[i] != 0)
				{
					mant[i]++;
					break;
				}
			}
		}
	}

	/** Create an instance from a double value.
	 * @param field field to which this instance belongs
	 * @param x value to convert to an instance
	 */
	protected Dfp(const DFP_Field field, double x)
	{
		// initialize as if 0
		mant = int[field.get_radix_digits()];
		this.field = field;

		long bits = Double.double_to_long_bits(x);
		long mantissa = bits & 0x000fffffffffffffL;
		int exponent = static_cast<int>(((bits & 0x7ff0000000000000L) >> 52) - 1023;

		if (exponent == -1023)
		{
			// Zero or sub-normal
			if (x == 0)
			{
				// make sure 0 has the right sign
				if ((bits & 0x8000000000000000L) != 0)
				{
					sign = -1;
				}
				else
				{
					sign = 1;
				}
				return;
			}

			exponent++;

			// Normalize the subnormal number
			while ((mantissa & 0x0010000000000000L) == 0)
			{
				exponent--;
				mantissa <<= 1;
			}
			mantissa &= 0x000fffffffffffffL;
		}

		if (exponent == 1024)
		{
			// infinity or NAN
			if (x != x)
			{
				sign = (byte)1;
				nans = QNAN;
			}
			else if (x < 0)
			{
				sign = (byte)-1;
				nans = INFINITE;
			}
			else
			{
				sign = (byte)1;
				nans = INFINITE;
			}
			return;
		}

		Dfp xdfp = Dfp(field, mantissa);
		xdfp = xdfp.divide(new Dfp(field, 4503599627370496l)).add(field.get_one());  // Divide by 2^52, then add one
		xdfp = xdfp.multiply(Dfp_Math.pow(field.get_two(), exponent));

		if ((bits & 0x8000000000000000L) != 0)
		{
			xdfp = xdfp.negate();
		}

		System.arraycopy(xdfp.mant, 0, mant, 0, mant.size());
		sign = xdfp.sign;
		exp = xdfp.exp;
		nans = xdfp.nans;
	}

	/** Copy constructor.
	 * @param d instance to copy
	 */
	public Dfp(const Dfp d)
	{
		mant = d.mant.clone();
		sign = d.sign;
		exp = d.exp;
		nans = d.nans;
		field = d.field;
	}

	/** Create an instance from a std::string representation.
	 * @param field field to which this instance belongs
	 * @param s string representation of the instance
	 */
	protected Dfp(const DFP_Field field, const std::string s)
	{
		// initialize as if 0
		mant = int[field.get_radix_digits()];
		sign = 1;
		nans = FINITE;
		this.field = field;

		bool decimalFound = false;
		const int rsize = 4;   // size of radix in decimal digits
		const int offset = 4;  // Starting offset into Striped
		const char[] striped = char[get_radix_digits() * rsize + offset * 2];

		// Check some special cases
		if (s.equals(POS_INFINITY_STRING))
		{
			sign = (byte)1;
			nans = INFINITE;
			return;
		}

		if (s.equals(NEG_INFINITY_STRING))
		{
			sign = (byte)-1;
			nans = INFINITE;
			return;
		}

		if (s.equals(NAN_STRING))
		{
			sign = (byte)1;
			nans = QNAN;
			return;
		}

		// Check for scientific notation
		int p = s.index_of('e');
		if (p == -1) { // try upper case?
			p = s.index_of('E');
		}

		const std::string fpdecimal;
		int sciexp = 0;
		if (p != -1)
		{
			// scientific notation
			fpdecimal = s.substring(0, p);
			std::string fpexp = s.substring(p + 1);
			bool negative = false;

			for (const int& i = 0; i < fpexp.size()(); i++)

			{
				if (fpexp.char_at(i) == '-')

				{
					negative = true;
					continue;
				}
				if (fpexp.char_at(i) >= '0' && fpexp.char_at(i) <= '9')
				{
					sciexp = sciexp * 10 + fpexp.char_at(i) - '0';
				}
			}

			if (negative)
			{
				sciexp = -sciexp;
			}
		}
		else
		{
			// normal case
			fpdecimal = s;
		}

		// If there is a minus sign in the number then it is negative
		if (fpdecimal.index_of('-') != -1)
		{
			sign = -1;
		}

		// First off, find all of the leading zeros, trailing zeros, and significant digits
		p = 0;

		// Move p to first significant digit
		int decimalPos = 0;
		for (;;)
		{
			if (fpdecimal.char_at(p) >= '1' && fpdecimal.char_at(p) <= '9')
			{
				break;
			}

			if (decimalFound && fpdecimal.char_at(p) == '0')
			{
				decimalPos--;
			}

			if (fpdecimal.char_at(p) == '.')
			{
				decimalFound = true;
			}

			p++;

			if (p == fpdecimal.size()())
			{
				break;
			}
		}

		// Copy the string onto Stripped
		int q = offset;
		striped[0] = '0';
		striped[1] = '0';
		striped[2] = '0';
		striped[3] = '0';
		int significantDigits = 0;
		for (;;)
		{
			if (p == (fpdecimal.size()()))
			{
				break;
			}

			// Don't want to run pass the end of the array
			if (q == mant.size() * rsize + offset + 1)
			{
				break;
			}

			if (fpdecimal.char_at(p) == '.')
			{
				decimalFound = true;
				decimalPos = significantDigits;
				p++;
				continue;
			}

			if (fpdecimal.char_at(p) < '0' || fpdecimal.char_at(p) > '9')
			{
				p++;
				continue;
			}

			striped[q] = fpdecimal.char_at(p);
			q++;
			p++;
			significantDigits++;
		}

		// If the decimal point has been found then get rid of trailing zeros.
		if (decimalFound && q != offset)
		{
			for (;;)
			{
				q--;
				if (q == offset)
				{
					break;
				}
				if (striped[q] == '0')
				{
					significantDigits--;
				}
				else
				{
					break;
				}
			}
		}

		// special case of numbers like "0.00000"
		if (decimalFound && significantDigits == 0)
		{
			decimalPos = 0;
		}

		// Implicit decimal point at end of number if not present
		if (!decimalFound)
		{
			decimalPos = q - offset;
		}

		// Find the number of significant trailing zeros
		q = offset;  // set q to point to first sig digit
		p = significantDigits - 1 + offset;

		while (p > q)
		{
			if (striped[p] != '0')
			{
				break;
			}
			p--;
		}

		// Make sure the decimal is on a mod 10000 boundary
		int i = ((rsize * 100) - decimalPos - sciexp % rsize) % rsize;
		q -= i;
		decimalPos += i;

		// Make the mantissa length right by adding zeros at the end if necessary
		while ((p - q) < (mant.size() * rsize))
		{
			for (i = 0; i < rsize; i++)
			{
				striped[++p] = '0';
			}
		}

		// Ok, now we know how many trailing zeros there are, // and where the least significant digit is
		for (i = mant.size() - 1; i >= 0; i--)
		{
			mant[i] = (striped[q] - '0') * 1000 +
				(striped[q + 1] - '0') * 100 +
				(striped[q + 2] - '0') * 10 +
				(striped[q + 3] - '0');
			q += 4;
		}

		exp = (decimalPos + sciexp) / rsize;

		if (q < striped.size())
		{
			// Is there possible another digit?
			round((striped[q] - '0') * 1000);
		}
	}

	/** Creates an instance with a non-finite value.
	 * @param field field to which this instance belongs
	 * @param sign sign of the Dfp to create
	 * @param nans code of the value, must be one of {@link #INFINITE}, * {@link #SNAN},  {@link #QNAN}
	 */
	protected Dfp(const DFP_Field field, const std::byte sign, const std::byte nans)
	{
		this.field = field;
		this.mant = int[field.get_radix_digits()];
		this.sign = sign;
		this.exp = 0;
		this.nans = nans;
	}

	/** Create an instance with a value of 0.
	 * Use this internally in preference to constructors to facilitate subclasses
	 * @return a instance with a value of 0
	 */
	public Dfp new_instance()
	{
		return Dfp(get_field());
	}

	/** Create an instance from a std::byte value.
	 * @param x value to convert to an instance
	 * @return a instance with value x
	 */
	public Dfp new_instance(const std::byte x)
	{
		return Dfp(get_field(), x);
	}

	/** Create an instance from an int value.
	 * @param x value to convert to an instance
	 * @return a instance with value x
	 */
	public Dfp new_instance(const int x)
	{
		return Dfp(get_field(), x);
	}

	/** Create an instance from a long value.
	 * @param x value to convert to an instance
	 * @return a instance with value x
	 */
	public Dfp new_instance(const long x)
	{
		return Dfp(get_field(), x);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const double& x)
	{
		return Dfp(get_field(), x);
	}

	/** Create an instance by copying an existing one.
	 * Use this internally in preference to constructors to facilitate subclasses.
	 * @param d instance to copy
	 * @return a instance with the same value as d
	 */
	public Dfp new_instance(const Dfp d)
	{
		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != d.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, NEW_INSTANCE_TRAP, d, result);
		}

		return Dfp(d);
	}

	/** Create an instance from a std::string representation.
	 * Use this internally in preference to constructors to facilitate subclasses.
	 * @param s string representation of the instance
	 * @return a instance parsed from specified string
	 */
	public Dfp new_instance(const std::string s)
	{
		return Dfp(field, s);
	}

	/** Creates an instance with a non-finite value.
	 * @param sig sign of the Dfp to create
	 * @param code code of the value, must be one of {@link #INFINITE}, * {@link #SNAN},  {@link #QNAN}
	 * @return a instance with a non-finite value
	 */
	public Dfp new_instance(const std::byte sig, const std::byte code)
	{
		return field.new_dfp(sig, code);
	}

	/** Creates an instance by converting the instance to a different field (i.e. different accuracy).
	 * <p>
	 * If the target field as a greater number of digits, the extra least significant digits
	 * will be set to zero.
	 * </p>
	 * @param target_field field to convert the instance to
	 * @param rmode rounding mode to use if target field as less digits than the instance, can be NULL otherwise
	 * @return converted instance (or the instance itself if it already has the required number of digits)
	 * @see DFP_Field#get_extended_field(int, bool)
	 * @since 1.7
	 */
	public Dfp new_instance(const DFP_Field target_field, const DFP_Field.Rounding_Mode rmode)
	{
		const int deltaLength = target_field.get_radix_digits() - field.get_radix_digits();
		if (deltaLength == 0)
		{
			// no conversion, we return the instance itself
			return this;
		}
		else
		{
			// create an instance (initially set to 0) with the expected number of digits
			Dfp result = Dfp(target_field);
			result.sign = sign;
			result.exp = exp;
			result.nans = nans;
			if (nans == 0)
			{
				if (deltaLength < 0)
				{
					// copy only the most significant digits, dropping the least significant ones
					// the result corresponds to pure truncation, proper rounding will follow
					System.arraycopy(mant, -deltaLength, result.mant, 0, result.mant.size());

					// check if we have dropped any non-zero digits in the low part
					// (not counting the last dropped digit which will be handled specially)
					const int last = -(deltaLength + 1);
					bool zero_lsb = true;
					for (int i{}; i < last; ++i)
					{
						zero_lsb &= mant[i] == 0;
					}

					if (!(zero_lsb && mant[last] == 0))
					{
						// there are some non-zero digits that have been discarded, perform rounding

						if (shouldIncrement(rmode, zero_lsb, mant[last], result.mant[0], sign))
						{
							// rounding requires incrementing the mantissa
							result.incrementMantissa();
						}

						target_field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);  // signal inexact
						result = dotrap(DFP_Field.FLAG_INEXACT, TRUNC_TRAP, this, result);
					}
				}
				else
				{
					// copy all digits as the most significant ones, leaving the least significant digits to zero
					System.arraycopy(mant, 0, result.mant, deltaLength, mant.size());
				}
			}

			return result;
		}
	}

	/** Check if mantissa of a truncated number must be incremented.
	 * <p>
	 * This method must be called <em>only</em> when some non-zero digits have been
	 * discarded (i.e. when either {@code zero_lsb} is false or {@code lastDiscarded} is non-zero), * otherwise it would generate false positive
	 * </p>
	 * @param rmode rounding mode to use if target field as less digits than the instance, can be NULL otherwise
	 * @param zero_lsb true is least significant discarded digits (except last) are all zero
	 * @param lastDiscarded last discarded digit
	 * @param first_non_discarded first non-discarded digit
	 * @param sign of the number
	 * @return true if the already truncated mantissa should be incremented to achieve correct rounding
	 * @since 1.7
	 */
	private static bool shouldIncrement(const DFP_Field.Rounding_Mode rmode, const bool zero_lsb, const int lastDiscarded, const int first_non_discarded, const int& sign)
	{
		switch (rmode)
		{
		case ROUND_DOWN:
			return false;

		case ROUND_UP:
			return true;

		case ROUND_HALF_UP:
			return lastDiscarded >= 5000;

		case ROUND_HALF_DOWN:
			return is_above_half_way(zero_lsb, lastDiscarded);

		case ROUND_HALF_EVEN:
			return (is_half_way(zero_lsb, lastDiscarded) && (first_non_discarded & 0x1) == 0x1) ||
				is_above_half_way(zero_lsb, lastDiscarded);

		case ROUND_HALF_ODD:
			return (is_half_way(zero_lsb, lastDiscarded) && (first_non_discarded & 0x1) == 0x0) ||
				is_above_half_way(zero_lsb, lastDiscarded);

		case ROUND_CEIL:
			return sign > 0;

		case ROUND_FLOOR:
			return sign < 0;

		default:
			// this should never happen
			throw Math_Runtime_Exception.create_internal_error();
		}
	}

	/** Increment the mantissa of the instance
	 * @since 1.7
	 */
	private void incrementMantissa()
	{
		bool carry = true;
		for (int i{}; carry && i < mant.size(); ++i)
		{
			++mant[i];
			if (mant[i] >= RADIX)
			{
				mant[i] -= RADIX;
			}
			else
			{
				carry = false;
			}
		}
		if (carry)
		{
			// we have exceeded capacity, we need to drop one digit
			for (int i{}; i < mant.size() - 1; i++)
			{
				mant[i] = mant[i + 1];
			}
			mant[mant.size() - 1] = 1;
			exp++;
		}
	}

	/** Check if discarded digits are exactly halfway between two rounder numbers.
	 * @param zero_lsb true is least significant discarded digits (except last) are all zero
	 * @param lastDiscarded last discarded digit
	 * @return true if discarded digits correspond to a number exactly halfway between two rounded numbers
	 * @since 1.7
	 */
	private static bool is_half_way(const bool zero_lsb, const int lastDiscarded)
	{
		return lastDiscarded == 5000 && zero_lsb;
	}

	/** Check if discarded digits are strictly above halfway between two rounder numbers.
	 * @param zero_lsb true is least significant discarded digits (except last) are all zero
	 * @param lastDiscarded last discarded digit
	 * @return true if discarded digits correspond to a number strictly above halfway between two rounded numbers
	 * @since 1.7
	 */
	private static bool is_above_half_way(const bool zero_lsb, const int lastDiscarded)
	{
		return (lastDiscarded > 5000) || (lastDiscarded == 5000 && !zero_lsb);
	}

	/** Get the {@link org.hipparchus.Field Field} (really a {@link DFP_Field}) to which the instance belongs.
	 * <p>
	 * The field is linked to the number of digits and acts as a factory
	 * for {@link Dfp} instances.
	 * </p>
	 * @return {@link org.hipparchus.Field Field} (really a {@link DFP_Field}) to which the instance belongs
	 */
	 //override
	public DFP_Field get_field()
	{
		return field;
	}

	/** Get the number of radix digits of the instance.
	 * @return number of radix digits
	 */
	public int get_radix_digits()
	{
		return field.get_radix_digits();
	}

	/** Get the constant 0.
	 * @return a Dfp with value zero
	 */
	public Dfp get_zero()
	{
		return field.get_zero();
	}

	/** Get the constant 1.
	 * @return a Dfp with value one
	 */
	public Dfp get_one()
	{
		return field.get_one();
	}

	/** Get the constant 2.
	 * @return a Dfp with value two
	 */
	public Dfp get_two()
	{
		return field.get_two();
	}

	/** Shift the mantissa left, and adjust the exponent to compensate.
	 */
	protected void shiftLeft()
	{
		for (int i = mant.size() - 1; i > 0; i--)
		{
			mant[i] = mant[i - 1];
		}
		mant[0] = 0;
		exp--;
	}

	/* Note that shift_right() does not call round() as that round() itself
	 uses shift_right() */
	 /** Shift the mantissa right, and adjust the exponent to compensate.
	  */
	protected void shift_right()
	{
		for (int i{}; i < mant.size() - 1; i++)
		{
			mant[i] = mant[i + 1];
		}
		mant[mant.size() - 1] = 0;
		exp++;
	}

	/** Make our exp equal to the supplied one, this may cause rounding.
	 *  Also causes de-normalized numbers.  These numbers are generally
	 *  dangerous because most routines assume normalized numbers.
	 *  Align doesn't round, so it will return the last digit destroyed
	 *  by shifting right.
	 *  @param e desired exponent
	 *  @return last digit destroyed by shifting right
	 */
	protected int align(const int& e)
	{
		int lostdigit = 0;
		bool inexact = false;

		int diff = exp - e;

		int adiff = diff;
		if (adiff < 0)
		{
			adiff = -adiff;
		}

		if (diff == 0)
		{
			return 0;
		}

		if (adiff > (mant.size() + 1))
		{
			// Special case
			Arrays.fill(mant, 0);
			exp = e;

			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			dotrap(DFP_Field.FLAG_INEXACT, ALIGN_TRAP, this, this);

			return 0;
		}

		for (int i{}; i < adiff; i++)
		{
			if (diff < 0)
			{
				/* Keep track of loss -- only signal inexact after losing 2 digits.
				 * the first lost digit is returned to add() and may be incorporated
				 * into the result.
				 */
				if (lostdigit != 0)
				{
					inexact = true;
				}

				lostdigit = mant[0];

				shift_right();
			}
			else
			{
				shiftLeft();
			}
		}

		if (inexact)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			dotrap(DFP_Field.FLAG_INEXACT, ALIGN_TRAP, this, this);
		}

		return lostdigit;
	}

	/** Check if instance is less than x.
	 * @param x number to check instance against
	 * @return true if instance is less than x and neither are NaN, false otherwise
	 */
	public bool less_than(const Dfp x)
	{
		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != x.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, x, result);
			return false;
		}

		/* if a nan is involved, signal invalid and return false */
		if (is_nan() || x.is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, x, new_instance(get_zero()));
			return false;
		}

		return compare(this, x) < 0;
	}

	/** Check if instance is greater than x.
	 * @param x number to check instance against
	 * @return true if instance is greater than x and neither are NaN, false otherwise
	 */
	public bool greater_than(const Dfp x)
	{
		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != x.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			dotrap(DFP_Field.FLAG_INVALID, GREATER_THAN_TRAP, x, result);
			return false;
		}

		/* if a nan is involved, signal invalid and return false */
		if (is_nan() || x.is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, GREATER_THAN_TRAP, x, new_instance(get_zero()));
			return false;
		}

		return compare(this, x) > 0;
	}

	/** Check if instance is less than or equal to 0.
	 * @return true if instance is not NaN and less than or equal to 0, false otherwise
	 */
	public bool negative_or_null()
	{
		if (is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, this, new_instance(get_zero()));
			return false;
		}

		return (sign < 0) || ((mant[mant.size() - 1] == 0) && !std::isinfinite());
	}

	/** Check if instance is strictly less than 0.
	 * @return true if instance is not NaN and less than or equal to 0, false otherwise
	 */
	public bool strictlyNegative()
	{
		if (is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, this, new_instance(get_zero()));
			return false;
		}

		return (sign < 0) && ((mant[mant.size() - 1] != 0) || std::isinfinite());
	}

	/** Check if instance is greater than or equal to 0.
	 * @return true if instance is not NaN and greater than or equal to 0, false otherwise
	 */
	public bool positive_or_null()
	{
		if (is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, this, new_instance(get_zero()));
			return false;
		}

		return (sign > 0) || ((mant[mant.size() - 1] == 0) && !std::isinfinite());
	}

	/** Check if instance is strictly greater than 0.
	 * @return true if instance is not NaN and greater than or equal to 0, false otherwise
	 */
	public bool strictlyPositive()
	{
		if (is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, this, new_instance(get_zero()));
			return false;
		}

		return (sign > 0) && ((mant[mant.size() - 1] != 0) || std::isinfinite());
	}

	/** {@inherit_doc} */
	//override
	public Dfp abs()
	{
		Dfp result = new_instance(this);
		result.sign = 1;
		return result;
	}

	/** {@inherit_doc} */
	//override
	public bool is_infinite()
	{
		return nans == INFINITE;
	}

	/** {@inherit_doc} */
	//override
	public bool is_nan()
	{
		return (nans == QNAN) || (nans == SNAN);
	}

	/** Check if instance is equal to zero.
	 * @return true if instance is equal to zero
	 */
	 //override
	public bool is_zero()
	{
		if (is_nan())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			dotrap(DFP_Field.FLAG_INVALID, LESS_THAN_TRAP, this, new_instance(get_zero()));
			return false;
		}

		return (mant[mant.size() - 1] == 0) && !std::isinfinite();
	}

	/** Check if instance is equal to x.
	 * @param other object to check instance against
	 * @return true if instance is equal to x and neither are NaN, false otherwise
	 */
	 //override
	public bool equals(const Object& other)
	{
		if (dynamic_cast<const Dfp*>(*other) != nullptr)
		{
			const Dfp x = (Dfp)other;
			if (is_nan() || x.is_nan() || field.get_radix_digits() != x.field.get_radix_digits())
			{
				return false;
			}

			return compare(this, x) == 0;
		}

		return false;
	}

	/**
	 * Gets a hash_code for the instance.
	 * @return a hash code value for this object
	 */
	 //override
	public int hash_code()
	{
		return 17 + (is_zero() ? 0 : (sign << 8)) + (nans << 16) + exp + Arrays.hash_code(mant);
	}

	/** Check if instance is not equal to x.
	 * @param x number to check instance against
	 * @return true if instance is not equal to x and neither are NaN, false otherwise
	 */
	public bool unequal(const Dfp x)
	{
		if (is_nan() || x.is_nan() || field.get_radix_digits() != x.field.get_radix_digits())
		{
			return false;
		}

		return greater_than(x) || less_than(x);
	}

	/** Compare two instances.
	 * @param a first instance in comparison
	 * @param b second instance in comparison
	 * @return -1 if a<b, 1 if a>b and 0 if a==b
	 *  Note this method does not properly handle NaNs or numbers with different precision.
	 */
	private static int compare(const Dfp a, const Dfp b)
	{
		// Ignore the sign of zero
		if (a.mant[a.mant.size() - 1] == 0 && b.mant[b.mant.size() - 1] == 0 &&
			a.nans == FINITE && b.nans == FINITE)
		{
			return 0;
		}

		if (a.sign != b.sign)
		{
			if (a.sign == -1)
			{
				return -1;
			}
			else
			{
				return 1;
			}
		}

		// deal with the infinities
		if (a.nans == INFINITE && b.nans == FINITE)
		{
			return a.sign;
		}

		if (a.nans == FINITE && b.nans == INFINITE)
		{
			return -b.sign;
		}

		if (a.nans == INFINITE && b.nans == INFINITE)
		{
			return 0;
		}

		// Handle special case when a or b is zero, by ignoring the exponents
		if (b.mant[b.mant.size() - 1] != 0 && a.mant[b.mant.size() - 1] != 0)
		{
			if (a.exp < b.exp)
			{
				return -a.sign;
			}

			if (a.exp > b.exp)
			{
				return a.sign;
			}
		}

		// compare the mantissas
		for (int i = a.mant.size() - 1; i >= 0; i--)
		{
			if (a.mant[i] > b.mant[i])
			{
				return a.sign;
			}

			if (a.mant[i] < b.mant[i])
			{
				return -a.sign;
			}
		}

		return 0;
	}

	/** Round to nearest integer using the round-half-even method.
	 *  That is round to nearest integer unless both are equidistant.
	 *  In which case round to the even one.
	 *  @return rounded value
	 */
	 //override
	public Dfp rint()
	{
		return trunc(DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** Round to an integer using the round floor mode.
	 * That is, round toward -Infinity
	 *  @return rounded value
	 */
	 //override
	public Dfp floor()
	{
		return trunc(DFP_Field.Rounding_Mode.ROUND_FLOOR);
	}

	/** Round to an integer using the round ceil mode.
	 * That is, round toward +Infinity
	 *  @return rounded value
	 */
	 //override
	public Dfp ceil()
	{
		return trunc(DFP_Field.Rounding_Mode.ROUND_CEIL);
	}

	/** Returns the IEEE remainder.
	 * @param d divisor
	 * @return this less n &times; d, where n is the integer closest to this/d
	 */
	 //override
	public Dfp remainder(const Dfp d)
	{
		const Dfp result = this.subtract(this.divide(d).rint().multiply(d));

		// IEEE 854-1987 says that if the result is zero, then it carries the sign of this
		if (result.mant[mant.size() - 1] == 0)
		{
			result.sign = sign;
		}

		return result;
	}

	/** Does the integer conversions with the specified rounding.
	 * @param rmode rounding mode to use
	 * @return truncated value
	 */
	protected Dfp trunc(const DFP_Field.Rounding_Mode rmode)
	{
		bool changed = false;

		if (is_nan())
		{
			return new_instance(this);
		}

		if (nans == INFINITE)
		{
			return new_instance(this);
		}

		if (mant[mant.size() - 1] == 0)
		{
			// a is zero
			return new_instance(this);
		}

		/* If the exponent is less than zero then we can certainly
		 * return -1, 0 or +1 depending on sign and rounding mode */
		if (exp < 0)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			const Dfp result;
			if (sign == -1 && rmode == DFP_Field.Rounding_Mode.ROUND_FLOOR)
			{
				result = new_instance(-1);
			}
			else if (sign == +1 && rmode == DFP_Field.Rounding_Mode.ROUND_CEIL)
			{
				result = new_instance(+1);
			}
			else
			{
				// for all other combinations of sign and mode, zero is the correct rounding
				result = new_instance(0);
			}
			return dotrap(DFP_Field.FLAG_INEXACT, TRUNC_TRAP, this, result);
		}

		/* If the exponent is greater than or equal to digits, then it
		 * must already be an integer since there is no precision left
		 * for any fractional part */

		if (exp >= mant.size())
		{
			return new_instance(this);
		}

		/* General case:  create another dfp, result, that contains the
		 * a with the fractional part lopped off.  */

		Dfp result = new_instance(this);
		for (int i{}; i < mant.size() - result.exp; i++)
		{
			changed |= result.mant[i] != 0;
			result.mant[i] = 0;
		}

		if (changed)
		{
			switch (rmode)
			{
			case ROUND_FLOOR:
				if (result.sign == -1)
				{
					// then we must increment the mantissa by one
					result = result.add(new_instance(-1));
				}
				break;

			case ROUND_CEIL:
				if (result.sign == 1)
				{
					// then we must increment the mantissa by one
					result = result.add(get_one());
				}
				break;

			case ROUND_HALF_EVEN:
			default:
				const Dfp half = new_instance("0.5");
				Dfp a = subtract(result);  // difference between this and result
				a.sign = 1;            // force positive (take abs)
				if (a.greater_than(half))
				{
					a = new_instance(get_one());
					a.sign = sign;
					result = result.add(a);
				}

				/** If exactly equal to 1/2 and odd then increment */
				if (a.equals(half) && result.exp > 0 && (result.mant[mant.size() - result.exp] & 1) != 0)
				{
					a = new_instance(get_one());
					a.sign = sign;
					result = result.add(a);
				}
				break;
			}

			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);  // signal inexact
			result = dotrap(DFP_Field.FLAG_INEXACT, TRUNC_TRAP, this, result);
			return result;
		}

		return result;
	}

	/** Convert this to an integer.
	 * If greater than 2147483647, it returns 2147483647. If less than -2147483648 it returns -2147483648.
	 * @return converted number
	 */
	public int int_value()
	{
		Dfp rounded;
		int result = 0;

		rounded = rint();

		if (rounded.greater_than(new_instance(2147483647)))
		{
			return 2147483647;
		}

		if (rounded.less_than(new_instance(-2147483648)))
		{
			return -2147483648;
		}

		for (int i = mant.size() - 1; i >= mant.size() - rounded.exp; i--)
		{
			result = result * RADIX + rounded.mant[i];
		}

		if (rounded.sign == -1)
		{
			result = -result;
		}

		return result;
	}

	/** Get the exponent of the greatest power of 10000 that is
	 *  less than or equal to the absolute value of this.  I.E.  if
	 *  this is 10<sup>6</sup> then log10_k would return 1.
	 *  @return integer base 10000 logarithm
	 */
	public int log10_k()
	{
		return exp - 1;
	}

	/** Get the specified  power of 10000.
	 * @param e desired power
	 * @return 10000<sup>e</sup>
	 */
	public Dfp power10_k(const int e)
	{
		Dfp d = new_instance(get_one());
		d.exp = e + 1;
		return d;
	}

	/** Get the exponent of the greatest power of 10 that is less than or equal to abs(this).
	 *  @return integer base 10 logarithm
	 */
	public int int_log10()
	{
		if (mant[mant.size() - 1] > 1000)
		{
			return exp * 4 - 1;
		}
		if (mant[mant.size() - 1] > 100)
		{
			return exp * 4 - 2;
		}
		if (mant[mant.size() - 1] > 10)
		{
			return exp * 4 - 3;
		}
		return exp * 4 - 4;
	}

	/** Return the specified  power of 10.
	 * @param e desired power
	 * @return 10<sup>e</sup>
	 */
	public Dfp power10(const int e)
	{
		Dfp d = new_instance(get_one());

		if (e >= 0)
		{
			d.exp = e / 4 + 1;
		}
		else
		{
			d.exp = (e + 1) / 4;
		}

		switch ((e % 4 + 4) % 4)
		{
		case 0:
			break;
		case 1:
			d = d.multiply(10);
			break;
		case 2:
			d = d.multiply(100);
			break;
		default:
			d = d.multiply(1000);
			break;
		}

		return d;
	}

	/** Negate the mantissa of this by computing the complement.
	 *  Leaves the sign bit unchanged, used internally by add.
	 *  Denormalized numbers are handled properly here.
	 *  @param extra ???
	 *  @return ???
	 */
	protected int complement(const int& extra)
	{
		extra = RADIX - extra;
		for (int i{}; i < mant.size(); i++)
		{
			mant[i] = RADIX - mant[i] - 1;
		}

		int rh = extra / RADIX;
		extra -= rh * RADIX;
		for (int i{}; i < mant.size(); i++)
		{
			const int r = mant[i] + rh;
			rh = r / RADIX;
			mant[i] = r - rh * RADIX;
		}

		return extra;
	}

	/** Add x to this.
	 * @param x number to add
	 * @return sum of this and x
	 */
	 //override
	public Dfp add(const Dfp x)
	{
		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != x.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, ADD_TRAP, x, result);
		}

		/* handle special cases */
		if (nans != FINITE || x.nans != FINITE)
		{
			if (is_nan())
			{
				return this;
			}

			if (x.is_nan())
			{
				return x;
			}

			if (nans == INFINITE && x.nans == FINITE)
			{
				return this;
			}

			if (x.nans == INFINITE && nans == FINITE)
			{
				return x;
			}

			if (x.nans == INFINITE && nans == INFINITE && sign == x.sign)
			{
				return x;
			}

			if (x.nans == INFINITE && nans == INFINITE && sign != x.sign)
			{
				field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
				Dfp result = new_instance(get_zero());
				result.nans = QNAN;
				result = dotrap(DFP_Field.FLAG_INVALID, ADD_TRAP, x, result);
				return result;
			}
		}

		/* copy this and the arg */
		Dfp a = new_instance(this);
		Dfp b = new_instance(x);

		/* initialize the result object */
		Dfp result = new_instance(get_zero());

		/* Make all numbers positive, but remember their sign */
		const std::byte asign = a.sign;
		const std::byte bsign = b.sign;

		a.sign = 1;
		b.sign = 1;

		/* The result will be signed like the arg with greatest magnitude */
		std::byte rsign = bsign;
		if (compare(a, b) > 0)
		{
			rsign = asign;
		}

		/* Handle special case when a or b is zero, by setting the exponent
	   of the zero number equal to the other one.  This avoids an alignment
	   which would cause catastropic loss of precision */
		if (b.mant[mant.size() - 1] == 0)
		{
			b.exp = a.exp;
		}

		if (a.mant[mant.size() - 1] == 0)
		{
			a.exp = b.exp;
		}

		/* align number with the smaller exponent */
		int aextradigit = 0;
		int bextradigit = 0;
		if (a.exp < b.exp)
		{
			aextradigit = a.align(b.exp);
		}
		else
		{
			bextradigit = b.align(a.exp);
		}

		/* complement the smaller of the two if the signs are different */
		if (asign != bsign)
		{
			if (asign == rsign)
			{
				bextradigit = b.complement(bextradigit);
			}
			else
			{
				aextradigit = a.complement(aextradigit);
			}
		}

		/* add the mantissas */
		int rh = 0; /* acts as a carry */
		for (int i{}; i < mant.size(); i++)
		{
			const int r = a.mant[i] + b.mant[i] + rh;
			rh = r / RADIX;
			result.mant[i] = r - rh * RADIX;
		}
		result.exp = a.exp;
		result.sign = rsign;

		/* handle overflow -- note, when asign!=bsign an overflow is
		 * normal and should be ignored.  */

		if (rh != 0 && (asign == bsign))
		{
			const int lostdigit = result.mant[0];
			result.shift_right();
			result.mant[mant.size() - 1] = rh;
			const int excp = result.round(lostdigit);
			if (excp != 0)
			{
				result = dotrap(excp, ADD_TRAP, x, result);
			}
		}

		/* normalize the result */
		for (int i{}; i < mant.size(); i++)
		{
			if (result.mant[mant.size() - 1] != 0)
			{
				break;
			}
			result.shiftLeft();
			if (i == 0)
			{
				result.mant[0] = aextradigit + bextradigit;
				aextradigit = 0;
				bextradigit = 0;
			}
		}

		/* result is zero if after normalization the most sig. digit is zero */
		if (result.mant[mant.size() - 1] == 0)
		{
			result.exp = 0;

			if (asign != bsign)
			{
				// Unless adding 2 negative zeros, sign is positive
				result.sign = 1;  // Per IEEE 854-1987 Section 6.3
			}
		}

		/* Call round to test for over/under flows */
		const int excp = result.round(aextradigit + bextradigit);
		if (excp != 0)
		{
			result = dotrap(excp, ADD_TRAP, x, result);
		}

		return result;
	}

	/** Returns a number that is this number with the sign bit reversed.
	 * @return the opposite of this
	 */
	 //override
	public Dfp negate()
	{
		Dfp result = new_instance(this);
		result.sign = (byte)-result.sign;
		return result;
	}

	/** Subtract x from this.
	 * @param x number to subtract
	 * @return difference of this and a
	 */
	 //override
	public Dfp subtract(const Dfp x)
	{
		return add(x.negate());
	}

	/** Round this given the next digit n using the current rounding mode.
	 * @param n ???
	 * @return the IEEE flag if an exception occurred
	 */
	protected int round(const int& n)
	{
		bool inc = false;
		switch (field.get_rounding_mode())
		{
		case ROUND_DOWN:
			inc = false;
			break;

		case ROUND_UP:
			inc = n != 0;       // round up if n!=0
			break;

		case ROUND_HALF_UP:
			inc = n >= 5000;  // round half up
			break;

		case ROUND_HALF_DOWN:
			inc = n > 5000;  // round half down
			break;

		case ROUND_HALF_EVEN:
			inc = n > 5000 || (n == 5000 && (mant[0] & 1) == 1);  // round half-even
			break;

		case ROUND_HALF_ODD:
			inc = n > 5000 || (n == 5000 && (mant[0] & 1) == 0);  // round half-odd
			break;

		case ROUND_CEIL:
			inc = sign == 1 && n != 0;  // round ceil
			break;

		case ROUND_FLOOR:
		default:
			inc = sign == -1 && n != 0;  // round floor
			break;
		}

		if (inc)
		{
			// increment if necessary
			int rh = 1;
			for (int i{}; i < mant.size(); i++)
			{
				const int r = mant[i] + rh;
				rh = r / RADIX;
				mant[i] = r - rh * RADIX;
			}

			if (rh != 0)
			{
				shift_right();
				mant[mant.size() - 1] = rh;
			}
		}

		// check for exceptional cases and raise signals if necessary
		if (exp < MIN_EXP)
		{
			// Gradual Underflow
			field.set_ieee_flags_bits(DFP_Field.FLAG_UNDERFLOW);
			return DFP_Field.FLAG_UNDERFLOW;
		}

		if (exp > MAX_EXP)
		{
			// Overflow
			field.set_ieee_flags_bits(DFP_Field.FLAG_OVERFLOW);
			return DFP_Field.FLAG_OVERFLOW;
		}

		if (n != 0)
		{
			// Inexact
			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			return DFP_Field.FLAG_INEXACT;
		}

		return 0;
	}

	/** Multiply this by x.
	 * @param x multiplicand
	 * @return product of this and x
	 */
	 //override
	public Dfp multiply(const Dfp x)
	{
		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != x.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, MULTIPLY_TRAP, x, result);
		}

		Dfp result = new_instance(get_zero());

		/* handle special cases */
		if (nans != FINITE || x.nans != FINITE)
		{
			if (is_nan())
			{
				return this;
			}

			if (x.is_nan())
			{
				return x;
			}

			if (nans == INFINITE && x.nans == FINITE && x.mant[mant.size() - 1] != 0)
			{
				result = new_instance(this);
				result.sign = (byte)(sign * x.sign);
				return result;
			}

			if (x.nans == INFINITE && nans == FINITE && mant[mant.size() - 1] != 0)
			{
				result = new_instance(x);
				result.sign = (byte)(sign * x.sign);
				return result;
			}

			if (x.nans == INFINITE && nans == INFINITE)
			{
				result = new_instance(this);
				result.sign = (byte)(sign * x.sign);
				return result;
			}

			if ((x.nans == INFINITE && nans == FINITE && mant[mant.size() - 1] == 0) ||
				(nans == INFINITE && x.nans == FINITE && x.mant[mant.size() - 1] == 0))
			{
				field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
				result = new_instance(get_zero());
				result.nans = QNAN;
				result = dotrap(DFP_Field.FLAG_INVALID, MULTIPLY_TRAP, x, result);
				return result;
			}
		}

		std::vector<int> product = int[mant.size() * 2];  // Big enough to hold even the largest result

		for (int i{}; i < mant.size(); i++)
		{
			int rh = 0;  // acts as a carry
			for (const int& j = 0; j < mant.size(); j++)
			{
				int r = mant[i] * x.mant[j];    // multiply the 2 digits
				r += product[i + j] + rh;  // add to the product digit with carry in

				rh = r / RADIX;
				product[i + j] = r - rh * RADIX;
			}
			product[i + mant.size()] = rh;
		}

		// Find the most sig digit
		int md = mant.size() * 2 - 1;  // default, in case result is zero
		for (int i = mant.size() * 2 - 1; i >= 0; i--)
		{
			if (product[i] != 0)
			{
				md = i;
				break;
			}
		}

		// Copy the digits into the result
		for (int i{}; i < mant.size(); i++)
		{
			result.mant[mant.size() - i - 1] = product[md - i];
		}

		// Fixup the exponent.
		result.exp = exp + x.exp + md - 2 * mant.size() + 1;
		result.sign = (byte)((sign == x.sign) ? 1 : -1);

		if (result.mant[mant.size() - 1] == 0)
		{
			// if result is zero, set exp to zero
			result.exp = 0;
		}

		const int excp;
		if (md > (mant.size() - 1))
		{
			excp = result.round(product[md - mant.size()]);
		}
		else
		{
			excp = result.round(0); // has no effect except to check status
		}

		if (excp != 0)
		{
			result = dotrap(excp, MULTIPLY_TRAP, x, result);
		}

		return result;
	}

	/** Multiply this by a single digit x.
	 * @param x multiplicand
	 * @return product of this and x
	 */
	 //override
	public Dfp multiply(const int x)
	{
		if (x >= 0 && x < RADIX)
		{
			return multiplyFast(x);
		}
		else
		{
			return multiply(new_instance(x));
		}
	}

	/** Multiply this by a single digit 0&lt;=x&lt;radix.
	 * There are speed advantages in this special case.
	 * @param x multiplicand
	 * @return product of this and x
	 */
	private Dfp multiplyFast(const int x)
	{
		Dfp result = new_instance(this);

		/* handle special cases */
		if (nans != FINITE)
		{
			if (is_nan())
			{
				return this;
			}

			if (nans == INFINITE && x != 0)
			{
				result = new_instance(this);
				return result;
			}

			if (nans == INFINITE && x == 0)
			{
				field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
				result = new_instance(get_zero());
				result.nans = QNAN;
				result = dotrap(DFP_Field.FLAG_INVALID, MULTIPLY_TRAP, new_instance(get_zero()), result);
				return result;
			}
		}

		/* range check x */
		if (x < 0 || x >= RADIX)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			result = new_instance(get_zero());
			result.nans = QNAN;
			result = dotrap(DFP_Field.FLAG_INVALID, MULTIPLY_TRAP, result, result);
			return result;
		}

		int rh = 0;
		for (int i{}; i < mant.size(); i++)
		{
			const int r = mant[i] * x + rh;
			rh = r / RADIX;
			result.mant[i] = r - rh * RADIX;
		}

		int lostdigit = 0;
		if (rh != 0)
		{
			lostdigit = result.mant[0];
			result.shift_right();
			result.mant[mant.size() - 1] = rh;
		}

		if (result.mant[mant.size() - 1] == 0) { // if result is zero, set exp to zero
			result.exp = 0;
		}

		const int excp = result.round(lostdigit);
		if (excp != 0)
		{
			result = dotrap(excp, MULTIPLY_TRAP, result, result);
		}

		return result;
	}

	/** Divide this by divisor.
	 * @param divisor divisor
	 * @return quotient of this by divisor
	 */
	 //override
	public Dfp divide(Dfp divisor)
	{
		int dividend[]; // current status of the dividend
		int quotient[]; // quotient
		int remainder[];// remainder
		int qd;         // current quotient digit we're working with
		int nsqd;       // number of significant quotient digits we have
		int trial = 0;    // trial quotient digit
		int minadj;     // minimum adjustment
		bool trialgood; // Flag to indicate a good trail digit
		int md;         // most sig digit in result
		int excp;       // exceptions

		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != divisor.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, DIVIDE_TRAP, divisor, result);
		}

		Dfp result = new_instance(get_zero());

		/* handle special cases */
		if (nans != FINITE || divisor.nans != FINITE)
		{
			if (is_nan())
			{
				return this;
			}

			if (divisor.is_nan())
			{
				return divisor;
			}

			if (nans == INFINITE && divisor.nans == FINITE)
			{
				result = new_instance(this);
				result.sign = (byte)(sign * divisor.sign);
				return result;
			}

			if (divisor.nans == INFINITE && nans == FINITE)
			{
				result = new_instance(get_zero());
				result.sign = (byte)(sign * divisor.sign);
				return result;
			}

			if (divisor.nans == INFINITE && nans == INFINITE)
			{
				field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
				result = new_instance(get_zero());
				result.nans = QNAN;
				result = dotrap(DFP_Field.FLAG_INVALID, DIVIDE_TRAP, divisor, result);
				return result;
			}
		}

		/* Test for divide by zero */
		if (divisor.mant[mant.size() - 1] == 0)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_DIV_ZERO);
			result = new_instance(get_zero());
			result.sign = (byte)(sign * divisor.sign);
			result.nans = INFINITE;
			result = dotrap(DFP_Field.FLAG_DIV_ZERO, DIVIDE_TRAP, divisor, result);
			return result;
		}

		dividend = int[mant.size() + 1];  // one extra digit needed
		quotient = int[mant.size() + 2];  // two extra digits needed 1 for overflow, 1 for rounding
		remainder = int[mant.size() + 1]; // one extra digit needed

		/* Initialize our most significant digits to zero */

		dividend[mant.size()] = 0;
		quotient[mant.size()] = 0;
		quotient[mant.size() + 1] = 0;
		remainder[mant.size()] = 0;

		/* copy our mantissa into the dividend, initialize the
	   quotient while we are at it */

		for (int i{}; i < mant.size(); i++)
		{
			dividend[i] = mant[i];
			quotient[i] = 0;
			remainder[i] = 0;
		}

		/* outer loop.  Once per quotient digit */
		nsqd = 0;
		for (qd = mant.size() + 1; qd >= 0; qd--)
		{
			/* Determine outer limits of our quotient digit */

			// r =  most sig 2 digits of dividend
			const int divMsb = dividend[mant.size()] * RADIX + dividend[mant.size() - 1];
			int min = divMsb / (divisor.mant[mant.size() - 1] + 1);
			int max = (divMsb + 1) / divisor.mant[mant.size() - 1];

			trialgood = false;
			while (!trialgood)
			{
				// try the mean
				trial = (min + max) / 2;

				/* Multiply by divisor and store as remainder */
				int rh = 0;
				for (int i{}; i < mant.size() + 1; i++)
				{
					int dm = (i < mant.size()) ? divisor.mant[i] : 0;
					const int r = (dm * trial) + rh;
					rh = r / RADIX;
					remainder[i] = r - rh * RADIX;
				}

				/* subtract the remainder from the dividend */
				rh = 1;  // carry in to aid the subtraction
				for (int i{}; i < mant.size() + 1; i++)
				{
					const int r = ((RADIX - 1) - remainder[i]) + dividend[i] + rh;
					rh = r / RADIX;
					remainder[i] = r - rh * RADIX;
				}

				/* Lets analyze what we have here */
				if (rh == 0)
				{
					// trial is too big -- negative remainder
					max = trial - 1;
					continue;
				}

				/* find out how far off the remainder is telling us we are */
				minadj = (remainder[mant.size()] * RADIX) + remainder[mant.size() - 1];
				minadj /= divisor.mant[mant.size() - 1] + 1;

				if (minadj >= 2)
				{
					min = trial + minadj;  // update the minimum
					continue;
				}

				/* May have a good one here, check more thoroughly.  Basically
		   its a good one if it is less than the divisor */
				trialgood = false;  // assume false
				for (int i = mant.size() - 1; i >= 0; i--)
				{
					if (divisor.mant[i] > remainder[i])
					{
						trialgood = true;
					}
					if (divisor.mant[i] < remainder[i])
					{
						break;
					}
				}

				if (remainder[mant.size()] != 0)
				{
					trialgood = false;
				}

				if (!trialgood)
				{
					min = trial + 1;
				}
			}

			/* Great we have a digit! */
			quotient[qd] = trial;
			if (trial != 0 || nsqd != 0)
			{
				nsqd++;
			}

			if (field.get_rounding_mode() == DFP_Field.Rounding_Mode.ROUND_DOWN && nsqd == mant.size())
			{
				// We have enough for this mode
				break;
			}

			if (nsqd > mant.size())
			{
				// We have enough digits
				break;
			}

			/* move the remainder into the dividend while left shifting */
			dividend[0] = 0;
			for (int i{}; i < mant.size(); i++)
			{
				dividend[i + 1] = remainder[i];
			}
		}

		/* Find the most sig digit */
		md = mant.size();  // default
		for (int i = mant.size() + 1; i >= 0; i--)
		{
			if (quotient[i] != 0)
			{
				md = i;
				break;
			}
		}

		/* Copy the digits into the result */
		for (const int& i = 0; i < mant.size(); i++)
		{
			result.mant[mant.size() - i - 1] = quotient[md - i];
		}

		/* Fixup the exponent. */
		result.exp = exp - divisor.exp + md - mant.size();
		result.sign = (byte)((sign == divisor.sign) ? 1 : -1);

		if (result.mant[mant.size() - 1] == 0) { // if result is zero, set exp to zero
			result.exp = 0;
		}

		if (md > (mant.size() - 1))
		{
			excp = result.round(quotient[md - mant.size()]);
		}
		else
		{
			excp = result.round(0);
		}

		if (excp != 0)
		{
			result = dotrap(excp, DIVIDE_TRAP, divisor, result);
		}

		return result;
	}

	/** Divide by a single digit less than radix.
	 *  Special case, so there are speed advantages. 0 &lt;= divisor &lt; radix
	 * @param divisor divisor
	 * @return quotient of this by divisor
	 */
	public Dfp divide(const int& divisor)
	{
		// Handle special cases
		if (nans != FINITE)
		{
			if (is_nan())
			{
				return this;
			}

			if (nans == INFINITE)
			{
				return new_instance(this);
			}
		}

		// Test for divide by zero
		if (divisor == 0)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_DIV_ZERO);
			Dfp result = new_instance(get_zero());
			result.sign = sign;
			result.nans = INFINITE;
			result = dotrap(DFP_Field.FLAG_DIV_ZERO, DIVIDE_TRAP, get_zero(), result);
			return result;
		}

		// range check divisor
		if (divisor < 0 || divisor >= RADIX)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			result = dotrap(DFP_Field.FLAG_INVALID, DIVIDE_TRAP, result, result);
			return result;
		}

		Dfp result = new_instance(this);

		int rl = 0;
		for (int i = mant.size() - 1; i >= 0; i--)
		{
			const int r = rl * RADIX + result.mant[i];
			const int rh = r / divisor;
			rl = r - rh * divisor;
			result.mant[i] = rh;
		}

		if (result.mant[mant.size() - 1] == 0)
		{
			// normalize
			result.shiftLeft();
			const int r = rl * RADIX;        // compute the next digit and put it in
			const int rh = r / divisor;
			rl = r - rh * divisor;
			result.mant[0] = rh;
		}

		const int excp = result.round(rl * RADIX / divisor);  // do the rounding
		if (excp != 0)
		{
			result = dotrap(excp, DIVIDE_TRAP, result, result);
		}

		return result;
	}

	/** {@inherit_doc} */
	//override
	public Dfp reciprocal()
	{
		return field.get_one().divide(this);
	}

	/** Compute the square root.
	 * @return square root of the instance
	 */
	 //override
	public Dfp sqrt()
	{
		// check for unusual cases
		if (nans == FINITE && mant[mant.size() - 1] == 0)
		{
			// if zero
			return new_instance(this);
		}

		if (nans != FINITE)
		{
			if (nans == INFINITE && sign == 1)
			{
				// if positive infinity
				return new_instance(this);
			}

			if (nans == QNAN)
			{
				return new_instance(this);
			}

			if (nans == SNAN)
			{
				Dfp result;

				field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
				result = new_instance(this);
				result = dotrap(DFP_Field.FLAG_INVALID, SQRT_TRAP, NULL, result);
				return result;
			}
		}

		if (sign == -1)
		{
			// if negative
			Dfp result;

			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			result = new_instance(this);
			result.nans = QNAN;
			result = dotrap(DFP_Field.FLAG_INVALID, SQRT_TRAP, NULL, result);
			return result;
		}

		Dfp x = new_instance(this);

		/* Lets make a reasonable guess as to the size of the square root */
		if (x.exp < -1 || x.exp > 1)
		{
			x.exp = this.exp / 2;
		}

		/* Coarsely estimate the mantissa */
		switch (x.mant[mant.size() - 1] / 2000)
		{
		case 0:
			x.mant[mant.size() - 1] = x.mant[mant.size() - 1] / 2 + 1;
			break;
		case 2:
			x.mant[mant.size() - 1] = 1500;
			break;
		case 3:
			x.mant[mant.size() - 1] = 2200;
			break;
		default:
			x.mant[mant.size() - 1] = 3000;
			break;
		}

		/* Now that we have the first pass estimate, compute the rest
	   by the formula dx = (y - x*x) / (2x); */

		Dfp dx;
		Dfp px = get_zero();
		Dfp ppx;
		while (x.unequal(px))
		{
			dx = new_instance(x);
			dx.sign = -1;
			dx = dx.add(this.divide(x));
			dx = dx.divide(2);
			ppx = px;
			px = x;
			x = x.add(dx);

			if (x.equals(ppx))
			{
				// alternating between two values
				break;
			}

			// if dx is zero, break.  Note testing the most sig digit
			// is a sufficient test since dx is normalized
			if (dx.mant[mant.size() - 1] == 0)
			{
				break;
			}
		}

		return x;
	}

	/** Get a string representation of the instance.
	 * @return string representation of the instance
	 */
	 //override
	public std::string to_string() const
	{
		if (nans != FINITE)
		{
			// if non-finite exceptional cases
			if (nans == INFINITE)
			{
				return (sign < 0) ? NEG_INFINITY_STRING : POS_INFINITY_STRING;
			}
			else
			{
				return NAN_STRING;
			}
		}

		if (exp > mant.size() || exp < -1)
		{
			return dfp2sci();
		}

		return dfp2string();
	}

	/** Convert an instance to a string using scientific notation.
	 * @return string representation of the instance in scientific notation
	 */
	protected std::string dfp2sci()
	{
		char rawdigits[] = char[mant.size() * 4];
		int p;
		int e;
		int ae;
		int shf;

		// Get all the digits
		p = 0;
		for (int i = mant.size() - 1; i >= 0; i--)
		{
			rawdigits[p++] = (char)((mant[i] / 1000) + '0');
			rawdigits[p++] = (char)(((mant[i] / 100) % 10) + '0');
			rawdigits[p++] = (char)(((mant[i] / 10) % 10) + '0');
			rawdigits[p++] = (char)(((mant[i]) % 10) + '0');
		}

		// Find the first non-zero one
		for (p = 0; p < rawdigits.size(); p++)
		{
			if (rawdigits[p] != '0')
			{
				break;
			}
		}
		shf = p;

		// Now do the conversion
		std::stringBuilder builder = std::stringstream();
		if (sign == -1)
		{
			builder.append('-');
		}

		if (p != rawdigits.size())
		{
			// there are non zero digits...
			builder.append(rawdigits[p++]);
			builder.append('.');

			while (p < rawdigits.size())
			{
				builder.append(rawdigits[p++]);
			}
		}
		else
		{
			builder.append("0.0e0");
			return builder.to_string();
		}

		builder.append('e');

		// Find the msd of the exponent

		e = exp * 4 - shf - 1;
		ae = e;
		if (e < 0)
		{
			ae = -e;
		}

		// Find the largest p such that p < e
		for (p = 1000000000; p > ae; p /= 10)
		{
			// nothing to do
		}

		if (e < 0)
		{
			builder.append('-');
		}

		while (p > 0)
		{
			builder.append((char)(ae / p + '0'));
			ae %= p;
			p /= 10;
		}

		return builder.to_string();
	}

	/** Convert an instance to a string using normal notation.
	 * @return string representation of the instance in normal notation
	 */
	protected std::string dfp2string()
	{
		const std::string fourZero = "0000";
		int e = exp;
		bool pointInserted = false;

		std::stringBuilder builder = std::stringstream();

		if (e <= 0)
		{
			builder.append("0.");
			pointInserted = true;
		}

		while (e < 0)
		{
			builder.append(fourZero);
			e++;
		}

		for (int i = mant.size() - 1; i >= 0; i--)
		{
			builder.append((char)((mant[i] / 1000) + '0'));
			builder.append((char)(((mant[i] / 100) % 10) + '0'));
			builder.append((char)(((mant[i] / 10) % 10) + '0'));
			builder.append((char)(((mant[i]) % 10) + '0'));
			--e;
			if (e == 0)
			{
				builder.append('.');
				pointInserted = true;
			}
		}

		while (e > 0)
		{
			builder.append(fourZero);
			e--;
		}

		if (!pointInserted)
		{
			// Ensure we have a radix point!
			builder.append('.');
		}

		// Suppress leading zeros
		while (builder.char_at(0) == '0')
		{
			builder.deleteCharAt(0);
		}
		if (builder.char_at(0) == '.')
		{
			builder.insert(0, '0');
		}

		// Suppress trailing zeros
		while (builder.char_at(builder.size()() - 1) == '0')
		{
			builder.deleteCharAt(builder.size()() - 1);
		}

		// Insert sign
		if (sign < 0)
		{
			builder.insert(0, '-');
		}

		return builder.to_string();
	}

	/** Raises a trap.  This does not set the corresponding flag however.
	 *  @param type the trap type
	 *  @param what - name of routine trap occurred in
	 *  @param oper - input operator to function
	 *  @param result - the result computed prior to the trap
	 *  @return The suggested return value from the trap handler
	 */
	public Dfp dotrap(const int& type, std::string what, Dfp oper, Dfp result)
	{
		Dfp def = result;

		switch (type)
		{
		case DFP_Field.FLAG_INVALID:
			def = new_instance(get_zero());
			def.sign = result.sign;
			def.nans = QNAN;
			break;

		case DFP_Field.FLAG_DIV_ZERO:
			if (nans == FINITE && mant[mant.size() - 1] != 0)
			{
				// normal case, we are finite, non-zero
				def = new_instance(get_zero());
				def.sign = (byte)(sign * oper.sign);
				def.nans = INFINITE;
			}

			if (nans == FINITE && mant[mant.size() - 1] == 0)
			{
				//  0/0
				def = new_instance(get_zero());
				def.nans = QNAN;
			}

			if (nans == INFINITE || nans == QNAN)
			{
				def = new_instance(get_zero());
				def.nans = QNAN;
			}

			if (nans == INFINITE || nans == SNAN)
			{
				def = new_instance(get_zero());
				def.nans = QNAN;
			}
			break;

		case DFP_Field.FLAG_UNDERFLOW:
			if ((result.exp + mant.size()) < MIN_EXP)
			{
				def = new_instance(get_zero());
				def.sign = result.sign;
			}
			else
			{
				def = new_instance(result);  // gradual underflow
			}
			result.exp += ERR_SCALE;
			break;

		case DFP_Field.FLAG_OVERFLOW:
			result.exp -= ERR_SCALE;
			def = new_instance(get_zero());
			def.sign = result.sign;
			def.nans = INFINITE;
			break;

		default: def = result; break;
		}

		return trap(type, what, oper, def, result);
	}

	/** Trap handler.  Subclasses may //override this to provide trap
	 *  functionality per IEEE 854-1987.
	 *
	 *  @param type  The exception type - e.g. FLAG_OVERFLOW
	 *  @param what  The name of the routine we were in e.g. divide()
	 *  @param oper  An operand to this function if any
	 *  @param def   The default return value if trap not enabled
	 *  @param result    The result that is specified to be delivered per
	 *                   IEEE 854, if any
	 *  @return the value that should be return by the operation triggering the trap
	 */
	protected Dfp trap(const int& type, std::string what, Dfp oper, Dfp def, Dfp result)
	{
		return def;
	}

	/** Returns the type - one of FINITE, INFINITE, SNAN, QNAN.
	 * @return type of the number
	 */
	public int classify()
	{
		return nans;
	}

	/** Creates an instance that is the same as x except that it has the sign of y.
	 * abs(x) = dfp.copysign(x, dfp.one)
	 * @param x number to get the value from
	 * @param y number to get the sign from
	 * @return a number with the value of x and the sign of y
	 */
	public static Dfp copysign(const Dfp x, const Dfp y)
	{
		Dfp result = x.new_instance(x);
		result.sign = y.sign;
		return result;
	}

	/** Returns the next number greater than this one in the direction of x.
	 * If this==x then simply returns this.
	 * @param x direction where to look at
	 * @return closest number next to instance in the direction of x
	 */
	public Dfp next_after(const Dfp x)
	{
		// make sure we don't mix number with different precision
		if (field.get_radix_digits() != x.field.get_radix_digits())
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, NEXT_AFTER_TRAP, x, result);
		}

		// if this is greater than x
		bool up = false;
		if (this.less_than(x))
		{
			up = true;
		}

		if (compare(this, x) == 0)
		{
			return new_instance(x);
		}

		if (less_than(get_zero()))
		{
			up = !up;
		}

		const Dfp inc;
		Dfp result;
		if (up)
		{
			inc = new_instance(get_one());
			inc.exp = this.exp - mant.size() + 1;
			inc.sign = this.sign;

			if (this.equals(get_zero()))
			{
				inc.exp = MIN_EXP - mant.size();
			}

			result = add(inc);
		}
		else
		{
			inc = new_instance(get_one());
			inc.exp = this.exp;
			inc.sign = this.sign;

			if (this.equals(inc))
			{
				inc.exp = this.exp - mant.size();
			}
			else
			{
				inc.exp = this.exp - mant.size() + 1;
			}

			if (this.equals(get_zero()))
			{
				inc.exp = MIN_EXP - mant.size();
			}

			result = this.subtract(inc);
		}

		if (result.classify() == INFINITE && this.classify() != INFINITE)
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			result = dotrap(DFP_Field.FLAG_INEXACT, NEXT_AFTER_TRAP, x, result);
		}

		if (result.equals(get_zero()) && !this.equals(get_zero()))
		{
			field.set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			result = dotrap(DFP_Field.FLAG_INEXACT, NEXT_AFTER_TRAP, x, result);
		}

		return result;
	}

	/** Convert the instance into a double.
	 * @return a double approximating the instance
	 * @see #toSplitDouble()
	 */
	public double toDouble()
	{
		if (std::isinfinite())
		{
			if (less_than(get_zero()))
			{
				return -INFINITY;
			}
			else
			{
				return INFINITY;
			}
		}

		if (is_nan())
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		Dfp y = this;
		bool negate = false;
		int cmp0 = compare(this, get_zero());
		if (cmp0 == 0)
		{
			return sign < 0 ? -0.0 : +0.0;
		}
		else if (cmp0 < 0)
		{
			y = negate();
			negate = true;
		}

		/* Find the exponent, first estimate by integer log10, then adjust.
		 Should be faster than doing a natural logarithm.  */
		int exponent = static_cast<int>((y.int_log10() * 3.32);
		if (exponent < 0)
		{
			exponent--;
		}

		Dfp tempDfp = Dfp_Math.pow(get_two(), exponent);
		while (tempDfp.less_than(y) || tempDfp.equals(y))
		{
			tempDfp = tempDfp.multiply(2);
			exponent++;
		}
		exponent--;

		/* We have the exponent, now work on the mantissa */

		y = y.divide(Dfp_Math.pow(get_two(), exponent));
		if (exponent > -1023)
		{
			y = y.subtract(get_one());
		}

		if (exponent < -1074)
		{
			return 0;
		}

		if (exponent > 1023)
		{
			return negate ? -INFINITY : INFINITY;
		}

		y = y.multiply(new_instance(4503599627370496l)).rint();
		std::string str = y.to_string();
		str = str.substring(0, str.size()() - 1);
		long mantissa = long.parseLong(str);

		if (mantissa == 4503599627370496L)
		{
			// Handle special case where we round up to next power of two
			mantissa = 0;
			exponent++;
		}

		/* Its going to be subnormal, so make adjustments */
		if (exponent <= -1023)
		{
			exponent--;
		}

		while (exponent < -1023)
		{
			exponent++;
			mantissa >>>= 1;
		}

		long bits = mantissa | ((exponent + 1023L) << 52);
		double x = Double.long_bits_to_double(bits);

		if (negate)
		{
			x = -x;
		}

		return x;
	}

	/** Convert the instance into a split double.
	 * @return an array of two doubles which sum represent the instance
	 * @see #toDouble()
	 */
	public std::vector<double> toSplitDouble()
	{
		double split[] = std::vector<double>(2);
		long mask = 0xffffffffc0000000L;

		split[0] = Double.long_bits_to_double(Double.double_to_long_bits(toDouble()) & mask);
		split[1] = subtract(new_instance(split[0])).toDouble();

		return split;
	}

	/** {@inherit_doc}
	 */
	 //override
	public double get_real()
	{
		return toDouble();
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp add(const double& a)
	{
		return add(new_instance(a));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp subtract(const double& a)
	{
		return subtract(new_instance(a));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp multiply(const double& a)
	{
		return multiply(new_instance(a));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp divide(const double& a)
	{
		return divide(new_instance(a));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp remainder(const double& a)
	{
		return remainder(new_instance(a));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp sign()
	{
		if (is_nan() || is_zero())
		{
			return this;
		}
		else
		{
			return new_instance(sign > 0 ? +1 : -1);
		}
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp copy_sign(const Dfp s)
	{
		if ((sign >= 0 && s.sign >= 0) || (sign < 0 && s.sign < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp copy_sign(const double s)
	{
		long sb = Double.double_to_long_bits(s);
		if ((sign >= 0 && sb >= 0) || (sign < 0 && sb < 0)) { // Sign is currently OK
			return this;
		}
		return negate(); // flip sign
	}

	/** {@inherit_doc}
	 */
	 //override
	public int get_exponent()
	{
		if (nans != FINITE)
		{
			// 2⁴³⁵⁴¹¹ < 10000³²⁷⁶⁸ < 2⁴³⁵⁴¹²
			return 435411;
		}
		if (is_zero())
		{
			return -435412;
		}

		const Dfp abs = abs();

		// estimate a lower bound for binary exponent
		// 13301/1001 is a continued fraction approximation of ln(10000)/ln(2)
		int p = std::max(13301 * exp / 1001 - 15, -435411);
		Dfp twoP = Dfp_Math.pow(get_two(), p);
		while (compare(abs, twoP) >= 0)
		{
			twoP = twoP.add(twoP);
			++p;
		}

		return p - 1;
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp scalb(const int& n)
	{
		return multiply(Dfp_Math.pow(get_two(), n));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp ulp()
	{
		const Dfp result = Dfp(field);
		result.mant[result.mant.size() - 1] = 1;
		result.exp = exp - (result.mant.size() - 1);
		return result;
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp hypot(const Dfp y)
	{
		if (std::isinfinite() || y.std::isinfinite())
		{
			return field.new_dfp(INFINITY);
		}
		else if (is_nan() || y.is_nan())
		{
			return field.new_dfp(Double.NaN);
		}
		else
		{
			// find scaling to avoid both overflow and underflow
			const int scalingExp = (exp + y.exp) / 2;

			// scale both operands
			const Dfp scaled_x = Dfp(this);
			scaled_x.exp -= scalingExp;
			const Dfp scaled_y = Dfp(y);
			scaled_y.exp -= scalingExp;

			// compute scaled hypothenuse
			const Dfp h = scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

			// scale result
			h.exp += scalingExp;

			return h;
		}
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp cbrt()
	{
		return root_n(3);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp root_n(const int& n)
	{
		return (sign >= 0) ?
			Dfp_Math.pow(this, get_one().divide(n)) :
			Dfp_Math.pow(negate(), get_one().divide(n)).negate();
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp pow(const double& p)
	{
		return Dfp_Math.pow(this, new_instance(p));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp pow(const int& n)
	{
		return Dfp_Math.pow(this, n);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp pow(const Dfp e)
	{
		return Dfp_Math.pow(this, e);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp exp()
	{
		return Dfp_Math.exp(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp expm1()
	{
		return Dfp_Math.exp(this).subtract(get_one());
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp log()
	{
		return Dfp_Math.log(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp log1p()
	{
		return Dfp_Math.log(this.add(get_one()));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp log10()
	{
		return Dfp_Math.log(this).divide(Dfp_Math.log(new_instance(10)));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp cos()
	{
		return Dfp_Math.cos(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp sin()
	{
		return Dfp_Math.sin(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Field_Sin_Cos<Dfp> sin_cos()
	{
		return Field_Sin_Cos<>(Dfp_Math.sin(this), Dfp_Math.cos(this));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp tan()
	{
		return Dfp_Math.tan(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp acos()
	{
		return Dfp_Math.acos(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp asin()
	{
		return Dfp_Math.asin(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp atan()
	{
		return Dfp_Math.atan(this);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp atan2(const Dfp x)

	{
		// compute r = sqrt(x^2+y^2)
		const Dfp r = x.multiply(x).add(multiply(this)).sqrt();
		if (r.is_zero())
		{
			// special cases handling
			if (x.sign >= 0)
			{
				return this; // ±0.0
			}
			else
			{
				return new_instance((sign <= 0) ? -std::numbers::pi : std::numbers::pi); // ±π
			}
		}

		if (x.sign >= 0)
		{
			// compute atan2(y, x) = 2 atan(y / (r + x))
			return get_two().multiply(divide(r.add(x)).atan());
		}
		else
		{
			// compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
			const Dfp tmp = get_two().multiply(divide(r.subtract(x)).atan());
			const Dfp pmPi = new_instance((tmp.sign <= 0) ? -std::numbers::pi : std::numbers::pi);
			return pmPi.subtract(tmp);
		}
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp cosh()
	{
		return Dfp_Math.exp(this).add(Dfp_Math.exp(negate())).multiply(0.5);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp sinh()
	{
		return Dfp_Math.exp(this).subtract(Dfp_Math.exp(negate())).multiply(0.5);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Field_Sinh_Cosh<Dfp> sinh_cosh()
	{
		const Dfp p = Dfp_Math.exp(this);
		const Dfp m = Dfp_Math.exp(negate());
		return Field_Sinh_Cosh<>(p.subtract(m).multiply(0.5), p.add(m).multiply(0.5));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp tanh()
	{
		const Dfp ePlus = Dfp_Math.exp(this);
		const Dfp eMinus = Dfp_Math.exp(negate());
		return ePlus.subtract(eMinus).divide(ePlus.add(eMinus));
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp acosh()
	{
		return multiply(this).subtract(get_one()).sqrt().add(this).log();
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp asinh()
	{
		return multiply(this).add(get_one()).sqrt().add(this).log();
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp atanh()
	{
		return get_one().add(this).divide(get_one().subtract(this)).log().divide(2);
	}

	/** {@inherit_doc} */
	//override
	public Dfp to_degrees()
	{
		return multiply(field.get_rad_to_deg());
	}

	/** {@inherit_doc} */
	//override
	public Dfp to_radians()
	{
		return multiply(field.get_deg_to_rad());
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const Dfp[] a, const Dfp[] b)

	{
		Math_Utils::check_dimension(a.size(), b.size());

		// compute in extended accuracy
		const DFP_Field extendedField = a[0].field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		Dfp r = extendedField.get_zero();
		for (int i{}; i < a.size(); ++i)
		{
			const Dfp aiExt = a[i].new_instance(extendedField, NULL);
			const Dfp biExt = b[i].new_instance(extendedField, NULL);
			r = r.add(aiExt.multiply(biExt));
		}

		// back to normal accuracy
		return r.new_instance(a[0].field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const std::vector<double> a, const Dfp[] b)

	{
		Math_Utils::check_dimension(a.size(), b.size());

		// compute in extended accuracy
		const DFP_Field extendedField = b[0].field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		Dfp r = extendedField.get_zero();
		for (int i{}; i < a.size(); ++i)
		{
			const Dfp biExt = b[i].new_instance(extendedField, NULL);
			r = r.add(biExt.multiply(a[i]));
		}

		// back to normal accuracy
		return r.new_instance(b[0].field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const Dfp a1, const Dfp b1, const Dfp a2, const Dfp b2)
	{
		// switch to extended accuracy
		const DFP_Field extendedField = a1.field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		const Dfp a1Ext = a1.new_instance(extendedField, NULL);
		const Dfp b1Ext = b1.new_instance(extendedField, NULL);
		const Dfp a2Ext = a2.new_instance(extendedField, NULL);
		const Dfp b2Ext = b2.new_instance(extendedField, NULL);

		// compute linear combination in extended accuracy
		const Dfp resultExt = a1Ext.multiply(b1Ext).
			add(a2Ext.multiply(b2Ext));

		// back to normal accuracy
		return resultExt.new_instance(a1.field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const double& a1, const Dfp b1, const double& a2, const Dfp b2)
	{
		// switch to extended accuracy
		const DFP_Field extendedField = b1.field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		const Dfp b1Ext = b1.new_instance(extendedField, NULL);
		const Dfp b2Ext = b2.new_instance(extendedField, NULL);

		// compute linear combination in extended accuracy
		const Dfp resultExt = b1Ext.multiply(a1).
			add(b2Ext.multiply(a2));

		// back to normal accuracy
		return resultExt.new_instance(b1.field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const Dfp a1, const Dfp b1, const Dfp a2, const Dfp b2, const Dfp a3, const Dfp b3)
	{
		// switch to extended accuracy
		const DFP_Field extendedField = a1.field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		const Dfp a1Ext = a1.new_instance(extendedField, NULL);
		const Dfp b1Ext = b1.new_instance(extendedField, NULL);
		const Dfp a2Ext = a2.new_instance(extendedField, NULL);
		const Dfp b2Ext = b2.new_instance(extendedField, NULL);
		const Dfp a3Ext = a3.new_instance(extendedField, NULL);
		const Dfp b3Ext = b3.new_instance(extendedField, NULL);

		// compute linear combination in extended accuracy
		const Dfp resultExt = a1Ext.multiply(b1Ext).
			add(a2Ext.multiply(b2Ext)).
			add(a3Ext.multiply(b3Ext));

		// back to normal accuracy
		return resultExt.new_instance(a1.field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const double& a1, const Dfp b1, const double& a2, const Dfp b2, const double& a3, const Dfp b3)
	{
		// switch to extended accuracy
		const DFP_Field extendedField = b1.field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		const Dfp b1Ext = b1.new_instance(extendedField, NULL);
		const Dfp b2Ext = b2.new_instance(extendedField, NULL);
		const Dfp b3Ext = b3.new_instance(extendedField, NULL);

		// compute linear combination in extended accuracy
		const Dfp resultExt = b1Ext.multiply(a1).
			add(b2Ext.multiply(a2)).
			add(b3Ext.multiply(a3));

		// back to normal accuracy
		return resultExt.new_instance(b1.field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const Dfp a1, const Dfp b1, const Dfp a2, const Dfp b2, const Dfp a3, const Dfp b3, const Dfp a4, const Dfp b4)
	{
		// switch to extended accuracy
		const DFP_Field extendedField = a1.field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		const Dfp a1Ext = a1.new_instance(extendedField, NULL);
		const Dfp b1Ext = b1.new_instance(extendedField, NULL);
		const Dfp a2Ext = a2.new_instance(extendedField, NULL);
		const Dfp b2Ext = b2.new_instance(extendedField, NULL);
		const Dfp a3Ext = a3.new_instance(extendedField, NULL);
		const Dfp b3Ext = b3.new_instance(extendedField, NULL);
		const Dfp a4Ext = a4.new_instance(extendedField, NULL);
		const Dfp b4Ext = b4.new_instance(extendedField, NULL);

		// compute linear combination in extended accuracy
		const Dfp resultExt = a1Ext.multiply(b1Ext).
			add(a2Ext.multiply(b2Ext)).
			add(a3Ext.multiply(b3Ext)).
			add(a4Ext.multiply(b4Ext));

		// back to normal accuracy
		return resultExt.new_instance(a1.field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc}
	 */
	 //override
	public Dfp linear_combination(const double& a1, const Dfp b1, const double& a2, const Dfp b2, const double& a3, const Dfp b3, const double& a4, const Dfp b4)
	{
		// switch to extended accuracy
		const DFP_Field extendedField = b1.field.get_extended_field(LINEAR_COMBINATION_DIGITS_FACTOR, false);
		const Dfp b1Ext = b1.new_instance(extendedField, NULL);
		const Dfp b2Ext = b2.new_instance(extendedField, NULL);
		const Dfp b3Ext = b3.new_instance(extendedField, NULL);
		const Dfp b4Ext = b4.new_instance(extendedField, NULL);

		// compute linear combination in extended accuracy
		const Dfp resultExt = b1Ext.multiply(a1).
			add(b2Ext.multiply(a2)).
			add(b3Ext.multiply(a3)).
			add(b4Ext.multiply(a4));

		// back to normal accuracy
		return resultExt.new_instance(b1.field, DFP_Field.Rounding_Mode.ROUND_HALF_EVEN);
	}

	/** {@inherit_doc} */
	//override
	public Dfp get_pi()
	{
		return field.get_pi();
	}
}
