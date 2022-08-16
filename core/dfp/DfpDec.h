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

  /** Subclass of {@link Dfp} which hides the radix-10000 artifacts of the superclass.
   * This should give outward appearances of being a decimal number with DIGITS*4-3
   * decimal digits. This class can be subclassed to appear to be an arbitrary number
   * of decimal digits less than DIGITS*4-3.
   */
class Dfp_Dec extends Dfp
{
	/** Makes an instance with a value of zero.
	 * @param factory factory linked to this instance
	 */
	protected Dfp_Dec(const DFP_Field factory)
	{
		super(factory);
	}

	/** Create an instance from a std::byte value.
	 * @param factory factory linked to this instance
	 * @param x value to convert to an instance
	 */
	protected Dfp_Dec(const DFP_Field factory, std::byte x)
	{
		super(factory, x);
	}

	/** Create an instance from an int value.
	 * @param factory factory linked to this instance
	 * @param x value to convert to an instance
	 */
	protected Dfp_Dec(const DFP_Field factory, int x)
	{
		super(factory, x);
	}

	/** Create an instance from a long value.
	 * @param factory factory linked to this instance
	 * @param x value to convert to an instance
	 */
	protected Dfp_Dec(const DFP_Field factory, long x)
	{
		super(factory, x);
	}

	/** Create an instance from a double value.
	 * @param factory factory linked to this instance
	 * @param x value to convert to an instance
	 */
	protected Dfp_Dec(const DFP_Field factory, double x)
	{
		super(factory, x);
		round(0);
	}

	/** Copy constructor.
	 * @param d instance to copy
	 */
	public Dfp_Dec(const Dfp d)
	{
		super(d);
		round(0);
	}

	/** Create an instance from a std::string representation.
	 * @param factory factory linked to this instance
	 * @param s string representation of the instance
	 */
	protected Dfp_Dec(const DFP_Field factory, const std::string s)
	{
		super(factory, s);
		round(0);
	}

	/** Creates an instance with a non-finite value.
	 * @param factory factory linked to this instance
	 * @param sign sign of the Dfp to create
	 * @param nans code of the value, must be one of {@link #INFINITE}, * {@link #SNAN},  {@link #QNAN}
	 */
	protected Dfp_Dec(const DFP_Field factory, const std::byte sign, const std::byte nans)
	{
		super(factory, sign, nans);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance()
	{
		return Dfp_Dec(get_field());
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const std::byte x)
	{
		return Dfp_Dec(get_field(), x);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const int x)
	{
		return Dfp_Dec(get_field(), x);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const long x)
	{
		return Dfp_Dec(get_field(), x);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const double& x)
	{
		return Dfp_Dec(get_field(), x);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const Dfp d)
	{
		// make sure we don't mix number with different precision
		if (get_field().get_radix_digits() != d.get_field().get_radix_digits())
		{
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, "new_instance", d, result);
		}

		return Dfp_Dec(d);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const std::string s)
	{
		return Dfp_Dec(get_field(), s);
	}

	/** {@inherit_doc} */
	//override
	public Dfp new_instance(const std::byte sign, const std::byte nans)
	{
		return Dfp_Dec(get_field(), sign, nans);
	}

	/** Get the number of decimal digits this class is going to represent.
	 * Default implementation returns {@link #get_radix_digits()}*4-3. Subclasses can
	 * //override this to return something less.
	 * @return number of decimal digits this class is going to represent
	 */
	protected int get_decimal_digits()
	{
		return get_radix_digits() * 4 - 3;
	}

	/** {@inherit_doc} */
	//override
	protected int round(const int& in)
	{
		int msb = mant[mant.size() - 1];
		if (msb == 0)
		{
			// special case -- this == zero
			return 0;
		}

		int cmaxdigits = mant.size() * 4;
		int lsbthreshold = 1000;
		while (lsbthreshold > msb)
		{
			lsbthreshold /= 10;
			cmaxdigits--;
		}

		const int digits = get_decimal_digits();
		const int lsbshift = cmaxdigits - digits;
		const int lsd = lsbshift / 4;

		lsbthreshold = 1;
		for (int i{}; i < lsbshift % 4; i++)
		{
			lsbthreshold *= 10;
		}

		const int lsb = mant[lsd];

		if (lsbthreshold <= 1 && digits == 4 * mant.size() - 3)
		{
			return super.round(in);
		}

		int discarded = in;  // not looking at this after this point
		const int& n;
		if (lsbthreshold == 1)
		{
			// look to the next digit for rounding
			n = (mant[lsd - 1] / 1000) % 10;
			mant[lsd - 1] %= 1000;
			discarded |= mant[lsd - 1];
		}
		else
		{
			n = (lsb * 10 / lsbthreshold) % 10;
			discarded |= lsb % (lsbthreshold / 10);
		}

		for (int i{}; i < lsd; i++)
		{
			discarded |= mant[i];    // need to know if there are any discarded bits
			mant[i] = 0;
		}

		mant[lsd] = lsb / lsbthreshold * lsbthreshold;

		const bool inc;
		switch (get_field().get_rounding_mode())
		{
		case ROUND_DOWN:
			inc = false;
			break;

		case ROUND_UP:
			inc = (n != 0) || (discarded != 0); // round up if n!=0
			break;

		case ROUND_HALF_UP:
			inc = n >= 5;  // round half up
			break;

		case ROUND_HALF_DOWN:
			inc = n > 5;  // round half down
			break;

		case ROUND_HALF_EVEN:
			inc = (n > 5) ||
				(n == 5 && discarded != 0) ||
				(n == 5 && discarded == 0 && ((lsb / lsbthreshold) & 1) == 1);  // round half-even
			break;

		case ROUND_HALF_ODD:
			inc = (n > 5) ||
				(n == 5 && discarded != 0) ||
				(n == 5 && discarded == 0 && ((lsb / lsbthreshold) & 1) == 0);  // round half-odd
			break;

		case ROUND_CEIL:
			inc = (sign == 1) && (n != 0 || discarded != 0);  // round ceil
			break;

		case ROUND_FLOOR:
		default:
			inc = (sign == -1) && (n != 0 || discarded != 0);  // round floor
			break;
		}

		if (inc)
		{
			// increment if necessary
			int rh = lsbthreshold;
			for (int i = lsd; i < mant.size(); i++)
			{
				const int r = mant[i] + rh;
				rh = r / RADIX;
				mant[i] = r % RADIX;
			}

			if (rh != 0)
			{
				shift_right();
				mant[mant.size() - 1] = rh;
			}
		}

		// Check for exceptional cases and raise signals if necessary
		if (exp < MIN_EXP)
		{
			// Gradual Underflow
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_UNDERFLOW);
			return DFP_Field.FLAG_UNDERFLOW;
		}

		if (exp > MAX_EXP)
		{
			// Overflow
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_OVERFLOW);
			return DFP_Field.FLAG_OVERFLOW;
		}

		if (n != 0 || discarded != 0)
		{
			// Inexact
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			return DFP_Field.FLAG_INEXACT;
		}
		return 0;
	}

	/** {@inherit_doc} */
	//override
	public Dfp next_after(Dfp x)
	{
		const std::string trap_name = "next_after";

		// make sure we don't mix number with different precision
		if (get_field().get_radix_digits() != x.get_field().get_radix_digits())
		{
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_INVALID);
			const Dfp result = new_instance(get_zero());
			result.nans = QNAN;
			return dotrap(DFP_Field.FLAG_INVALID, trap_name, x, result);
		}

		bool up = false;
		Dfp result;
		Dfp inc;

		// if this is greater than x
		if (this.less_than(x))
		{
			up = true;
		}

		if (equals(x))
		{
			return new_instance(x);
		}

		if (less_than(get_zero()))
		{
			up = !up;
		}

		if (up)
		{
			inc = power10(int_log10() - get_decimal_digits() + 1);
			inc = copysign(inc, this);

			if (this.equals(get_zero()))
			{
				inc = power10_k(MIN_EXP - mant.size() - 1);
			}

			if (inc.equals(get_zero()))
			{
				result = copysign(new_instance(get_zero()), this);
			}
			else
			{
				result = add(inc);
			}
		}
		else
		{
			inc = power10(int_log10());
			inc = copysign(inc, this);

			if (this.equals(inc))
			{
				inc = inc.divide(power10(get_decimal_digits()));
			}
			else
			{
				inc = inc.divide(power10(get_decimal_digits() - 1));
			}

			if (this.equals(get_zero()))
			{
				inc = power10_k(MIN_EXP - mant.size() - 1);
			}

			if (inc.equals(get_zero()))
			{
				result = copysign(new_instance(get_zero()), this);
			}
			else
			{
				result = subtract(inc);
			}
		}

		if (result.classify() == INFINITE && this.classify() != INFINITE)
		{
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			result = dotrap(DFP_Field.FLAG_INEXACT, trap_name, x, result);
		}

		if (result.equals(get_zero()) && !this.equals(get_zero()))
		{
			get_field().set_ieee_flags_bits(DFP_Field.FLAG_INEXACT);
			result = dotrap(DFP_Field.FLAG_INEXACT, trap_name, x, result);
		}

		return result;
	}
}
