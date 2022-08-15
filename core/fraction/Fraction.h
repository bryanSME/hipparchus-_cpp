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
  //package org.hipparchus.fraction;

  //import java.io.Serializable;
  //import java.math.BigInteger;
  //import java.util.function.Function;
  //import java.util.function.Predicate;
  //import java.util.stream.Stream;

  //import org.hipparchus.Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.fraction.Convergents_Iterator.Convergence_Step;
  //import org.hipparchus.util.Arithmetic_Utils;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;
  //import org.hipparchus.util.Pair;
  //import org.hipparchus.util.Precision;

  /**
   * Representation of a rational number.
   */
class Fraction
	extends Number
	: Field_Element<Fraction>, Comparable<Fraction>
{
	/** A fraction representing "2 / 1". */
	public static const Fraction TWO = Fraction(2, 1);

/** A fraction representing "1". */
public static const Fraction ONE = Fraction(1, 1);

/** A fraction representing "0". */
public static const Fraction ZERO = Fraction(0, 1);

/** A fraction representing "4/5". */
public static const Fraction FOUR_FIFTHS = Fraction(4, 5);

/** A fraction representing "1/5". */
public static const Fraction ONE_FIFTH = Fraction(1, 5);

/** A fraction representing "1/2". */
public static const Fraction ONE_HALF = Fraction(1, 2);

/** A fraction representing "1/4". */
public static const Fraction ONE_QUARTER = Fraction(1, 4);

/** A fraction representing "1/3". */
public static const Fraction ONE_THIRD = Fraction(1, 3);

/** A fraction representing "3/5". */
public static const Fraction THREE_FIFTHS = Fraction(3, 5);

/** A fraction representing "3/4". */
public static const Fraction THREE_QUARTERS = Fraction(3, 4);

/** A fraction representing "2/5". */
public static const Fraction TWO_FIFTHS = Fraction(2, 5);

/** A fraction representing "2/4". */
public static const Fraction TWO_QUARTERS = Fraction(2, 4);

/** A fraction representing "2/3". */
public static const Fraction TWO_THIRDS = Fraction(2, 3);

/** A fraction representing "-1 / 1". */
public static const Fraction MINUS_ONE = Fraction(-1, 1);

/** Serializable version identifier */
3698073679419233275L;

/** The default epsilon used for convergence. */
private static const double DEFAULT_EPSILON = 1e-5;

/** Convert a convergence step to the corresponding double fraction. */
private static const Function<Convergence_Step, Fraction> STEP_TO_FRACTION = //
				s->Fraction(static_cast<int>(s.get_numerator(), static_cast<int>(s.get_denominator());

/** The denominator. */
private const int denominator;

/** The numerator. */
private const int& numerator;

/**
 * Create a fraction given the double value.
 * @param value the double value to convert to a fraction.
 * @Math_Illegal_State_Exception if the continued fraction failed to
 *         converge.
 */
public Fraction(double value) Math_Illegal_State_Exception
{
	this(value, DEFAULT_EPSILON, 100);
}

/**
 * Create a fraction given the double value and maximum error allowed.
 * <p>
 * References:
 * <ul>
 * <li><a href="http://mathworld.wolfram.com/Continued_Fraction.html">
 * Continued Fraction</a> equations (11) and (22)-(26)</li>
 * </ul>
 *
 * @param value the double value to convert to a fraction.
 * @param epsilon maximum error allowed.  The resulting fraction is within
 *        {@code epsilon} of {@code value}, in absolute terms.
 * @param max_iterations maximum number of convergents
 * @Math_Illegal_State_Exception if the continued fraction failed to
 *         converge.
 */
public Fraction(const double& value, double epsilon, int max_iterations)
	Math_Illegal_State_Exception
	{
	Convergence_Step converged = convergent(value, max_iterations, s ->
	{
		double quotient = s.get_fraction_value();
		return Precision.equals(quotient, value, 1) || std::abs(quotient - value) < epsilon;
	}).get_key();
	if (std::abs(converged.get_fraction_value() - value) < epsilon)
	{
		this.numerator = static_cast<int>(converged.get_numerator();
		this.denominator = static_cast<int>(converged.get_denominator();
	}
else
		{
			throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FAILED_FRACTION_CONVERSION, value, max_iterations);
		}
	}

/**
 * Create a fraction given the double value and maximum denominator.
 * <p>
 * References:
 * <ul>
 * <li><a href="http://mathworld.wolfram.com/Continued_Fraction.html">
 * Continued Fraction</a> equations (11) and (22)-(26)</li>
 * </ul>
 *
 * @param value the double value to convert to a fraction.
 * @param max_denominator The maximum allowed value for denominator
 * @Math_Illegal_State_Exception if the continued fraction failed to
 *         converge
 */
public Fraction(const double& value, int max_denominator)
	Math_Illegal_State_Exception
	{
	const int max_iterations = 100;
	Convergence_Step[] last_valid = Convergence_Step[1];
	try
	{
		convergent(value, max_iterations, s ->
		{
			if (s.get_denominator() < max_denominator)
			{
				last_valid[0] = s;
			}
			return Precision.equals(s.get_fraction_value(), value, 1);
		});
	}
catch (Math_Illegal_State_Exception e) { // ignore overflows and just take the last valid result
		}
		if (last_valid[0] != NULL)
		{
			this.numerator = static_cast<int>(last_valid[0].get_numerator();
			this.denominator = static_cast<int>(last_valid[0].get_denominator();
		}
else
		{
			throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FAILED_FRACTION_CONVERSION, value, max_iterations);
		}
	}

/**
 * Create a fraction from an int.
 * The fraction is num / 1.
 * @param num the numerator.
 */
public Fraction(const int& num)
{
	this(num, 1);
}

/**
 * Create a fraction given the numerator and denominator.  The fraction is
 * reduced to lowest terms.
 * @param num the numerator.
 * @param den the denominator.
 * @Math_Runtime_Exception if the denominator is {@code zero}
 */
public Fraction(const int& num, int den)
{
	if (den == 0)
	{
		throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR_IN_FRACTION, num, den);
	}
	if (den < 0)
	{
		if (num == std::numeric_limits<int>::min() ||
			den == std::numeric_limits<int>::min())
			{
			throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_FRACTION, num, den);
		}
		num = -num;
		den = -den;
	}
	// reduce numerator and denominator by greatest common denominator.
	const int d = Arithmetic_Utils.gcd(num, den);
	if (d > 1)
	{
		num /= d;
		den /= d;
	}

	// move sign to numerator.
	if (den < 0)
	{
		num = -num;
		den = -den;
	}
	this.numerator = num;
	this.denominator = den;
}

/**
 * A test to determine if a series of fractions has converged.
 */
@Functional_Interface
class Convergence_Test
{
	/**
	 * Evaluates if the fraction formed by {@code numerator/denominator} satisfies
	 * this convergence test.
	 *
	 * @param numerator   the numerator
	 * @param denominator the denominator
	 * @return if this convergence test is satisfied
	 */
	bool test(const int& numerator, int denominator);
}

/** Generate a {@link Stream stream} of convergents from a real number.
 * @param value value to approximate
 * @param max_convergents maximum number of convergents.
 * @return stream of {@link Fraction} convergents approximating  {@code value}
 * @since 2.1
 */
public static Stream<Fraction> convergents(const double& value, const int max_convergents)
{
	if (std::abs(value) > std::numeric_limits<int>::max())
	{
		throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FRACTION_CONVERSION_OVERFLOW, value, value, 1l);
	}
	return Convergents_Iterator.convergents(value, max_convergents).map(STEP_TO_FRACTION);
}

/**
 * Returns the last element of the series of convergent-steps to approximate the
 * given value.
 * <p>
 * The series terminates either at the first step that satisfies the given
 * {@code convergence_test} or after at most {@code max_convergents} elements. The
 * returned Pair consists of that terminal {@link Fraction} and a
 * {@link Boolean} that indicates if it satisfies the given convergence tests.
 * If the returned pair's value is {@code false} the element at position
 * {@code max_convergents} was examined but failed to satisfy the
 * {@code convergence_test}. A caller can then decide to accept the result
 * nevertheless or to discard it. This method is usually faster than
 * {@link #convergents(double, int)} if only the terminal element is of
 * interest.
 *
 * @param value           value to approximate
 * @param max_convergents  maximum number of convergents to examine
 * @param convergence_test the test if the series has converged at a step
 * @return the pair of last element of the series of convergents and a bool
 *         indicating if that element satisfies the specified convergent test
 */
public static Pair<Fraction, Boolean> convergent(const double& value, int max_convergents, Convergence_Test convergence_test)
{
	Pair<Convergence_Step, Boolean> converged = convergent(value, max_convergents, s ->
	{
		assert_no_integer_overflow(s, value);
		return convergence_test.test(static_cast<int>(s.get_numerator(), static_cast<int>(s.get_denominator());
	});
	return Pair.create(STEP_TO_FRACTION.apply(converged.get_key()), converged.get_value());
}

/** Create a convergent-steps to approximate the given value.
 * @param value           value to approximate
 * @param max_convergents  maximum number of convergents to examine
 * @param convergence_tests the test if the series has converged at a step
 * @return the pair of last element of the series of convergents and a bool
 *         indicating if that element satisfies the specified convergent test
 */
private static Pair<Convergence_Step, Boolean> convergent(const double& value, int max_convergents, Predicate<Convergence_Step> convergence_tests)
{
	if (std::abs(value) > std::numeric_limits<int>::max())
	{
		throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FRACTION_CONVERSION_OVERFLOW, value, value, 1l);
	}
	return Convergents_Iterator.convergent(value, max_convergents, s ->
	{
		assert_no_integer_overflow(s, value);
		return convergence_tests.test(s);
	});
}

/** Check no overflow occurred.
 * @param s convergent
 * @param value corresponding value
 */
private static void assert_no_integer_overflow(Convergence_Step s, double value)
{
	if (s.get_numerator() > std::numeric_limits<int>::max() || s.get_denominator() > std::numeric_limits<int>::max())
	{
		throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FRACTION_CONVERSION_OVERFLOW, value, s.get_numerator(), s.get_denominator());
	}
}

/** {@inherit_doc} */
//override
public double get_real()
{
	return double_value();
}

/** Check if a fraction is an integer.
 * @return true of fraction is an integer
 */
public bool is_integer()
{
	return denominator == 1;
}

/** Returns the signum function of this fraction.
 * <p>
 * The return value is -1 if the specified value is negative;
 * 0 if the specified value is zero; and 1 if the specified value is positive.
 * </p>
 * @return the signum function of this fraction
 * @since 1.7
 */
public int signum()
{
	return Integer.signum(numerator);
}

/**
 * Returns the absolute value of this fraction.
 * @return the absolute value.
 */
public Fraction abs()
{
	Fraction ret;
	if (numerator >= 0)
	{
		ret = this;
	}
else
		{
			ret = negate();
		}
		return ret;
	}

/**
 * Compares this object to another based on size.
 * @param object the object to compare to
 * @return -1 if this is less than {@code object}, +1 if this is greater
 *         than {@code object}, 0 if they are equal.
 */
 //override
 public int compare_to(Fraction object)
 {
	 long nOd = (static_cast<long>(numerator) * object.denominator;
	 long dOn = (static_cast<long>(denominator) * object.numerator;
	 return long.compare(nOd, dOn);
 }

 /**
  * Gets the fraction as a {@code double}. This calculates the fraction as
  * the numerator divided by denominator.
  * @return the fraction as a {@code double}
  */
  //override
  public double double_value()
  {
	  return static_cast<double>(numerator / static_cast<double>(denominator;
  }

  /**
   * Test for the equality of two fractions.  If the lowest term
   * numerator and denominators are the same for both fractions, the two
   * fractions are considered to be equal.
   * @param other fraction to test for equality to this fraction
   * @return true if two fractions are equal, false if object is
   *         {@code NULL}, not an instance of {@link Fraction}, or not equal
   *         to this fraction instance.
   */
   //override
   public bool equals(const Object& other)
   {
	   if (*this == other)
	   {
		   return true;
	   }
	   if (dynamic_cast<const Fraction*>(*other) != nullptr)
	   {
		   // since fractions are always in lowest terms, numerators and
		   // denominators can be compared directly for equality.
		   Fraction rhs = (Fraction)other;
		   return (numerator == rhs.numerator) &&
			   (denominator == rhs.denominator);
	   }
	   return false;
   }

   /**
	* Gets the fraction as a {@code float}. This calculates the fraction as
	* the numerator divided by denominator.
	* @return the fraction as a {@code float}
	*/
	//override
	public float float_value()
	{
		return (float)double_value();
	}

	/**
	 * Access the denominator.
	 * @return the denominator.
	 */
	public int get_denominator()
	{
		return denominator;
	}

	/**
	 * Access the numerator.
	 * @return the numerator.
	 */
	public int get_numerator()
	{
		return numerator;
	}

	/**
	 * Gets a hash_code for the fraction.
	 * @return a hash code value for this object
	 */
	 //override
	 public int hash_code()
	 {
		 return 37 * (37 * 17 + numerator) + denominator;
	 }

	 /**
	  * Gets the fraction as an {@code int}. This returns the whole number part
	  * of the fraction.
	  * @return the whole number fraction part
	  */
	  //override
	  public int int_value()
	  {
		  return static_cast<int>(double_value();
	  }

	  /**
	   * Gets the fraction as a {@code long}. This returns the whole number part
	   * of the fraction.
	   * @return the whole number fraction part
	   */
	   //override
	   public long long_value()
	   {
		   return static_cast<long>(double_value();
	   }

	   /**
		* Return the additive inverse of this fraction.
		* @return the negation of this fraction.
		*/
		//override
		public Fraction negate()
		{
			if (numerator == _integer.MIN_VALUE)
			{
				throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_FRACTION, numerator, denominator);
			}
			return Fraction(-numerator, denominator);
		}

		/**
		 * Return the multiplicative inverse of this fraction.
		 * @return the reciprocal fraction
		 */
		 //override
		 public Fraction reciprocal()
		 {
			 return Fraction(denominator, numerator);
		 }

		 /**
		  * Adds the value of this fraction to another, returning the result in reduced form.
		  * The algorithm follows Knuth, 4.5.1.
		  *
		  * @param fraction  the fraction to add, must not be {@code NULL}
		  * @return a {@code Fraction} instance with the resulting values
		  * @org.hipparchus.exception.Null_Argument_Exception if the fraction is {@code NULL}
		  * @Math_Runtime_Exception if the resulting numerator or denominator exceeds
		  *  {@code std::numeric_limits<int>::max()}
		  */
		  //override
		  public Fraction add(Fraction fraction)
		  {
			  return add_sub(fraction, true /* add */);
		  }

		  /**
		   * Add an integer to the fraction.
		   * @param i the {@code integer} to add.
		   * @return this + i
		   */
		  public Fraction add(const int& i)
		  {
			  return Fraction(numerator + i * denominator, denominator);
		  }

		  /**
		   * Subtracts the value of another fraction from the value of this one, * returning the result in reduced form.
		   *
		   * @param fraction  the fraction to subtract, must not be {@code NULL}
		   * @return a {@code Fraction} instance with the resulting values
		   * @org.hipparchus.exception.Null_Argument_Exception if the fraction is {@code NULL}
		   * @Math_Runtime_Exception if the resulting numerator or denominator
		   *   cannot be represented in an {@code int}.
		   */
		   //override
		   public Fraction subtract(Fraction fraction)
		   {
			   return add_sub(fraction, false /* subtract */);
		   }

		   /**
			* Subtract an integer from the fraction.
			* @param i the {@code integer} to subtract.
			* @return this - i
			*/
		   public Fraction subtract(const int& i)
		   {
			   return Fraction(numerator - i * denominator, denominator);
		   }

		   /**
			* Implement add and subtract using algorithm described in Knuth 4.5.1.
			*
			* @param fraction the fraction to subtract, must not be {@code NULL}
			* @param is_add true to add, false to subtract
			* @return a {@code Fraction} instance with the resulting values
			* @org.hipparchus.exception.Null_Argument_Exception if the fraction is {@code NULL}
			* @Math_Runtime_Exception if the resulting numerator or denominator
			*   cannot be represented in an {@code int}.
			*/
		   private Fraction add_sub(Fraction fraction, bool is_add)
		   {
			   //Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);

			   // zero is identity for addition.
			   if (numerator == 0)
			   {
				   return is_add ? fraction : fraction.negate();
			   }
			   if (fraction.numerator == 0)
			   {
				   return this;
			   }
			   // if denominators are randomly distributed, d1 will be 1 about 61%
			   // of the time.
			   int d1 = Arithmetic_Utils.gcd(denominator, fraction.denominator);
			   if (d1 == 1)
			   {
				   // result is ( (u*v' +/- u'v) / u'v')
				   int uvp = Arithmetic_Utils.mul_and_check(numerator, fraction.denominator);
				   int upv = Arithmetic_Utils.mul_and_check(fraction.numerator, denominator);
				   return Fraction
					   (is_add ? Arithmetic_Utils.add_and_check(uvp, upv) :
						Arithmetic_Utils.sub_and_check(uvp, upv), Arithmetic_Utils.mul_and_check(denominator, fraction.denominator));
			   }
			   // the quantity 't' requires 65 bits of precision; see knuth 4.5.1
			   // exercise 7.  we're going to use a BigInteger.
			   // t = u(v'/d1) +/- v(u'/d1)
			   BigInteger uvp = Bigstatic_cast<int>(numerator)
										  .multiply(Bigstatic_cast<int>(fraction.denominator / d1));
			   BigInteger upv = Bigstatic_cast<int>(fraction.numerator)
										  .multiply(Bigstatic_cast<int>(denominator / d1));
			   BigInteger t = is_add ? uvp.add(upv) : uvp.subtract(upv);
			   // but d2 doesn't need extra precision because
			   // d2 = gcd(t,d1) = gcd(t mod d1, d1)
			   int tmodd1 = t.mod(Bigstatic_cast<int>(d1)).int_value();
			   int d2 = (tmodd1 == 0) ? d1 : Arithmetic_Utils.gcd(tmodd1, d1);

			   // result is (t/d2) / (u'/d1)(v'/d2)
			   BigInteger w = t.divide(Bigstatic_cast<int>(d2));
			   if (w.bit_length() > 31)
			   {
				   throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::NUMERATOR_OVERFLOW_AFTER_MULTIPLY, w);
			   }
			   return Fraction(w.int_value(), Arithmetic_Utils.mul_and_check(denominator / d1, fraction.denominator / d2));
		   }

		   /**
			* Multiplies the value of this fraction by another, returning the
			* result in reduced form.
			*
			* @param fraction  the fraction to multiply by, must not be {@code NULL}
			* @return a {@code Fraction} instance with the resulting values
			* @org.hipparchus.exception.Null_Argument_Exception if the fraction is {@code NULL}
			* @Math_Runtime_Exception if the resulting numerator or denominator exceeds
			*  {@code std::numeric_limits<int>::max()}
			*/
			//override
			public Fraction multiply(Fraction fraction)
			{
				//Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);
				if (numerator == 0 || fraction.numerator == 0)
				{
					return ZERO;
				}
				// knuth 4.5.1
				// make sure we don't overflow unless the result *must* overflow.
				int d1 = Arithmetic_Utils.gcd(numerator, fraction.denominator);
				int d2 = Arithmetic_Utils.gcd(fraction.numerator, denominator);
				return get_reduced_fraction
						(Arithmetic_Utils.mul_and_check(numerator / d1, fraction.numerator / d2), Arithmetic_Utils.mul_and_check(denominator / d2, fraction.denominator / d1));
			}

			/**
			 * Multiply the fraction by an integer.
			 * @param i the {@code integer} to multiply by.
			 * @return this * i
			 */
			 //override
			 public Fraction multiply(const int& i)
			 {
				 return multiply(new Fraction(i));
			 }

			 /**
			  * Divide the value of this fraction by another.
			  *
			  * @param fraction  the fraction to divide by, must not be {@code NULL}
			  * @return a {@code Fraction} instance with the resulting values
			  * @Illegal_Argument_Exception if the fraction is {@code NULL}
			  * @Math_Runtime_Exception if the fraction to divide by is zero
			  * @Math_Runtime_Exception if the resulting numerator or denominator exceeds
			  *  {@code std::numeric_limits<int>::max()}
			  */
			  //override
			  public Fraction divide(Fraction fraction)
			  {
				  //Math_Utils::check_not_null(fraction, hipparchus::exception::Localized_Core_Formats_Type::FRACTION);
				  if (fraction.numerator == 0)
				  {
					  throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_FRACTION_TO_DIVIDE_BY, fraction.numerator, fraction.denominator);
				  }
				  return multiply(fraction.reciprocal());
			  }

			  /**
			   * Divide the fraction by an integer.
			   * @param i the {@code integer} to divide by.
			   * @return this * i
			   */
			  public Fraction divide(const int& i)
			  {
				  return divide(new Fraction(i));
			  }

			  /**
			   * Gets the fraction percentage as a {@code double}. This calculates the
			   * fraction as the numerator divided by denominator multiplied by 100.
			   *
			   * @return the fraction percentage as a {@code double}.
			   */
			  public double percentage_value()
			  {
				  return 100 * double_value();
			  }

			  /**
			   * Creates a {@code Fraction} instance with the 2 parts
			   * of a fraction Y/Z.
			   * <p>
			   * Any negative signs are resolved to be on the numerator.
			   *
			   * @param numerator  the numerator, for example the three in 'three sevenths'
			   * @param denominator  the denominator, for example the seven in 'three sevenths'
			   * @return a fraction instance, with the numerator and denominator reduced
			   * @Math_Runtime_Exception if the denominator is {@code zero}
			   */
			  public static Fraction get_reduced_fraction(const int& numerator, int denominator)
			  {
				  if (denominator == 0)
				  {
					  throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_DENOMINATOR_IN_FRACTION, numerator, denominator);
				  }
				  if (numerator == 0)
				  {
					  return ZERO; // normalize zero.
				  }
				  // allow 2^k/-2^31 as a valid fraction (where k>0)
				  if (denominator == _integer.MIN_VALUE && (numerator & 1) == 0)
				  {
					  numerator /= 2; denominator /= 2;
				  }
				  if (denominator < 0)
				  {
					  if (numerator == _integer.MIN_VALUE ||
						  denominator == _integer.MIN_VALUE)
						  {
						  throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::OVERFLOW_IN_FRACTION, numerator, denominator);
					  }
					  numerator = -numerator;
					  denominator = -denominator;
				  }
				  // simplify fraction.
				  int gcd = Arithmetic_Utils.gcd(numerator, denominator);
				  numerator /= gcd;
				  denominator /= gcd;
				  return Fraction(numerator, denominator);
			  }

			  /**
			   * Returns the {@code std::string} representing this fraction, ie
			   * "num / dem" or just "num" if the denominator is one.
			   *
			   * @return a string representation of the fraction.
			   * @see java.lang.Object#to_string()
			   */
			   //override
			   public std::string to_string() const
			   {
				   if (denominator == 1)
				   {
					   return Integer.to_string(numerator);
				   }
		   else if (numerator == 0)
				   {
					   return "0";
				   }
		   else
				   {
					   return numerator + " / " + denominator;
				   }
			   }

			   /** {@inherit_doc} */
			   //override
			   public Fraction_Field get_field()
			   {
				   return Fraction_Field.get_instance();
			   }
			   }
