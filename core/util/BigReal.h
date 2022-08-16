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

  //import java.io.Serializable;
  //import java.math.BigDecimal;
  //import java.math.BigInteger;
  //import java.math.Math_Context;
  //import java.math.Rounding_Mode;

  //import org.hipparchus.Field;
  //import org.hipparchus.Field_Element;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Runtime_Exception;
#include "../FieldElement.h"
#include "Comparable.h"
#include "../Field.h"

/**
 * Arbitrary precision decimal number.
 * <p>
 * This class is a simple wrapper around the standard <code>BigDecimal</code>
 * in order to implement the {@link Field_Element} interface.
 */
class BigReal : Field_Element<BigReal>, Comparable<BigReal>
{
private:

	/** Underlying BigDecimal. */
	const BigDecimal my_d;

	/** Rounding mode for divisions. **/
	Rounding_Mode my_rounding_mode = Rounding_Mode.HALF_UP;

	/*** BigDecimal scale ***/
	int my_scale{ 64 };

public:
	/** A big real representing 0. */
	static const BigReal ZERO = BigReal(BigDecimal.ZERO);

	/** A big real representing 1. */
	static const BigReal ONE = BigReal(BigDecimal.ONE);

	/** Build an instance from a BigDecimal.
	 * @param val value of the instance
	 */
	BigReal(const BigDecimal& val) : my_d{ val } {};

	/** Build an instance from a BigInteger.
	 * @param val value of the instance
	 */
	BigReal(const BigInteger& val) : my_d{ BigDecimal(val) } {};

	/** Build an instance from an unscaled BigInteger.
	 * @param unscaled_val unscaled value
	 * @param scale scale to use
	 */
	BigReal(const BigInteger& unscaled_val, const int& scale) : my_d{ BigDecimal(unscaled_val, scale) } {};

	/** Build an instance from an unscaled BigInteger.
	 * @param unscaled_val unscaled value
	 * @param scale scale to use
	 * @param mc to used
	 */
	BigReal(const BigInteger& unscaled_val, const int& scale, const Math_Context& mc) : my_d{ BigDecimal(unscaled_val, scale, mc) } {};

	/** Build an instance from a BigInteger.
	 * @param val value of the instance
	 * @param mc context to use
	 */
	BigReal(const BigInteger& val, const Math_Context& mc) : my_d{ BigDecimal(val, mc) } {};

	/** Build an instance from a characters representation.
	 * @param in character representation of the value
	 */
	BigReal(const char* in) : my_d{ BigDecimal(in) } {};

	/** Build an instance from a characters representation.
	 * @param in character representation of the value
	 * @param offset offset of the first character to analyze
	 * @param len length of the array slice to analyze
	 */
	BigReal(const char* in, const int& offset, const int& len) : my_d{ BigDecimal(in, offset, len) } {};

	/** Build an instance from a characters representation.
	 * @param in character representation of the value
	 * @param offset offset of the first character to analyze
	 * @param len length of the array slice to analyze
	 * @param mc context to use
	 */
	BigReal(const char* in, const int& offset, const int& len, const Math_Context& mc) : my_d{ BigDecimal(in, offset, len, mc) } {};

	/** Build an instance from a characters representation.
	 * @param in character representation of the value
	 * @param mc context to use
	 */
	BigReal(const char* in, const Math_Context& mc) : my_d{ BigDecimal(in, mc) } {};

	/** Build an instance from a double.
	 * @param val value of the instance
	 */
	BigReal(const double& val) : my_d{ BigDecimal(val) } {};
	// NOPMD - we really want double conversion here

	/** Build an instance from a double.
	 * @param val value of the instance
	 * @param mc context to use
	 */
	BigReal(const double& val, const Math_Context& mc) : my_d{ BigDecimal(val, mc) } {};
	// NOPMD - we really want double conversion here

	/** Build an instance from an int.
	 * @param val value of the instance
	 */
	BigReal(const int& val) : my_d{ BigDecimal(val) } {};

	/** Build an instance from an int.
	 * @param val value of the instance
	 * @param mc context to use
	 */
	BigReal(const int& val, const Math_Context& mc) : my_d{ BigDecimal(val, mc) } {};

	/** Build an instance from a long.
	 * @param val value of the instance
	 */
	BigReal(long val) : my_d{ BigDecimal(val) } {};

	/** Build an instance from a long.
	 * @param val value of the instance
	 * @param mc context to use
	 */
	BigReal(const long& val, const Math_Context& mc) : my_d{ BigDecimal(val, mc) } {};

	/** Build an instance from a std::string representation.
	 * @param val character representation of the value
	 */
	BigReal(std::string val) : my_d{ BigDecimal(val) } {};

	/** Build an instance from a std::string representation.
	 * @param val character representation of the value
	 * @param mc context to use
	 */
	BigReal(std::string val, const Math_Context& mc) : my_d{ BigDecimal(val, mc) } {};

	/** {@inherit_doc} */
	//override
	double get_real()
	{
		return double_value();
	}

	/***
	 * Gets the rounding mode for division operations
	 * The default is {@code Rounding_Mode.HALF_UP}
	 * @return the rounding mode.
	 */
	Rounding_Mode get_rounding_mode() const
	{
		return my_rounding_mode;
	}

	/***
	 * Sets the rounding mode for decimal divisions.
	 * @param rounding_mode rounding mode for decimal divisions
	 */
	void set_rounding_mode(const Rounding_Mode& rounding_mode)
	{
		my_rounding_mode = rounding_mode;
	}

	/***
	 * Sets the scale for division operations.
	 * The default is 64
	 * @return the scale
	 */
	int get_scale()
	{
		return my_scale;
	}

	/***
	 * Sets the scale for division operations.
	 * @param scale scale for division operations
	 */
	void set_scale(const int& scale)
	{
		this.scale = scale;
	}

	/** {@inherit_doc} */
	//override
	BigReal add(const BigReal& a)
	{
		return BigReal(d.add(a.d));
	}

	/** {@inherit_doc} */
	//override
	BigReal subtract(const BigReal& a)
	{
		return BigReal(d.subtract(a.d));
	}

	/** {@inherit_doc} */
	//override
	BigReal negate()
	{
		return BigReal(d.negate());
	}

	/**
	 * {@inherit_doc}
	 *
	 * @Math_Runtime_Exception if {@code a} is zero
	 */
	 //override
	BigReal divide(const BigReal& a)
	{
		try
		{
			return BigReal(d.divide(a.d, scale, rounding_mode));
		}
		catch (Arithmetic_Exception e)
		{
			throw std::exception("not implemented");
			// Division by zero has occurred
			//throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
		}
	}

	/**
	 * {@inherit_doc}
	 *
	 * @Math_Runtime_Exception if {@code this} is zero
	 */
	 //override
	BigReal reciprocal()
	{
		throw std::exception("not implemented");
		//try
		//{
		//    return BigReal(BigDecimal.ONE.divide(d, scale, rounding_mode));
		//}
		//catch (Arithmetic_Exception e)
		//{
		//    // Division by zero has occurred
		//    //throw Math_Runtime_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::ZERO_NOT_ALLOWED);
		//}
	}

	/** {@inherit_doc} */
	//override
	BigReal multiply(const BigReal& a)
	{
		return BigReal(d.multiply(a.d));
	}

	/** {@inherit_doc} */
	//override
	BigReal multiply(const int& n)
	{
		return BigReal(d.multiply(new BigDecimal(n)));
	}

	/** {@inherit_doc} */
	//override
	int compare_to(const BigReal& a)
	{
		return d.compare_to(a.d);
	}

	/** Get the double value corresponding to the instance.
	 * @return double value corresponding to the instance
	 */
	double double_value()
	{
		return d.double_value();
	}

	/** Get the BigDecimal value corresponding to the instance.
	 * @return BigDecimal value corresponding to the instance
	 */
	BigDecimal big_decimal_value()
	{
		return my_d;
	}

	/** {@inherit_doc} */
	//override
	bool equals(Object other)
	{
		if (this == other)
		{
			return true;
		}

		if (dynamic_cast<const BigReal*>(*other) != nullptr)
		{
			return d.equals(((BigReal)other).d);
		}
		return false;
	}

	/** {@inherit_doc} */
	//override
	int hash_code()
	{
		return d.hash_code();
	}

	/** {@inherit_doc} */
	//override
	Field<BigReal> get_field()
	{
		return Big_Real_Field.get_instance();
	}
};