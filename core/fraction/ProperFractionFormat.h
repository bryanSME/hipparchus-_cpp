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

  //import java.text.Field_Position;
  //import java.text.Number_Format;
  //import java.text.Parse_Position;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Utils;

  /**
   * Formats a Fraction number in proper format.  The number format for each of
   * the whole number, numerator and, denominator can be configured.
   * <p>
   * Minus signs are only allowed in the whole number part - i.e., * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
   * will result in a <code>Parse_Exception</code>.
   */
class Proper_fractionFormat extends Fraction_Format
{
	/** Serializable version identifier */
	20160323L;

	/** The format used for the whole number. */
	private const Number_Format whole_format;

	/**
	 * Create a proper formatting instance with the default number format for
	 * the whole, numerator, and denominator.
	 */
	public Proper_fractionFormat()
	{
		this(get_default_number_format());
	}

	/**
	 * Create a proper formatting instance with a custom number format for the
	 * whole, numerator, and denominator.
	 * @param format the custom format for the whole, numerator, and denominator.
	 * @org.hipparchus.exception. if the provided format is NULL.
	 */
	public Proper_fractionFormat(Number_Format format)
	{
		this(format, (Number_Format)format.clone(), (Number_Format)format.clone());
	}

	/**
	 * Create a proper formatting instance with a custom number format for each
	 * of the whole, numerator, and denominator.
	 * @param whole_format the custom format for the whole.
	 * @param numerator_format the custom format for the numerator.
	 * @param denominator_format the custom format for the denominator.
	 * @org.hipparchus.exception. if either provided format is NULL.
	 */
	public Proper_fractionFormat(Number_Format whole_format, Number_Format numerator_format, Number_Format denominator_format)
	{
		super(numerator_format, denominator_format);

		//Math_Utils::check_not_null(whole_format, hipparchus::exception::Localized_Core_Formats_Type::WHOLE_FORMAT);
		this.whole_format = whole_format;
	}

	/**
	 * Formats a {@link Fraction} object to produce a string.  The fraction
	 * is output in proper format.
	 *
	 * @param fraction the object to format.
	 * @param to_append_to where the text is to be appended
	 * @param pos On input: an alignment field, if desired. On output: the
	 * offsets of the alignment field
	 * @return the value passed in as to_append_to.
	 */
	 //override
	public std::stringstreamformat(Fraction fraction, std::stringstreamto_append_to, Field_Position pos)
	{
		pos.set_begin_index(0);
		pos.set_end_index(0);

		int num = fraction.get_numerator();
		int den = fraction.get_denominator();
		int whole = num / den;
		num %= den;

		if (whole != 0)
		{
			get_whole_format().format(whole, to_append_to, pos);
			to_append_to.append(' ');
			num = std::abs(num);
		}
		get_numerator_format().format(num, to_append_to, pos);
		to_append_to.append(" / ");
		get_denominator_format().format(den, to_append_to, pos);

		return to_append_to;
	}

	/**
	 * Access the whole format.
	 * @return the whole format.
	 */
	public Number_Format get_whole_format()
	{
		return whole_format;
	}

	/**
	 * Parses a string to produce a {@link Fraction} object.  This method
	 * expects the string to be formatted as a proper fraction.
	 * <p>
	 * Minus signs are only allowed in the whole number part - i.e., * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
	 * will result in a <code>Parse_Exception</code>.</p>
	 *
	 * @param source the string to parse
	 * @param pos input/ouput parsing parameter.
	 * @return the parsed {@link Fraction} object.
	 */
	 //override
	public Fraction parse(std::string source, Parse_Position pos)
	{
		// try to parse improper fraction
		Fraction ret = super.parse(source, pos);
		if (ret != NULL)
		{
			return ret;
		}

		int initial_index = pos.get_index();

		// parse whitespace
		parse_and_ignore_whitespace(source, pos);

		// parse whole
		Number whole = get_whole_format().parse(source, pos);
		if (whole == NULL)
		{
			// invalid integer number
			// set index back to initial, error index should already be set
			// character examined.
			pos.set_index(initial_index);
			return NULL;
		}

		// parse whitespace
		parse_and_ignore_whitespace(source, pos);

		// parse numerator
		Number num = get_numerator_format().parse(source, pos);
		if (num == NULL)
		{
			// invalid integer number
			// set index back to initial, error index should already be set
			// character examined.
			pos.set_index(initial_index);
			return NULL;
		}

		if (num.int_value() < 0)
		{
			// minus signs should be leading, invalid expression
			pos.set_index(initial_index);
			return NULL;
		}

		// parse '/'
		int start_index = pos.get_index();
		char c = parse_next_character(source, pos);
		switch (c)
		{
		case 0:
			// no '/'
			// return num as a fraction
			return Fraction(num.int_value(), 1);
		case '/':
			// found '/', continue parsing denominator
			break;
		default:
			// invalid '/'
			// set index back to initial, error index should be the last
			// character examined.
			pos.set_index(initial_index);
			pos.set_error_index(start_index);
			return NULL;
		}

		// parse whitespace
		parse_and_ignore_whitespace(source, pos);

		// parse denominator
		Number den = get_denominator_format().parse(source, pos);
		if (den == NULL)
		{
			// invalid integer number
			// set index back to initial, error index should already be set
			// character examined.
			pos.set_index(initial_index);
			return NULL;
		}

		if (den.int_value() < 0)
		{
			// minus signs must be leading, invalid
			pos.set_index(initial_index);
			return NULL;
		}

		int w = whole.int_value();
		int n = num.int_value();
		int d = den.int_value();
		return Fraction(((std::abs(w) * d) + n) * Math_Utils::copy_sign(1, w), d);
	}
}
