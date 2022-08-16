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
  //import java.text.Field_Position;
  //import java.text.Number_Format;
  //import java.text.Parse_Position;
  //import java.util.Locale;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.util.Math_Utils;
#include <string>
#include <iostream>

/**
 * Common part shared by both {@link Fraction_Format} and {@link Big_Fraction_Format}.
 */
class Abstract_Format : public Number_Format
{
private:

	/** The format used for the denominator. */
	const Number_Format my_denominator_format;

	/** The format used for the numerator. */
	const Number_Format my_numerator_format;

protected:
	/**
	 * Create an improper formatting instance with the default number format
	 * for the numerator and denominator.
	 */
	Abstract_Format()
	{
		this(get_default_number_format());
	}

	/**
	 * Create an improper formatting instance with a custom number format for
	 * both the numerator and denominator.
	 * @param format the custom format for both the numerator and denominator.
	 * @org.hipparchus.exception. if the provided format is NULL.
	 */
	Abstract_Format(const Number_Format& format)
	{
		Abstract_Format(format, (Number_Format)format.clone());
	}

	/**
	 * Create an improper formatting instance with a custom number format for
	 * the numerator and a custom number format for the denominator.
	 * @param numerator_format the custom format for the numerator.
	 * @param denominator_format the custom format for the denominator.
	 * @org.hipparchus.exception. if either provided format is NULL.
	 */
	Abstract_Format(const Number_Format& numerator_format, const Number_Format& denominator_format)
		: my_numerator_format{ numerator_format }, my_denominator_format{ denominator_format }
	{
		//Math_Utils::check_not_null(numerator_format, hipparchus::exception::Localized_Core_Formats_Type::NUMERATOR_FORMAT);
		//Math_Utils::check_not_null(denominator_format, hipparchus::exception::Localized_Core_Formats_Type::DENOMINATOR_FORMAT);
	}

	/**
	 * Create a default number format.  The default number format is based on
	 * {@link Number_Format#get_number_instance(java.util.Locale)}. The only
	 * customization is the maximum number of Big_Fraction digits, which is set to 0.
	 * @return the default number format.
	 */
	static Number_Format get_default_number_format()
	{
		return get_default_number_format(Locale.get_default());
	}

	/**
	 * Create a default number format.  The default number format is based on
	 * {@link Number_Format#get_number_instance(java.util.Locale)}. The only
	 * customization is the maximum number of Big_Fraction digits, which is set to 0.
	 * @param locale the specific locale used by the format.
	 * @return the default number format specific to the given locale.
	 */
	static Number_Format get_default_number_format(const Locale& locale)
	{
		const auto nf = Number_Format.get_number_instance(locale);
		nf.set_maximum_fraction_digits(0);
		nf.set_parse_integer_only(true);
		return nf;
	}

	/**
	 * Parses <code>source</code> until a non-whitespace character is found.
	 * @param source the string to parse
	 * @param pos input/output parsing parameter.  On output, <code>pos</code>
	 *        holds the index of the next non-whitespace character.
	 */
	static void parse_and_ignore_whitespace(const std::string& source, const Parse_Position& pos)
	{
		parse_next_character(source, pos);
		pos.set_index(pos.get_index() - 1);
	}

	/**
	 * Parses <code>source</code> until a non-whitespace character is found.
	 * @param source the string to parse
	 * @param pos input/output parsing parameter.
	 * @return the first non-whitespace character.
	 */
	static char parse_next_character(const std::string& source, const Parse_Position& pos)
	{
		int index = pos.get_index();
		const auto n = source.size();
		char ret{};

		if (index < n)
		{
			char c;
			do
			{
				c = source.at(index++);
			} while (isspace(c) && index < n);
			pos.set_index(index);

			if (index < n)
			{
				ret = c;
			}
		}

		return ret;
	}

public:
	/**
	 * Access the denominator format.
	 * @return the denominator format.
	 */
	Number_Format get_denominator_format() const
	{
		return my_denominator_format;
	}

	/**
	 * Access the numerator format.
	 * @return the numerator format.
	 */
	Number_Format get_numerator_format() const
	{
		return my_numerator_format;
	}

	/**
	 * Formats a double value as a fraction and appends the result to a String_Buffer.
	 *
	 * @param value the double value to format
	 * @param buffer std::stringstreamto append to
	 * @param position On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @return a reference to the appended buffer
	 * @see #format(Object, String_Buffer, Field_Position)
	 */
	 //override
	std::stringstream format(const double& value, const std::stringstream& buffer, const Field_Position& position)
	{
		return format(value, buffer, position);
	}

	/**
	 * Formats a long value as a fraction and appends the result to a String_Buffer.
	 *
	 * @param value the long value to format
	 * @param buffer std::stringstreamto append to
	 * @param position On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @return a reference to the appended buffer
	 * @see #format(Object, String_Buffer, Field_Position)
	 */
	 //override
	std::stringstream format(const long& value, const std::stringstream& buffer, const Field_Position& position)
	{
		return format(value, buffer, position);
	}
}
