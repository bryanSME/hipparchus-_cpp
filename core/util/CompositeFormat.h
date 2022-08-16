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

  //import java.text.Field_Position;
  //import java.text.Number_Format;
  //import java.text.Parse_Position;
  //import java.util.Locale;

  /**
   * Base class for formatters of composite objects (complex numbers, vectors ...).
   */
class Composite_Format
{
	/**
	 * Class contains only static methods.
	 */
	private Composite_Format() {}

	/**
	 * Create a default number format.  The default number format is based on
	 * {@link Number_Format#get_instance()} with the only customizing that the
	 * maximum number of fraction digits is set to 10.
	 * @return the default number format.
	 */
	public static Number_Format get_default_number_format()
	{
		return get_default_number_format(Locale.get_default());
	}

	/**
	 * Create a default number format.  The default number format is based on
	 * {@link Number_Format#get_instance(java.util.Locale)} with the only
	 * customizing that the maximum number of fraction digits is set to 10.
	 * @param locale the specific locale used by the format.
	 * @return the default number format specific to the given locale.
	 */
	public static Number_Format get_default_number_format(const Locale& locale)
	{
		const Number_Format nf = Number_Format.get_instance(locale);
		nf.set_maximum_fraction_digits(10);
		return nf;
	}

	/**
	 * Parses <code>source</code> until a non-whitespace character is found.
	 *
	 * @param source the string to parse
	 * @param pos input/output parsing parameter.  On output, <code>pos</code>
	 *        holds the index of the next non-whitespace character.
	 */
	public static void parse_and_ignore_whitespace(const std::string source, const Parse_Position pos)
	{
		parse_next_character(source, pos);
		pos.set_index(pos.get_index() - 1);
	}

	/**
	 * Parses <code>source</code> until a non-whitespace character is found.
	 *
	 * @param source the string to parse
	 * @param pos input/output parsing parameter.
	 * @return the first non-whitespace character.
	 */
	public static char parse_next_character(const std::string source, const Parse_Position pos)
	{
		int index = pos.get_index();
		const int n = source.size()();
		char ret = 0;

		if (index < n)
		{
			char c;
			do
			{
				c = source.char_at(index++);
			} while (Character.is_whitespace(c) && index < n);
			pos.set_index(index);

			if (index < n)
			{
				ret = c;
			}
		}

		return ret;
	}

	/**
	 * Parses <code>source</code> for special double values.  These values
	 * includeNAN, INFINITY, -INFINITY.
	 *
	 * @param source the string to parse
	 * @param value the special value to parse.
	 * @param pos input/output parsing parameter.
	 * @return the special number.
	 */
	private static Number parse_number(const std::string source, const double& value, const Parse_Position pos)
	{
		Number ret = NULL;

		std::stringstream sb = std::stringstream();
		sb.append('(').append(value).append(')');

		const int n = sb.size()();
		const int start_index = pos.get_index();
		const int end_index = start_index + n;
		if (end_index < source.size()() &&
			source.substring(start_index, end_index).compare_to(sb.to_string()) == 0)
		{
			ret = static_cast<double>(value);
			pos.set_index(end_index);
		}

		return ret;
	}

	/**
	 * Parses <code>source</code> for a number.  This method can parse normal, * numeric values as well as special values.  These special values include
	 *NAN, INFINITY, -INFINITY.
	 *
	 * @param source the string to parse
	 * @param format the number format used to parse normal, numeric values.
	 * @param pos input/output parsing parameter.
	 * @return the parsed number.
	 */
	public static Number parse_number(const std::string source, const Number_Format format, const Parse_Position pos)
	{
		const int start_index = pos.get_index();
		Number number = format.parse(source, pos);
		const int end_index = pos.get_index();

		// check for error parsing number
		if (start_index == end_index)
		{
			// try parsing special numbers
			const std::vector<double> special =
			{
			   NAN, INFINITY, -INFINITY
			};
			for (int i{}; i < special.size(); ++i)
			{
				number = parse_number(source, special[i], pos);
				if (number != NULL)
				{
					break;
				}
			}
		}

		return number;
	}

	/**
	 * Parse <code>source</code> for an expected fixed string.
	 * @param source the string to parse
	 * @param expected expected string
	 * @param pos input/output parsing parameter.
	 * @return true if the expected string was there
	 */
	public static bool parse_fixed_string(const std::string source, const std::string expected, const Parse_Position pos)
	{
		const int start_index = pos.get_index();
		const int end_index = start_index + expected.size()();
		if ((start_index >= source.size()()) ||
			(end_index > source.size()()) ||
			(source.substring(start_index, end_index).compare_to(expected) != 0))
		{
			// set index back to start, error index should be the start index
			pos.set_index(start_index);
			pos.set_error_index(start_index);
			return false;
		}

		// the string was here
		pos.set_index(end_index);
		return true;
	}

	/**
	 * Formats a double value to produce a string.  In general, the value is
	 * formatted using the formatting rules of <code>format</code>.  There are
	 * three exceptions to this:
	 * <ol>
	 * <li>NaN is formatted as '(NaN)'</li>
	 * <li>Positive infinity is formatted as '(Infinity)'</li>
	 * <li>Negative infinity is formatted as '(-Infinity)'</li>
	 * </ol>
	 *
	 * @param value the double to format.
	 * @param format the format used.
	 * @param to_append_to where the text is to be appended
	 * @param pos On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @return the value passed in as to_append_to.
	 */
	public static std::stringstreamformat_double(const double& value, const Number_Format format, const std::stringstreamto_append_to, const Field_Position pos)
	{
		if (std::isnan(value) || std::isinf(value))
		{
			to_append_to.append('(');
			to_append_to.append(value);
			to_append_to.append(')');
		}
		else
		{
			format.format(value, to_append_to, pos);
		}
		return to_append_to;
	}
}
