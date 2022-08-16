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

  //package org.hipparchus.geometry;

  //import java.text.Field_Position;
  //import java.text.Number_Format;
  //import java.text.Parse_Position;
  //import java.util.Locale;

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.Composite_Format;

#include <type_traits>
#include <vector>
#include <string>

#include "Space.h"

  /**
   * Formats a vector in components list format "{x; y; ...}".
   * <p>The prefix and suffix "{" and "}" and the separator "; " can be replaced by
   * any user-defined strings. The number format for components can be configured.</p>
   * <p>White space is ignored at parse time, even if it is in the prefix, suffix
   * or separator specifications. So even if the default separator does include a space
   * character that is used at format time, both input string "{1;1;1}" and
   * " { 1 ; 1 ; 1 } " will be parsed without error and the same vector will be
   * returned. In the second case, however, the parse position after parsing will be
   * just after the closing curly brace, i.e. just before the trailing space.</p>
   * <p><b>Note:</b> using "," as a separator may interfere with the grouping separator
   * of the default {@link Number_Format} for the current locale. Thus it is advised
   * to use a {@link Number_Format} instance with disabled grouping in such a case.</p>
   *
   * @param <S> Type of the space.
   */
template<
	typename S,
	typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Vector_Format
{
public:
	/** The default prefix: "{". */
	static const std::string DEFAULT_PREFIX = "{";

	/** The default suffix: "}". */
	static const std::string DEFAULT_SUFFIX = "}";

	/** The default separator: ", ". */
	static const std::string DEFAULT_SEPARATOR = "; ";

	/**
 * Get the set of locales for which point/vector formats are available.
 * <p>This is the same set as the {@link Number_Format} set.</p>
 * @return available point/vector format locales.
 */
	static Locale[] get_available_locales()
	{
		return Number_Format.get_available_locales();
	}

	/**
	 * Get the format prefix.
	 * @return format prefix.
	 */
	std::string get_prefix()
	{
		return prefix;
	}

	/**
	 * Get the format suffix.
	 * @return format suffix.
	 */
	std::string get_suffix()
	{
		return suffix;
	}

	/**
	 * Get the format separator between components.
	 * @return format separator.
	 */
	std::string get_separator()
	{
		return separator;
	}

	/**
	 * Get the components format.
	 * @return components format.
	 */
	Number_Format get_format()
	{
		return format;
	}

	/**
	 * Formats a {@link Vector} object to produce a string.
	 * @param vector the object to format.
	 * @return a formatted string.
	 */
	std::string format(Vector<S> vector)
	{
		return format(vector, String_Buffer(), Field_Position(0)).to_string();
	}

	/**
	 * Formats a {@link Vector} object to produce a string.
	 * @param vector the object to format.
	 * @param to_append_to where the text is to be appended
	 * @param pos On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @return the value passed in as to_append_to.
	 */
	virtual std::stringstreamformat(Vector<S> vector, std::stringstreamto_append_to, Field_Position pos);

	/**
	 * Parses a string to produce a {@link Vector} object.
	 * @param source the string to parse
	 * @return the parsed {@link Vector} object.
	 * @Math_Illegal_State_Exception if the beginning of the specified string
	 * cannot be parsed.
	 */
	virtual Vector<S> parse(std::string source) Math_Illegal_State_Exception;

	/**
	 * Parses a string to produce a {@link Vector} object.
	 * @param source the string to parse
	 * @param pos input/output parsing parameter.
	 * @return the parsed {@link Vector} object.
	 */
	virtual Vector<S> parse(std::string source, Parse_Position pos);

private:
	/** Prefix. */
	const std::string my_prefix;

	/** Suffix. */
	const std::string suffix;

	/** Separator. */
	const std::string my_separator;

	/** Trimmed prefix. */
	const std::string my_trimmed_prefix;

	/** Trimmed suffix. */
	const std::string my_trimmed_suffix;

	/** Trimmed separator. */
	const std::string my_trimmed_separator;

	/** The format used for components. */
	const Number_Format my_format;

protected:
	/**
	 * Create an instance with default settings.
	 * <p>The instance uses the default prefix, suffix and separator:
	 * "{", "}", and "; " and the default number format for components.</p>
	 */
	Vector_Format()
	{
		this(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, Composite_Format.get_default_number_format());
	}

	/**
	 * Create an instance with a custom number format for components.
	 * @param format the custom format for components.
	 */
	Vector_Format(const Number_Format format)
	{
		this(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, format);
	}

	/**
	 * Create an instance with custom prefix, suffix and separator.
	 * @param prefix prefix to use instead of the default "{"
	 * @param suffix suffix to use instead of the default "}"
	 * @param separator separator to use instead of the default "; "
	 */
	Vector_Format(const std::string prefix, const std::string suffix, const std::string separator)
	{
		this(prefix, suffix, separator, Composite_Format.get_default_number_format());
	}

	/**
	 * Create an instance with custom prefix, suffix, separator and format
	 * for components.
	 * @param prefix prefix to use instead of the default "{"
	 * @param suffix suffix to use instead of the default "}"
	 * @param separator separator to use instead of the default "; "
	 * @param format the custom format for components.
	 */
	Vector_Format(const std::string prefix, const std::string suffix, const std::string separator, const Number_Format format)
	{
		this.prefix = prefix;
		this.suffix = suffix;
		this.separator = separator;
		trimmed_prefix = prefix.trim();
		trimmed_suffix = suffix.trim();
		trimmed_separator = separator.trim();
		this.format = format;
	}

	/**
	 * Parses a string to produce an array of coordinates.
	 * @param dimension dimension of the space
	 * @param source the string to parse
	 * @param pos input/output parsing parameter.
	 * @return coordinates array.
	 */
	std::vector<double> parse_coordinates(const int& dimension, std::string source, Parse_Position pos)
	{
		int initial_index = pos.get_index();
		std::vector<double> coordinates = std::vector<double>(dimension];

		// parse prefix
		Composite_Format.parse_and_ignore_whitespace(source, pos);
		if (!Composite_Format.parse_fixed_string(source, trimmed_prefix, pos))
		{
			return NULL;
		}

		for (int i{}; i < dimension; ++i)
		{
			// skip whitespace
			Composite_Format.parse_and_ignore_whitespace(source, pos);

			// parse separator
			if (i > 0 && !Composite_Format.parse_fixed_string(source, trimmed_separator, pos))
			{
				return NULL;
			}

			// skip whitespace
			Composite_Format.parse_and_ignore_whitespace(source, pos);

			// parse coordinate
			Number c = Composite_Format.parse_number(source, format, pos);
			if (c == NULL)
			{
				// invalid coordinate
				// set index back to initial, error index should already be set
				pos.set_index(initial_index);
				return NULL;
			}

			// store coordinate
			coordinates[i] = c.double_value();
		}

		// parse suffix
		Composite_Format.parse_and_ignore_whitespace(source, pos);
		if (!Composite_Format.parse_fixed_string(source, trimmed_suffix, pos))
		{
			return NULL;
		}

		return coordinates;
	}

	/**
	 * Formats the coordinates of a {@link Vector} to produce a string.
	 * @param to_append_to where the text is to be appended
	 * @param pos On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @param coordinates coordinates of the object to format.
	 * @return the value passed in as to_append_to.
	 */
	std::stringstreamformat(std::stringstreamto_append_to, Field_Position pos, double ... coordinates)
	{
		pos.set_begin_index(0);
		pos.set_end_index(0);

		// format prefix
		to_append_to.append(prefix);

		// format components
		for (int i{}; i < coordinates.size(); ++i)
		{
			if (i > 0)
			{
				to_append_to.append(separator);
			}
			Composite_Format.format_double(coordinates[i], format, to_append_to, pos);
		}

		// format suffix
		to_append_to.append(suffix);

		return to_append_to;
	}
};