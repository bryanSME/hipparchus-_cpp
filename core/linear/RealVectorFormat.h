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

  //package org.hipparchus.linear;

  //import java.text.Field_Position;
  //import java.text.Number_Format;
  //import java.text.Parse_Position;
  //import java.util.Array_list;
  //import java.util.List;
  //import java.util.Locale;
#include <string>
#include <vector>

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.Composite_Format;

  /**
   * Formats a vector in components list format "{v0; v1; ...; vk-1}".
   * <p>The prefix and suffix "{" and "}" and the separator "; " can be replaced by
   * any user-defined strings. The number format for components can be configured.</p>
   * <p>White space is ignored at parse time, even if it is in the prefix, suffix
   * or separator specifications. So even if the default separator does include a space
   * character that is used at format time, both input string "{1;1;1}" and
   * " { 1 ; 1 ; 1 } " will be parsed without error and the same vector will be
   * returned. In the second case, however, the parse position after parsing will be
   * just after the closing curly brace, i.e. just before the trailing space.</p>
   *
   */
class Real_Vector_Format
{
private:
	/** The default prefix: "{". */
	static constexpr char DEFAULT_PREFIX = '{';
	/** The default suffix: "}". */
	static constexpr char DEFAULT_SUFFIX = "}";
	/** The default separator: ", ". */
	static constexpr char DEFAULT_SEPARATOR = "; ";
	/** Prefix. */
	const std::string my_prefix;
	/** Suffix. */
	const std::string my_suffix;
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

public:
	/**
	 * Create an instance with default settings.
	 * <p>The instance uses the default prefix, suffix and separator:
	 * "{", "}", and "; " and the default number format for components.</p>
	 */
	Real_Vector_Format()
	{
		Real_Vector_Format(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, Composite_Format::get_default_number_format());
	}

	/**
	 * Create an instance with a custom number format for components.
	 * @param format the custom format for components.
	 */
	Real_Vector_Format(const Number_Format& format)
	{
		Real_Vector_Format(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, format);
	}

	/**
	 * Create an instance with custom prefix, suffix and separator.
	 * @param prefix prefix to use instead of the default "{"
	 * @param suffix suffix to use instead of the default "}"
	 * @param separator separator to use instead of the default "; "
	 */
	Real_Vector_Format(const std::string& prefix, const std::string& suffix, const std::string& separator)
	{
		Real_Vector_Format(prefix, suffix, separator, Composite_Format::get_default_number_format());
	}

	/**
	 * Create an instance with custom prefix, suffix, separator and format
	 * for components.
	 * @param prefix prefix to use instead of the default "{"
	 * @param suffix suffix to use instead of the default "}"
	 * @param separator separator to use instead of the default "; "
	 * @param format the custom format for components.
	 */
	Real_Vector_Format(const std::string& prefix, const std::string& suffix, const std::string& separator, const Number_Format& format)
		:
		my_prefix{ prefix },
		my_suffix{ suffix },
		my_separator{ separator },
		my_trimmed_prefix{ prefix.trim() },
		my_trimmed_suffix{ suffix.trim() },
		my_trimmed_separator{ separator.trim() },
		my_format{ format }
	{
	}

	/**
	 * Get the set of locales for which real vectors formats are available.
	 * <p>This is the same set as the {@link Number_Format} set.</p>
	 * @return available real vector format locales.
	 */
	static std::vector<Locale> get_available_locales()
	{
		return Number_Format::get_available_locales();
	}

	/**
	 * Get the format prefix.
	 * @return format prefix.
	 */
	std::string get_prefix() const
	{
		return my_prefix;
	}

	/**
	 * Get the format suffix.
	 * @return format suffix.
	 */
	std::string get_suffix() const
	{
		return my_suffix;
	}

	/**
	 * Get the format separator between components.
	 * @return format separator.
	 */
	std::string get_separator() const
	{
		return my_separator;
	}

	/**
	 * Get the components format.
	 * @return components format.
	 */
	Number_Format get_format() const
	{
		return my_format;
	}

	/**
	 * Returns the default real vector format for the current locale.
	 * @return the default real vector format.
	 * @since 1.4
	 */
	static Real_Vector_Format get_real_vector_format()
	{
		return get_real_vector_format(Locale.get_default());
	}

	/**
	 * Returns the default real vector format for the given locale.
	 * @param locale the specific locale used by the format.
	 * @return the real vector format specific to the given locale.
	 * @since 1.4
	 */
	static Real_Vector_Format get_real_vector_format(const Locale& locale)
	{
		return Real_Vector_Format(Composite_Format::get_default_number_format(locale));
	}

	/**
	 * This method calls {@link #format(Real_Vector,String_Buffer,Field_Position)}.
	 *
	 * @param v Real_Vector object to format.
	 * @return a formatted vector.
	 */
	std::string format(const Real_Vector& v)
	{
		return format(v, String_Buffer(), Field_Position(0)).to_string();
	}

	/**
	 * Formats a {@link Real_Vector} object to produce a string.
	 * @param vector the object to format.
	 * @param to_append_to where the text is to be appended
	 * @param pos On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @return the value passed in as to_append_to.
	 */
	std::stringstream format(const Real_Vector& vector, std::stringstream& to_append_to, const Field_Position& pos)
	{
		pos.set_begin_index(0);
		pos.set_end_index(0);

		// format prefix
		to_append_to.append(prefix);

		// format components
		for (int i{}; i < vector.get_dimension(); ++i)
		{
			if (i > 0)
			{
				to_append_to.append(separator);
			}
			Composite_Format.format_double(vector.get_entry(i), format, to_append_to, pos);
		}

		// format suffix
		to_append_to.append(suffix);

		return to_append_to;
	}

	/**
	 * Parse a string to produce a {@link Real_Vector} object.
	 *
	 * @param source std::string to parse.
	 * @return the parsed {@link Real_Vector} object.
	 * @Math_Illegal_State_Exception if the beginning of the specified string
	 * cannot be parsed.
	 */
	Array_Real_Vector parse(const std::string& source)
	{
		const auto parse_position = Parse_Position(0);
		const Array_Real_Vector result = parse(source, parse_position);
		if (parse_position.get_index() == 0)
		{
			throw std::exception("not implemented");
			//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE_AS_TYPE, source, parse_position.get_error_index(), Array_Real_Vector.class);
		}
		return result;
	}

	/**
	 * Parse a string to produce a {@link Real_Vector} object.
	 *
	 * @param source std::string to parse.
	 * @param pos input/ouput parsing parameter.
	 * @return the parsed {@link Real_Vector} object.
	 */
	Array_Real_Vector parse(const std::string& source, const Parse_Position& pos)
	{
		int initial_index = pos.get_index();

		// parse prefix
		Composite_Format.parse_and_ignore_whitespace(source, pos);
		if (!Composite_Format.parse_fixed_string(source, trimmed_prefix, pos))
		{
			return NULL;
		}

		// parse components
		List<Number> components = Array_list<>();
		for (bool loop{ true }; loop;)
		{
			if (!components.is_empty())
			{
				Composite_Format.parse_and_ignore_whitespace(source, pos);
				if (!Composite_Format.parse_fixed_string(source, trimmed_separator, pos))
				{
					loop = false;
				}
			}

			if (loop)
			{
				Composite_Format.parse_and_ignore_whitespace(source, pos);
				Number component = Composite_Format.parse_number(source, format, pos);
				if (component != NULL)
				{
					components.add(component);
				}
				else
				{
					// invalid component
					// set index back to initial, error index should already be set
					pos.set_index(initial_index);
					return NULL;
				}
			}
		}

		// parse suffix
		Composite_Format::parse_and_ignore_whitespace(source, pos);
		if (!Composite_Format.parse_fixed_string(source, trimmed_suffix, pos))
		{
			return NULL;
		}

		// build vector
		auto data = std::vector<double>(components.size());
		for (int i{}; i < data.size(); ++i)
		{
			data[i] = components.get(i).double_value();
		}
		return Array_Real_Vector(data, false);
	}
};