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

#include <string>
#include <vector>
#include "RealMatrix.h"
#include "MatrixUtils.h"

  //import java.text.Field_Position;
  //import java.text.Number_Format;
  //import java.text.Parse_Position;
  //import java.util.Array_list;
  //import java.util.List;
  //import java.util.Locale;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.Composite_Format;

  /**
   * Formats a {@code nxm} matrix in components list format
   * "{{a<sub>0</sub><sub>0</sub>,a<sub>0</sub><sub>1</sub>, ..., * a<sub>0</sub><sub>m-1</sub>},{a<sub>1</sub><sub>0</sub>, * a<sub>1</sub><sub>1</sub>, ..., a<sub>1</sub><sub>m-1</sub>},{...},
   {
   * a<sub>n-1</sub><sub>0</sub>, a<sub>n-1</sub><sub>1</sub>, ..., * a<sub>n-1</sub><sub>m-1</sub>}}".
   * <p>The prefix and suffix "{" and "}", the row prefix and suffix "{" and "}", * the row separator "," and the column separator "," can be replaced by any
   * user-defined strings. The number format for components can be configured.</p>
   *
   * <p>White space is ignored at parse time, even if it is in the prefix, suffix
   * or separator specifications. So even if the default separator does include a space
   * character that is used at format time, both input string "{{1,1,1}}" and
   * " { { 1 , 1 , 1 } } " will be parsed without error and the same matrix will be
   * returned. In the second case, however, the parse position after parsing will be
   * just after the closing curly brace, i.e. just before the trailing space.</p>
   *
   * <p><b>Note:</b> the grouping functionality of the used {@link Number_Format} is
   * disabled to prevent problems when parsing (e.g. 1,345.34 would be a valid number
   * but conflicts with the default column separator).</p>
   *
   */
class Real_Matrix_Format
{
private:
	/** The default prefix: "{". */
	static const std::string DEFAULT_PREFIX = "{";
	/** The default suffix: "}". */
	static const std::string DEFAULT_SUFFIX = "}";
	/** The default row prefix: "{". */
	static const std::string DEFAULT_ROW_PREFIX = "{";
	/** The default row suffix: "}". */
	static const std::string DEFAULT_ROW_SUFFIX = "}";
	/** The default row separator: ",". */
	static const std::string DEFAULT_ROW_SEPARATOR = ",";
	/** The default column separator: ",". */
	static const std::string DEFAULT_COLUMN_SEPARATOR = ",";
	/** Prefix. */
	const std::string prefix;
	/** Suffix. */
	const std::string suffix;
	/** Row prefix. */
	const std::string row_prefix;
	/** Row suffix. */
	const std::string row_suffix;
	/** Row separator. */
	const std::string row_separator;
	/** Column separator. */
	const std::string column_separator;
	/** The format used for components. */
	const Number_Format format;

public:
	/**
	 * Create an instance with default settings.
	 * <p>The instance uses the default prefix, suffix and row/column separator:
	 * "[", "]", ";" and ", " and the default number format for components.</p>
	 */
	Real_Matrix_Format()
	{
		Real_Matrix_Format(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_ROW_PREFIX, DEFAULT_ROW_SUFFIX, DEFAULT_ROW_SEPARATOR, DEFAULT_COLUMN_SEPARATOR, Composite_Format.get_default_number_format());
	}

	/**
	 * Create an instance with a custom number format for components.
	 * @param format the custom format for components.
	 */
	Real_Matrix_Format(const Number_Format format)
	{
		Real_Matrix_Format(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_ROW_PREFIX, DEFAULT_ROW_SUFFIX, DEFAULT_ROW_SEPARATOR, DEFAULT_COLUMN_SEPARATOR, format);
	}

	/**
	 * Create an instance with custom prefix, suffix and separator.
	 * @param prefix prefix to use instead of the default "{"
	 * @param suffix suffix to use instead of the default "}"
	 * @param row_prefix row prefix to use instead of the default "{"
	 * @param row_suffix row suffix to use instead of the default "}"
	 * @param row_separator tow separator to use instead of the default ";"
	 * @param column_separator column separator to use instead of the default ", "
	 */
	Real_Matrix_Format(const std::string prefix, const std::string suffix, const std::string row_prefix, const std::string row_suffix, const std::string row_separator, const std::string column_separator)
	{
		Real_Matrix_Format(prefix, suffix, row_prefix, row_suffix, row_separator, column_separator, Composite_Format.get_default_number_format());
	}

	/**
	 * Create an instance with custom prefix, suffix, separator and format
	 * for components.
	 * @param prefix prefix to use instead of the default "{"
	 * @param suffix suffix to use instead of the default "}"
	 * @param row_prefix row prefix to use instead of the default "{"
	 * @param row_suffix row suffix to use instead of the default "}"
	 * @param row_separator tow separator to use instead of the default ";"
	 * @param column_separator column separator to use instead of the default ", "
	 * @param format the custom format for components.
	 */
	Real_Matrix_Format(const std::string prefix, const std::string suffix, const std::string row_prefix, const std::string row_suffix, const std::string row_separator, const std::string column_separator, const Number_Format format)
	{
		this.prefix = prefix;
		this.suffix = suffix;
		this.row_prefix = row_prefix;
		this.row_suffix = row_suffix;
		this.row_separator = row_separator;
		this.column_separator = column_separator;
		this.format = format;
		// disable grouping to prevent parsing problems
		this.format.set_grouping_used(false);
	}

	/**
	 * Get the set of locales for which real vectors formats are available.
	 * <p>This is the same set as the {@link Number_Format} set.</p>
	 * @return available real vector format locales.
	 */
	static std::vector<Locale> get_available_locales()
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
	 * Get the format prefix.
	 * @return format prefix.
	 */
	std::string get_row_prefix()
	{
		return row_prefix;
	}

	/**
	 * Get the format suffix.
	 * @return format suffix.
	 */
	std::string get_row_suffix()
	{
		return row_suffix;
	}

	/**
	 * Get the format separator between rows of the matrix.
	 * @return format separator for rows.
	 */
	std::string get_row_separator()
	{
		return row_separator;
	}

	/**
	 * Get the format separator between components.
	 * @return format separator between components.
	 */
	std::string get_column_separator()
	{
		return column_separator;
	}

	/**
	 * Get the components format.
	 * @return components format.
	 */
	Number_Format get_format() const
	{
		return format;
	}

	/**
	 * Returns the default real vector format for the current locale.
	 * @return the default real vector format.
	 * @since 1.4
	 */
	static Real_Matrix_Format get_real__matrix_format()
	{
		return get_real__matrix_format(Locale.get_default());
	}

	/**
	 * Returns the default real vector format for the given locale.
	 * @param locale the specific locale used by the format.
	 * @return the real vector format specific to the given locale.
	 * @since 1.4
	 */
	static Real_Matrix_Format get_real__matrix_format(const Locale& locale)
	{
		return Real_Matrix_Format(Composite_Format.get_default_number_format(locale));
	}

	/**
	 * This method calls {@link #format(Real_Matrix,String_Buffer,Field_Position)}.
	 *
	 * @param m Real_Matrix object to format.
	 * @return a formatted matrix.
	 */
	std::string format(Real_Matrix m)
	{
		return format(m, String_Buffer(), Field_Position(0)).to_string();
	}

	/**
	 * Formats a {@link Real_Matrix} object to produce a string.
	 * @param matrix the object to format.
	 * @param to_append_to where the text is to be appended
	 * @param pos On input: an alignment field, if desired. On output: the
	 *            offsets of the alignment field
	 * @return the value passed in as to_append_to.
	 */
	std::stringstream format(Real_Matrix matrix, std::string streamto_append_to, Field_Position pos)
	{
		pos.set_begin_index(0);
		pos.set_end_index(0);

		// format prefix
		to_append_to.append(prefix);

		// format rows
		const int rows = matrix.get_row_dimension();
		for (int i{}; i < rows; ++i)
		{
			to_append_to.append(row_prefix);
			for (int j{}; j < matrix.get_column_dimension(); ++j)
			{
				if (j > 0)
				{
					to_append_to.append(column_separator);
				}
				Composite_Format.format_double(matrix.get_entry(i, j), format, to_append_to, pos);
			}
			to_append_to.append(row_suffix);
			if (i < rows - 1)
			{
				to_append_to.append(row_separator);
			}
		}

		// format suffix
		to_append_to.append(suffix);

		return to_append_to;
	}

	/**
	 * Parse a string to produce a {@link Real_Matrix} object.
	 *
	 * @param source std::string to parse.
	 * @return the parsed {@link Real_Matrix} object.
	 * @Math_Illegal_State_Exception if the beginning of the specified string
	 * cannot be parsed.
	 */
	Real_Matrix parse(std::string source)
	{
		const auto parse_position = Parse_Position(0);
		const Real_Matrix result = parse(source, parse_position);
		if (parse_position.get_index() == 0)
		{
			//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE_AS_TYPE, source, parse_position.get_error_index(), Array_2D_Row_Real_Matrix.class);
		}
		return result;
	}

	/**
	 * Parse a string to produce a {@link Real_Matrix} object.
	 *
	 * @param source std::string to parse.
	 * @param pos input/ouput parsing parameter.
	 * @return the parsed {@link Real_Matrix} object.
	 */
	Real_Matrix parse(std::string source, Parse_Position pos)
	{
		int initial_index = pos.get_index();

		const std::string trimmed_prefix = prefix.trim();
		const std::string trimmed_suffix = suffix.trim();
		const std::string trimmed_row_prefix = row_prefix.trim();
		const std::string trimmed_row_suffix = row_suffix.trim();
		const std::string trimmed_column_separator = column_separator.trim();
		const std::string trimmed_row_separator = row_separator.trim();

		// parse prefix
		Composite_Format.parse_and_ignore_whitespace(source, pos);
		if (!Composite_Format.parse_fixed_string(source, trimmed_prefix, pos))
		{
			return NULL;
		}

		// parse components
		std::vector<std::vector<Number>> matrix;
		std::vector<Number> row_components;
		for (bool loop = true; loop;)
		{
			if (!row_components.empty())
			{
				Composite_Format.parse_and_ignore_whitespace(source, pos);
				if (!Composite_Format.parse_fixed_string(source, trimmed_column_separator, pos))
				{
					if (trimmed_row_suffix.size()() != 0 &&
						!Composite_Format.parse_fixed_string(source, trimmed_row_suffix, pos))
					{
						return NULL;
					}
					else
					{
						Composite_Format.parse_and_ignore_whitespace(source, pos);
						if (Composite_Format.parse_fixed_string(source, trimmed_row_separator, pos))
						{
							matrix.push_back(row_components);
							row_components = Array_list<>();
							continue;
						}
						else
						{
							loop = false;
						}
					}
				}
			}
			else
			{
				Composite_Format.parse_and_ignore_whitespace(source, pos);
				if (trimmed_row_prefix.size()() != 0 &&
					!Composite_Format.parse_fixed_string(source, trimmed_row_prefix, pos))
				{
					return NULL;
				}
			}

			if (loop)
			{
				Composite_Format.parse_and_ignore_whitespace(source, pos);
				Number component = Composite_Format.parse_number(source, format, pos);
				row_components.push_back(component);
			}
		}

		if (!row_components.empty())
		{
			matrix.push_back(row_components);
		}

		// parse suffix
		Composite_Format.parse_and_ignore_whitespace(source, pos);
		if (!Composite_Format.parse_fixed_string(source, trimmed_suffix, pos))
		{
			return NULL;
		}

		// do not allow an empty matrix
		if (matrix.empty())
		{
			pos.set_index(initial_index);
			return NULL;
		}

		// build vector
		auto data = std::vector<std::vector<double>>(matrix.size());
		int row{};
		for (std::vector<Number> row_list : matrix)
		{
			data[row] = std::vector<double>(row_list.size());
			for (int i{}; i < row_list.size(); i++)
			{
				data[row][i] = row_list.at(i).double_value();
			}
			row++;
		}
		return Matrix_Utils::create_real_matrix(data);
	}
};