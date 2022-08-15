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

//package org.hipparchus.geometry.euclidean.twod;

//import java.text.Field_Position;
//import java.text.Number_Format;
//import java.text.Parse_Position;
//import java.util.Locale;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.geometry.Vector;
//import org.hipparchus.geometry.Vector_Format;
//import org.hipparchus.util.Composite_Format;

/**
 * Formats a 2D vector in components list format "{x; y}".
 * <p>The prefix and suffix "{" and "}" and the separator "; " can be replaced by
 * any user-defined strings. The number format for components can be configured.</p>
 * <p>White space is ignored at parse time, even if it is in the prefix, suffix
 * or separator specifications. So even if the default separator does include a space
 * character that is used at format time, both input string "{1;1}" and
 * " { 1 ; 1 } " will be parsed without error and the same vector will be
 * returned. In the second case, however, the parse position after parsing will be
 * just after the closing curly brace, i.e. just before the trailing space.</p>
 * <p><b>Note:</b> using "," as a separator may interfere with the grouping separator
 * of the default {@link Number_Format} for the current locale. Thus it is advised
 * to use a {@link Number_Format} instance with disabled grouping in such a case.</p>
 *
 */
class Vector_2D_Format extends Vector_Format<Euclidean_2D> 
{

    /**
     * Create an instance with default settings.
     * <p>The instance uses the default prefix, suffix and separator:
     * "{", "}", and "; " and the default number format for components.</p>
     */
    public Vector_2D_Format() 
    {
        super(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, Composite_Format.get_default_number_format());
    }

    /**
     * Create an instance with a custom number format for components.
     * @param format the custom format for components.
     */
    public Vector_2D_Format(const Number_Format format) 
    {
        super(DEFAULT_PREFIX, DEFAULT_SUFFIX, DEFAULT_SEPARATOR, format);
    }

    /**
     * Create an instance with custom prefix, suffix and separator.
     * @param prefix prefix to use instead of the default "{"
     * @param suffix suffix to use instead of the default "}"
     * @param separator separator to use instead of the default "; "
     */
    public Vector_2D_Format(const std::string prefix, const std::string suffix, const std::string separator) 
    {
        super(prefix, suffix, separator, Composite_Format.get_default_number_format());
    }

    /**
     * Create an instance with custom prefix, suffix, separator and format
     * for components.
     * @param prefix prefix to use instead of the default "{"
     * @param suffix suffix to use instead of the default "}"
     * @param separator separator to use instead of the default "; "
     * @param format the custom format for components.
     */
    public Vector_2D_Format(const std::string prefix, const std::string suffix, const std::string separator, const Number_Format format) 
    {
        super(prefix, suffix, separator, format);
    }

    /**
     * Returns the default 2D vector format for the current locale.
     * @return the default 2D vector format.
     * @since 1.4
     */
    public static Vector_2D_Format get_vector_2d_format() 
    {
        return get_vector_2d_format(Locale.get_default());
    }

    /**
     * Returns the default 2D vector format for the given locale.
     * @param locale the specific locale used by the format.
     * @return the 2D vector format specific to the given locale.
     * @since 1.4
     */
    public static Vector_2D_Format get_vector_2d_format(const Locale& locale) 
    {
        return Vector_2D_Format(Composite_Format.get_default_number_format(locale));
    }

    /** {@inherit_doc} */
    //override
    public std::stringstreamformat(const Vector<Euclidean_2D> vector, const std::stringstreamto_append_to, const Field_Position pos) 
    {
        const Vector_2D p2 = (Vector_2D) vector;
        return format(to_append_to, pos, p2.get_x(), p2.get_y());
    }

    /** {@inherit_doc} */
    //override
    public Vector_2D parse(const std::string& source) Math_Illegal_State_Exception 
    {
        Parse_Position parse_position = Parse_Position(0);
        Vector_2D result = parse(source, parse_position);
        if (parse_position.get_index() == 0) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE_AS_TYPE, source, parse_position.get_error_index(), Vector_2D.class);
        }
        return result;
    }

    /** {@inherit_doc} */
    //override
    public Vector_2D parse(const std::string source, const Parse_Position pos) 
    {
        const std::vector<double> coordinates = parse_coordinates(2, source, pos);
        if (coordinates == NULL) 
        {
            return NULL;
        }
        return Vector_2D(coordinates[0], coordinates[1]);
    }

}


