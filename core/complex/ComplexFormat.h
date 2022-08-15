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

//package org.hipparchus.complex;

//import java.text.Field_Position;
//import java.text.Number_Format;
//import java.text.Parse_Position;
//import java.util.Locale;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.util.Composite_Format;
//import org.hipparchus.util.Math_Utils;
#include <string>
#include "../exception/NullArgumentException.h"

/**
 * Formats a std::complex<double> number in cartesian format "Re(c) + Im(c)i".  'i' can
 * be replaced with 'j' (or anything else), and the number format for both real
 * and imaginary parts can be configured.
 */
class std::complex<double>_Format 
{
private:
     /** The default imaginary character. */
    inline static const std::string DEFAULT_IMAGINARY_CHARACTER = "i";
    /** The notation used to signify the imaginary part of the complex number. */
    const std::string my_imaginary_character;
    /** The format used for the imaginary part. */
    const Number_Format my_imaginary_format;
    /** The format used for the real part. */
    const Number_Format my_real_format;

    /**
     * Format the absolute value of the imaginary part.
     *
     * @param abs_im Absolute value of the imaginary part of a complex number.
     * @param to_append_to where the text is to be appended.
     * @param pos On input: an alignment field, if desired. On output: the
     * offsets of the alignment field.
     * @return the value passed in as to_append_to.
     */
    std::stringstream format_imaginary(double abs_im, std::stringstream to_append_to, Field_Position pos)
    {
        pos.set_begin_index(0);
        pos.set_end_index(0);

        Composite_Format.format_double(abs_im, get_imaginary_format(), to_append_to, pos);
        if ("1".equals(to_append_to.to_string()))
        {
            // Remove the character "1" if it is the only one.
            to_append_to.set_length(0);
        }

        return to_append_to;
    }

public:
    /**
     * Create an instance with the default imaginary character, 'i', and the
     * default number format for both real and imaginary parts.
     */
    std::complex<double>_Format() 
    {
        my_imaginary_character = DEFAULT_IMAGINARY_CHARACTER;
        my_imaginary_format = Composite_Format.get_default_number_format();
        my_real_format = my_imaginary_format;
    }

    /**
     * Create an instance with a custom number format for both real and
     * imaginary parts.
     * @param format the custom format for both real and imaginary parts.
     * @Null_Argument_Exception if {@code real_format} is {@code NULL}.
     */
    std::complex<double>_Format(const Number_Format& format) : my_imaginary_character{ DEFAULT_IMAGINARY_CHARACTER }, my_imaginary_format{ format }, my_real_format{ format } {};

    /**
     * Create an instance with a custom number format for the real part and a
     * custom number format for the imaginary part.
     * @param real_format the custom format for the real part.
     * @param imaginary_format the custom format for the imaginary part.
     * @Null_Argument_Exception if {@code imaginary_format} is {@code NULL}.
     * @Null_Argument_Exception if {@code real_format} is {@code NULL}.
      */
    std::complex<double>_Format(const Number_Format& real_format, const Number_Format& imaginary_format) : my_imaginary_character{ DEFAULT_IMAGINARY_CHARACTER }, my_imaginary_format{ imaginary_format }, my_real_format{ real_format } {};
    {
        //Math_Utils::check_not_null(imaginary_format, hipparchus::exception::Localized_Core_Formats_Type::IMAGINARY_FORMAT);
        //Math_Utils::check_not_null(real_format, hipparchus::exception::Localized_Core_Formats_Type::REAL_FORMAT);
    }

    /**
     * Create an instance with a custom imaginary character, and the default
     * number format for both real and imaginary parts.
     * @param imaginary_character The custom imaginary character.
     * @Null_Argument_Exception if {@code imaginary_character} is
     * {@code NULL}.
     * @ if {@code imaginary_character} is an
     * empty string.
     */
    std::complex<double>_Format(std::string imaginary_character)
    {
        this(imaginary_character, Composite_Format.get_default_number_format());
    }

    /**
     * Create an instance with a custom imaginary character, and a custom number
     * format for both real and imaginary parts.
     * @param imaginary_character The custom imaginary character.
     * @param format the custom format for both real and imaginary parts.
     * @Null_Argument_Exception if {@code imaginary_character} is
     * {@code NULL}.
     * @ if {@code imaginary_character} is an
     * empty string.
     * @Null_Argument_Exception if {@code format} is {@code NULL}.
     */
    std::complex<double>_Format(std::string imaginary_character, Number_Format format)
    {
        this(imaginary_character, format, format);
    }

    /**
     * Create an instance with a custom imaginary character, a custom number
     * format for the real part, and a custom number format for the imaginary
     * part.
     *
     * @param imaginary_character The custom imaginary character.
     * @param real_format the custom format for the real part.
     * @param imaginary_format the custom format for the imaginary part.
     * @Null_Argument_Exception if {@code imaginary_character} is
     * {@code NULL}.
     * @ if {@code imaginary_character} is an
     * empty string.
     * @Null_Argument_Exception if {@code imaginary_format} is {@code NULL}.
     * @Null_Argument_Exception if {@code real_format} is {@code NULL}.
     */
    std::complex<double>_Format(std::string imaginary_character, Number_Format real_format, Number_Format imaginary_format)
    {
        if (imaginary_character == NULL) 
        {
            throw Null_Argument_Exception();
        }
        if (imaginary_character.size() == 0) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
        }
        //Math_Utils::check_not_null(imaginary_format, hipparchus::exception::Localized_Core_Formats_Type::IMAGINARY_FORMAT);
        //Math_Utils::check_not_null(real_format, hipparchus::exception::Localized_Core_Formats_Type::REAL_FORMAT);

        this.imaginary_character = imaginary_character;
        this.imaginary_format = imaginary_format;
        this.real_format = real_format;
    }

    /**
     * Get the set of locales for which complex formats are available.
     * <p>This is the same set as the {@link Number_Format} set.</p>
     * @return available complex format locales.
     */
    static Locale[] get_available_locales() 
    {
        return Number_Format.get_available_locales();
    }

    /**
     * This method calls {@link #format(Object,String_Buffer,Field_Position)}.
     *
     * @param c std::complex<double> object to format.
     * @return A formatted number in the form "Re(c) + Im(c)i".
     */
    std::string format(std::complex<double> c) 
    {
        return format(c, String_Buffer(), Field_Position(0)).to_string();
    }

    /**
     * This method calls {@link #format(Object,String_Buffer,Field_Position)}.
     *
     * @param c Double object to format.
     * @return A formatted number.
     */
    std::string format(Double c) 
    {
        return format(new std::complex<double>(c, 0), String_Buffer(), Field_Position(0)).to_string();
    }

    /**
     * Formats a {@link std::complex<double>} object to produce a string.
     *
     * @param complex the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     *            offsets of the alignment field
     * @return the value passed in as to_append_to.
     */
    std::stringstream format(std::complex<double> complex, std::stringstream to_append_to, Field_Position pos) 
    {
        pos.set_begin_index(0);
        pos.set_end_index(0);

        // format real
        double re = complex.get_real();
        Composite_Format.format_double(re, get_real_format(), to_append_to, pos);

        // format sign and imaginary
        double im = complex.get_imaginary();
        std::stringstreamim_append_to;
        if (im < 0.0) 
        {
            to_append_to.append(" - ");
            im_append_to = format_imaginary(-im, String_Buffer(), pos);
            to_append_to.append(im_append_to);
            to_append_to.append(get_imaginary_character());
        }
        else if (im > 0.0 || std::isnan(im)) 
        {
            to_append_to.append(" + ");
            im_append_to = format_imaginary(im, String_Buffer(), pos);
            to_append_to.append(im_append_to);
            to_append_to.append(get_imaginary_character());
        }

        return to_append_to;
    }

    

    /**
     * Formats a object to produce a string.  {@code obj} must be either a
     * {@link std::complex<double>} object or a {@link Number} object.  Any other type of
     * object will result in an {@link Illegal_Argument_Exception} being thrown.
     *
     * @param obj the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     *            offsets of the alignment field
     * @return the value passed in as to_append_to.
     * @see java.text.Format#format(java.lang.Object, java.lang.String_Buffer, java.text.Field_Position)
     * @ is {@code obj} is not a valid type.
     */
    std::stringstreamformat(Object obj, std::stringstreamto_append_to, Field_Position pos)
    {
        if (dynamic_cast<const std::complex<double>*>(*obj) != nullptr)
        {
            return format( (std::complex<double>)obj, to_append_to, pos);
        }
        if (dynamic_cast<const Number*>(*obj) != nullptr)
        {
            return format(new std::complex<double>(((Number)obj).double_value(), 0.0), to_append_to, pos);
        }

        throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_FORMAT_INSTANCE_AS_COMPLEX, obj.get_class().get_name());
    }

    /**
     * Access the imaginary_character.
     * @return the imaginary_character.
     */
    std::string get_imaginary_character() const 
    {
        return imaginary_character;
    }

    /**
     * Access the imaginary_format.
     * @return the imaginary_format.
     */
    Number_Format get_imaginary_format() const 
    {
        return imaginary_format;
    }

    /**
     * Returns the default complex format for the current locale.
     * @return the default complex format.
     * @since 1.4
     */
    static std::complex<double>_Format get_complex_format() 
    {
        return get_complex_format(Locale.get_default());
    }

    /**
     * Returns the default complex format for the given locale.
     * @param locale the specific locale used by the format.
     * @return the complex format specific to the given locale.
     * @since 1.4
     */
    static std::complex<double>_Format get_complex_format(Locale locale) 
    {
        Number_Format f = Composite_Format.get_default_number_format(locale);
        return std::complex<double>_Format(f);
    }

    /**
     * Returns the default complex format for the given locale.
     * @param locale the specific locale used by the format.
     * @param imaginary_character Imaginary character.
     * @return the complex format specific to the given locale.
     * @Null_Argument_Exception if {@code imaginary_character} is
     * {@code NULL}.
     * @ if {@code imaginary_character} is an
     * empty string.
     * @since 1.4
     */
    static std::complex<double>_Format get_complex_format(std::string imaginary_character, Locale locale)
    {
        Number_Format f = Composite_Format.get_default_number_format(locale);
        return std::complex<double>_Format(imaginary_character, f);
    }

    /**
     * Access the real_format.
     * @return the real_format.
     */
    Number_Format get_real_format() const 
    {
        return real_format;
    }

    /**
     * Parses a string to produce a {@link std::complex<double>} object.
     *
     * @param source the string to parse.
     * @return the parsed {@link std::complex<double>} object.
     * @Math_Illegal_State_Exception if the beginning of the specified string
     * cannot be parsed.
     */
    std::complex<double> parse(std::string source) Math_Illegal_State_Exception 
    {
        Parse_Position parse_position = Parse_Position(0);
        std::complex<double> result = parse(source, parse_position);
        if (parse_position.get_index() == 0) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE_AS_TYPE, source, parse_position.get_error_index(), std::complex<double>.class);
        }
        return result;
    }

    /**
     * Parses a string to produce a {@link std::complex<double>} object.
     *
     * @param source the string to parse
     * @param pos input/ouput parsing parameter.
     * @return the parsed {@link std::complex<double>} object.
     */
    std::complex<double> parse(std::string source, Parse_Position pos) 
    {
        int initial_index = pos.get_index();

        // parse whitespace
        Composite_Format.parse_and_ignore_whitespace(source, pos);

        // parse real
        Number re = Composite_Format.parse_number(source, get_real_format(), pos);
        if (re == NULL) 
        {
            // invalid real number
            // set index back to initial, error index should already be set
            pos.set_index(initial_index);
            return NULL;
        }

        // parse sign
        int start_index = pos.get_index();
        char c = Composite_Format.parse_next_character(source, pos);
        int sign = 0;
        switch (c) 
        {
        case 0 :
            // no sign
            // return real only complex number
            return std::complex<double>(re.double_value(), 0.0);
        case '-' :
            sign = -1;
            break;
        case '+' :
            sign = 1;
            break;
        default :
            // invalid sign
            // set index back to initial, error index should be the last
            // character examined.
            pos.set_index(initial_index);
            pos.set_error_index(start_index);
            return NULL;
        }

        // parse whitespace
        Composite_Format.parse_and_ignore_whitespace(source, pos);

        // parse imaginary
        Number im = Composite_Format.parse_number(source, get_real_format(), pos);
        if (im == NULL) 
        {
            // invalid imaginary number
            // set index back to initial, error index should already be set
            pos.set_index(initial_index);
            return NULL;
        }

        // parse imaginary character
        if (!Composite_Format.parse_fixed_string(source, get_imaginary_character(), pos)) 
        {
            return NULL;
        }

        return std::complex<double>(re.double_value(), im.double_value() * sign);
    }
};