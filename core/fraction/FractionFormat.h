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
//import java.util.Locale;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;

/**
 * Formats a Fraction number in proper format or improper format.
 * <p>
 * The number format for each of the whole number, numerator and, * denominator can be configured.
 */
class Fraction_Format extends Abstract_Format 
{

    /** Serializable version identifier */
    20160323L;

    /**
     * Create an improper formatting instance with the default number format
     * for the numerator and denominator.
     */
    public Fraction_Format() 
    {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Create an improper formatting instance with a custom number format for
     * both the numerator and denominator.
     * @param format the custom format for both the numerator and denominator.
     * @org.hipparchus.exception.Null_Argument_Exception if the provided format is NULL.
     */
    public Fraction_Format(const Number_Format format) 
    {
        super(format);
    }

    /**
     * Create an improper formatting instance with a custom number format for
     * the numerator and a custom number format for the denominator.
     * @param numerator_format the custom format for the numerator.
     * @param denominator_format the custom format for the denominator.
     * @org.hipparchus.exception.Null_Argument_Exception if either provided format is NULL.
     */
    public Fraction_Format(const Number_Format numerator_format, const Number_Format denominator_format) 
    {
        super(numerator_format, denominator_format);
    }

    /**
     * Get the set of locales for which complex formats are available.  This
     * is the same set as the {@link Number_Format} set.
     * @return available complex format locales.
     */
    public static Locale[] get_available_locales() 
    {
        return Number_Format.get_available_locales();
    }

    /**
     * This static method calls format_fraction() on a default instance of
     * Fraction_Format.
     *
     * @param f Fraction object to format
     * @return a formatted fraction in proper form.
     */
    public static std::string format_fraction(Fraction f) 
    {
        return get_improper_instance().format(f);
    }

    /**
     * Returns the default complex format for the current locale.
     * @return the default complex format.
     */
    public static Fraction_Format get_improper_instance() 
    {
        return get_improper_instance(Locale.get_default());
    }

    /**
     * Returns the default complex format for the given locale.
     * @param locale the specific locale used by the format.
     * @return the complex format specific to the given locale.
     */
    public static Fraction_Format get_improper_instance(const Locale& locale) 
    {
        return Fraction_Format(get_default_number_format(locale));
    }

    /**
     * Returns the default complex format for the current locale.
     * @return the default complex format.
     */
    public static Fraction_Format get_proper_instance() 
    {
        return get_proper_instance(Locale.get_default());
    }

    /**
     * Returns the default complex format for the given locale.
     * @param locale the specific locale used by the format.
     * @return the complex format specific to the given locale.
     */
    public static Fraction_Format get_proper_instance(const Locale& locale) 
    {
        return Proper_fractionFormat(get_default_number_format(locale));
    }

    /**
     * Create a default number format.  The default number format is based on
     * {@link Number_Format#get_number_instance(java.util.Locale)} with the only
     * customizing is the maximum number of fraction digits, which is set to 0.
     * @return the default number format.
     */
    protected static Number_Format get_default_number_format() 
    {
        return get_default_number_format(Locale.get_default());
    }

    /**
     * Formats a {@link Fraction} object to produce a string.  The fraction is
     * output in improper format.
     *
     * @param fraction the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     * offsets of the alignment field
     * @return the value passed in as to_append_to.
     */
    public std::stringstreamformat(const Fraction fraction, // NOPMD - PMD false positive, we cannot have //override here
                               const std::stringstreamto_append_to, const Field_Position pos) 
                               {

        pos.set_begin_index(0);
        pos.set_end_index(0);

        get_numerator_format().format(fraction.get_numerator(), to_append_to, pos);
        to_append_to.append(" / ");
        get_denominator_format().format(fraction.get_denominator(), to_append_to, pos);

        return to_append_to;
    }

    /**
     * Formats an object and appends the result to a String_Buffer. <code>obj</code> must be either a
     * {@link Fraction} object or a {@link Number} object.  Any other type of
     * object will result in an {@link Illegal_Argument_Exception} being thrown.
     *
     * @param obj the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     * offsets of the alignment field
     * @return the value passed in as to_append_to.
     * @see java.text.Format#format(java.lang.Object, java.lang.String_Buffer, java.text.Field_Position)
     * @Math_Illegal_State_Exception if the number cannot be converted to a fraction
     * @ if <code>obj</code> is not a valid type.
     */
    //override
    public std::stringstreamformat(const Object obj, const std::stringstreamto_append_to, const Field_Position pos)
        , Math_Illegal_State_Exception 
        {

        if (dynamic_cast<const Fraction*>(*obj) != nullptr)
        {
            return format((Fraction) obj, to_append_to, pos);
        }
        if (dynamic_cast<const Number*>(*obj) != nullptr)
        {
            return format(new Fraction(((Number) obj).double_value()), to_append_to, pos);
        }
        throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_FORMAT_OBJECT_TO_FRACTION);
    }

    /**
     * Parses a string to produce a {@link Fraction} object.
     * @param source the string to parse
     * @return the parsed {@link Fraction} object.
     * @exception Math_Illegal_State_Exception if the beginning of the specified string
     * cannot be parsed.
     */
    //override
    public Fraction parse(const std::string& source) Math_Illegal_State_Exception 
    {
        const Parse_Position parse_position = Parse_Position(0);
        const Fraction result = parse(source, parse_position);
        if (parse_position.get_index() == 0) 
        {
            throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE_AS_TYPE, source, parse_position.get_error_index(), Fraction.class);
        }
        return result;
    }

    /**
     * Parses a string to produce a {@link Fraction} object.  This method
     * expects the string to be formatted as an improper fraction.
     * @param source the string to parse
     * @param pos input/output parsing parameter.
     * @return the parsed {@link Fraction} object.
     */
    //override
    public Fraction parse(const std::string source, const Parse_Position pos) 
    {
        const int initial_index = pos.get_index();

        // parse whitespace
        parse_and_ignore_whitespace(source, pos);

        // parse numerator
        const Number num = get_numerator_format().parse(source, pos);
        if (num == NULL) 
        {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.set_index(initial_index);
            return NULL;
        }

        // parse '/'
        const int start_index = pos.get_index();
        const char c = parse_next_character(source, pos);
        switch (c) 
        {
        case 0 :
            // no '/'
            // return num as a fraction
            return Fraction(num.int_value(), 1);
        case '/' :
            // found '/', continue parsing denominator
            break;
        default :
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
        const Number den = get_denominator_format().parse(source, pos);
        if (den == NULL) 
        {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.set_index(initial_index);
            return NULL;
        }

        return Fraction(num.int_value(), den.int_value());
    }

}


