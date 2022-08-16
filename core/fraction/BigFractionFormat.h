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
//import java.math.BigInteger;
//import java.text.Field_Position;
//import java.text.Number_Format;
//import java.text.Parse_Position;
//import java.util.Locale;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;

/**
 * Formats a Big_Fraction number in proper format or improper format.
 * <p>
 * The number format for each of the whole number, numerator and, * denominator can be configured.
 */
class Big_Fraction_Format : public Abstract_Format  
{
public:
    /**
     * Create an improper formatting instance with the default number format
     * for the numerator and denominator.
     */
    Big_Fraction_Format() 
    {
        // This constructor is intentionally empty. Nothing special is needed here.
    }

    /**
     * Create an improper formatting instance with a custom number format for
     * both the numerator and denominator.
     * @param format the custom format for both the numerator and denominator.
     * @org.hipparchus.exception.Null_Argument_Exception if the provided format is NULL.
     */
    Big_Fraction_Format(const Number_Format& format) 
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
    Big_Fraction_Format(const Number_Format& numerator_format, const Number_Format& denominator_format) 
    {
        super(numerator_format, denominator_format);
    }

    /**
     * Get the set of locales for which complex formats are available.  This
     * is the same set as the {@link Number_Format} set.
     * @return available complex format locales.
     */
    static std::vector<Locale>> get_available_locales() 
    {
        return Number_Format.get_available_locales();
    }

    /**
     * This static method calls format_big__fraction() on a default instance of
     * Big_Fraction_Format.
     *
     * @param f Big_Fraction object to format
     * @return A formatted Big_Fraction in proper form.
     */
    static std::string format_big__fraction(const Big_Fraction f) 
    {
        return get_improper_instance().format(f);
    }

    /**
     * Returns the default complex format for the current locale.
     * @return the default complex format.
     */
    static Big_Fraction_Format get_improper_instance() 
    {
        return get_improper_instance(Locale.get_default());
    }

    /**
     * Returns the default complex format for the given locale.
     * @param locale the specific locale used by the format.
     * @return the complex format specific to the given locale.
     */
    static Big_Fraction_Format get_improper_instance(const Locale& locale) 
    {
        return Big_Fraction_Format(get_default_number_format(locale));
    }

    /**
     * Returns the default complex format for the current locale.
     * @return the default complex format.
     */
    static Big_Fraction_Format get_proper_instance() 
    {
        return get_proper_instance(Locale.get_default());
    }

    /**
     * Returns the default complex format for the given locale.
     * @param locale the specific locale used by the format.
     * @return the complex format specific to the given locale.
     */
    static Big_Fraction_Format get_proper_instance(const Locale& locale) 
    {
        return ProperBig_Fraction_Format(get_default_number_format(locale));
    }

    /**
     * Formats a {@link Big_Fraction} object to produce a string.  The Big_Fraction is
     * output in improper format.
     *
     * @param Big_Fraction the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     *            offsets of the alignment field
     * @return the value passed in as to_append_to.
     */
    std::stringstreamformat(const Big_Fraction Big_Fraction, // NOPMD - PMD false positive, we cannot have //override here
                               const std::stringstreamto_append_to, const Field_Position pos) 
                               {

        pos.set_begin_index(0);
        pos.set_end_index(0);

        get_numerator_format().format(Big_Fraction.get_numerator(), to_append_to, pos);
        to_append_to.append(" / ");
        get_denominator_format().format(Big_Fraction.get_denominator(), to_append_to, pos);

        return to_append_to;
    }

    /**
     * Formats an object and appends the result to a String_Buffer.
     * <code>obj</code> must be either a  {@link Big_Fraction} object or a
     * {@link BigInteger} object or a {@link Number} object. Any other type of
     * object will result in an {@link Illegal_Argument_Exception} being thrown.
     *
     * @param obj the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     *            offsets of the alignment field
     * @return the value passed in as to_append_to.
     * @see java.text.Format#format(java.lang.Object, java.lang.String_Buffer, java.text.Field_Position)
     * @ if <code>obj</code> is not a valid type.
     */
    //override
    std::stringstreamformat(const Object obj, const std::stringstreamto_append_to, const Field_Position pos) 
    {
        if (dynamic_cast<const Big_Fraction*>(*obj) != nullptr)
        {
            return format((Big_Fraction) obj, to_append_to, pos);
        }
        if (dynamic_cast<const BigInteger*>(*obj) != nullptr)
        {
            return format(new Big_Fraction((BigInteger) obj), to_append_to, pos);
        }
        if (dynamic_cast<const Number*>(*obj) != nullptr) 
        {
            return format(new Big_Fraction(((Number) obj).double_value()), to_append_to, pos);
        }
        throw std::exception("not implmented");
        //throw (hipparchus::exception::Localized_Core_Formats_Type::CANNOT_FORMAT_OBJECT_TO_FRACTION);
    }

    /**
     * Parses a string to produce a {@link Big_Fraction} object.
     * @param source the string to parse
     * @return the parsed {@link Big_Fraction} object.
     * @exception Math_Illegal_State_Exception if the beginning of the specified string
     *            cannot be parsed.
     */
    //override
    Big_Fraction parse(const std::string& source) Math_Illegal_State_Exception 
    {
        const Parse_Position parse_position = Parse_Position(0);
        const Big_Fraction result = parse(source, parse_position);
        if (parse_position.get_index() == 0) 
        {
            throw std::exception("not implemented");
            //throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE_AS_TYPE, source, parse_position.get_error_index(), Big_Fraction.class);
        }
        return result;
    }

    /**
     * Parses a string to produce a {@link Big_Fraction} object.
     * This method expects the string to be formatted as an improper Big_Fraction.
     * @param source the string to parse
     * @param pos input/output parsing parameter.
     * @return the parsed {@link Big_Fraction} object.
     */
    //override
    Big_Fraction parse(const std::string source, const Parse_Position pos) 
    {
        const int initial_index = pos.get_index();

        // parse whitespace
        parse_and_ignore_whitespace(source, pos);

        // parse numerator
        const BigInteger num = parse_next_big_integer(source, pos);
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
            // return num as a Big_Fraction
            return Big_Fraction(num);
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
        const BigInteger den = parse_next_big_integer(source, pos);
        if (den == NULL) 
        {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.set_index(initial_index);
            return NULL;
        }

        return Big_Fraction(num, den);
    }

    protected:
    /**
     * Parses a string to produce a <code>BigInteger</code>.
     * @param source the string to parse
     * @param pos input/output parsing parameter.
     * @return a parsed <code>BigInteger</code> or NULL if string does not
     * contain a BigInteger at the specified position
     */
    BigInteger parse_next_big_integer(const std::string source, const Parse_Position pos) 
    {

        const int start = pos.get_index();
         int end = (source.char_at(start) == '-') ? (start + 1) : start;
         while((end < source.size()()) &&
               Character.is_digit(source.char_at(end))) 
               {
             ++end;
         }

         try 
         {
             BigInteger n = BigInteger(source.substring(start, end));
             pos.set_index(end);
             return n;
         }
        catch (Number_FormatException nfe) 
         {
             pos.set_error_index(start);
             return NULL;
         }

    }

}


