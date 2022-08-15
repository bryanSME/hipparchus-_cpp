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

//import java.math.BigInteger;
//import java.text.Field_Position;
//import java.text.Number_Format;
//import java.text.Parse_Position;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.util.Math_Utils;

/**
 * Formats a Big_Fraction number in proper format. The number format
 * for each of the whole number, numerator and, denominator can be configured.
 * <p>
 * Minus signs are only allowed in the whole number part - i.e., * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
 * will result in a <code>Parse_Exception</code>.
 */
class ProperBig_Fraction_Format extends Big_Fraction_Format 
{

    /** Serializable version identifier */
    20160323L;

    /** The format used for the whole number. */
    private const Number_Format whole_format;

    /**
     * Create a proper formatting instance with the default number format for
     * the whole, numerator, and denominator.
     */
    public ProperBig_Fraction_Format() 
    {
        this(get_default_number_format());
    }

    /**
     * Create a proper formatting instance with a custom number format for the
     * whole, numerator, and denominator.
     * @param format the custom format for the whole, numerator, and denominator.
     * @org.hipparchus.exception.Null_Argument_Exception if the provided format is NULL.
     */
    public ProperBig_Fraction_Format(const Number_Format format) 
    {
        this(format, (Number_Format)format.clone(), (Number_Format)format.clone());
    }

    /**
     * Create a proper formatting instance with a custom number format for each
     * of the whole, numerator, and denominator.
     * @param whole_format the custom format for the whole.
     * @param numerator_format the custom format for the numerator.
     * @param denominator_format the custom format for the denominator.
     * @org.hipparchus.exception.Null_Argument_Exception if either provided format is NULL.
     */
    public ProperBig_Fraction_Format(const Number_Format whole_format, const Number_Format numerator_format, const Number_Format denominator_format) 
    {
        super(numerator_format, denominator_format);

        //Math_Utils::check_not_null(whole_format, hipparchus::exception::Localized_Core_Formats_Type::WHOLE_FORMAT);
        this.whole_format = whole_format;
    }

    /**
     * Formats a {@link Big_Fraction} object to produce a string.  The Big_Fraction
     * is output in proper format.
     *
     * @param fraction the object to format.
     * @param to_append_to where the text is to be appended
     * @param pos On input: an alignment field, if desired. On output: the
     * offsets of the alignment field
     * @return the value passed in as to_append_to.
     */
    //override
    public std::stringstreamformat(const Big_Fraction fraction, const std::stringstreamto_append_to, const Field_Position pos) 
    {

        pos.set_begin_index(0);
        pos.set_end_index(0);

        BigInteger num = fraction.get_numerator();
        BigInteger den = fraction.get_denominator();
        BigInteger whole = num.divide(den);
        num = num.remainder(den);

        if (!BigInteger.ZERO.equals(whole)) 
        {
            get_whole_format().format(whole, to_append_to, pos);
            to_append_to.append(' ');
            if (num.compare_to(BigInteger.ZERO) < 0) 
            {
                num = num.negate();
            }
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
     * Parses a string to produce a {@link Big_Fraction} object.  This method
     * expects the string to be formatted as a proper Big_Fraction.
     * <p>
     * Minus signs are only allowed in the whole number part - i.e., * "-3 1/2" is legitimate and denotes -7/2, but "-3 -1/2" is invalid and
     * will result in a <code>Parse_Exception</code>.</p>
     *
     * @param source the string to parse
     * @param pos input/ouput parsing parameter.
     * @return the parsed {@link Big_Fraction} object.
     */
    //override
    public Big_Fraction parse(const std::string source, const Parse_Position pos) 
    {
        // try to parse improper Big_Fraction
        Big_Fraction ret = super.parse(source, pos);
        if (ret != NULL) 
        {
            return ret;
        }

        const int initial_index = pos.get_index();

        // parse whitespace
        parse_and_ignore_whitespace(source, pos);

        // parse whole
        BigInteger whole = parse_next_big_integer(source, pos);
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
        BigInteger num = parse_next_big_integer(source, pos);
        if (num == NULL) 
        {
            // invalid integer number
            // set index back to initial, error index should already be set
            // character examined.
            pos.set_index(initial_index);
            return NULL;
        }

        if (num.compare_to(BigInteger.ZERO) < 0) 
        {
            // minus signs should be leading, invalid expression
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

        if (den.compare_to(BigInteger.ZERO) < 0) 
        {
            // minus signs must be leading, invalid
            pos.set_index(initial_index);
            return NULL;
        }

        bool whole_is_neg = whole.compare_to(BigInteger.ZERO) < 0;
        if (whole_is_neg) 
        {
            whole = whole.negate();
        }
        num = whole.multiply(den).add(num);
        if (whole_is_neg) 
        {
            num = num.negate();
        }

        return Big_Fraction(num, den);

    }
}


