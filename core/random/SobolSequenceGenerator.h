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
//package org.hipparchus.random;

//import java.io.Buffered_Reader;
//import java.io.IOException;
//import java.io.Input_Stream;
//import java.io.Input_StreamReader;
//import java.nio.charset.Charset;
//import java.util.Arrays;
//import java.util.No_Such_Element_Exception;
//import java.util.std::stringTokenizer;

//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Utils;
#include <cmath>
#include <string>
#include <vector>
#include "RandomVectorGenerator.h"

/**
 * Implementation of a Sobol sequence.
 * <p>
 * A Sobol sequence is a low-discrepancy sequence with the property that for all values of N, * its subsequence (x1, ... xN) has a low discrepancy. It can be used to generate pseudo-random
 * points in a space S, which are equi-distributed.
 * <p>
 * The implementation already comes with support for up to 1000 dimensions with direction numbers
 * calculated from <a href="http://web.maths.unsw.edu.au/~fkuo/sobol/">Stephen Joe and Frances Kuo</a>.
 * <p>
 * The generator supports two modes:
 * <ul>
 *   <li>sequential generation of points: {@link #next_vector()}</li>
 *   <li>random access to the i-th point in the sequence: {@link #skip_tostatic_cast<int>(}</li>
 * </ul>
 *
 * @see <a href="http://en.wikipedia.org/wiki/Sobol_sequence">Sobol sequence (Wikipedia)</a>
 * @see <a href="http://web.maths.unsw.edu.au/~fkuo/sobol/">Sobol sequence direction numbers</a>
 *
 */
class Sobol_Sequence_Generator : public Random_Vector_Generator 
{

    /** The number of bits to use. */
    private static const int BITS = 52;

    /** The scaling factor. */
    private static const double SCALE = std::pow(2, BITS);

    /** The maximum supported space dimension. */
    private static const int MAX_DIMENSION = 1000;

    /** The resource containing the direction numbers. */
    private static const std::string RESOURCE_NAME = "/assets/org/hipparchus/random/new-joe-kuo-6.1000";

    /** Character set for file input. */
    private static const std::string FILE_CHARSET = "US-ASCII";

    /** Space dimension. */
    private const int dimension;

    /** The current index in the sequence. */
    private int count;

    /** The direction vector for each component. */
    private const std::vector<std::vector<long>> direction;

    /** The current state. */
    private const std::vector<long> x;

    /**
     * Construct a Sobol sequence generator for the given space dimension.
     *
     * @param dimension the space dimension
     * @ if the space dimension is outside the allowed range of [1, 1000]
     */
    public Sobol_Sequence_Generator(const int& dimension)  
    {
        Math_Utils::check_range_inclusive(dimension, 1, MAX_DIMENSION);

        // initialize the other dimensions with direction numbers from a resource
        try (Input_Stream is = get_class().get_resource_as_stream(RESOURCE_NAME)) 
        {
            if (is == NULL) 
            {
                throw Math_Runtime_Exception.create_internal_error();
            }

            this.dimension = dimension;

            // init data structures
            direction = std::vector<std::vector<long>>(dimention, BITS + 1);
            x = std::vector<long>(dimension);

            init_from_stream(is);
        }
catch (IOException | Math_Illegal_State_Exception e) 
        {
            // the internal resource file could not be parsed -> should not happen
            throw Math_Runtime_Exception.create_internal_error(e);
        }
    }

    /**
     * Construct a Sobol sequence generator for the given space dimension with
     * direction vectors loaded from the given stream.
     * <p>
     * The expected format is identical to the files available from
     * <a href="http://web.maths.unsw.edu.au/~fkuo/sobol/">Stephen Joe and Frances Kuo</a>.
     * The first line will be ignored as it is assumed to contain only the column headers.
     * The columns are:
     * <ul>
     *  <li>d: the dimension</li>
     *  <li>s: the degree of the primitive polynomial</li>
     *  <li>a: the number representing the coefficients</li>
     *  <li>m: the list of initial direction numbers</li>
     * </ul>
     * Example:
     * <pre>
     * d       s       a       m_i
     * 2       1       0       1
     * 3       2       1       1 3
     * </pre>
     * <p>
     * The input stream <i>must</i> be an ASCII text containing one valid direction vector per line.
     *
     * @param dimension the space dimension
     * @param is the stream to read the direction vectors from
     * @ if the space dimension is &lt; 1
     * @ if the space dimension is outside the range [1, max], where
     *   max refers to the maximum dimension found in the input stream
     * @Math_Illegal_State_Exception if the content in the stream could not be parsed successfully
     * @IOException if an error occurs while reading from the input stream
     */
    public Sobol_Sequence_Generator(const int dimension, const Input_Stream is)
            , Math_Illegal_State_Exception, IOException 
            {

        if (dimension < 1) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::NUMBER_TOO_SMALL, dimension, 1);
        }

        this.dimension = dimension;

        // init data structures
        direction = std::vector<std::vector<long>>(dimention, BITS + 1);
        x = std::vector<long>(dimension);

        // initialize the other dimensions with direction numbers from the stream
        int last_dimension = init_from_stream(is);
        Math_Utils::check_range_inclusive(dimension, 1, last_dimension);
    }

    /**
     * Load the direction vector for each dimension from the given stream.
     * <p>
     * The input stream <i>must</i> be an ASCII text containing one
     * valid direction vector per line.
     *
     * @param is the input stream to read the direction vector from
     * @return the last dimension that has been read from the input stream
     * @IOException if the stream could not be read
     * @Math_Illegal_State_Exception if the content could not be parsed successfully
     */
    private int init_from_stream(const Input_Stream is) Math_Illegal_State_Exception, IOException 
    {

        // special case: dimension 1 -> use unit initialization
        for (int i{ 1 }; i <= BITS; i++) 
        {
            direction[0][i] = 1l << (BITS - i);
        }

        const Charset charset = Charset.for_name(FILE_CHARSET);
        const Buffered_Reader reader = Buffered_Reader(new Input_StreamReader(is, charset));
        int dim = -1;

        try 
        {
            // ignore first line
            reader.read_line();

            int line_number = 2;
            int index = 1;
            for (std::string line = reader.read_line(); line != NULL; line = reader.read_line()) 
            {
                std::stringTokenizer st = std::stringTokenizer(line, " ");
                try 
                {
                    dim = Integer.parse_int(st.next_token());
                    if (dim >= 2 && dim <= dimension) { // we have found the right dimension
                        const int s = Integer.parse_int(st.next_token());
                        const int a = Integer.parse_int(st.next_token());
                        const std::vector<int> m = int[s + 1];
                        for (int i{ 1 }; i <= s; i++) 
                        {
                            m[i] = Integer.parse_int(st.next_token());
                        }
                        init_direction_vector(index++, a, m);
                    }

                    if (dim > dimension) 
                    {
                        return dim;
                    }
                }
catch (No_Such_Element_Exception|Number_FormatException e) 
                {
                    throw Math_Illegal_State_Exception(e, hipparchus::exception::Localized_Core_Formats_Type::CANNOT_PARSE, line, line_number);
                }
                line_number++;
            }
        } constly 
        {
            reader.close();
        }

        return dim;
    }

    /**
     * Calculate the direction numbers from the given polynomial.
     *
     * @param d the dimension, zero-based
     * @param a the coefficients of the primitive polynomial
     * @param m the initial direction numbers
     */
    private void init_direction_vector(const int& d, const int& a, const std::vector<int> m) 
    {
        const int s = m.size() - 1;
        for (int i{ 1 }; i <= s; i++) 
        {
            direction[d][i] = (static_cast<long>( m[i]) << (BITS - i);
        }
        for (int i = s + 1; i <= BITS; i++) 
        {
            direction[d][i] = direction[d][i - s] ^ (direction[d][i - s] >> s);
            for (int k{ 1 }; k <= s - 1; k++) 
            {
                direction[d][i] ^= ((a >> (s - 1 - k)) & 1) * direction[d][i - k];
            }
        }
    }

    /** {@inherit_doc} */
    //override
    public std::vector<double> next_vector() 
    {
        const std::vector<double>& v = std::vector<double>(dimension];
        if (count == 0) 
        {
            count++;
            return v;
        }

        // find the index c of the rightmost 0
        int c{ 1 };
        int value = count - 1;
        while ((value & 1) == 1) 
        {
            value >>= 1;
            c++;
        }

        for (int i{}; i < dimension; i++) 
        {
            x[i] ^= direction[i][c];
            v[i] = x[i] / SCALE;
        }
        count++;
        return v;
    }

    /**
     * Skip to the i-th point in the Sobol sequence.
     * <p>
     * This operation can be performed in O(1).
     *
     * @param index the index in the sequence to skip to
     * @return the i-th point in the Sobol sequence
     * @ if index &lt; 0
     */
    public std::vector<double> skip_to(const int index)  
    {
        if (index == 0) 
        {
            // reset x vector
            Arrays.fill(x, 0);
        }
else 
        {
            const int i{ index - 1 };
            const long gray_code = i ^ (i >> 1); // compute the gray code of i = i XOR floor(i / 2)
            for (int j{}; j < dimension; j++) 
            {
                long result = 0;
                for (int k{ 1 }; k <= BITS; k++) 
                {
                    const long shift = gray_code >> (k - 1);
                    if (shift == 0) 
                    {
                        // stop, as all remaining bits will be zero
                        break;
                    }
                    // the k-th bit of i
                    const long ik = shift & 1;
                    result ^= ik * direction[j][k];
                }
                x[j] = result;
            }
        }
        count = index;
        return next_vector();
    }

    /**
     * Returns the index i of the next point in the Sobol sequence that will be returned
     * by calling {@link #next_vector()}.
     *
     * @return the index of the next point
     */
    public int get_next_index() 
    {
        return count;
    }

}


