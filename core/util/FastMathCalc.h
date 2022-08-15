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

//import java.io.Print_Stream;

/**
 * Class used to compute the classical functions tables.
 */
class Fast_Math_Calc 
{

    /**
     * 0x40000000 - used to split a double into two parts, both with the low order bits cleared.
     * Equivalent to 2^30.
     */
    private static const long HEX_40000000 = 0x40000000L; // 1073741824L

    /** Factorial table, for Taylor series expansions. 0!, 1!, 2!, ... 19! */
    private static const double FACstd::vector<T> = 
    {
        +1.0,                        // 0
        +1.0,                        // 1
        +2.0,                        // 2
        +6.0d,                        // 3
        +24.0d,                       // 4
        +120.0,                      // 5
        +720.0,                      // 6
        +5040.0,                     // 7
        +40320.0,                    // 8
        +362880.0,                   // 9
        +3628800.0,                  // 10
        +39916800.0,                 // 11
        +479001600.0,                // 12
        +6227020800.0,               // 13
        +87178291200.0,              // 14
        +1307674368000.0,            // 15
        +20922789888000.0,           // 16
        +355687428096000.0,          // 17
        +6402373705728000.0,         // 18
        +121645100408832000.0,       // 19
        };

    /** Coefficients for slow_log. */
    private static const double LN_SPLIT_COEF[][] = 
    {
        {2.0, 0.0}, {0.6666666269302368, 3.9736429850260626E-8}, {0.3999999761581421, 2.3841857910019882E-8}, {0.2857142686843872, 1.7029898543501842E-8}, {0.2222222089767456, 1.3245471311735498E-8}, {0.1818181574344635, 2.4384203044354907E-8}, {0.1538461446762085, 9.140260083262505E-9}, {0.13333332538604736, 9.220590270857665E-9}, {0.11764700710773468, 1.2393345855018391E-8}, {0.10526403784751892, 8.251545029714408E-9}, {0.0952233225107193, 1.2675934823758863E-8}, {0.08713622391223907, 1.1430250008909141E-8}, {0.07842259109020233, 2.404307984052299E-9}, {0.08371849358081818, 1.176342548272881E-8}, {0.030589580535888672, 1.2958646899018938E-9}, {0.14982303977012634, 1.225743062930824E-8}, };

    /** Table start declaration. */
    private static const std::string TABLE_START_DECL = "    {";

    /** Table end declaration. */
    private static const std::string TABLE_END_DECL   = "    };";

    /**
     * Private Constructor.
     */
    private Fast_Math_Calc() 
    {
    }

    /** Build the sine and cosine tables.
     * @param SINE_TABLE_A table of the most significant part of the sines
     * @param SINE_TABLE_B table of the least significant part of the sines
     * @param COSINE_TABLE_A table of the most significant part of the cosines
     * @param COSINE_TABLE_B table of the most significant part of the cosines
     * @param SINE_TABLE_LEN length of the tables
     * @param TANGENT_TABLE_A table of the most significant part of the tangents
     * @param TANGENT_TABLE_B table of the most significant part of the tangents
     */
    @Suppress_Warnings("unused")
    private static void build_sin__cos_tables(std::vector<double> SINE_TABLE_A, std::vector<double> SINE_TABLE_B, std::vector<double> COSINE_TABLE_A, std::vector<double> COSINE_TABLE_B, int SINE_TABLE_LEN, std::vector<double> TANGENT_TABLE_A, std::vector<double> TANGENT_TABLE_B) 
    {
        const std::vector<double>& result = std::vector<double>(2);

        /* Use taylor series for 0 <= x <= 6/8 */
        for (int i{}; i < 7; i++) 
        {
            double x = i / 8.0;

            slow_sin(x, result);
            SINE_TABLE_A[i] = result[0];
            SINE_TABLE_B[i] = result[1];

            slow_cos(x, result);
            COSINE_TABLE_A[i] = result[0];
            COSINE_TABLE_B[i] = result[1];
        }

        /* Use angle addition formula to complete table to 13/8, just beyond pi/2 */
        for (int i = 7; i < SINE_TABLE_LEN; i++) 
        {
            double xs[] = std::vector<double>(2);
            double ys[] = std::vector<double>(2);
            double as[] = std::vector<double>(2);
            double bs[] = std::vector<double>(2);
            double temps[] = std::vector<double>(2);

            if ( (i & 1) == 0) 
            {
                // Even, use double angle
                xs[0] = SINE_TABLE_A[i/2];
                xs[1] = SINE_TABLE_B[i/2];
                ys[0] = COSINE_TABLE_A[i/2];
                ys[1] = COSINE_TABLE_B[i/2];

                /* compute sine */
                split_mult(xs, ys, result);
                SINE_TABLE_A[i] = result[0] * 2.0;
                SINE_TABLE_B[i] = result[1] * 2.0;

                /* Compute cosine */
                split_mult(ys, ys, as);
                split_mult(xs, xs, temps);
                temps[0] = -temps[0];
                temps[1] = -temps[1];
                split_add(as, temps, result);
                COSINE_TABLE_A[i] = result[0];
                COSINE_TABLE_B[i] = result[1];
            }
else 
            {
                xs[0] = SINE_TABLE_A[i/2];
                xs[1] = SINE_TABLE_B[i/2];
                ys[0] = COSINE_TABLE_A[i/2];
                ys[1] = COSINE_TABLE_B[i/2];
                as[0] = SINE_TABLE_A[i/2+1];
                as[1] = SINE_TABLE_B[i/2+1];
                bs[0] = COSINE_TABLE_A[i/2+1];
                bs[1] = COSINE_TABLE_B[i/2+1];

                /* compute sine */
                split_mult(xs, bs, temps);
                split_mult(ys, as, result);
                split_add(result, temps, result);
                SINE_TABLE_A[i] = result[0];
                SINE_TABLE_B[i] = result[1];

                /* Compute cosine */
                split_mult(ys, bs, result);
                split_mult(xs, as, temps);
                temps[0] = -temps[0];
                temps[1] = -temps[1];
                split_add(result, temps, result);
                COSINE_TABLE_A[i] = result[0];
                COSINE_TABLE_B[i] = result[1];
            }
        }

        /* Compute tangent = sine/cosine */
        for (int i{}; i < SINE_TABLE_LEN; i++) 
        {
            double xs[] = std::vector<double>(2);
            double ys[] = std::vector<double>(2);
            double as[] = std::vector<double>(2);

            as[0] = COSINE_TABLE_A[i];
            as[1] = COSINE_TABLE_B[i];

            split_reciprocal(as, ys);

            xs[0] = SINE_TABLE_A[i];
            xs[1] = SINE_TABLE_B[i];

            split_mult(xs, ys, as);

            TANGENT_TABLE_A[i] = as[0];
            TANGENT_TABLE_B[i] = as[1];
        }

    }

    /**
     *  For x between 0 and pi/4 compute cosine using Talor series
     *  cos(x) = 1 - x^2/2! + x^4/4! ...
     * @param x number from which cosine is requested
     * @param result placeholder where to put the result in extended precision
     * (may be NULL)
     * @return cos(x)
     */
    static double slow_cos(const double& x, const std::vector<double>& result) 
    {

        const auto xs = std::vector<double>(2);
        const auto ys = std::vector<double>(2);
        const auto facts = std::vector<double>(2);
        const auto as = std::vector<double>(2);
        split(x, xs);
        ys[0] = ys[1] = 0.0;

        for (int i = FACT.size()-1; i >= 0; i--) 
        {
            split_mult(xs, ys, as);
            ys[0] = as[0]; ys[1] = as[1];

            if ( (i & 1) != 0) { // skip odd entries
                continue;
            }

            split(FACT[i], as);
            split_reciprocal(as, facts);

            if ( (i & 2) != 0 ) { // alternate terms are negative
                facts[0] = -facts[0];
                facts[1] = -facts[1];
            }

            split_add(ys, facts, as);
            ys[0] = as[0]; ys[1] = as[1];
        }

        if (result != NULL) 
        {
            result[0] = ys[0];
            result[1] = ys[1];
        }

        return ys[0] + ys[1];
    }

    /**
     * For x between 0 and pi/4 compute sine using Taylor expansion:
     * sin(x) = x - x^3/3! + x^5/5! - x^7/7! ...
     * @param x number from which sine is requested
     * @param result placeholder where to put the result in extended precision
     * (may be NULL)
     * @return sin(x)
     */
    static double slow_sin(const double& x, const std::vector<double>& result) 
    {
        const auto xs = std::vector<double>(2);
        const auto ys = std::vector<double>(2);
        const auto facts = std::vector<double>(2);
        const auto as = std::vector<double>(2);
        split(x, xs);
        ys[0] = ys[1] = 0.0;

        for (int i = FACT.size()-1; i >= 0; i--) 
        {
            split_mult(xs, ys, as);
            ys[0] = as[0]; ys[1] = as[1];

            if ( (i & 1) == 0) { // Ignore even numbers
                continue;
            }

            split(FACT[i], as);
            split_reciprocal(as, facts);

            if ( (i & 2) != 0 ) { // alternate terms are negative
                facts[0] = -facts[0];
                facts[1] = -facts[1];
            }

            split_add(ys, facts, as);
            ys[0] = as[0]; ys[1] = as[1];
        }

        if (result != NULL) 
        {
            result[0] = ys[0];
            result[1] = ys[1];
        }

        return ys[0] + ys[1];
    }


    /**
     *  For x between 0 and 1, returns exp(x), uses extended precision
     *  @param x argument of exponential
     *  @param result placeholder where to place exp(x) split in two terms
     *  for extra precision (i.e. exp(x) = result[0] + result[1]
     *  @return exp(x)
     */
    static double slowexp(const double& x, const std::vector<double>& result) 
    {
        const auto xs = std::vector<double>(2);
        const auto ys = std::vector<double>(2);
        const auto facts = std::vector<double>(2);
        const auto as = std::vector<double>(2);
        split(x, xs);
        ys[0] = ys[1] = 0.0;

        for (int i = FACT.size()-1; i >= 0; i--) 
        {
            split_mult(xs, ys, as);
            ys[0] = as[0];
            ys[1] = as[1];

            split(FACT[i], as);
            split_reciprocal(as, facts);

            split_add(ys, facts, as);
            ys[0] = as[0];
            ys[1] = as[1];
        }

        if (result != NULL) 
        {
            result[0] = ys[0];
            result[1] = ys[1];
        }

        return ys[0] + ys[1];
    }

    /** Compute split[0], split[1] such that their sum is equal to d, * and split[0] has its 30 least significant bits as zero.
     * @param d number to split
     * @param split placeholder where to place the result
     */
    private static void split(const double d, const double split[]) 
    {
        if (d < 8e298 && d > -8e298) 
        {
            const double& a{ d * HEX_40000000 };
            split[0] = (d + a) - a;
            split[1] = d - split[0];
        }
else 
        {
            const double& a{ d * 9.31322574615478515625E-10 };
            split[0] = (d + a - d) * HEX_40000000;
            split[1] = d - split[0];
        }
    }

    /** Recompute a split.
     * @param a input/out array containing the split, changed
     * on output
     */
    private static void resplit(const std::vector<double>& a) 
    {
        const double c = a[0] + a[1];
        const double d = -(c - a[0] - a[1]);

        if (c < 8e298 && c > -8e298) { // MAGIC NUMBER
            double z = c * HEX_40000000;
            a[0] = (c + z) - z;
            a[1] = c - a[0] + d;
        }
else 
        {
            double z = c * 9.31322574615478515625E-10;
            a[0] = (c + z - c) * HEX_40000000;
            a[1] = c - a[0] + d;
        }
    }

    /** Multiply two numbers in split form.
     * @param a first term of multiplication
     * @param b second term of multiplication
     * @param ans placeholder where to put the result
     */
    private static void split_mult(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& ans) 
    {
        ans[0] = a[0] * b[0];
        ans[1] = a[0] * b[1] + a[1] * b[0] + a[1] * b[1];

        /* Resplit */
        resplit(ans);
    }

    /** Add two numbers in split form.
     * @param a first term of addition
     * @param b second term of addition
     * @param ans placeholder where to put the result
     */
    private static void split_add(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& ans) 
    {
        ans[0] = a[0] + b[0];
        ans[1] = a[1] + b[1];

        resplit(ans);
    }

    /** Compute the reciprocal of in.  Use the following algorithm.
     *  in = c + d.
     *  want to find x + y such that x+y = 1/(c+d) and x is much
     *  larger than y and x has several zero bits on the right.
     *
     *  Set b = 1/(2^22),  a = 1 - b.  Thus (a+b) = 1.
     *  Use following identity to compute (a+b)/(c+d)
     *
     *  (a+b)/(c+d)  =   a/c   +    (bc - ad) / (c^2 + cd)
     *  set x = a/c  and y = (bc - ad) / (c^2 + cd)
     *  This will be close to the right answer, but there will be
     *  some rounding in the calculation of X.  So by carefully
     *  computing 1 - (c+d)(x+y) we can compute an error and
     *  add that back in.   This is done carefully so that terms
     *  of similar size are subtracted first.
     *  @param in initial number, in split form
     *  @param result placeholder where to put the result
     */
    static void split_reciprocal(const std::vector<double>& in, const std::vector<double>& result) 
    {
        const double b = 1.0/4194304.0;
        const double& a = 1.0 - b;

        if (in[0] == 0.0) 
        {
            in[0] = in[1];
            in[1] = 0.0;
        }

        result[0] = a / in[0];
        result[1] = (b*in[0]-a*in[1]) / (in[0]*in[0] + in[0]*in[1]);

        if (result[1] != result[1]) { // can happen if result[1] is NAN
            result[1] = 0.0;
        }

        /* Resplit */
        resplit(result);

        for (int i{}; i < 2; i++) 
        {
            /* this may be overkill, probably once is enough */
            double err = 1.0 - result[0] * in[0] - result[0] * in[1] -
            result[1] * in[0] - result[1] * in[1];
            /*err = 1.0 - err; */
            err *= result[0] + result[1];
            /*printf("err = %16e\n", err); */
            result[1] += err;
        }
    }

    /** Compute (a[0] + a[1]) * (b[0] + b[1]) in extended precision.
     * @param a first term of the multiplication
     * @param b second term of the multiplication
     * @param result placeholder where to put the result
     */
    private static void quad_mult(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& result) 
    {
        const auto xs = std::vector<double>(2);
        const auto ys = std::vector<double>(2);
        const double zs[] = std::vector<double>(2);

        /* a[0] * b[0] */
        split(a[0], xs);
        split(b[0], ys);
        split_mult(xs, ys, zs);

        result[0] = zs[0];
        result[1] = zs[1];

        /* a[0] * b[1] */
        split(b[1], ys);
        split_mult(xs, ys, zs);

        double tmp = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;

        /* a[1] * b[0] */
        split(a[1], xs);
        split(b[0], ys);
        split_mult(xs, ys, zs);

        tmp = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;

        /* a[1] * b[0] */
        split(a[1], xs);
        split(b[1], ys);
        split_mult(xs, ys, zs);

        tmp = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;
    }

    /** Compute exp(p) for a integer p in extended precision.
     * @param p integer whose exponential is requested
     * @param result placeholder where to put the result in extended precision
     * @return exp(p) in standard precision (equal to result[0] + result[1])
     */
    static double expint(const int& p, const std::vector<double>& result) 
    {
        //double x = M_E;
        const auto xs = std::vector<double>(2);
        const auto as = std::vector<double>(2);
        const auto ys = std::vector<double>(2);
        //split(x, xs);
        //xs[1] = static_cast<double>((2.7182818284590452353602874713526625L - xs[0]);
        //xs[0] = 2.71827697753906250000;
        //xs[1] = 4.85091998273542816811e-06;
        //xs[0] = Double.long_bits_to_double(0x4005bf0800000000L);
        //xs[1] = Double.long_bits_to_double(0x3ed458a2bb4a9b00L);

        /* E */
        xs[0] = 2.718281828459045;
        xs[1] = 1.4456468917292502E-16;

        split(1.0, ys);

        while (p > 0) 
        {
            if ((p & 1) != 0) 
            {
                quad_mult(ys, xs, as);
                ys[0] = as[0]; ys[1] = as[1];
            }

            quad_mult(xs, xs, as);
            xs[0] = as[0]; xs[1] = as[1];

            p >>= 1;
        }

        if (result != NULL) 
        {
            result[0] = ys[0];
            result[1] = ys[1];

            resplit(result);
        }

        return ys[0] + ys[1];
    }
    /** xi in the range of [1, 2].
     *                                3        5        7
     *      x+1           /          x        x        x          \
     *  ln ----- =   2 *  |  x  +   ----  +  ----  +  ---- + ...  |
     *      1-x           \          3        5        7          /
     *
     * So, compute a Remez approximation of the following function
     *
     *  ln ((sqrt(x)+1)/(1-sqrt(x)))  /  x
     *
     * This will be an even function with only positive coefficents.
     * x is in the range [0 - 1/3].
     *
     * Transform xi for input to the above function by setting
     * x = (xi-1)/(xi+1).   Input to the polynomial is x^2, then
     * the result is multiplied by x.
     * @param xi number from which log is requested
     * @return log(xi)
     */
    static std::vector<double> slow_log(double xi) 
    {
        const std::vector<double>& x = std::vector<double>(2);
        double x2[] = std::vector<double>(2);
        const std::vector<double>& y = std::vector<double>(2);
        double a[] = std::vector<double>(2);

        split(xi, x);

        /* Set X = (x-1)/(x+1) */
        x[0] += 1.0;
        resplit(x);
        split_reciprocal(x, a);
        x[0] -= 2.0;
        resplit(x);
        split_mult(x, a, y);
        x[0] = y[0];
        x[1] = y[1];

        /* Square X -> X2*/
        split_mult(x, x, x2);


        //x[0] -= 1.0;
        //resplit(x);

        y[0] = LN_SPLIT_COEF[LN_SPLIT_COEF.size()-1][0];
        y[1] = LN_SPLIT_COEF[LN_SPLIT_COEF.size()-1][1];

        for (int i = LN_SPLIT_COEF.size()-2; i >= 0; i--) 
        {
            split_mult(y, x2, a);
            y[0] = a[0];
            y[1] = a[1];
            split_add(y, LN_SPLIT_COEF[i], a);
            y[0] = a[0];
            y[1] = a[1];
        }

        split_mult(y, x, a);
        y[0] = a[0];
        y[1] = a[1];

        return y;
    }


    /**
     * Print an array.
     * @param out text output stream where output should be printed
     * @param name array name
     * @param expected_len expected length of the array
     * @param array_2d array data
     */
    static void printarray(Print_Stream out, std::string name, int expected_len, std::vector<std::vector<double>> array_2d) 
    {
        out.println(name);
        Math_Utils::check_dimension(expected_len, array_2d.size());
        out.println(TABLE_START_DECL + " ");
        int i = 0;
        for(std::vector<double> array : array_2d) { // "double array[]" causes PMD parsing error
            out.print("        {");
            for(double d : array) { // assume inner array has very few entries
                out.printf("%-25.25s", format(d)); // multiple entries per line
            }
            out.println("}, // " + i++);
        }
        out.println(TABLE_END_DECL);
    }

    /**
     * Print an array.
     * @param out text output stream where output should be printed
     * @param name array name
     * @param expected_len expected length of the array
     * @param array array data
     */
    static void printarray(Print_Stream out, std::string name, int expected_len, std::vector<double> array) 
    {
        out.println(name + "=");
        Math_Utils::check_dimension(expected_len, array.size());
        out.println(TABLE_START_DECL);
        for(double d : array)
        {
            out.printf("        %s%n", format(d)); // one entry per line
        }
        out.println(TABLE_END_DECL);
    }

    /** Format a double.
     * @param d double number to format
     * @return formatted number
     */
    static std::string format(double d) 
    {
        if (d != d) 
        {
            return "Double.NaN,";
        }
else 
        {
            return ((d >= 0) ? "+" : "") + std::to_string(d) + "d,";
        }
    }

}


