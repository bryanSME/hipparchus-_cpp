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

#include <cmath>
//import java.io.Serializable;
//import java.util.Collections;
//import java.util.Hash_Map;
//import java.util.Map;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Precision;
//import org.hipparchus.util.Sin_Cos;
//import org.hipparchus.util.Sinh_Cosh;
#include <unordered_map>
#include "../../CalculusFieldElement.hpp"

/**
 * First derivative computation with large number of variables.
 * <p>
 * This class plays a similar role to {@link Derivative_Structure}, with
 * a focus on efficiency when dealing with large number of independent variables
 * and most computation depend only on a few of them, and when only first derivative
 * is desired. When these conditions are met, this class should be much faster than
 * {@link Derivative_Structure} and use less memory.
 * </p>
 *
 */
class Sparse_Gradient : public Calculus_Field_Element<Sparse_Gradient>
{
private:

    /** Value of the calculation. */
    double my_value;

    /** Stored derivative, each key representing a different independent variable. */
    std::unordered_map<int, double> my_derivatives;

    /** Internal constructor.
     * @param value value of the function
     * @param derivatives derivatives map, a deep copy will be performed, * so the map given here will remain safe from changes in the instance, * may be null to create an empty derivatives map, i.e. a constant value
     */
    Sparse_Gradient(const double& value, const std::unordered_map<int, double>& derivatives) 
    {
        my_value = value;
        my_derivatives.put_all(derivatives);
    }

    /** Internal constructor.
     * @param value value of the function
     * @param scale scaling factor to apply to all derivatives
     * @param derivatives derivatives map, a deep copy will be performed, * so the map given here will remain safe from changes in the instance, * may be null to create an empty derivatives map, i.e. a constant value
     */
    Sparse_Gradient(const double& value, const double scale, const std::unordered_map<int, double>& derivatives) 
    {
        my_value = value;
        for (const Map.Entry<Integer, Double> entry : derivatives.entry_set()) 
        {
            my_derivatives.put(entry.get_key(), scale * entry.get_value());
        }
    }

    template<typename T>
    static int sgn(const T& val)
    {
        return (T(0) < val) - (val < T(0));
    }

public:
    /** {@inherit_doc} */
    //override
    Sparse_Gradient new_instance(const double& v) 
    {
        return create_constant(v);
    }

    /** Factory method creating a constant.
     * @param value value of the constant
     * @return a instance
     */
    static Sparse_Gradient create_constant(const double value) 
    {
        return Sparse_Gradient(value, Collections.<Integer, Double> empty_map());
    }

    /** Factory method creating an independent variable.
     * @param idx index of the variable
     * @param value value of the variable
     * @return a instance
     */
    static Sparse_Gradient create_variable(const int& idx, const double& value) 
    {
        return Sparse_Gradient(value, Collections.singleton_map(idx, 1.0));
    }

    /**
     * Find the number of variables.
     * @return number of variables
     */
    int num_vars() 
    {
        return my_derivatives.size();
    }

    /**
     * Get the derivative with respect to a particular index variable.
     *
     * @param index index to differentiate with.
     * @return derivative with respect to a particular index variable
     */
    double get_derivative(const int& index) const
    {
        return my_derivatives.contains(index)
            ? my_derivatives.at(index)
            : 0;
    }

    /**
     * Get the value of the function.
     * @return value of the function.
     */
    double get_value() const
    {
        return my_value;
    }

    /** {@inherit_doc} */
    //override
    double get_real() const
    {
        return my_value;
    }

    /** {@inherit_doc} */
    //override
    Sparse_Gradient add(const Sparse_Gradient& a) 
    {
        const Sparse_Gradient out = Sparse_Gradient(value + a.value, derivatives);
        for (Map.Entry<Integer, Double> entry : a.derivatives.entry_set()) 
        {
            const int id = entry.get_key();
            const Double old = out.derivatives.get(id);
            if (old == null) 
            {
                out.derivatives.put(id, entry.get_value());
            }
            else 
            {
                out.derivatives.put(id, old + entry.get_value());
            }
        }

        return out;
    }

    /**
     * Add in place.
     * <p>
     * This method is designed to be faster when used multiple times in a loop.
     * </p>
     * <p>
     * The instance is changed here, in order to not change the
     * instance the {@link #add(Sparse_Gradient)} method should
     * be used.
     * </p>
     * @param a instance to add
     */
    public void add_in_place(const Sparse_Gradient a) 
    {
        value += a.value;
        for (const Map.Entry<Integer, Double> entry : a.derivatives.entry_set()) 
        {
            const int id = entry.get_key();
            const Double old = derivatives.get(id);
            if (old == null) 
            {
                derivatives.put(id, entry.get_value());
            }
else 
            {
                derivatives.put(id, old + entry.get_value());
            }
        }
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient add(const double& c) 
    {
        return Sparse_Gradient(value + c, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient subtract(const Sparse_Gradient a) 
    {
        const Sparse_Gradient out = Sparse_Gradient(value - a.value, derivatives);
        for (Map.Entry<Integer, Double> entry : a.derivatives.entry_set()) 
        {
            const int id = entry.get_key();
            const Double old = out.derivatives.get(id);
            if (old == null) 
            {
                out.derivatives.put(id, -entry.get_value());
            }
else 
            {
                out.derivatives.put(id, old - entry.get_value());
            }
        }
        return out;
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient subtract(double c) 
    {
        return Sparse_Gradient(value - c, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient multiply(const Sparse_Gradient a) 
    {
        const Sparse_Gradient out =
            Sparse_Gradient(value * a.value, Collections.<Integer, Double> empty_map());

        // Derivatives.
        for (Map.Entry<Integer, Double> entry : derivatives.entry_set()) 
        {
            out.derivatives.put(entry.get_key(), a.value * entry.get_value());
        }
        for (Map.Entry<Integer, Double> entry : a.derivatives.entry_set()) 
        {
            const int id = entry.get_key();
            const Double old = out.derivatives.get(id);
            if (old == null) 
            {
                out.derivatives.put(id, value * entry.get_value());
            }
else 
            {
                out.derivatives.put(id, old + value * entry.get_value());
            }
        }
        return out;
    }

    /**
     * Multiply in place.
     * <p>
     * This method is designed to be faster when used multiple times in a loop.
     * </p>
     * <p>
     * The instance is changed here, in order to not change the
     * instance the {@link #add(Sparse_Gradient)} method should
     * be used.
     * </p>
     * @param a instance to multiply
     */
    public void multiply_in_place(const Sparse_Gradient a) 
    {
        // Derivatives.
        for (Map.Entry<Integer, Double> entry : derivatives.entry_set()) 
        {
            derivatives.put(entry.get_key(), a.value * entry.get_value());
        }
        for (Map.Entry<Integer, Double> entry : a.derivatives.entry_set()) 
        {
            const int id = entry.get_key();
            const Double old = derivatives.get(id);
            if (old == null) 
            {
                derivatives.put(id, value * entry.get_value());
            }
else 
            {
                derivatives.put(id, old + value * entry.get_value());
            }
        }
        value *= a.value;
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient multiply(const double& c) 
    {
        return Sparse_Gradient(value * c, c, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient multiply(const int& n) 
    {
        return Sparse_Gradient(value * n, n, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient divide(const Sparse_Gradient a) 
    {
        const Sparse_Gradient out = Sparse_Gradient(value / a.value, Collections.<Integer, Double> empty_map());

        // Derivatives.
        for (Map.Entry<Integer, Double> entry : derivatives.entry_set()) 
        {
            out.derivatives.put(entry.get_key(), entry.get_value() / a.value);
        }
        for (Map.Entry<Integer, Double> entry : a.derivatives.entry_set()) 
        {
            const int id = entry.get_key();
            const Double old = out.derivatives.get(id);
            if (old == null) 
            {
                out.derivatives.put(id, -out.value / a.value * entry.get_value());
            }
else 
            {
                out.derivatives.put(id, old - out.value / a.value * entry.get_value());
            }
        }
        return out;
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient divide(const double& c) 
    {
        return Sparse_Gradient(value / c, 1.0 / c, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient negate() 
    {
        return Sparse_Gradient(-value, -1.0, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Field<Sparse_Gradient> get_field() 
    {
        return Sparse_Gradient_Field.get_instance();
    }

    /** Local field for sparse gradient. */
    private static class Sparse_Gradient_Field : Field<Sparse_Gradient> 
    {

        /** Zero constant. */
        private const Sparse_Gradient zero;

        /** One constant. */
        private const Sparse_Gradient one;

        /** Private constructor for the singleton.
         */
        private Sparse_Gradient_Field() 
        {
            zero = create_constant(0);
            one  = create_constant(1);
        }

        /** Get the unique instance.
         * @return the unique instance
         */
        public static Sparse_Gradient_Field get_instance() 
        {
            return Lazy_Holder.INSTANCE;
        }

        /** {@inherit_doc} */
        //override
        public Sparse_Gradient get_zero() 
        {
            return zero;
        }

        /** {@inherit_doc} */
        //override
        public Sparse_Gradient get_one() 
        {
            return one;
        }

        /** {@inherit_doc} */
        //override
        public Class<Sparse_Gradient> get_runtime_class() 
        {
            return Sparse_Gradient.class;
        }

        /** {@inherit_doc} */
        //override
        public bool equals(const Object& other) 
        {
            return this == other;
        }

        /** {@inherit_doc} */
        //override
        public int hash_code() 
        {
            return 0x142aeff7;
        }

        // CHECKSTYLE: stop Hide_Utility_Class_Constructor
        /** Holder for the instance.
         * <p>We use here the Initialization On Demand Holder Idiom.</p>
         */
        private static class Lazy_Holder 
        {
            /** Cached field instance. */
            private static const Sparse_Gradient_Field INSTANCE = Sparse_Gradient_Field();
        }
        // CHECKSTYLE: resume Hide_Utility_Class_Constructor

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient remainder(const double& a) 
    {
        return Sparse_Gradient(std::remainder(value, a), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient remainder(const Sparse_Gradient a) 
    {

        // compute k such that lhs % rhs = lhs - k rhs
        const double rem = std::remainder(value, a.value);
        const double k   = std::rint((value - rem) / a.value);

        return subtract(a.multiply(k));

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient abs() 
    {
        if (Double.double_to_long_bits(value) < 0) 
        {
            // we use the bits representation to also handle -0.0
            return negate();
        }
else 
        {
            return this;
        }
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient ceil() 
    {
        return create_constant(std::ceil(value));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient floor() 
    {
        return create_constant(std::floor(value));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient rint() 
    {
        return create_constant(std::rint(value));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient sign() 
    {
        return create_constant(sgn(value));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient copy_sign(const Sparse_Gradient sign) 
    {
        const long m = Double.double_to_long_bits(value);
        const long s = Double.double_to_long_bits(sign.value);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient copy_sign(const double sign) 
    {
        const long m = Double.double_to_long_bits(value);
        const long s = Double.double_to_long_bits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient scalb(const int& n) 
    {
        const Sparse_Gradient out = Sparse_Gradient(std::scalbn(value, n), Collections.<Integer, Double> empty_map());
        for (Map.Entry<Integer, Double> entry : derivatives.entry_set()) 
        {
            out.derivatives.put(entry.get_key(), std::scalbn(entry.get_value(), n));
        }
        return out;
    }

    /** {@inherit_doc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    //override
    public Sparse_Gradient ulp() 
    {
        return new_instance(FastMath.ulp(value));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient hypot(const Sparse_Gradient y) 
    {
        if (std::isinf(value) || std::isinf(y.value)) 
        {
            return create_constant(INFINITY);
        }
else if (std::isnan(value) || std::isnan(y.value)) 
        {
            return create_constant(Double.NaN);
        }
else 
        {

            const int exp_x = FastMath.get_exponent(value);
            const int exp_y = FastMath.get_exponent(y.value);
            if (exp_x > exp_y + 27) 
            {
                // y is negligible with respect to x
                return abs();
            }
else if (exp_y > exp_x + 27) 
            {
                // x is negligible with respect to y
                return y.abs();
            }
else 
            {

                // find an intermediate scale to avoid both overflow and underflow
                const int middle_exp = (exp_x + exp_y) / 2;

                // scale parameters without losing precision
                const Sparse_Gradient scaled_x = scalb(-middle_exp);
                const Sparse_Gradient scaled_y = y.scalb(-middle_exp);

                // compute scaled hypotenuse
                const Sparse_Gradient scaled_h =
                        scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

                // remove scaling
                return scaled_h.scalb(middle_exp);

            }

        }
    }

    /**
     * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
     * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * @param x a value
     * @param y a value
     * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     */
    public static Sparse_Gradient hypot(const Sparse_Gradient x, const Sparse_Gradient y) 
    {
        return x.hypot(y);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient reciprocal() 
    {
        return Sparse_Gradient(1.0 / value, -1.0 / (value * value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient sqrt() 
    {
        const double sqrt = std::sqrt(value);
        return Sparse_Gradient(sqrt, 0.5 / sqrt, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient cbrt() 
    {
        const double cbrt = std::cbrt(value);
        return Sparse_Gradient(cbrt, 1.0 / (3 * cbrt * cbrt), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient root_n(const int& n) 
    {
        if (n == 2) 
        {
            return sqrt();
        }
else if (n == 3) 
        {
            return cbrt();
        }
else 
        {
            const double root = std::pow(value, 1.0 / n);
            return Sparse_Gradient(root, 1.0 / (n * std::pow(root, n - 1)), derivatives);
        }
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient pow(const double& p) 
    {
        return Sparse_Gradient(std::pow(value,  p), p * std::pow(value,  p - 1), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient pow(const int& n) 
    {
        if (n == 0) 
        {
            return get_field().get_one();
        }
else 
        {
            const double value_nm1 = std::pow(value,  n - 1);
            return Sparse_Gradient(value * value_nm1, n * value_nm1, derivatives);
        }
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient pow(const Sparse_Gradient e) 
    {
        return log().multiply(e).exp();
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link Sparse_Gradient}
     * @param a number to exponentiate
     * @param x power to apply
     * @return a<sup>x</sup>
     */
    public static Sparse_Gradient pow(const double& a, const Sparse_Gradient x) 
    {
        if (a == 0) 
        {
            if (x.value == 0) 
            {
                return x.compose(1.0, -INFINITY);
            }
else if (x.value < 0) 
            {
                return x.compose(Double.NaN,NAN);
            }
else 
            {
                return x.get_field().get_zero();
            }
        }
else 
        {
            const double& ax = std::pow(a, x.value);
            return Sparse_Gradient(ax, ax * std::log(a), x.derivatives);
        }
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient exp() 
    {
        const double e = std::exp(value);
        return Sparse_Gradient(e, e, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient expm1() 
    {
        return Sparse_Gradient(std::expm1(value), std::exp(value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient log() 
    {
        return Sparse_Gradient(std::log(value), 1.0 / value, derivatives);
    }

    /** Base 10 logarithm.
     * @return base 10 logarithm of the instance
     */
    //override
    public Sparse_Gradient log10() 
    {
        return Sparse_Gradient(std::log10(value), 1.0 / (std::log(10.0) * value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient log1p() 
    {
        return Sparse_Gradient(std::log1p(value), 1.0 / (1.0 + value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient cos() 
    {
        const Sin_Cos sc = Sin_Cos(value);
        return Sparse_Gradient(sc.cos(), -sc.sin(), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient sin() 
    {
        const Sin_Cos sc = Sin_Cos(value);
        return Sparse_Gradient(sc.sin(), sc.cos(), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Field_Sin_Cos<Sparse_Gradient> sin_cos() 
    {
        const Sin_Cos sc = Sin_Cos(value);
        return Field_Sin_Cos<>(new Sparse_Gradient(sc.sin(),  sc.cos(), derivatives), Sparse_Gradient(sc.cos(), -sc.sin(), derivatives));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient tan() 
    {
        const double t = std::tan(value);
        return Sparse_Gradient(t, 1 + t * t, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient acos() 
    {
        return Sparse_Gradient(std::acos(value), -1.0 / std::sqrt(1 - value * value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient asin() 
    {
        return Sparse_Gradient(std::asin(value), 1.0 / std::sqrt(1 - value * value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient atan() 
    {
        return Sparse_Gradient(std::atan(value), 1.0 / (1 + value * value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient atan2(const Sparse_Gradient x) 
    {

        // compute r = sqrt(x^2+y^2)
        const Sparse_Gradient r = multiply(this).add(x.multiply(x)).sqrt();

        const Sparse_Gradient a;
        if (x.value >= 0) 
        {

            // compute atan2(y, x) = 2 atan(y / (r + x))
            a = divide(r.add(x)).atan().multiply(2);

        }
else 
        {

            // compute atan2(y, x) = +/- pi - 2 atan(y / (r - x))
            const Sparse_Gradient tmp = divide(r.subtract(x)).atan().multiply(-2);
            a = tmp.add(tmp.value <= 0 ? -std::numbers::pi : std::numbers::pi);

        }

        // fix value to take special cases (+0/+0, +0/-0, -0/+0, -0/-0, +/-infinity) correctly
        a.value = std::atan2(value, x.value);

        return a;

    }

    /** Two arguments arc tangent operation.
     * @param y first argument of the arc tangent
     * @param x second argument of the arc tangent
     * @return atan2(y, x)
     */
    public static Sparse_Gradient atan2(const Sparse_Gradient y, const Sparse_Gradient x) 
    {
        return y.atan2(x);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient cosh() 
    {
        return Sparse_Gradient(std::cosh(value), std::sinh(value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient sinh() 
    {
        return Sparse_Gradient(std::sinh(value), std::cosh(value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Field_Sinh_Cosh<Sparse_Gradient> sinh_cosh() 
    {
        const Sinh_Cosh sch = std::sinh_cosh(value);
        return Field_Sinh_Cosh<>(new Sparse_Gradient(sch.sinh(), sch.cosh(), derivatives), Sparse_Gradient(sch.cosh(), sch.sinh(), derivatives));
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient tanh() 
    {
        const double t = std::tanh(value);
        return Sparse_Gradient(t, 1 - t * t, derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient acosh() 
    {
        return Sparse_Gradient(std::acosh(value), 1.0 / std::sqrt(value * value - 1.0), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient asinh() 
    {
        return Sparse_Gradient(std::asinh(value), 1.0 / std::sqrt(value * value + 1.0), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient atanh() 
    {
        return Sparse_Gradient(std::atanh(value), 1.0 / (1.0 - value * value), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient to_degrees() 
    {
        return Sparse_Gradient(FastMath.to_degrees(value), FastMath.to_degrees(1.0), derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient to_radians() 
    {
        return Sparse_Gradient(FastMath.to_radians(value), FastMath.to_radians(1.0), derivatives);
    }

    /** Evaluate Taylor expansion of a sparse gradient.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     */
    public double taylor(const double ... delta) 
    {
        double y = value;
        for (int i{}; i < delta.size(); ++i) 
        {
            y += delta[i] * get_derivative(i);
        }
        return y;
    }

    /** Compute composition of the instance by a univariate function.
     * @param f array of value and derivatives of the function at
     * the current point (i.e. [f({@link #get_value()}), * f'({@link #get_value()}), f''({@link #get_value()})...]).
     * @return f(this)
     * @exception  if the number of elements
     * in the array is not equal to 2 (i.e. value and first derivative)
     */
    public Sparse_Gradient compose(const double... f) 
    {
        Math_Utils::check_dimension(f.size(), 2);
        return Sparse_Gradient(f[0], f[1], derivatives);
    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const Sparse_Gradient[] a, const Sparse_Gradient[] b)
         
        {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = a[0].get_field().get_zero();
        for (int i{}; i < a.size(); ++i) 
        {
            out = out.add(a[i].multiply(b[i]));
        }

        // recompute an accurate value, taking care of cancellations
        const std::vector<double> a_double = std::vector<double>(a.size()];
        for (int i{}; i < a.size(); ++i) 
        {
            a_double[i] = a[i].get_value();
        }
        const std::vector<double> b_double = std::vector<double>(b.size()];
        for (int i{}; i < b.size(); ++i) 
        {
            b_double[i] = b[i].get_value();
        }
        out.value = Math_Arrays::linear_combination(a_double, b_double);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const std::vector<double> a, const Sparse_Gradient[] b) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = b[0].get_field().get_zero();
        for (int i{}; i < a.size(); ++i) 
        {
            out = out.add(b[i].multiply(a[i]));
        }

        // recompute an accurate value, taking care of cancellations
        const std::vector<double> b_double = std::vector<double>(b.size()];
        for (int i{}; i < b.size(); ++i) 
        {
            b_double[i] = b[i].get_value();
        }
        out.value = Math_Arrays::linear_combination(a, b_double);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const Sparse_Gradient a1, const Sparse_Gradient b1, const Sparse_Gradient a2, const Sparse_Gradient b2) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = a1.multiply(b1).add(a2.multiply(b2));

        // recompute an accurate value, taking care of cancellations
        out.value = Math_Arrays::linear_combination(a1.value, b1.value, a2.value, b2.value);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const double& a1, const Sparse_Gradient b1, const double& a2, const Sparse_Gradient b2) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = b1.multiply(a1).add(b2.multiply(a2));

        // recompute an accurate value, taking care of cancellations
        out.value = Math_Arrays::linear_combination(a1, b1.value, a2, b2.value);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const Sparse_Gradient a1, const Sparse_Gradient b1, const Sparse_Gradient a2, const Sparse_Gradient b2, const Sparse_Gradient a3, const Sparse_Gradient b3) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3));

        // recompute an accurate value, taking care of cancellations
        out.value = Math_Arrays::linear_combination(a1.value, b1.value, a2.value, b2.value, a3.value, b3.value);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const double& a1, const Sparse_Gradient b1, const double& a2, const Sparse_Gradient b2, const double& a3, const Sparse_Gradient b3) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = b1.multiply(a1).add(b2.multiply(a2)).add(b3.multiply(a3));

        // recompute an accurate value, taking care of cancellations
        out.value = Math_Arrays::linear_combination(a1, b1.value, a2, b2.value, a3, b3.value);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const Sparse_Gradient a1, const Sparse_Gradient b1, const Sparse_Gradient a2, const Sparse_Gradient b2, const Sparse_Gradient a3, const Sparse_Gradient b3, const Sparse_Gradient a4, const Sparse_Gradient b4) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3)).add(a4.multiply(b4));

        // recompute an accurate value, taking care of cancellations
        out.value = Math_Arrays::linear_combination(a1.value, b1.value, a2.value, b2.value, a3.value, b3.value, a4.value, b4.value);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient linear_combination(const double& a1, const Sparse_Gradient b1, const double& a2, const Sparse_Gradient b2, const double& a3, const Sparse_Gradient b3, const double& a4, const Sparse_Gradient b4) 
    {

        // compute a simple value, with all partial derivatives
        Sparse_Gradient out = b1.multiply(a1).add(b2.multiply(a2)).add(b3.multiply(a3)).add(b4.multiply(a4));

        // recompute an accurate value, taking care of cancellations
        out.value = Math_Arrays::linear_combination(a1, b1.value, a2, b2.value, a3, b3.value, a4, b4.value);

        return out;

    }

    /** {@inherit_doc} */
    //override
    public Sparse_Gradient get_pi() 
    {
        return Sparse_Gradient(std::numbers::pi, null);
    }

    /**
     * Test for the equality of two sparse gradients.
     * <p>
     * Sparse gradients are considered equal if they have the same value
     * and the same derivatives.
     * </p>
     * @param other Object to test for equality to this
     * @return true if two sparse gradients are equal
     */
    //override
    public bool equals(Object other) 
    {

        if (this == other) 
        {
            return true;
        }

        if (dynamic_cast<const Sparse_Gradient*>(*other) != nullptr)
        {
            const Sparse_Gradient rhs = (Sparse_Gradient)other;
            if (!Precision.equals(value, rhs.value, 1)) 
            {
                return false;
            }
            if (derivatives.size() != rhs.derivatives.size()) 
            {
                return false;
            }
            for (const Map.Entry<Integer, Double> entry : derivatives.entry_set()) 
            {
                if (!rhs.derivatives.contains_key(entry.get_key())) 
                {
                    return false;
                }
                if (!Precision.equals(entry.get_value(), rhs.derivatives.get(entry.get_key()), 1)) 
                {
                    return false;
                }
            }
            return true;
        }

        return false;

    }

    /**
     * Get a hash_code for the derivative structure.
     * @return a hash code value for this object
     */
    //override
    public int hash_code() 
    {
        return 743 + 809 * Math_Utils::hash(value) + 167 * derivatives.hash_code();
    }

}


