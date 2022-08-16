#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
#include <cmath>
#include <numbers>
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Field_Sinh_Cosh;
//import org.hipparchus.util.Math_Arrays;
//import org.hipparchus.util.Math_Utils;
//import org.hipparchus.util.Sin_Cos;
//import org.hipparchus.util.Sinh_Cosh;
#include "UnivariateDerivative.hpp"


/** Class representing both the value and the differentials of a function.
 * <p>This class is a stripped-down version of {@link Derivative_Structure}
 * with only one {@link Derivative_Structure#get_free_parameters() free parameter}
 * and {@link Derivative_Structure#get_order() derivation order} also limited to one.
 * It should have less overhead than {@link Derivative_Structure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link Univariate_Derivative_1} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see Derivative_Structure
 * @see Univariate_Derivative_2
 * @see Gradient
 * @see Field_Derivative_Structure
 * @see Field_Univariate_Derivative_1
 * @see Field_Univariate_Derivative_2
 * @see Field_Gradient
 * @since 1.7
 */
class Univariate_Derivative_1 : public Univariate_Derivative<Univariate_Derivative_1> 
{

    /** The constant value of π as a {@code Univariate_Derivative_1}.
     * @since 2.0
     */
    public static const Univariate_Derivative_1 PI = Univariate_Derivative_1(std::numbers::pi, 0.0);

    /** Value of the function. */
    private const double f0;

    /** First derivative of the function. */
    private const double f1;

    /** Build an instance with values and derivative.
     * @param f0 value of the function
     * @param f1 first derivative of the function
     */
    public Univariate_Derivative_1(const double f0, const double f1) 
    {
        this.f0 = f0;
        this.f1 = f1;
    }

    /** Build an instance from a {@link Derivative_Structure}.
     * @param ds derivative structure
     * @exception  if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 1
     */
    public Univariate_Derivative_1(const Derivative_Structure ds)  
    {
        Math_Utils::check_dimension(ds.get_free_parameters(), 1);
        Math_Utils::check_dimension(ds.get_order(), 1);
        this.f0 = ds.get_value();
        this.f1 = ds.get_partial_derivative(1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 new_instance(const double value) 
    {
        return Univariate_Derivative_1(value, 0.0);
    }

    /** {@inherit_doc} */
    //override
    public double get_real() 
    {
        return get_value();
    }

    /** {@inherit_doc} */
    //override
    public double get_value() 
    {
        return f0;
    }

    /** {@inherit_doc} */
    //override
    public double get_derivative(const int& n) 
    {
        switch (n) 
        {
            case 0 :
                return f0;
            case 1 :
                return f1;
            default :
                throw std::exception("not implmented");
                //throw (hipparchus::exception::Localized_Core_Formats_Type::DERIVATION_ORDER_NOT_ALLOWED, n);
        }
    }

    /** {@inherit_doc} */
    //override
    public int get_order() 
    {
        return 1;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #get_value()
     */
    public double get_first_derivative() 
    {
        return f1;
    }

    /** {@inherit_doc} */
    //override
    public Derivative_Structure to_derivative_structure() 
    {
        return get_field().get_conversion_factory().build(f0, f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 add(const double& a) 
    {
        return Univariate_Derivative_1(f0 + a, f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 add(const Univariate_Derivative_1 a) 
    {
        return Univariate_Derivative_1(f0 + a.f0, f1 + a.f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 subtract(const double& a) 
    {
        return Univariate_Derivative_1(f0 - a, f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 subtract(const Univariate_Derivative_1 a) 
    {
        return Univariate_Derivative_1(f0 - a.f0, f1 - a.f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 multiply(const int& n) 
    {
        return Univariate_Derivative_1(f0 * n, f1 * n);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 multiply(const double& a) 
    {
        return Univariate_Derivative_1(f0 * a, f1 * a);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 multiply(const Univariate_Derivative_1 a) 
    {
        return Univariate_Derivative_1(f0 * a.f0, Math_Arrays::linear_combination(f1, a.f0, f0, a.f1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 divide(const double& a) 
    {
        const double inv1 = 1.0 / a;
        return Univariate_Derivative_1(f0 * inv1, f1 * inv1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 divide(const Univariate_Derivative_1 a) 
    {
        const double inv1 = 1.0 / a.f0;
        const double inv2 = inv1 * inv1;
        return Univariate_Derivative_1(f0 * inv1, Math_Arrays::linear_combination(f1, a.f0, -f0, a.f1) * inv2);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 remainder(const double& a) 
    {
        return Univariate_Derivative_1(std::remainder(f0, a), f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 remainder(const Univariate_Derivative_1 a) 
    {

        // compute k such that lhs % rhs = lhs - k rhs
        const double rem = std::remainder(f0, a.f0);
        const double k   = std::rint((f0 - rem) / a.f0);

        return Univariate_Derivative_1(rem, f1 - k * a.f1);

    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 negate() 
    {
        return Univariate_Derivative_1(-f0, -f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 abs() 
    {
        if (Double.double_to_long_bits(f0) < 0) 
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
    public Univariate_Derivative_1 ceil() 
    {
        return Univariate_Derivative_1(std::ceil(f0), 0.0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 floor() 
    {
        return Univariate_Derivative_1(std::floor(f0), 0.0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 rint() 
    {
        return Univariate_Derivative_1(std::rint(f0), 0.0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 sign() 
    {
        return Univariate_Derivative_1(FastMath.signum(f0), 0.0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 copy_sign(const Univariate_Derivative_1 sign) 
    {
        long m = Double.double_to_long_bits(f0);
        long s = Double.double_to_long_bits(sign.f0);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 copy_sign(const double sign) 
    {
        long m = Double.double_to_long_bits(f0);
        long s = Double.double_to_long_bits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inherit_doc} */
    //override
    public int get_exponent() 
    {
        return FastMath.get_exponent(f0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 scalb(const int& n) 
    {
        return Univariate_Derivative_1(std::scalbn(f0, n), std::scalbn(f1, n));
    }

    /** {@inherit_doc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    //override
    public Univariate_Derivative_1 ulp() 
    {
        return Univariate_Derivative_1(FastMath.ulp(f0), 0.0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 hypot(const Univariate_Derivative_1 y) 
    {

        if (std::isinf(f0) || std::isinf(y.f0)) 
        {
            return Univariate_Derivative_1(INFINITY, 0.0);
        }
else if (std::isnan(f0) || std::isnan(y.f0)) 
        {
            return Univariate_Derivative_1(Double.NaN, 0.0);
        }
else 
        {

            const int exp_x = get_exponent();
            const int exp_y = y.get_exponent();
            if (exp_x > exp_y + 27) 
            {
                // y is neglectible with respect to x
                return abs();
            }
else if (exp_y > exp_x + 27) 
            {
                // x is neglectible with respect to y
                return y.abs();
            }
else 
            {

                // find an intermediate scale to avoid both overflow and underflow
                const int middle_exp = (exp_x + exp_y) / 2;

                // scale parameters without losing precision
                const Univariate_Derivative_1 scaled_x = scalb(-middle_exp);
                const Univariate_Derivative_1 scaled_y = y.scalb(-middle_exp);

                // compute scaled hypotenuse
                const Univariate_Derivative_1 scaled_h =
                        scaled_x.multiply(scaled_x).add(scaled_y.multiply(scaled_y)).sqrt();

                // remove scaling
                return scaled_h.scalb(middle_exp);

            }

        }
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 reciprocal() 
    {
        const double inv1 = 1.0 / f0;
        const double inv2 = inv1 * inv1;
        return Univariate_Derivative_1(inv1, -f1 * inv2);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 compose(const double... f) 
    {
        Math_Utils::check_dimension(f.size(), get_order() + 1);
        return Univariate_Derivative_1(f[0], f[1] * f1);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 sqrt() 
    {
        const double s = std::sqrt(f0);
        return compose(s, 1 / (2 * s));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 cbrt() 
    {
        const double c = std::cbrt(f0);
        return compose(c, 1 / (3 * c * c));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 root_n(const int& n) 
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
            const double r = std::pow(f0, 1.0 / n);
            return compose(r, 1 / (n * std::pow(r, n - 1)));
        }
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1Field get_field() 
    {
        return Univariate_Derivative_1Field.get_instance();
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link Univariate_Derivative_1}
     * @param a number to exponentiate
     * @param x power to apply
     * @return a<sup>x</sup>
     */
    public static Univariate_Derivative_1 pow(const double& a, const Univariate_Derivative_1 x) 
    {
        if (a == 0) 
        {
            return x.get_field().get_zero();
        }
else 
        {
            const double& a_x = std::pow(a, x.f0);
            return Univariate_Derivative_1(a_x, std::log(a) * a_x * x.f1);
        }
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 pow(const double& p) 
    {
        if (p == 0) 
        {
            return get_field().get_one();
        }
else 
        {
            const double f0_pm1 = std::pow(f0, p - 1);
            return compose(f0_pm1 * f0, p * f0_pm1);
        }
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 pow(const int& n) 
    {
        if (n == 0) 
        {
            return get_field().get_one();
        }
else 
        {
            const double f0_nm1 = std::pow(f0, n - 1);
            return compose(f0_nm1 * f0, n * f0_nm1);
        }
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 pow(const Univariate_Derivative_1 e) 
    {
        return log().multiply(e).exp();
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 exp() 
    {
        const double exp = std::exp(f0);
        return compose(exp, exp);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 expm1() 
    {
        const double exp   = std::exp(f0);
        const double exp_m1 = std::expm1(f0);
        return compose(exp_m1, exp);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 log() 
    {
        return compose(std::log(f0), 1 / f0);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 log1p() 
    {
        return compose(std::log1p(f0), 1 / (1 + f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 log10() 
    {
        return compose(std::log10(f0), 1 / (f0 * std::log(10.0)));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 cos() 
    {
        const Sin_Cos sin_cos = Sin_Cos(f0);
        return compose(sin_cos.cos(), -sin_cos.sin());
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 sin() 
    {
        const Sin_Cos sin_cos = Sin_Cos(f0);
        return compose(sin_cos.sin(), sin_cos.cos());
    }

    /** {@inherit_doc} */
    //override
    public Field_Sin_Cos<Univariate_Derivative_1> sin_cos() 
    {
        const Sin_Cos sin_cos = Sin_Cos(f0);
        return Field_Sin_Cos<>(new Univariate_Derivative_1(sin_cos.sin(),  f1 * sin_cos.cos()), Univariate_Derivative_1(sin_cos.cos(), -f1 * sin_cos.sin()));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 tan() 
    {
        const double tan = std::tan(f0);
        return compose(tan, 1 + tan * tan);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 acos() 
    {
        return compose(std::acos(f0), -1 / std::sqrt(1 - f0 * f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 asin() 
    {
        return compose(std::asin(f0), 1 / std::sqrt(1 - f0 * f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 atan() 
    {
        return compose(std::atan(f0), 1 / (1 + f0 * f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 atan2(const Univariate_Derivative_1 x) 
    {
        const double inv = 1.0 / (f0 * f0 + x.f0 * x.f0);
        return Univariate_Derivative_1(std::atan2(f0, x.f0), Math_Arrays::linear_combination(x.f0, f1, -x.f1, f0) * inv);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 cosh() 
    {
        return compose(std::cosh(f0), std::sinh(f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 sinh() 
    {
        return compose(std::sinh(f0), std::cosh(f0));
    }

    /** {@inherit_doc} */
    //override
    public Field_Sinh_Cosh<Univariate_Derivative_1> sinh_cosh() 
    {
        const Sinh_Cosh sinh_cosh = std::sinh_cosh(f0);
        return Field_Sinh_Cosh<>(new Univariate_Derivative_1(sinh_cosh.sinh(), f1 * sinh_cosh.cosh()), Univariate_Derivative_1(sinh_cosh.cosh(), f1 * sinh_cosh.sinh()));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 tanh() 
    {
        const double tanh = std::tanh(f0);
        return compose(tanh, 1 - tanh * tanh);
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 acosh() 
    {
        return compose(std::acosh(f0), 1 / std::sqrt(f0 * f0 - 1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 asinh() 
    {
        return compose(std::asinh(f0), 1 / std::sqrt(f0 * f0 + 1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 atanh() 
    {
        return compose(std::atanh(f0), 1 / (1 - f0 * f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 to_degrees() 
    {
        return Univariate_Derivative_1(FastMath.to_degrees(f0), FastMath.to_degrees(f1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 to_radians() 
    {
        return Univariate_Derivative_1(FastMath.to_radians(f0), FastMath.to_radians(f1));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public double taylor(const double delta) 
    {
        return f0 + delta * f1;
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const Univariate_Derivative_1[] a, const Univariate_Derivative_1[] b) 
    {

        // extract values and first derivatives
        const int      n  = a.size();
        const std::vector<double> a0 = std::vector<double>(n];
        const std::vector<double> b0 = std::vector<double>(n];
        const std::vector<double> a1 = std::vector<double>(2 * n];
        const std::vector<double> b1 = std::vector<double>(2 * n];
        for (int i{}; i < n; ++i) 
        {
            const Univariate_Derivative_1 ai = a[i];
            const Univariate_Derivative_1 bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
        }

        return Univariate_Derivative_1(Math_Arrays::linear_combination(a0, b0), Math_Arrays::linear_combination(a1, b1));

    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const std::vector<double> a, const Univariate_Derivative_1[] b) 
    {

        // extract values and first derivatives
        const int      n  = b.size();
        const std::vector<double> b0 = std::vector<double>(n];
        const std::vector<double> b1 = std::vector<double>(n];
        for (int i{}; i < n; ++i) 
        {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
        }

        return Univariate_Derivative_1(Math_Arrays::linear_combination(a, b0), Math_Arrays::linear_combination(a, b1));

    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const Univariate_Derivative_1 a1, const Univariate_Derivative_1 b1, const Univariate_Derivative_1 a2, const Univariate_Derivative_1 b2) 
    {
        return Univariate_Derivative_1(Math_Arrays::linear_combination(a1.f0, b1.f0, a2.f0, b2.f0), Math_Arrays::linear_combination(a1.f0, b1.f1, a1.f1, b1.f0, a2.f0, b2.f1, a2.f1, b2.f0));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const double& a1, const Univariate_Derivative_1 b1, const double& a2, const Univariate_Derivative_1 b2) 
    {
        return Univariate_Derivative_1(Math_Arrays::linear_combination(a1, b1.f0, a2, b2.f0), Math_Arrays::linear_combination(a1, b1.f1, a2, b2.f1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const Univariate_Derivative_1 a1, const Univariate_Derivative_1 b1, const Univariate_Derivative_1 a2, const Univariate_Derivative_1 b2, const Univariate_Derivative_1 a3, const Univariate_Derivative_1 b3) 
    {
        return Univariate_Derivative_1(Math_Arrays::linear_combination(a1.f0, b1.f0, a2.f0, b2.f0, a3.f0, b3.f0), Math_Arrays::linear_combination(std::vector<double> 
        {
                                                                          a1.f0, a1.f1, a2.f0, a2.f1, a3.f0, a3.f1
                                                                      }, std::vector<double> 
                                                                      {
                                                                          b1.f1, b1.f0, b2.f1, b2.f0, b3.f1, b3.f0
                                                                      }));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const double& a1, const Univariate_Derivative_1 b1, const double& a2, const Univariate_Derivative_1 b2, const double& a3, const Univariate_Derivative_1 b3) 
    {
        return Univariate_Derivative_1(Math_Arrays::linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0), Math_Arrays::linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const Univariate_Derivative_1 a1, const Univariate_Derivative_1 b1, const Univariate_Derivative_1 a2, const Univariate_Derivative_1 b2, const Univariate_Derivative_1 a3, const Univariate_Derivative_1 b3, const Univariate_Derivative_1 a4, const Univariate_Derivative_1 b4) 
    {
        return Univariate_Derivative_1(Math_Arrays::linear_combination(a1.f0, b1.f0, a2.f0, b2.f0, a3.f0, b3.f0, a4.f0, b4.f0), Math_Arrays::linear_combination(std::vector<double> 
        {
                                                                          a1.f0, a1.f1, a2.f0, a2.f1, a3.f0, a3.f1, a4.f0, a4.f1
                                                                      }, std::vector<double> 
                                                                      {
                                                                          b1.f1, b1.f0, b2.f1, b2.f0, b3.f1, b3.f0, b4.f1, b4.f0
                                                                      }));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 linear_combination(const double& a1, const Univariate_Derivative_1 b1, const double& a2, const Univariate_Derivative_1 b2, const double& a3, const Univariate_Derivative_1 b3, const double& a4, const Univariate_Derivative_1 b4) 
    {
        return Univariate_Derivative_1(Math_Arrays::linear_combination(a1, b1.f0, a2, b2.f0, a3, b3.f0, a4, b4.f0), Math_Arrays::linear_combination(a1, b1.f1, a2, b2.f1, a3, b3.f1, a4, b4.f1));
    }

    /** {@inherit_doc} */
    //override
    public Univariate_Derivative_1 get_pi() 
    {
        return PI;
    }

    /** Test for the equality of two univariate derivatives.
     * <p>
     * univariate derivatives are considered equal if they have the same derivatives.
     * </p>
     * @param other Object to test for equality to this
     * @return true if two univariate derivatives are equal
     */
    //override
    public bool equals(Object other) 
    {

        if (this == other) 
        {
            return true;
        }

        if (dynamic_cast<const Univariate_Derivative_1*>(*other) != nullptr)
        {
            const Univariate_Derivative_1 rhs = (Univariate_Derivative_1) other;
            return f0 == rhs.f0 && f1 == rhs.f1;
        }

        return false;

    }

    /** Get a hash_code for the univariate derivative.
     * @return a hash code value for this object
     */
    //override
    public int hash_code() 
    {
        return 453 - 19 * Double.hash_code(f0) + 37 * Double.hash_code(f1);
    }

};