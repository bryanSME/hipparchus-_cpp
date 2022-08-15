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

#include <complex>
#include "RcRealDuplication.h"
#include "RdRealDuplication.h"
#include "RfRealDuplication.h"
#include "RcFieldDuplication.hpp"
#include "RfFieldDuplication.hpp"
#include "RjRealDuplication.h"
#include "RjFieldDuplication.hpp"
#include <type_traits>
#include "../CalculusFieldElement.hpp"

/** Elliptic integrals in Carlson symmetric form.
 * <p>
 * This utility class computes the various symmetric elliptic
 * integrals defined as:
 * \[
 *   \left\{\begin{align}
 *   R_F(x,y,z)   &= \frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)}\\
 *   R_J(x,y,z,p) &= \frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)(t+p)}\\
 *   R_G(x,y,z)   &= \frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
                     \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t\\
 *   R_D(x,y,z)   &= R_J(x,y,z,z)\\
 *   R_C(x,y)     &= R_F(x,y,y)
 *   \end{align}\right.
 * \]
 * </p>
 * <p>
 * where
 * \[
 *   s(t) = \sqrt{t+x}\sqrt{t+y}\sqrt{t+z}
 * \]
 * </p>
 * <p>
 * The algorithms used are based on the duplication method as described in
 * B. C. Carlson 1995 paper "Numerical computation of real or complex
 * elliptic integrals", with the improvements described in the appendix
 * of B. C. Carlson and James FitzSimons 2000 paper "Reduction theorems
 * for elliptic integrands with the square root of two quadratic factors".
 * They are also described in <a href="https://dlmf.nist.gov/19.36#i">section 19.36(i)</a>
 * of Digital Library of Mathematical Functions.
 * </p>
 * @since 2.0
 */
class Carlson_Elliptic_Integral 
{

private:
    /** Private constructor for a utility class.
     */
    Carlson_Elliptic_Integral() = default;

    /** Compute Carlson elliptic integral R<sub>G</sub> in the general case.
 * @param x first symmetric variable of the integral
 * @param y second symmetric variable of the integral
 * @param z third symmetric variable of the integral
 * @return Carlson elliptic integral R<sub>G</sub>
 */
    static double general_compute_rg(const double& x, const double& y, const double& z)
    {
        // permute parameters if needed to avoid cancellations
        if (x <= y)
        {
            if (y <= z)
            {
                // x ≤ y ≤ z
                return permuted_compute_rg(x, z, y);
            }
            if (x <= z)
            {
                // x ≤ z < y
                return permuted_compute_rg(x, y, z);
            }
            // z < x ≤ y
            return permuted_compute_rg(z, y, x);
        }
        if (x <= z)
        {
            // y < x ≤ z
            return permuted_compute_rg(y, z, x);
        }
        if (y <= z)
        {
            // y ≤ z < x
            return permuted_compute_rg(y, x, z);
        }
        // z < y < x
        return permuted_compute_rg(z, x, y);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> in the general case.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T general_compute_rg(const T& x, const T y, const T z)
    {
        // permute parameters if needed to avoid cancellations
        const double x_r = x.get_real();
        const double y_r = y.get_real();
        const double zR = z.get_real();
        if (x_r <= y_r)
        {
            if (y_r <= zR)
            {
                // x ≤ y ≤ z
                return permuted_compute_rg(x, z, y);
            }
            if (x_r <= zR)
            {
                // x ≤ z < y
                return permuted_compute_rg(x, y, z);
            }
            // z < x ≤ y
            return permuted_compute_rg(z, y, x);
        }
        if (x_r <= zR)
        {
            // y < x ≤ z
            return permuted_compute_rg(y, z, x);
        }
        if (y_r <= zR)
        {
            // y ≤ z < x
            return permuted_compute_rg(y, x, z);
        }
        // z < y < x
        return permuted_compute_rg(z, x, y);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with already permuted variables to avoid cancellations.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    static double permuted_compute_rg(const double& x, const double& y, const double& z)
    {
        // permute parameters if needed to avoid divisions by zero
        if (z == 0)
        {
            return x == 0
                ? safe_compute_rg(z, x, y)
                : safe_compute_rg(y, z, x);
        }
        return safe_compute_rg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with already permuted variables to avoid cancellations.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T permuted_compute_rg(const T& x, const T y, const T z)
    {
        // permute parameters if needed to avoid divisions by zero
        if (z.is_zero())
        {
            return x.is_zero()
                ? safe_compute_rg(z, x, y)
                : safe_compute_rg(y, z, x);
        }
        return safe_compute_rg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with non-zero third variable.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @see <a href="https://dlmf.nist.gov/19.21#E10">Digital Library of Mathematical Functions, equation 19.21.10</a>
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    static double safe_compute_rg(const double& x, const double& y, const double& z)
    {

        // contribution of the R_F integral
        const double term_f = Rf_Real_Duplication(x, y, z).integral() * z;

        // contribution of the R_D integral
        const double term_d = (x - z) * (y - z) * Rd_Real_Duplication(x, y, z).integral() / 3;

        // contribution of the square roots
        const double term_s = std::sqrt(x * y / z);

        // equation 19.21.10
        return (term_f - term_d + term_s) * 0.5;

    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with non-zero third variable.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link std::complex<double>} or {@link Field_Complex<double>})
     * @see <a href="https://dlmf.nist.gov/19.21#E10">Digital Library of Mathematical Functions, equation 19.21.10</a>
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T safe_compute_rg(const T& x, const T& y, const T& z)
    {
        // contribution of the R_F integral
        const T term_f = Rf_Field_Duplication<>(x, y, z).integral().multiply(z);

        // contribution of the R_D integral
        const T term_d = x.subtract(z).multiply(y.subtract(z)).multiply(new Rd_Field_Duplication<>(x, y, z).integral()).divide(3);

        // contribution of the square roots
        // BEWARE: this term MUST be computed as √x√y/√z with all square roots selected with positive real part
        // and NOT as √(xy/z), otherwise sign errors may occur
        const T term_s = x.sqrt().multiply(y.sqrt()).divide(z.sqrt());

        // equation 19.21.10
        return term_f.subtract(term_d).add(term_s).multiply(0.5);

    }

public:
    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    static double rC(const double& x, const double& y) 
    {
        if (y < 0) 
        {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            const double x_my = x - y;
            return std::sqrt(x / x_my) * Rc_Real_Duplication(x_my, -y).integral();
        }

        return Rc_Real_Duplication(x, y).integral();
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T rC(const T& x, const T& y) 
    {
        if (y.get_real() < 0) 
        {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            const T x_my = x.subtract(y);
            return std::sqrt(x.divide(x_my)).multiply(new Rc_Field_Duplication<>(x_my, y.negate()).integral());
        }
        return Rc_Field_Duplication<>(x, y).integral();
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    static std::complex<double> rC(const std::complex<double> x, const std::complex<double> y) 
    {
        if (y.get_imaginary_part() == 0 && y.get_real_part() < 0) 
        {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            const std::complex<double> x_my = x.subtract(y);
            return std::sqrt(x.divide(x_my)).multiply(new Rc_Field_Duplication<>(x_my, y.negate()).integral());
        }
        return Rc_Field_Duplication<>(x, y).integral();
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Complex<T> rC(const Field_Complex<T>& x, const Field_Complex<T>& y) 
    {
        if (y.get_imaginary_part().is_zero() && y.get_real_part().get_real() < 0) 
        {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            const Field_Complex<T> x_my = x.subtract(y);
            return std::sqrt(x.divide(x_my)).multiply(new Rc_Field_Duplication<>(x_my, y.negate()).integral());
        }
        return Rc_Field_Duplication<>(x, y).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    static double r_f(const double& x, const double& y, const double& z) 
    {
        return Rf_Real_Duplication(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T r_f(const T& x, const T y, const T z) 
    {
        return Rf_Field_Duplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    static std::complex<double> r_f(const std::complex<double>& x, const std::complex<double>& y, const std::complex<double>& z) 
    {
        return Rf_Field_Duplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Complex<T> r_f(const Field_Complex<T> x, const Field_Complex<T> y, const Field_Complex<T> z) 
    {
        return Rf_Field_Duplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    static double rJ(const double& x, const double& y, const double z, const double p) 
    {
        const double delta = (p - x) * (p - y) * (p - z);
        return rJ(x, y, z, p, delta);
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    static double rJ(const double& x, const double& y, const double& z, const double& p, const double& delta) 
    {
        return Rj_Real_Duplication(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T rJ(const T& x, const T y, const T z, const T p) 
    {
        const T delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        return Rj_Field_Duplication<>(x, y, z, p, delta). integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T rJ(const T& x, const T y, const T z, const T p, const T delta) 
    {
        return Rj_Field_Duplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    static std::complex<double> rJ(const std::complex<double>& x, const std::complex<double>& y, const std::complex<double>& z, const std::complex<double> p) 
    {
        const std::complex<double> delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        return Rj_Field_Duplication(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    static std::complex<double> rJ(const std::complex<double> x, const std::complex<double> y, const std::complex<double> z, const std::complex<double> p, const std::complex<double> delta) 
    {
        return Rj_Field_Duplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Complex<T> rJ(const Field_Complex<T>& x, const Field_Complex<T>& y, const Field_Complex<T>& z, const Field_Complex<T>& p) 
    {
        const Field_Complex<T> delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        return Rj_Field_Duplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Complex<T> rJ(const Field_Complex<T>& x, const Field_Complex<T>& y, const Field_Complex<T>& z, const Field_Complex<T>& p, const Field_Complex<T> delta) 
    {
        return Rj_Field_Duplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    static double rD(const double& x, const double& y, const double& z) 
    {
        return Rd_Real_Duplication(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T rD(const T& x, const T y, const T z) 
    {
        return Rd_Field_Duplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    static std::complex<double> rD(const std::complex<double>& x, const std::complex<double>& y, const std::complex<double>& z) 
    {
        return Rd_Field_Duplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Complex<T> rD(const Field_Complex<T>& x, const Field_Complex<T>& y, const Field_Complex<T>& z) 
    {
        return Rd_Field_Duplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    static double rG(const double& x, const double& y, const double& z) 
    {
        return general_compute_rg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static T rG(const T& x, const T& y, const T& z) 
    {
        return general_compute_rg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    static std::complex<double> rG(const std::complex<double>& x, const std::complex<double>& y, const std::complex<double>& z) 
    {
        return general_compute_rg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static Field_Complex<T> rG(const Field_Complex<T>& x, const Field_Complex<T>& y, const Field_Complex<T>& z) 
    {
        return general_compute_rg(x, y, z);
    }
};