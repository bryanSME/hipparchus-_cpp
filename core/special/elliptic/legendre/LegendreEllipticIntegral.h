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
 //package org.hipparchus.special.elliptic.legendre;

 //import java.util.function.Double_Function;
 //import java.util.function.Function;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.analysis.Calculus_Field_Univariate_Function;
 //import org.hipparchus.complex.std::complex<double>;
 //import org.hipparchus.complex.std::complex<double>_Univariate_Integrator;
 //import org.hipparchus.complex.Field_Complex<double>;
 //import org.hipparchus.complex.Field_Complex<double>_Univariate_Integrator;
 //import org.hipparchus.special.elliptic.carlson.Carlson_Elliptic_Integral;
 //import org.hipparchus.util.FastMath;
 //import org.hipparchus.util.Math_Utils;
#include <vector>
#include <numbers>
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"
#include "../../../util/MathUtils.h"
#include "../../../complex/FieldComplex.h"
#include <complex>

/** Complete and incomplete elliptic integrals in Legendre form.
 * <p>
 * The elliptic integrals are related to Jacobi elliptic functions.
 * </p>
 * <p>
 * <emph>
 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
 * are considered experimental for now, they have known issues:
 * <a href="https://github.com/Hipparchus-Math/hipparchus/issues/151">issue 151</a>
 * and <a href="https://github.com/Hipparchus-Math/hipparchus/issues/152">issue 152</a>.
 * </emph>
 * </p>
 * <p>
 * There are different conventions to interpret the arguments of
 * Legendre elliptic integrals. In mathematical texts, these conventions show
 * up using the separator between arguments. So for example for the incomplete
 * integral of the first kind F we have:
 * <ul>
 *   <li>F(\xcf\x86, k): the first argument \xcf\x86 is an angle and the second argument k
 *       is the elliptic modulus: this is the trigonometric form of the integral</li>
 *   <li>F(\xcf\x86; m): the first argument \xcf\x86 is an angle and the second argument m=k\xc2\xb2
 *       is the parameter: this is also a trigonometric form of the integral</li>
 *   <li>F(x|m): the first argument x=sin(\xcf\x86) is not an angle anymore and the
 *       second argument m=k\xc2\xb2 is the parameter: this is the Legendre form</li>
 *   <li>F(\xcf\x86\\\xce\xb1): the first argument \xcf\x86 is an angle and the second argument \xce\xb1 is the
 *       modular angle</li>
 * </ul>
 * As we have no separator in a method call, we have to adopt one convention
 * and stick to it. In Hipparchus, we adopted the Legendre form (i.e. F(x|m), * with x=sin(\xcf\x86) and m=k\xc2\xb2. These conventions are consistent with Wolfram Alpha
 * functions Elliptic_F, Elliptic_E, ElliptiPI\xe2\x80\xa6
 * </p>
 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
 * @since 2.0
 */
class Legendre_Elliptic_Integral
{
private:
	/** Private constructor for a utility class.
	 */
	Legendre_Elliptic_Integral()
	{
		// nothing to do
	}

public:
	/** Get the nome q.
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return nome q
	 */
	static double nome(const double& m)
	{
		if (m < 1.0e-16)
		{
			// first terms of infinite series in Abramowitz and Stegun 17.3.21
			const double m16{ m * 0.0625 };
			return m16 * (1 + 8 * m16);
		}
		return std::exp(-std::numbers::pi * big_k_prime(m) / big_k(m));
	}

	/** Get the nome q.
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return nome q
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T nome(const T& m)
	{
		const T one = m.get_field().get_one();
		if (m.norm() < 100 * one.ulp().get_real())
		{
			// first terms of infinite series in Abramowitz and Stegun 17.3.21
			const T m16 = m.multiply(0.0625);
			return m16.multiply(m16.multiply(8).add(1));
		}
		return std::exp(big_k_prime(m).divide(big_k(m)).multiply(one.get_pi().negate()));
	}

	/** Get the complete elliptic integral of the first kind K(m).
	 * <p>
	 * The complete elliptic integral of the first kind K(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * it corresponds to the real quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the first kind K(m)
	 * @see #big_k_primestatic_cast<double>(
	 * @see #big_f(double, double)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	static double big_k(const double& m)
	{
		return (m < 1.0e-8)
			?  (1 + 0.25 * m) * Math_Utils::SEMI_PI // first terms of infinite series in Abramowitz and Stegun 17.3.11
			: Carlson_Elliptic_Integral::r_f(0, 1.0 - m, 1);
	}

	/** Get the complete elliptic integral of the first kind K(m).
	 * <p>
	 * The complete elliptic integral of the first kind K(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * it corresponds to the real quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the first kind K(m)
	 * @see #big_k_prime(Calculus_Field_Element)
	 * @see #big_f(Calculus_Field_Element, Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T& big_k(const T m)
	{
		const T zero = m.get_field().get_zero();
		const T one = m.get_field().get_one();
		return (m.norm() < 1.0e7 * one.ulp().get_real())
			? one.add(m.multiply(0.25)).multiply(zero.get_pi().multiply(0.5)) // first terms of infinite series in Abramowitz and Stegun 17.3.11
			: Carlson_Elliptic_Integral::r_f(zero, one.subtract(m), one);
	}

	/** Get the complete elliptic integral of the first kind K(m).
	 * <p>
	 * The complete elliptic integral of the first kind K(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * it corresponds to the real quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the first kind K(m)
	 * @see #big_k_prime(std::complex<double>)
	 * @see #big_f(std::complex<double>, std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	static std::complex<double> big_k(const std::complex<double>& m)
	{
		if (m.norm() < 1.0e-8)
		{
			// first terms of infinite series in Abramowitz and Stegun 17.3.11
			return std::complex<double>.ONE.add(m.multiply(0.25)).multiply(Math_Utils::SEMI_PI);
		}
		return Carlson_Elliptic_Integral::r_f(std::complex<double>.ZERO, std::complex<double>.ONE.subtract(m), std::complex<double>.ONE);
	}

	/** Get the complete elliptic integral of the first kind K(m).
	 * <p>
	 * The complete elliptic integral of the first kind K(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * it corresponds to the real quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the first kind K(m)
	 * @see #big_k_prime(Field_Complex<double>)
	 * @see #big_f(Field_Complex<double>, Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_k(const Field_Complex<T> m)
	{
		const Field_Complex<T> zero = m.get_field().get_zero();
		const Field_Complex<T> one = m.get_field().get_one();
		if (m.norm() < 1.0e7 * one.ulp().get_real())
		{
			// first terms of infinite series in Abramowitz and Stegun 17.3.11
			return one.add(m.multiply(0.25)).multiply(zero.get_pi().multiply(0.5));
		}
		else
		{
			return Carlson_Elliptic_Integral::r_f(zero, one.subtract(m), one);
		}
	}

	/** Get the complete elliptic integral of the first kind K'(m).
	 * <p>
	 * The complete elliptic integral of the first kind K'(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-(1-m) \\sin^2	heta}}
	 * \\]
	 * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the first kind K'(m)
	 * @see #big_kstatic_cast<double>(
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static double big_k_prime(const double m)
	{
		return Carlson_Elliptic_Integral::r_f(0, m, 1);
	}

	/** Get the complete elliptic integral of the first kind K'(m).
	 * <p>
	 * The complete elliptic integral of the first kind K'(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-(1-m) \\sin^2	heta}}
	 * \\]
	 * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the first kind K'(m)
	 * @see #big_k(Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_k_prime(const T& m)
	{
		const T zero = m.get_field().get_zero();
		const T one = m.get_field().get_one();
		return Carlson_Elliptic_Integral::r_f(zero, m, one);
	}

	/** Get the complete elliptic integral of the first kind K'(m).
	 * <p>
	 * The complete elliptic integral of the first kind K'(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-(1-m) \\sin^2	heta}}
	 * \\]
	 * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the first kind K'(m)
	 * @see #big_k(std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_k_prime(const std::complex<double> m)
	{
		return Carlson_Elliptic_Integral::r_f(std::complex<double>.ZERO, m, std::complex<double>.ONE);
	}

	/** Get the complete elliptic integral of the first kind K'(m).
	 * <p>
	 * The complete elliptic integral of the first kind K'(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-(1-m) \\sin^2	heta}}
	 * \\]
	 * it corresponds to the imaginary quarter-period of Jacobi elliptic functions
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the first kind K'(m)
	 * @see #big_k(Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_complete_elliptic_integralofthe_first_kind.html">Complete Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_k_prime(const Field_Complex<T> m)
	{
		const Field_Complex<T> zero = m.get_field().get_zero();
		const Field_Complex<T> one = m.get_field().get_one();
		return Carlson_Elliptic_Integral::r_f(zero, m, one);
	}

	/** Get the complete elliptic integral of the second kind E(m).
	 * <p>
	 * The complete elliptic integral of the second kind E(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the second kind E(m)
	 * @see #big_e(double, double)
	 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static double big_e(const double m)
	{
		return Carlson_Elliptic_Integral.rG(0, 1 - m, 1) * 2;
	}

	/** Get the complete elliptic integral of the second kind E(m).
	 * <p>
	 * The complete elliptic integral of the second kind E(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the second kind E(m)
	 * @see #big_e(Calculus_Field_Element, Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_e(const T m)
	{
		const T zero = m.get_field().get_zero();
		const T one = m.get_field().get_one();
		return Carlson_Elliptic_Integral.rG(zero, one.subtract(m), one).multiply(2);
	}

	/** Get the complete elliptic integral of the second kind E(m).
	 * <p>
	 * The complete elliptic integral of the second kind E(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the second kind E(m)
	 * @see #big_e(std::complex<double>, std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_e(const std::complex<double> m)
	{
		return Carlson_Elliptic_Integral.rG(std::complex<double>.ZERO, std::complex<double>.ONE.subtract(m), std::complex<double>.ONE).multiply(2);
	}

	/** Get the complete elliptic integral of the second kind E(m).
	 * <p>
	 * The complete elliptic integral of the second kind E(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the second kind E(m)
	 * @see #big_e(Field_Complex<double>, Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html">Complete Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_e(const Field_Complex<T> m)
	{
		const Field_Complex<T> zero = m.get_field().get_zero();
		const Field_Complex<T> one = m.get_field().get_one();
		return Carlson_Elliptic_Integral.rG(zero, one.subtract(m), one).multiply(2);
	}

	/** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
	 * <p>
	 * The complete elliptic integral D(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral D(m)
	 * @see #big_d(double, double)
	 */
	public static double big_d(const double m)
	{
		return Carlson_Elliptic_Integral.rD(0, 1 - m, 1) / 3;
	}

	/** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
	 * <p>
	 * The complete elliptic integral D(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral D(m)
	 * @see #big_d(Calculus_Field_Element, Calculus_Field_Element)
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_d(const T m)
	{
		const T zero = m.get_field().get_zero();
		const T one = m.get_field().get_one();
		return Carlson_Elliptic_Integral.rD(zero, one.subtract(m), one).divide(3);
	}

	/** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
	 * <p>
	 * The complete elliptic integral D(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral D(m)
	 * @see #big_d(std::complex<double>, std::complex<double>)
	 */
	public static std::complex<double> big_d(const std::complex<double> m)
	{
		return Carlson_Elliptic_Integral.rD(std::complex<double>.ZERO, std::complex<double>.ONE.subtract(m), std::complex<double>.ONE).divide(3);
	}

	/** Get the complete elliptic integral D(m) = [K(m) - E(m)]/m.
	 * <p>
	 * The complete elliptic integral D(m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral D(m)
	 * @see #big_d(Field_Complex<double>, Field_Complex<double>)
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_d(const Field_Complex<T> m)
	{
		const Field_Complex<T> zero = m.get_field().get_zero();
		const Field_Complex<T> one = m.get_field().get_one();
		return Carlson_Elliptic_Integral.rD(zero, one.subtract(m), one).divide(3);
	}

	/** Get the complete elliptic integral of the third kind \xce\xa0(n, m).
	 * <p>
	 * The complete elliptic integral of the third kind \xce\xa0(n, m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the third kind \xce\xa0(n, m)
	 * @see #big_pi(double, double, double)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static double big_pi(const double n, const double m)
	{
		const double k_prime2 = 1 - m;
		const double delta = n * (m - n) * (n - 1);
		return Carlson_Elliptic_Integral::r_f(0, k_prime2, 1) +
			Carlson_Elliptic_Integral.rJ(0, k_prime2, 1, 1 - n, delta) * n / 3;
	}

	/** Get the complete elliptic integral of the third kind \xce\xa0(n, m).
	 * <p>
	 * The complete elliptic integral of the third kind \xce\xa0(n, m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the third kind \xce\xa0(n, m)
	 * @see #big_pi(Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_pi(const T n, const T m)
	{
		const T zero = m.get_field().get_zero();
		const T one = m.get_field().get_one();
		const T k_prime2 = one.subtract(m);
		const T delta = n.multiply(m.subtract(n)).multiply(n.subtract(1));
		return Carlson_Elliptic_Integral::r_f(zero, k_prime2, one).
			add(Carlson_Elliptic_Integral.rJ(zero, k_prime2, one, one.subtract(n), delta).multiply(n).divide(3));
	}

	/** Get the complete elliptic integral of the third kind \xce\xa0(n, m).
	 * <p>
	 * The complete elliptic integral of the third kind \xce\xa0(n, m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return complete elliptic integral of the third kind \xce\xa0(n, m)
	 * @see #big_pi(std::complex<double>, std::complex<double>, std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_pi(const std::complex<double> n, const std::complex<double> m)
	{
		const std::complex<double> k_prime2 = std::complex<double>.ONE.subtract(m);
		const std::complex<double> delta = n.multiply(m.subtract(n)).multiply(n.subtract(1));
		return Carlson_Elliptic_Integral::r_f(std::complex<double>.ZERO, k_prime2, std::complex<double>.ONE).
			add(Carlson_Elliptic_Integral.rJ(std::complex<double>.ZERO, k_prime2, std::complex<double>.ONE, std::complex<double>.ONE.subtract(n), delta).multiply(n).divide(3));
	}

	/** Get the complete elliptic integral of the third kind \xce\xa0(n, m).
	 * <p>
	 * The complete elliptic integral of the third kind \xce\xa0(n, m) is
	 * \\[
	 *    \\int_0^{\x0crac{\\pi}{2}} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return complete elliptic integral of the third kind \xce\xa0(n, m)
	 * @see #big_pi(Field_Complex<double>, Field_Complex<double>, Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_pi(const Field_Complex<T> n, const Field_Complex<T> m)
	{
		const Field_Complex<T> zero = m.get_field().get_zero();
		const Field_Complex<T> one = m.get_field().get_one();
		const Field_Complex<T> k_prime2 = one.subtract(m);
		const Field_Complex<T> delta = n.multiply(m.subtract(n)).multiply(n.subtract(1));
		return Carlson_Elliptic_Integral::r_f(zero, k_prime2, one).
			add(Carlson_Elliptic_Integral.rJ(zero, k_prime2, one, one.subtract(n), delta).multiply(n).divide(3));
	}

	/** Get the incomplete elliptic integral of the first kind F(\xcf\x86, m).
	 * <p>
	 * The incomplete elliptic integral of the first kind F(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral of the first kind F(\xcf\x86, m)
	 * @see #big_kstatic_cast<double>(
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static double big_f(const double& phi, const double m)
	{
		// argument reduction
		const Double_Argument_Reduction ar = Double_Argument_Reduction(phi, m, n->big_k(n));

		// integrate part between 0 and \xcf\x80/2
		const double c_m1 = ar.csc2 - 1.0;
		const double c_mm = ar.csc2 - m;
		const double incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;
	}

	/** Get the incomplete elliptic integral of the first kind F(\xcf\x86, m).
	 * <p>
	 * The incomplete elliptic integral of the first kind F(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the first kind F(\xcf\x86, m)
	 * @see #big_k(Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_f(const T phi, const T m)
	{
		// argument reduction
		const Field_Argument_Reduction<T> ar = Field_Argument_Reduction<>(phi, m, n->big_k(n));

		// integrate part between 0 and \xcf\x80/2
		const T c_m1 = ar.csc2.subtract(1);
		const T c_mm = ar.csc2.subtract(m);
		const T incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the first kind F(\xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the first kind F(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral of the first kind F(\xcf\x86, m)
	 * @see #big_k(std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_f(const std::complex<double> phi, const std::complex<double> m)
	{
		// argument reduction
		const Field_Argument_Reduction<std::complex<double>> ar = Field_Argument_Reduction<>(phi, m, n->big_k(n));

		// integrate part between 0 and \xcf\x80/2
		const std::complex<double> c_m1 = ar.csc2.subtract(1);
		const std::complex<double> c_mm = ar.csc2.subtract(m);
		const std::complex<double> incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the first kind F(\xcf\x86, m) using numerical integration.
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the first kind F(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on numerical integration.
	 * If integration path comes too close to a pole of the integrand, then integration will fail
	 * with a {@link org.hipparchus.exception.Math_Illegal_State_Exception Math_Illegal_State_Exception}
	 * even for very large {@code max_eval}. This is normal behavior.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param integrator integrator to use
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately, so up to twice this number may be used)
	 * @return incomplete elliptic integral of the first kind F(\xcf\x86, m)
	 * @see #big_k(std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_f(const std::complex<double> phi, const std::complex<double> m, const std::complex<double>_Univariate_Integrator integrator, const int max_eval)
	{
		return integrator.integrate(max_eval, First<>(m), phi.get_field().get_zero(), phi);
	}

	/** Get the incomplete elliptic integral of the first kind F(\xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the first kind F(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the first kind F(\xcf\x86, m)
	 * @see #big_k(Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_f(const Field_Complex<T> phi, const Field_Complex<T> m)
	{
		// argument reduction
		const Field_Argument_Reduction<Field_Complex<T>> ar = Field_Argument_Reduction<>(phi, m, n->big_k(n));

		// integrate part between 0 and \xcf\x80/2
		const Field_Complex<T> c_m1 = ar.csc2.subtract(1);
		const Field_Complex<T> c_mm = ar.csc2.subtract(m);
		const Field_Complex<T> incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the first kind F(\xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the first kind F(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on numerical integration.
	 * If integration path comes too close to a pole of the integrand, then integration will fail
	 * with a {@link org.hipparchus.exception.Math_Illegal_State_Exception Math_Illegal_State_Exception}
	 * even for very large {@code max_eval}. This is normal behavior.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param integrator integrator to use
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately, so up to twice this number may be used)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the first kind F(\xcf\x86, m)
	 * @see #big_k(Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html">Elliptic Integrals of the First Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_f(const Field_Complex<T> phi, const Field_Complex<T> m, const Field_Complex<double>_Univariate_Integrator<T> integrator, const int max_eval)
	{
		return integrator.integrate(max_eval, First<>(m), phi.get_field().get_zero(), phi);
	}

	/** Get the incomplete elliptic integral of the second kind E(\xcf\x86, m).
	 * <p>
	 * The incomplete elliptic integral of the second kind E(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral of the second kind E(\xcf\x86, m)
	 * @see #big_estatic_cast<double>(
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static double big_e(const double& phi, const double m)
	{
		// argument reduction
		const Double_Argument_Reduction ar = Double_Argument_Reduction(phi, m, n->big_e(n));

		// integrate part between 0 and \xcf\x80/2
		const double c_m1 = ar.csc2 - 1.0;
		const double c_mm = ar.csc2 - m;
		const double incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2) -
			Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2) * (m / 3);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;
	}

	/** Get the incomplete elliptic integral of the second kind E(\xcf\x86, m).
	 * <p>
	 * The incomplete elliptic integral of the second kind E(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the second kind E(\xcf\x86, m)
	 * @see #big_e(Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_e(const T phi, const T m)
	{
		// argument reduction
		const Field_Argument_Reduction<T> ar = Field_Argument_Reduction<>(phi, m, n->big_e(n));

		// integrate part between 0 and \xcf\x80/2
		const T c_m1 = ar.csc2.subtract(1);
		const T c_mm = ar.csc2.subtract(m);
		const T incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2).
			subtract(Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2).multiply(m.divide(3)));

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the second kind E(\xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the second kind E(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral of the second kind E(\xcf\x86, m)
	 * @see #big_e(std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_e(const std::complex<double> phi, const std::complex<double> m)
	{
		// argument reduction
		const Field_Argument_Reduction<std::complex<double>> ar = Field_Argument_Reduction<>(phi, m, n->big_e(n));

		// integrate part between 0 and \xcf\x80/2
		const std::complex<double> c_m1 = ar.csc2.subtract(1);
		const std::complex<double> c_mm = ar.csc2.subtract(m);
		const std::complex<double> incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2).
			subtract(Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2).multiply(m.divide(3)));

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the second kind E(\xcf\x86, m) using numerical integration.
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the second kind E(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on numerical integration.
	 * If integration path comes too close to a pole of the integrand, then integration will fail
	 * with a {@link org.hipparchus.exception.Math_Illegal_State_Exception Math_Illegal_State_Exception}
	 * even for very large {@code max_eval}. This is normal behavior.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param integrator integrator to use
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately, so up to twice this number may be used)
	 * @return incomplete elliptic integral of the second kind E(\xcf\x86, m)
	 * @see #big_e(std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_e(const std::complex<double> phi, const std::complex<double> m, const std::complex<double>_Univariate_Integrator integrator, const int max_eval)
	{
		return integrator.integrate(max_eval, Second<>(m), phi.get_field().get_zero(), phi);
	}

	/** Get the incomplete elliptic integral of the second kind E(\xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the second kind E(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the second kind E(\xcf\x86, m)
	 * @see #big_e(Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_e(const Field_Complex<T> phi, const Field_Complex<T> m)
	{
		// argument reduction
		const Field_Argument_Reduction<Field_Complex<T>> ar = Field_Argument_Reduction<>(phi, m, n->big_e(n));

		// integrate part between 0 and \xcf\x80/2
		const Field_Complex<T> c_m1 = ar.csc2.subtract(1);
		const Field_Complex<T> c_mm = ar.csc2.subtract(m);
		const Field_Complex<T> incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2).
			subtract(Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2).multiply(m.divide(3)));

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the second kind E(\xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the second kind E(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \\sqrt{1-m \\sin^2	heta} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on numerical integration.
	 * If integration path comes too close to a pole of the integrand, then integration will fail
	 * with a {@link org.hipparchus.exception.Math_Illegal_State_Exception Math_Illegal_State_Exception}
	 * even for very large {@code max_eval}. This is normal behavior.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param integrator integrator to use
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately, so up to twice this number may be used)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the second kind E(\xcf\x86, m)
	 * @see #big_e(Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html">Elliptic Integrals of the Second Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_e(const Field_Complex<T> phi, const Field_Complex<T> m, const Field_Complex<double>_Univariate_Integrator<T> integrator, const int max_eval)
	{
		return integrator.integrate(max_eval, Second<>(m), phi.get_field().get_zero(), phi);
	}

	/** Get the incomplete elliptic integral D(\xcf\x86, m) = [F(\xcf\x86, m) - E(\xcf\x86, m)]/m.
	 * <p>
	 * The incomplete elliptic integral D(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral D(\xcf\x86, m)
	 * @see #big_dstatic_cast<double>(
	 */
	public static double big_d(const double& phi, const double m)
	{
		// argument reduction
		const Double_Argument_Reduction ar = Double_Argument_Reduction(phi, m, n->big_d(n));

		// integrate part between 0 and \xcf\x80/2
		const double c_m1 = ar.csc2 - 1.0;
		const double c_mm = ar.csc2 - m;
		const double incomplete = Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2) / 3;

		// combine complete and incomplete parts
		return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;
	}

	/** Get the incomplete elliptic integral D(\xcf\x86, m) = [F(\xcf\x86, m) - E(\xcf\x86, m)]/m.
	 * <p>
	 * The incomplete elliptic integral D(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral D(\xcf\x86, m)
	 * @see #big_d(Calculus_Field_Element)
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_d(const T phi, const T m)
	{
		// argument reduction
		const Field_Argument_Reduction<T> ar = Field_Argument_Reduction<>(phi, m, n->big_d(n));

		// integrate part between 0 and \xcf\x80/2
		const T c_m1 = ar.csc2.subtract(1);
		const T c_mm = ar.csc2.subtract(m);
		const T incomplete = Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2).divide(3);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral D(\xcf\x86, m) = [F(\xcf\x86, m) - E(\xcf\x86, m)]/m.
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral D(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral D(\xcf\x86, m)
	 * @see #big_d(std::complex<double>)
	 */
	public static std::complex<double> big_d(const std::complex<double> phi, const std::complex<double> m)
	{
		// argument reduction
		const Field_Argument_Reduction<std::complex<double>> ar = Field_Argument_Reduction<>(phi, m, n->big_d(n));

		// integrate part between 0 and \xcf\x80/2
		const std::complex<double> c_m1 = ar.csc2.subtract(1);
		const std::complex<double> c_mm = ar.csc2.subtract(m);
		const std::complex<double> incomplete = Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2).divide(3);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral D(\xcf\x86, m) = [F(\xcf\x86, m) - E(\xcf\x86, m)]/m.
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral D(\xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{\\sin^2	heta}{\\sqrt{1-m \\sin^2	heta}} d	heta
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral D(\xcf\x86, m)
	 * @see #big_d(Calculus_Field_Element)
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_d(const Field_Complex<T> phi, const Field_Complex<T> m)
	{
		// argument reduction
		const Field_Argument_Reduction<Field_Complex<T>> ar = Field_Argument_Reduction<>(phi, m, n->big_d(n));

		// integrate part between 0 and \xcf\x80/2
		const Field_Complex<T> c_m1 = ar.csc2.subtract(1);
		const Field_Complex<T> c_mm = ar.csc2.subtract(m);
		const Field_Complex<T> incomplete = Carlson_Elliptic_Integral.rD(c_m1, c_mm, ar.csc2).divide(3);

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m).
	 * <p>
	 * The incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m)
	 * @see #big_pi(double, double)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static double big_pi(const double n, const double& phi, const double m)
	{
		// argument reduction
		const Double_Argument_Reduction ar = Double_Argument_Reduction(phi, m, parameter->big_pi(n, parameter));

		// integrate part between 0 and \xcf\x80/2
		const double c_m1 = ar.csc2 - 1.0;
		const double c_mm = ar.csc2 - m;
		const double c_mn = ar.csc2 - n;
		const double delta = n * (m - n) * (n - 1);
		const double incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2) +
			Carlson_Elliptic_Integral.rJ(c_m1, c_mm, ar.csc2, c_mn, delta) * n / 3;

		// combine complete and incomplete parts
		return ar.negate ? ar.complete - incomplete : ar.complete + incomplete;
	}

	/** Get the incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m).
	 * <p>
	 * The incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m)
	 * @see #big_pi(Calculus_Field_Element, Calculus_Field_Element)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static T big_pi(const T n, const T phi, const T m)
	{
		// argument reduction
		const Field_Argument_Reduction<T> ar = Field_Argument_Reduction<>(phi, m, parameter->big_pi(n, parameter));

		// integrate part between 0 and \xcf\x80/2
		const T c_m1 = ar.csc2.subtract(1);
		const T c_mm = ar.csc2.subtract(m);
		const T c_mn = ar.csc2.subtract(n);
		const T delta = n.multiply(m.subtract(n)).multiply(n.subtract(1));
		const T incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2).
			add(Carlson_Elliptic_Integral.rJ(c_m1, c_mm, ar.csc2, c_mn, delta).multiply(n).divide(3));

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @return incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m)
	 * @see #big_pi(std::complex<double>, std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_pi(const std::complex<double> n, const std::complex<double> phi, const std::complex<double> m)
	{
		// argument reduction
		const Field_Argument_Reduction<std::complex<double>> ar = Field_Argument_Reduction<>(phi, m, parameter->big_pi(n, parameter));

		// integrate part between 0 and \xcf\x80/2
		const std::complex<double> c_m1 = ar.csc2.subtract(1);
		const std::complex<double> c_mm = ar.csc2.subtract(m);
		const std::complex<double> c_mn = ar.csc2.subtract(n);
		const std::complex<double> delta = n.multiply(m.subtract(n)).multiply(n.subtract(1));
		const std::complex<double> incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2).
			add(Carlson_Elliptic_Integral.rJ(c_m1, c_mm, ar.csc2, c_mn, delta).multiply(n).divide(3));

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) using numerical integration.
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on numerical integration.
	 * If integration path comes too close to a pole of the integrand, then integration will fail
	 * with a {@link org.hipparchus.exception.Math_Illegal_State_Exception Math_Illegal_State_Exception}
	 * even for very large {@code max_eval}. This is normal behavior.
	 * </p>
	 * @param n elliptic characteristic
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param integrator integrator to use
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * @return incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m)
	 * @see #big_pi(std::complex<double>, std::complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
	public static std::complex<double> big_pi(const std::complex<double> n, const std::complex<double> phi, const std::complex<double> m, const std::complex<double>_Univariate_Integrator integrator, const int max_eval)
	{
		return integrator.integrate(max_eval, Third<>(n, m), phi.get_field().get_zero(), phi);
	}

	/** Get the incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on {@link Carlson_Elliptic_Integral
	 * Carlson elliptic integrals}.
	 * </p>
	 * @param n elliptic characteristic
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m)
	 * @see #big_pi(Field_Complex<double>, Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static Field_Complex<T> big_pi(const Field_Complex<T> n, const Field_Complex<T> phi, const Field_Complex<T> m)
	{
		// argument reduction
		const Field_Argument_Reduction<Field_Complex<T>> ar = Field_Argument_Reduction<>(phi, m, parameter->big_pi(n, parameter));

		// integrate part between 0 and \xcf\x80/2
		const Field_Complex<T> c_m1 = ar.csc2.subtract(1);
		const Field_Complex<T> c_mm = ar.csc2.subtract(m);
		const Field_Complex<T> c_mn = ar.csc2.subtract(n);
		const Field_Complex<T> delta = n.multiply(m.subtract(n)).multiply(n.subtract(1));
		const Field_Complex<T> incomplete = Carlson_Elliptic_Integral::r_f(c_m1, c_mm, ar.csc2).
			add(Carlson_Elliptic_Integral.rJ(c_m1, c_mm, ar.csc2, c_mn, delta).multiply(n).divide(3));

		// combine complete and incomplete parts
		return ar.negate ? ar.complete.subtract(incomplete) : ar.complete.add(incomplete);
	}

	/** Get the incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m).
	 * <p>
	 * <emph>
	 * BEWARE! Elliptic integrals for complex numbers in the incomplete case
	 * are considered experimental for now, they have known issues.
	 * </emph>
	 * </p>
	 * <p>
	 * The incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m) is
	 * \\[
	 *    \\int_0^{\\phi} \x0crac{d	heta}{\\sqrt{1-m \\sin^2	heta}(1-n \\sin^2	heta)}
	 * \\]
	 * </p>
	 * <p>
	 * The algorithm for evaluating the functions is based on numerical integration.
	 * If integration path comes too close to a pole of the integrand, then integration will fail
	 * with a {@link org.hipparchus.exception.Math_Illegal_State_Exception Math_Illegal_State_Exception}
	 * even for very large {@code max_eval}. This is normal behavior.
	 * </p>
	 * @param n elliptic characteristic
	 * @param phi amplitude (i.e. upper bound of the integral)
	 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
	 * @param integrator integrator to use
	 * @param max_eval maximum number of evaluations (real and imaginary
	 * parts are evaluated separately, so up to twice this number may be used)
	 * @param <T> the type of the field elements
	 * @return incomplete elliptic integral of the third kind \xce\xa0(n, \xcf\x86, m)
	 * @see #big_pi(Field_Complex<double>, Field_Complex<double>)
	 * @see <a href="https://mathworld.wolfram.com/_elliptic_integralofthe_third_kind.html">Elliptic Integrals of the Third Kind (MathWorld)</a>
	 * @see <a href="https://en.wikipedia.org/wiki/Elliptic_integral">Elliptic Integrals (Wikipedia)</a>
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	static  Field_Complex<T> big_pi(const Field_Complex<T> n, const Field_Complex<T> phi, const Field_Complex<T> m, const Field_Complex<double>_Univariate_Integrator<T> integrator, const int max_eval)
	{
		return integrator.integrate(max_eval, Third<>(n, m), phi.get_field().get_zero(), phi);
	}

	/** Argument reduction for an incomplete integral. */
	static class Double_Argument_Reduction
	{
	private:
		/** Complete part. */
		double my_complete;

		/** Squared cosecant of the Jacobi amplitude. */
		double my_csc2;

		/** Indicator for negated Jacobi amplitude. */
		bool my_negate;
	public:

		/** Simple constructor.
		 * @param phi amplitude (i.e. upper bound of the integral)
		 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
		 * @param integral provider for complete integral
		 */
		Double_Argument_Reduction(const double& phi, const double& m, const Double_Function<Double>& integral)
		{
			const double sin = std::sin(phi);
			const int    p = static_cast<int>(std::rint(phi / std::numbers::pi);
			my_complete = p == 0
				? 0
				: integral.apply(m) * 2 * p;
			my_negate = sin < 0 ^ (p & 0x1) == 1;
			my_csc2 = 1.0 / (sin * sin);
		}
	}

	/** Argument reduction for an incomplete integral.
	 * @param <T> type fo the field elements
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	class Field_Argument_Reduction
	{
	private:
		/** Complete part. */
		T my_complete;

		/** Squared cosecant of the Jacobi amplitude. */
		T my_csc2;

		/** Indicator for negated Jacobi amplitude. */
		bool my_negate;

	public:
		/** Simple constructor.
		 * @param phi amplitude (i.e. upper bound of the integral)
		 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
		 * @param integral provider for complete integral
		 */
		Field_Argument_Reduction(const T& phi, const T& m, const Function<T, T>& integral)
		{
			const T sin = std::sin(phi);
			const int p = static_cast<int>(std::rint(phi.get_real() / std::numbers::pi);
			my_complete = p == 0
				? phi.get_field().get_zero()
				: integral.apply(m).multiply(2 * p);
			my_negate = sin.get_real() < 0 ^ (p & 0x1) == 1;
			my_csc2 = sin.multiply(sin).reciprocal();
		}
	}

	/** Integrand for elliptic integrals of the first kind.
	 * @param <T> type of the field elements
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	class First : public Calculus_Field_Univariate_Function<T>
	{
	private:
		/** Parameter. */
		const T my_m;

	public:
		/** Simple constructor.
		 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
		 */
		First(const T m) : my_m{ m } {};

		/** {@inherit_doc} */
		//override
		T value(const T& theta)
		{
			const T sin = theta.sin();
			const T sin2 = sin.multiply(sin);
			return sin2.multiply(my_m).negate().add(1).sqrt().reciprocal();
		}
	}

	/** Integrand for elliptic integrals of the second kind.
	 * @param <T> type of the field elements
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	class Second : public Calculus_Field_Univariate_Function<T>
	{
		/** Parameter. */
		const T my_m;

	public:
		/** Simple constructor.
		 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
		 */
		Second(const T& m) : my_m{ m } {};

		/** {@inherit_doc} */
		//override
		T value(const T& theta)
		{
			const T sin = theta.sin();
			const T sin2 = sin.multiply(sin);
			return sin2.multiply(my_m).negate().add(1).sqrt();
		}
	}

	/** Integrand for elliptic integrals of the third kind.
	 * @param <T> type of the field elements
	 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
	class Third : Calculus_Field_Univariate_Function<T>
	{
		/** Elliptic characteristic. */
		const T my_n;

		/** Parameter. */
		const T my_m;

	public:
		/** Simple constructor.
		 * @param n elliptic characteristic
		 * @param m parameter (m=k\xc2\xb2 where k is the elliptic modulus)
		 */
		Third(const T& n, const T& m) : my_n{ n }, my_n{ n } {};

		/** {@inherit_doc} */
		//override
		T value(const T& theta)
		{
			const T sin = theta.sin();
			const T sin2 = sin.multiply(sin);
			const T d1 = sin2.multiply(my_m).negate().add(1).sqrt();
			const T da = sin2.multiply(my_n).negate().add(1);
			return d1.multiply(da).reciprocal();
		}
	}
};