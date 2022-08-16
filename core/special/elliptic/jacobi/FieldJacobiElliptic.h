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
 //package org.hipparchus.special.elliptic.jacobi;

 //import org.hipparchus.Calculus_Field_Element;
 //import org.hipparchus.special.elliptic.carlson.Carlson_Elliptic_Integral;
 //import org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral;
 //import org.hipparchus.util.FastMath;
#include <type_traits>
#include "../../../CalculusFieldElement.hpp"

/** Computation of Jacobi elliptic functions.
 * The Jacobi elliptic functions are related to elliptic integrals.
 * @param <T> the type of the field elements
 * @since 2.0
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Jacobi_Elliptic
{
private:
	/** Parameter of the function. */
	const T my_m;

	/** Evaluate inverse of Jacobi elliptic function ps.
	 * <p>
	 * Here p, q, r are any permutation of the letters c, d, n.
	 * </p>
	 * @param x value of Jacobi elliptic function {@code ps(u|m)}
	 * @param deltaQP \xce\x94\xe2\x81\xa1(q, p) = q\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - p\xe2\x81\xa3s\xc2\xb2(u|m) (equation 19.5.28 of DLMF)
	 * @param deltaRP \xce\x94\xe2\x81\xa1(r, p) = r\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - p\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) (equation 19.5.28 of DLMF)
	 * @return u such that {@code x=ps(u|m)}
	 * @since 2.1
	 */
	T arcps(const T& x, const T& deltaQP, const T& deltaRP)
	{
		// see equation 19.25.32 in Digital Library of Mathematical Functions
		// https://dlmf.nist.gov/19.25.E32
		const T x2 = x.multiply(x);
		const T rf = Carlson_Elliptic_Integral::r_f(x2, x2.add(deltaQP), x2.add(deltaRP));
		return std::copysign(1.0, rf.get_real()) * std::copysign(1.0, x.get_real()) < 0 ?
			rf.negate() : rf;
	}

	/** Evaluate inverse of Jacobi elliptic function sp.
	 * <p>
	 * Here p, q, r are any permutation of the letters c, d, n.
	 * </p>
	 * @param x value of Jacobi elliptic function {@code sp(u|m)}
	 * @param deltaQP \xce\x94\xe2\x81\xa1(q, p) = q\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - p\xe2\x81\xa3s\xc2\xb2(u|m) (equation 19.5.28 of DLMF)
	 * @param deltaRP \xce\x94\xe2\x81\xa1(r, p) = r\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - p\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) (equation 19.5.28 of DLMF)
	 * @return u such that {@code x=sp(u|m)}
	 * @since 2.1
	 */
	T arcsp(const T& x, const T& deltaQP, const T& deltaRP)
	{
		// see equation 19.25.33 in Digital Library of Mathematical Functions
		// https://dlmf.nist.gov/19.25.E33
		const T x2 = x.multiply(x);
		return x.multiply(Carlson_Elliptic_Integral::r_f(x.get_field().get_one(), deltaQP.multiply(x2).add(1), deltaRP.multiply(x2).add(1)));
	}

	/** Evaluate inverse of Jacobi elliptic function pq.
	 * <p>
	 * Here p, q, r are any permutation of the letters c, d, n.
	 * </p>
	 * @param x value of Jacobi elliptic function {@code pq(u|m)}
	 * @param deltaQP \xce\x94\xe2\x81\xa1(q, p) = q\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - p\xe2\x81\xa3s\xc2\xb2(u|m) (equation 19.5.28 of DLMF)
	 * @param deltaRQ \xce\x94\xe2\x81\xa1(r, q) = r\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - q\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) (equation 19.5.28 of DLMF)
	 * @return u such that {@code x=pq(u|m)}
	 * @since 2.1
	 */
	T arcpq(const T& x, const T& deltaQP, const T& deltaRQ)
	{
		// see equation 19.25.34 in Digital Library of Mathematical Functions
		// https://dlmf.nist.gov/19.25.E34
		const T x2 = x.multiply(x);
		const T w = x2.subtract(1).negate().divide(deltaQP);
		const T rf = Carlson_Elliptic_Integral::r_f(x2, x.get_field().get_one(), deltaRQ.multiply(w).add(1));
		const T positive = w.sqrt().multiply(rf);
		return x.get_real() < 0 ? Legendre_Elliptic_Integral.big_k(get_m()).multiply(2).subtract(positive) : positive;
	}

	/** Evaluate inverse of Jacobi elliptic function pq.
	 * <p>
	 * Here p, q, r are any permutation of the letters c, d, n.
	 * </p>
	 * <p>
	 * This computed the same thing as {@link #arcpq(Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element)}
	 * but uses the homogeneity property Rf(x, y, z) = Rf(ax, ay, az) / \xe2\x88\x9aa to get rid of the division
	 * by deltaRQ. This division induces problems in the complex case as it may lose the sign
	 * of zero for values exactly along the real or imaginary axis, hence perturbing branch cuts.
	 * </p>
	 * @param x value of Jacobi elliptic function {@code pq(u|m)}
	 * @param deltaQP \xce\x94\xe2\x81\xa1(q, p) = q\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - p\xe2\x81\xa3s\xc2\xb2(u|m) (equation 19.5.28 of DLMF)
	 * @param deltaRQ \xce\x94\xe2\x81\xa1(r, q) = r\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) - q\xe2\x81\xa3s\xc2\xb2\xe2\x81\xa1(u|m) (equation 19.5.28 of DLMF)
	 * @return u such that {@code x=pq(u|m)}
	 * @since 2.1
	 */
	T arcpqNoDivision(const T& x, const T& deltaQP, const T& deltaRQ)
	{
		// see equation 19.25.34 in Digital Library of Mathematical Functions
		// https://dlmf.nist.gov/19.25.E34
		const T x2 = x.multiply(x);
		const T wDeltaQP = x2.subtract(1).negate();
		const T rf = Carlson_Elliptic_Integral::r_f(x2.multiply(deltaQP), deltaQP, deltaRQ.multiply(wDeltaQP).add(deltaQP));
		const T positive = wDeltaQP.sqrt().multiply(rf);
		return std::copysign(1.0, x.get_real()) < 0 ?
			Legendre_Elliptic_Integral.big_k(get_m()).multiply(2).subtract(positive) :
			positive;
	}

protected:
	/** Simple constructor.
	 * @param m parameter of the function
	 */
	Field_Jacobi_Elliptic(const T& m) : my_m{ m } {};

public:
	/** Get the parameter of the function.
	 * @return parameter of the function
	 */
	T get_m()
	{
		return m;
	}

	/** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three principal Jacobi
	 * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
	 */
	virtual Field_Copolar_N<T> values_n(T u);

	/** Evaluate the three principal Jacobi elliptic functions with pole at point n in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three principal Jacobi
	 * elliptic functions {@code sn(u|m)}, {@code cn(u|m)}, and {@code dn(u|m)}.
	 */
	Field_Copolar_N<T> values_n(const double u)
	{
		return values_n(m.new_instance(u));
	}

	/** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three subsidiary Jacobi
	 * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
	 */
	FieldCopolar_S<T> values_s(const T u)
	{
		return FieldCopolar_S<>(values_n(u));
	}

	/** Evaluate the three subsidiary Jacobi elliptic functions with pole at point s in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three subsidiary Jacobi
	 * elliptic functions {@code cs(u|m)}, {@code ds(u|m)} and {@code ns(u|m)}.
	 */
	FieldCopolar_S<T> values_s(const double u)
	{
		return FieldCopolar_S<>(values_n(u));
	}

	/** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three subsidiary Jacobi
	 * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
	 */
	FieldCopolar_C<T> values_c(const T u)
	{
		return FieldCopolar_C<>(values_n(u));
	}

	/** Evaluate the three subsidiary Jacobi elliptic functions with pole at point c in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three subsidiary Jacobi
	 * elliptic functions {@code dc(u|m)}, {@code nc(u|m)}, and {@code sc(u|m)}.
	 */
	FieldCopolar_C<T> values_c(const double u)
	{
		return FieldCopolar_C<>(values_n(u));
	}

	/** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three subsidiary Jacobi
	 * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
	 */
	FieldCopolar_D<T> values_d(const T u)
	{
		return FieldCopolar_D<>(values_n(u));
	}

	/** Evaluate the three subsidiary Jacobi elliptic functions with pole at point d in Glaisher\xe2\x80\x99s Notation.
	 * @param u argument of the functions
	 * @return copolar trio containing the three subsidiary Jacobi
	 * elliptic functions {@code nd(u|m)}, {@code sd(u|m)}, and {@code cd(u|m)}.
	 */
	FieldCopolar_D<T> values_d(const double u)
	{
		return FieldCopolar_D<>(values_n(u));
	}

	/** Evaluate inverse of Jacobi elliptic function sn.
	 * @param x value of Jacobi elliptic function {@code sn(u|m)}
	 * @return u such that {@code x=sn(u|m)}
	 * @since 2.1
	 */
	T arcsn(const T& x)
	{
		// p = n, q = c, r = d, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, p)
		return arcsp(x, x.get_field().get_one().negate(), get_m().negate());
	}

	/** Evaluate inverse of Jacobi elliptic function sn.
	 * @param x value of Jacobi elliptic function {@code sn(u|m)}
	 * @return u such that {@code x=sn(u|m)}
	 * @since 2.1
	 */
	T arcsn(const double& x)
	{
		return arcsn(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function cn.
	 * @param x value of Jacobi elliptic function {@code cn(u|m)}
	 * @return u such that {@code x=cn(u|m)}
	 * @since 2.1
	 */
	T arccn(const T& x)
	{
		// p = c, q = n, r = d, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, q)
		return arcpqNoDivision(x, x.get_field().get_one(), get_m().negate());
	}

	/** Evaluate inverse of Jacobi elliptic function cn.
	 * @param x value of Jacobi elliptic function {@code cn(u|m)}
	 * @return u such that {@code x=cn(u|m)}
	 * @since 2.1
	 */
	T arccn(const double& x)
	{
		return arccn(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function dn.
	 * @param x value of Jacobi elliptic function {@code dn(u|m)}
	 * @return u such that {@code x=dn(u|m)}
	 * @since 2.1
	 */
	T arcdn(const T& x)
	{
		// p = d, q = n, r = c, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, q)
		return arcpqNoDivision(x, get_m(), x.get_field().get_one().negate());
	}

	/** Evaluate inverse of Jacobi elliptic function dn.
	 * @param x value of Jacobi elliptic function {@code dn(u|m)}
	 * @return u such that {@code x=dn(u|m)}
	 * @since 2.1
	 */
	T arcdn(const double& x)
	{
		return arcdn(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function cs.
	 * @param x value of Jacobi elliptic function {@code cs(u|m)}
	 * @return u such that {@code x=cs(u|m)}
	 * @since 2.1
	 */
	T arccs(const T& x)
	{
		// p = c, q = n, r = d, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, p)
		return arcps(x, x.get_field().get_one(), get_m().subtract(1).negate());
	}

	/** Evaluate inverse of Jacobi elliptic function cs.
	 * @param x value of Jacobi elliptic function {@code cs(u|m)}
	 * @return u such that {@code x=cs(u|m)}
	 * @since 2.1
	 */
	T arccs(const double& x)
	{
		return arccs(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function ds.
	 * @param x value of Jacobi elliptic function {@code ds(u|m)}
	 * @return u such that {@code x=ds(u|m)}
	 * @since 2.1
	 */
	T arcds(const T& x)
	{
		// p = d, q = c, r = n, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, p)
		return arcps(x, get_m().subtract(1), get_m());
	}

	/** Evaluate inverse of Jacobi elliptic function ds.
	 * @param x value of Jacobi elliptic function {@code ds(u|m)}
	 * @return u such that {@code x=ds(u|m)}
	 * @since 2.1
	 */
	T arcds(const double& x)
	{
		return arcds(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function ns.
	 * @param x value of Jacobi elliptic function {@code ns(u|m)}
	 * @return u such that {@code x=ns(u|m)}
	 * @since 2.1
	 */
	T arcns(const T& x)
	{
		// p = n, q = c, r = d, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, p)
		return arcps(x, x.get_field().get_one().negate(), get_m().negate());
	}

	/** Evaluate inverse of Jacobi elliptic function ns.
	 * @param x value of Jacobi elliptic function {@code ns(u|m)}
	 * @return u such that {@code x=ns(u|m)}
	 * @since 2.1
	 */
	T arcns(const double& x)
	{
		return arcns(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function dc.
	 * @param x value of Jacobi elliptic function {@code dc(u|m)}
	 * @return u such that {@code x=dc(u|m)}
	 * @since 2.1
	 */
	T arcdc(const T& x)
	{
		// p = d, q = c, r = n, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, q)
		return arcpq(x, get_m().subtract(1), x.get_field().get_one());
	}

	/** Evaluate inverse of Jacobi elliptic function dc.
	 * @param x value of Jacobi elliptic function {@code dc(u|m)}
	 * @return u such that {@code x=dc(u|m)}
	 * @since 2.1
	 */
	T arcdc(const double& x)
	{
		return arcdc(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function nc.
	 * @param x value of Jacobi elliptic function {@code nc(u|m)}
	 * @return u such that {@code x=nc(u|m)}
	 * @since 2.1
	 */
	T arcnc(const T& x)
	{
		// p = n, q = c, r = d, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, q)
		return arcpq(x, x.get_field().get_one().negate(), get_m().subtract(1).negate());
	}

	/** Evaluate inverse of Jacobi elliptic function nc.
	 * @param x value of Jacobi elliptic function {@code nc(u|m)}
	 * @return u such that {@code x=nc(u|m)}
	 * @since 2.1
	 */
	T arcnc(const double& x)
	{
		return arcnc(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function sc.
	 * @param x value of Jacobi elliptic function {@code sc(u|m)}
	 * @return u such that {@code x=sc(u|m)}
	 * @since 2.1
	 */
	T arcsc(const T& x)
	{
		// p = c, q = n, r = d, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, p)
		return arcsp(x, x.get_field().get_one(), get_m().subtract(1).negate());
	}

	/** Evaluate inverse of Jacobi elliptic function sc.
	 * @param x value of Jacobi elliptic function {@code sc(u|m)}
	 * @return u such that {@code x=sc(u|m)}
	 * @since 2.1
	 */
	T arcsc(const double& x)
	{
		return arcsc(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function nd.
	 * @param x value of Jacobi elliptic function {@code nd(u|m)}
	 * @return u such that {@code x=nd(u|m)}
	 * @since 2.1
	 */
	T arcnd(const T& x)
	{
		// p = n, q = d, r = c, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, q)
		return arcpq(x, get_m().negate(), get_m().subtract(1));
	}

	/** Evaluate inverse of Jacobi elliptic function nd.
	 * @param x value of Jacobi elliptic function {@code nd(u|m)}
	 * @return u such that {@code x=nd(u|m)}
	 * @since 2.1
	 */
	T arcnd(const double& x)
	{
		return arcnd(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function sd.
	 * @param x value of Jacobi elliptic function {@code sd(u|m)}
	 * @return u such that {@code x=sd(u|m)}
	 * @since 2.1
	 */
	T arcsd(const T& x)
	{
		// p = d, q = n, r = c, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, p)
		return arcsp(x, get_m(), get_m().subtract(1));
	}

	/** Evaluate inverse of Jacobi elliptic function sd.
	 * @param x value of Jacobi elliptic function {@code sd(u|m)}
	 * @return u such that {@code x=sd(u|m)}
	 * @since 2.1
	 */
	T arcsd(const double& x)
	{
		return arcsd(get_m().get_field().get_zero().new_instance(x));
	}

	/** Evaluate inverse of Jacobi elliptic function cd.
	 * @param x value of Jacobi elliptic function {@code cd(u|m)}
	 * @return u such that {@code x=cd(u|m)}
	 * @since 2.1
	 */
	T arccd(const T& x)
	{
		// p = c, q = d, r = n, see DLMF 19.25.29 for evaluating \xce\x94\xe2\x81\xa1(q, p) and \xce\x94\xe2\x81\xa1(r, q)
		return arcpq(x, get_m().subtract(1).negate(), get_m());
	}

	/** Evaluate inverse of Jacobi elliptic function cd.
	 * @param x value of Jacobi elliptic function {@code cd(u|m)}
	 * @return u such that {@code x=cd(u|m)}
	 * @since 2.1
	 */
	T arccd(const double& x)
	{
		return arccd(get_m().get_field().get_zero().new_instance(x));
	}
};