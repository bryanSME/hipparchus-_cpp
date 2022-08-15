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
//package org.hipparchus.geometry.euclidean.twod;

//import java.text.Number_Format;

//import org.hipparchus.Field;
//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Math_Arrays;
#include <type_traits>
#include "../../../core/CalculusFieldElement.hpp"

/**
 * This class is a re-implementation of {@link Vector_2D} using {@link Calculus_Field_Element}.
 * <p>Instance of this class are guaranteed to be immutable.</p>
 * @param <T> the type of the field elements
 * @since 1.6
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Vector_2D 
{

    /** Abscissa. */
    private const T x;

    /** Ordinate. */
    private const T y;

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param x abscissa
     * @param y ordinate
     * @see #get_x()
     * @see #get_y()
     */
    public Field_Vector_2D(const T& x, const T& y) 
    {
        this.x = x;
        this.y = y;
    }

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param v coordinates array
     * @exception  if array does not have 2 elements
     * @see #to_array()
     */
    public Field_Vector_2D(const std::vector<T> v)  
    {
        if (v.size() != 2) 
        {
            throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), 2);
        }
        this.x = v[0];
        this.y = v[1];
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public Field_Vector_2D(const T& a, const Field_Vector_2D<T> u) 
    {
        this.x = a.multiply(u.x);
        this.y = a.multiply(u.y);
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public Field_Vector_2D(const T& a, const Vector_2D u) 
    {
        this.x = a.multiply(u.get_x());
        this.y = a.multiply(u.get_y());
    }

    /** Multiplicative constructor
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    public Field_Vector_2D(const double& a, const Field_Vector_2D<T> u) 
    {
        this.x = u.x.multiply(a);
        this.y = u.y.multiply(a);
    }

    /** Linear constructor
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public Field_Vector_2D(const T a1, const Field_Vector_2D<T> u1, const T a2, const Field_Vector_2D<T> u2) 
    {
        const T prototype = a1;
        this.x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x());
        this.y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y());
    }

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public Field_Vector_2D(const T a1, const Vector_2D u1, const T a2, const Vector_2D u2) 
    {
        const T prototype = a1;
        this.x = prototype.linear_combination(u1.get_x(), a1, u2.get_x(), a2);
        this.y = prototype.linear_combination(u1.get_y(), a1, u2.get_y(), a2);
    }

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    public Field_Vector_2D(const double& a1, const Field_Vector_2D<T> u1, const double& a2, const Field_Vector_2D<T> u2) 
    {
        const T prototype = u1.get_x();
        this.x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x());
        this.y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y());
    }

    /** Linear constructor.
     * Build a vector from three other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     */
    public Field_Vector_2D(const T a1, const Field_Vector_2D<T> u1, const T a2, const Field_Vector_2D<T> u2, const T a3, const Field_Vector_2D<T> u3) 
    {
        const T prototype = a1;
        this.x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x());
        this.y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y());
    }

    /** Linear constructor.
     * Build a vector from three other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     */
    public Field_Vector_2D(const T a1, const Vector_2D u1, const T a2, const Vector_2D u2, const T a3, const Vector_2D u3) 
    {
        const T prototype = a1;
        this.x = prototype.linear_combination(u1.get_x(), a1, u2.get_x(), a2, u3.get_x(), a3);
        this.y = prototype.linear_combination(u1.get_y(), a1, u2.get_y(), a2, u3.get_y(), a3);
    }

    /** Linear constructor.
     * Build a vector from three other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     */
    public Field_Vector_2D(const double& a1, const Field_Vector_2D<T> u1, const double& a2, const Field_Vector_2D<T> u2, const double& a3, const Field_Vector_2D<T> u3) 
    {
        const T prototype = u1.get_x();
        this.x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x());
        this.y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y());
    }

    /** Linear constructor.
     * Build a vector from four other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     * @param a4 fourth scale factor
     * @param u4 fourth base (unscaled) vector
     */
    public Field_Vector_2D(const T a1, const Field_Vector_2D<T> u1, const T a2, const Field_Vector_2D<T> u2, const T a3, const Field_Vector_2D<T> u3, const T a4, const Field_Vector_2D<T> u4) 
    {
        const T prototype = a1;
        this.x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x(), a4, u4.get_x());
        this.y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y(), a4, u4.get_y());
    }

    /** Linear constructor.
     * Build a vector from four other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     * @param a4 fourth scale factor
     * @param u4 fourth base (unscaled) vector
     */
    public Field_Vector_2D(const T a1, const Vector_2D u1, const T a2, const Vector_2D u2, const T a3, const Vector_2D u3, const T a4, const Vector_2D u4) 
    {
        const T prototype = a1;
        this.x = prototype.linear_combination(u1.get_x(), a1, u2.get_x(), a2, u3.get_x(), a3, u4.get_x(), a4);
        this.y = prototype.linear_combination(u1.get_y(), a1, u2.get_y(), a2, u3.get_y(), a3, u4.get_y(), a4);
    }

    /** Linear constructor.
     * Build a vector from four other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2 + a3 * u3 + a4 * u4
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     * @param a3 third scale factor
     * @param u3 third base (unscaled) vector
     * @param a4 fourth scale factor
     * @param u4 fourth base (unscaled) vector
     */
    public Field_Vector_2D(const double& a1, const Field_Vector_2D<T> u1, const double& a2, const Field_Vector_2D<T> u2, const double& a3, const Field_Vector_2D<T> u3, const double& a4, const Field_Vector_2D<T> u4) 
    {
        const T prototype = u1.get_x();
        this.x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x(), a4, u4.get_x());
        this.y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y(), a4, u4.get_y());
    }

    /** Build a {@link Field_Vector_2D} from a {@link Vector_2D}.
     * @param field field for the components
     * @param v vector to convert
     */
    public Field_Vector_2D(const Field<T> field, const Vector_2D v) 
    {
        this.x = field.get_zero().add(v.get_x());
        this.y = field.get_zero().add(v.get_y());
    }

    /** Get NULL vector (coordinates: 0, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_zero(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.ZERO);
    }

    /** Get first canonical vector (coordinates: 1, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_plus_i(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.PLUS_I);
    }

    /** Get opposite of the first canonical vector (coordinates: -1).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_minus_i(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.MINUS_I);
    }

    /** Get second canonical vector (coordinates: 0, 1).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_plus_j(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.PLUS_J);
    }

    /** Get opposite of the second canonical vector (coordinates: 0, -1).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_minus_j(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.MINUS_J);
    }

    /** Get a vector with all coordinates set to NaN.
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_nan(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.NaN);
    }

    /** Get a vector with all coordinates set to positive infinity.
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_positive_infinity(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.POSITIVE_INFINITY);
    }

    /** Get a vector with all coordinates set to negative infinity.
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static Field_Vector_2D<T> get_negative_infinity(const Field<T> field) 
    {
        return Field_Vector_2D<>(field, Vector_2D.NEGATIVE_INFINITY);
    }

    /** Get the abscissa of the vector.
     * @return abscissa of the vector
     * @see #Field_Vector_2D(Calculus_Field_Element, Calculus_Field_Element)
     */
    public T get_x() 
    {
        return x;
    }

    /** Get the ordinate of the vector.
     * @return ordinate of the vector
    * @see #Field_Vector_2D(Calculus_Field_Element, Calculus_Field_Element)
     */
    public T get_y() 
    {
        return y;
    }

    /** Get the vector coordinates as a dimension 2 array.
     * @return vector coordinates
     * @see #Field_Vector_2D(Calculus_Field_Element[])
     */
    public std::vector<T> to_array() 
    {
        const std::vector<T> array = Math_Arrays::build_array(x.get_field(), 2);
        array[0] = x;
        array[1] = y;
        return array;
    }

    /** Convert to a constant vector without extra field parts.
     * @return a constant vector
     */
    public Vector_2D to_vector_2d() 
    {
        return Vector_2D(x.get_real(), y.get_real());
    }

    /** Get the L<sub>1</sub> norm for the vector.
     * @return L<sub>1</sub> norm for the vector
     */
    public T get_norm1() 
    {
        return x.abs().add(y.abs());
    }

    /** Get the L<sub>2</sub> norm for the vector.
     * @return Euclidean norm for the vector
     */
    public T get_norm() 
    {
        // there are no cancellation problems here, so we use the straightforward formula
        return x.multiply(x).add(y.multiply(y)).sqrt();
    }

    /** Get the square of the norm for the vector.
     * @return square of the Euclidean norm for the vector
     */
    public T get_norm_sq() 
    {
        // there are no cancellation problems here, so we use the straightforward formula
        return x.multiply(x).add(y.multiply(y));
    }

    /** Get the L<sub>&infin;</sub> norm for the vector.
     * @return L<sub>&infin;</sub> norm for the vector
     */
    public T get_norm_inf() 
    {
        return std::max(std::abs(x), std::abs(y));
    }

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a vector
     */
    public Field_Vector_2D<T> add(const Field_Vector_2D<T> v) 
    {
        return Field_Vector_2D<>(x.add(v.x), y.add(v.y));
    }

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a vector
     */
    public Field_Vector_2D<T> add(const Vector_2D v) 
    {
        return Field_Vector_2D<>(x.add(v.get_x()), y.add(v.get_y()));
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    public Field_Vector_2D<T> add(const T factor, const Field_Vector_2D<T> v) 
    {
        return Field_Vector_2D<>(x.get_field().get_one(), this, factor, v);
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    public Field_Vector_2D<T> add(const T factor, const Vector_2D v) 
    {
        return Field_Vector_2D<>(x.add(factor.multiply(v.get_x())), y.add(factor.multiply(v.get_y())));
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    public Field_Vector_2D<T> add(const double factor, const Field_Vector_2D<T> v) 
    {
        return Field_Vector_2D<>(1.0, this, factor, v);
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    public Field_Vector_2D<T> add(const double factor, const Vector_2D v) 
    {
        return Field_Vector_2D<>(x.add(factor * v.get_x()), y.add(factor * v.get_y()));
    }

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a vector
     */
    public Field_Vector_2D<T> subtract(const Field_Vector_2D<T> v) 
    {
        return Field_Vector_2D<>(x.subtract(v.x), y.subtract(v.y));
    }

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a vector
     */
    public Field_Vector_2D<T> subtract(const Vector_2D v) 
    {
        return Field_Vector_2D<>(x.subtract(v.get_x()), y.subtract(v.get_y()));
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    public Field_Vector_2D<T> subtract(const T factor, const Field_Vector_2D<T> v) 
    {
        return Field_Vector_2D<>(x.get_field().get_one(), this, factor.negate(), v);
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    public Field_Vector_2D<T> subtract(const T factor, const Vector_2D v) 
    {
        return Field_Vector_2D<>(x.subtract(factor.multiply(v.get_x())), y.subtract(factor.multiply(v.get_y())));
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    public Field_Vector_2D<T> subtract(const double factor, const Field_Vector_2D<T> v) 
    {
        return Field_Vector_2D<>(1.0, this, -factor, v);
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    public Field_Vector_2D<T> subtract(const double factor, const Vector_2D v) 
    {
        return Field_Vector_2D<>(x.subtract(factor * v.get_x()), y.subtract(factor * v.get_y()));
    }

    /** Get a normalized vector aligned with the instance.
     * @return a normalized vector
     * @exception Math_Runtime_Exception if the norm is zero
     */
    public Field_Vector_2D<T> normalize() Math_Runtime_Exception 
    {
        const T s = get_norm();
        if (s.get_real() == 0) 
        {
            throw Math_Runtime_Exception(Localized_Geometry_Formats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
        }
        return scalar_multiply(s.reciprocal());
    }

    /** Compute the angular separation between two vectors.
     * <p>This method computes the angular separation between two
     * vectors using the dot product for well separated vectors and the
     * cross product for almost aligned vectors. This allows to have a
     * good accuracy in all cases, even for vectors very close to each
     * other.</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return angular separation between v1 and v2
     * @exception Math_Runtime_Exception if either vector has a NULL norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T angle(const Field_Vector_2D<T> v1, const Field_Vector_2D<T> v2)
        Math_Runtime_Exception 
        {

        const T norm_product = v1.get_norm().multiply(v2.get_norm());
        if (norm_product.get_real() == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        const T dot = v1.dot_product(v2);
        const double threshold = norm_product.get_real() * 0.9999;
        if (std::abs(dot.get_real()) > threshold) 
        {
            // the vectors are almost aligned, compute using the sine
            const T n = std::abs(dot.linear_combination(v1.x, v2.y, v1.y.negate(), v2.x));
            if (dot.get_real() >= 0) 
            {
                return std::asin(n.divide(norm_product));
            }
            return std::asin(n.divide(norm_product)).negate().add(dot.get_pi());
        }

        // the vectors are sufficiently separated to use the cosine
        return std::acos(dot.divide(norm_product));

    }

    /** Compute the angular separation between two vectors.
     * <p>This method computes the angular separation between two
     * vectors using the dot product for well separated vectors and the
     * cross product for almost aligned vectors. This allows to have a
     * good accuracy in all cases, even for vectors very close to each
     * other.</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return angular separation between v1 and v2
     * @exception Math_Runtime_Exception if either vector has a NULL norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T angle(const Field_Vector_2D<T> v1, const Vector_2D v2)
        Math_Runtime_Exception 
        {

        const T norm_product = v1.get_norm().multiply(v2.get_norm());
        if (norm_product.get_real() == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        const T dot = v1.dot_product(v2);
        const double threshold = norm_product.get_real() * 0.9999;
        if (std::abs(dot.get_real()) > threshold) 
        {
            // the vectors are almost aligned, compute using the sine
            const T n = std::abs(dot.linear_combination(v2.get_y(), v1.x, v2.get_x(), v1.y.negate()));
            if (dot.get_real() >= 0) 
            {
                return std::asin(n.divide(norm_product));
            }
            return std::asin(n.divide(norm_product)).negate().add(dot.get_pi());
        }

        // the vectors are sufficiently separated to use the cosine
        return std::acos(dot.divide(norm_product));

    }

    /** Compute the angular separation between two vectors.
     * <p>This method computes the angular separation between two
     * vectors using the dot product for well separated vectors and the
     * cross product for almost aligned vectors. This allows to have a
     * good accuracy in all cases, even for vectors very close to each
     * other.</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return angular separation between v1 and v2
     * @exception Math_Runtime_Exception if either vector has a NULL norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T angle(const Vector_2D v1, const Field_Vector_2D<T> v2)
        Math_Runtime_Exception 
        {
        return angle(v2, v1);
    }

    /** Get the opposite of the instance.
     * @return a vector which is opposite to the instance
     */
    public Field_Vector_2D<T> negate() 
    {
        return Field_Vector_2D<>(x.negate(), y.negate());
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a vector
     */
    public Field_Vector_2D<T> scalar_multiply(const T a) 
    {
        return Field_Vector_2D<>(x.multiply(a), y.multiply(a));
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a vector
     */
    public Field_Vector_2D<T> scalar_multiply(const double& a) 
    {
        return Field_Vector_2D<>(x.multiply(a), y.multiply(a));
    }

    /**
     * Returns true if any coordinate of this vector is NaN; false otherwise
     * @return  true if any coordinate of this vector is NaN; false otherwise
     */
    public bool is_nan() 
    {
        return std::isnan(x.get_real()) || std::isnan(y.get_real());
    }

    /**
     * Returns true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     * @return  true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     */
    public bool is_infinite() 
    {
        return !is_nan() && (std::isinf(x.get_real()) || std::isinf(y.get_real()));
    }

    /**
     * Test for the equality of two 2D vectors.
     * <p>
     * If all coordinates of two 2D vectors are exactly the same, and none of their
     * {@link Calculus_Field_Element#get_real() real part} are <code>NaN</code>, the
     * two 2D vectors are considered to be equal.
     * </p>
     * <p>
     * <code>NaN</code> coordinates are considered to affect globally the vector
     * and be equals to each other - i.e, if either (or all) real part of the
     * coordinates of the 3D vector are <code>NaN</code>, the 2D vector is <code>NaN</code>.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 2D vector objects are equal, false if
     *         object is NULL, not an instance of Field_Vector_2D, or
     *         not equal to this Field_Vector_2D instance
     *
     */
    //override
    public bool equals(Object other) 
    {

        if (this == other) 
        {
            return true;
        }

        if (other instanceof Field_Vector_2D) 
        {
            //@Suppress_Warnings("unchecked")
            const Field_Vector_2D<T> rhs = (Field_Vector_2D<T>) other;
            if (rhs.is_nan()) 
            {
                return this.is_nan();
            }

            return x.equals(rhs.x) && y.equals(rhs.y);

        }
        return false;
    }

    /**
     * Get a hash_code for the 3D vector.
     * <p>
     * All NaN values have the same hash code.</p>
     *
     * @return a hash code value for this object
     */
    //override
    public int hash_code() 
    {
        if (is_nan()) 
        {
            return 542;
        }
        return 122 * (76 * x.hash_code() +  y.hash_code());
    }

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    public T distance1(const Field_Vector_2D<T> v) 
    {
        const T dx = v.x.subtract(x).abs();
        const T dy = v.y.subtract(y).abs();
        return dx.add(dy);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    public T distance1(const Vector_2D v) 
    {
        const T dx = x.subtract(v.get_x()).abs();
        const T dy = y.subtract(v.get_y()).abs();
        return dx.add(dy);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>2</sub> norm
     */
    public T distance(const Field_Vector_2D<T> v) 
    {
        const T dx = v.x.subtract(x);
        const T dy = v.y.subtract(y);
        return dx.multiply(dx).add(dy.multiply(dy)).sqrt();
    }

    /** Compute the distance between the instance and another vector according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>2</sub> norm
     */
    public T distance(const Vector_2D v) 
    {
        const T dx = x.subtract(v.get_x());
        const T dy = y.subtract(v.get_y());
        return dx.multiply(dx).add(dy.multiply(dy)).sqrt();
    }

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    public T distance_inf(const Field_Vector_2D<T> v) 
    {
        const T dx = std::abs(x.subtract(v.x));
        const T dy = std::abs(y.subtract(v.y));
        return std::max(dx, dy);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    public T distance_inf(const Vector_2D v) 
    {
        const T dx = std::abs(x.subtract(v.get_x()));
        const T dy = std::abs(y.subtract(v.get_y()));
        return std::max(dx, dy);
    }

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    public T distance_sq(const Field_Vector_2D<T> v) 
    {
        const T dx = v.x.subtract(x);
        const T dy = v.y.subtract(y);
        return dx.multiply(dx).add(dy.multiply(dy));
    }

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    public T distance_sq(const Vector_2D v) 
    {
        const T dx = x.subtract(v.get_x());
        const T dy = y.subtract(v.get_y());
        return dx.multiply(dx).add(dy.multiply(dy));
    }


    /** Compute the dot-product of the instance and another vector.
     * <p>
     * The implementation uses specific multiplication and addition
     * algorithms to preserve accuracy and reduce cancellation effects.
     * It should be very accurate even for nearly orthogonal vectors.
     * </p>
     * @see Math_Arrays#linear_combination(double, double, double, double, double, double)
     * @param v second vector
     * @return the dot product this.v
     */
    public T dot_product(const Field_Vector_2D<T> v) 
    {
        return x.linear_combination(x, v.get_x(), y, v.get_y());
    }

    /** Compute the dot-product of the instance and another vector.
     * <p>
     * The implementation uses specific multiplication and addition
     * algorithms to preserve accuracy and reduce cancellation effects.
     * It should be very accurate even for nearly orthogonal vectors.
     * </p>
     * @see Math_Arrays#linear_combination(double, double, double, double, double, double)
     * @param v second vector
     * @return the dot product this.v
     */
    public T dot_product(const Vector_2D v) 
    {
        return x.linear_combination(v.get_x(), x, v.get_y(), y);
    }

    /**
     * Compute the cross-product of the instance and the given points.
     * <p>
     * The cross product can be used to determine the location of a point
     * with regard to the line formed by (p1, p2) and is calculated as:
     * \\[
     *    P = (x_2 - x_1)(y_3 - y_1) - (y_2 - y_1)(x_3 - x_1)
     * \\]
     * with \\(p3 = (x_3, y_3)\\) being this instance.
     * <p>
     * If the result is 0, the points are collinear, i.e. lie on a single straight line L;
     * if it is positive, this point lies to the left, otherwise to the right of the line
     * formed by (p1, p2).
     *
     * @param p1 first point of the line
     * @param p2 second point of the line
     * @return the cross-product
     *
     * @see <a href="http://en.wikipedia.org/wiki/Cross_product">Cross product (Wikipedia)</a>
     */
    public T cross_product(const Field_Vector_2D<T> p1, const Field_Vector_2D<T> p2) 
    {
        const T x1  = p2.get_x().subtract(p1.get_x());
        const T y1  = get_y().subtract(p1.get_y());
        const T mx2 = p1.get_x().subtract(get_x());
        const T y2  = p2.get_y().subtract(p1.get_y());
        return x1.linear_combination(x1, y1, mx2, y2);
    }

    /**
     * Compute the cross-product of the instance and the given points.
     * <p>
     * The cross product can be used to determine the location of a point
     * with regard to the line formed by (p1, p2) and is calculated as:
     * \\[
     *    P = (x_2 - x_1)(y_3 - y_1) - (y_2 - y_1)(x_3 - x_1)
     * \\]
     * with \\(p3 = (x_3, y_3)\\) being this instance.
     * <p>
     * If the result is 0, the points are collinear, i.e. lie on a single straight line L;
     * if it is positive, this point lies to the left, otherwise to the right of the line
     * formed by (p1, p2).
     *
     * @param p1 first point of the line
     * @param p2 second point of the line
     * @return the cross-product
     *
     * @see <a href="http://en.wikipedia.org/wiki/Cross_product">Cross product (Wikipedia)</a>
     */
    public T cross_product(const Vector_2D& p1, const Vector_2D& p2) 
    {
        const double x1  = p2.get_x() - p1.get_x();
        const T      y1  = get_y().subtract(p1.get_y());
        const T      x2 = get_x().subtract(p1.get_x());
        const double y2  = p2.get_y() - p1.get_y();
        return y1.linear_combination(x1, y1, -y2, x2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance1(const Field_Vector_2D<T> p1, const Field_Vector_2D<T> p2) 
    {
        return p1.distance1(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T  distance1(const Field_Vector_2D<T> p1, const Vector_2D& p2) 
    {
        return p1.distance1(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T  distance1(const Vector_2D& p1, const Field_Vector_2D<T> p2) 
    {
        return p2.distance1(p1);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance(const Field_Vector_2D<T> p1, const Field_Vector_2D<T> p2) 
    {
        return p1.distance(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance(const Field_Vector_2D<T> p1, const Vector_2D& p2) 
    {
        return p1.distance(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance( const Vector_2D& p1, const Field_Vector_2D<T> p2) 
    {
        return p2.distance(p1);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance_inf(const Field_Vector_2D<T> p1, const Field_Vector_2D<T> p2) 
    {
        return p1.distance_inf(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance_inf(const Field_Vector_2D<T> p1, const Vector_2D& p2) 
    {
        return p1.distance_inf(p2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the distance between p1 and p2 according to the L<sub>&infin;</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance_inf(const Vector_2D& p1, const Field_Vector_2D<T> p2) 
    {
        return p2.distance_inf(p1);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between p1 and p2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance_sq(const Field_Vector_2D<T> p1, const Field_Vector_2D<T> p2) 
    {
        return p1.distance_sq(p2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between p1 and p2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance_sq(const Field_Vector_2D<T> p1, const Vector_2D& p2) 
    {
        return p1.distance_sq(p2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>p1.subtract(p2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param p1 first vector
     * @param p2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between p1 and p2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T distance_sq(const Vector_2D& p1, const Field_Vector_2D<T> p2) 
    {
        return p2.distance_sq(p1);
    }

    /** Compute the orientation of a triplet of points.
     * @param p first vector of the triplet
     * @param q second vector of the triplet
     * @param r third vector of the triplet
     * @param <T> the type of the field elements
     * @return a positive value if (p, q, r) defines a counterclockwise oriented
     * triangle, a negative value if (p, q, r) defines a clockwise oriented
     * triangle, and 0 if (p, q, r) are collinear or some points are equal
     * @since 1.2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    public static T orientation(const Field_Vector_2D<T> p, const Field_Vector_2D<T> q, const Field_Vector_2D<T> r) 
    {
        const T prototype = p.get_x();
        const std::vector<T> a = Math_Arrays::build_array(prototype.get_field(), 6);
        a[0] = p.get_x();
        a[1] = p.get_x().negate();
        a[2] = q.get_x();
        a[3] = q.get_x().negate();
        a[4] = r.get_x();
        a[5] = r.get_x().negate();
        const std::vector<T> b = Math_Arrays::build_array(prototype.get_field(), 6);
        b[0] = q.get_y();
        b[1] = r.get_y();
        b[2] = r.get_y();
        b[3] = p.get_y();
        b[4] = p.get_y();
        b[5] = q.get_y();
        return prototype.linear_combination(a, b);
    }

    /** Get a string representation of this vector.
     * @return a string representation of this vector
     */
    //override
    public std::string to_string() const 
    {
        return Vector_2D_Format.get_vector_2d_format().format(to_vector_2d());
    }

    /** Get a string representation of this vector.
     * @param format the custom format for components
     * @return a string representation of this vector
     */
    public std::string to_string(const Number_Format format) 
    {
        return Vector_2D_Format(format).format(to_vector_2d());
    }

};