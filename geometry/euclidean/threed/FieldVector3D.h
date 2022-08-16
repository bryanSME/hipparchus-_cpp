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

//package org.hipparchus.geometry.euclidean.threed;

//import java.io.Serializable;
//import java.text.Number_Format;

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.Localized_Core_Formats;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Math_Arrays;
#include <vector>
#include <string>
#include <type_traits>
#include "../../../core/CalculusFieldElement.hpp"
/**
 * This class is a re-implementation of {@link Vector_3D} using {@link Calculus_Field_Element}.
 * <p>Instance of this class are guaranteed to be immutable.</p>
 * @param <T> the type of the field elements
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Vector_3D
{
private:

    /** Abscissa. */
    const T my_x;

    /** Ordinate. */
    const T my_y;

    /** Height. */
    const T my_z;

public:
    /** Simple constructor.
     * Build a vector from its coordinates
     * @param x abscissa
     * @param y ordinate
     * @param z height
     * @see #get_x()
     * @see #get_y()
     * @see #get_z()
     */
    Field_Vector_3D(const T& x, const T y, const T z) : my_x{ x }, my_y{ y }, my_z{ z }{}

    /** Simple constructor.
     * Build a vector from its coordinates
     * @param v coordinates array
     * @exception  if array does not have 3 elements
     * @see #to_array()
     */
    Field_Vector_3D(const std::vector<T> v)  
    {
        if (v.size() != 3) 
        {
            throw std::exception("not implemented");
            //throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, v.size(), 3);
        }
        my_x = v[0];
        my_y = v[1];
        my_z = v[2];
    }

    /** Simple constructor.
     * Build a vector from its azimuthal coordinates
     * @param alpha azimuth (&alpha;) around Z
     *              (0 is +X, &pi;/2 is +Y, &pi; is -X and 3&pi;/2 is -Y)
     * @param delta elevation (&delta;) above (XY) plane, from -&pi;/2 to +&pi;/2
     * @see #get_alpha()
     * @see #get_delta()
     */
    Field_Vector_3D(const T alpha, const T delta) 
    {
        Field_Sin_Cos<T> sin_cos_alpha = Sin_Cos(alpha);
        Field_Sin_Cos<T> sin_cos_delta = Sin_Cos(delta);
        my_x = sin_cos_alpha.cos().multiply(sin_cos_delta.cos());
        my_y = sin_cos_alpha.sin().multiply(sin_cos_delta.cos());
        my_z = sin_cos_delta.sin();
    }

    /** Multiplicative constructor.
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    Field_Vector_3D(const T& a, const Field_Vector_3D<T>u)
        :
        my_x{ a.multiply(u.x) },
        my_y{ a.multiply(u.y) },
        my_z{ a.multiply(u.z) }
    {}

    /** Multiplicative constructor.
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    Field_Vector_3D(const T& a, const Vector_3D u) 
        :
        my_x{ a.multiply(u.get_x()) },
        my_y{ a.multiply(u.get_y()) },
        my_z{ a.multiply(u.get_z()) }
    {}

    /** Multiplicative constructor.
     * Build a vector from another one and a scale factor.
     * The vector built will be a * u
     * @param a scale factor
     * @param u base (unscaled) vector
     */
    Field_Vector_3D(const double& a, const Field_Vector_3D<T> u) 
        :
        my_x{ u.x.multiply(a) },
        my_y{ u.y.multiply(a) },
        my_z{ u.z.multiply(a) }
    {}

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    Field_Vector_3D(const T& a1, const Field_Vector_3D<T>& u1, const T a2, const Field_Vector_3D<T> u2) 
    {
        const T prototype = a1;
        my_x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x());
        my_y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y());
        my_z = prototype.linear_combination(a1, u1.get_z(), a2, u2.get_z());
    }

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    Field_Vector_3D(const T& a1, const Vector_3D& u1, const T& a2, const Vector_3D& u2) 
    {
        const T prototype = a1;
        my_x = prototype.linear_combination(u1.get_x(), a1, u2.get_x(), a2);
        my_y = prototype.linear_combination(u1.get_y(), a1, u2.get_y(), a2);
        my_z = prototype.linear_combination(u1.get_z(), a1, u2.get_z(), a2);
    }

    /** Linear constructor.
     * Build a vector from two other ones and corresponding scale factors.
     * The vector built will be a1 * u1 + a2 * u2
     * @param a1 first scale factor
     * @param u1 first base (unscaled) vector
     * @param a2 second scale factor
     * @param u2 second base (unscaled) vector
     */
    Field_Vector_3D(const double& a1, const Field_Vector_3D<T>& u1, const double& a2, const Field_Vector_3D<T> u2) 
    {
        const T prototype = u1.get_x();
        my_x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x());
        my_y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y());
        my_z = prototype.linear_combination(a1, u1.get_z(), a2, u2.get_z());
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
    Field_Vector_3D(const T& a1, const Field_Vector_3D<T>& u1, const T a2, const Field_Vector_3D<T>& u2, const T a3, const Field_Vector_3D<T> u3) 
    {
        const T prototype = a1;
        my_x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x());
        my_y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y());
        my_z = prototype.linear_combination(a1, u1.get_z(), a2, u2.get_z(), a3, u3.get_z());
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
    Field_Vector_3D(const T& a1, const Vector_3D u1, const T a2, const Vector_3D u2, const T a3, const Vector_3D u3) 
    {
        const T prototype = a1;
        my_x = prototype.linear_combination(u1.get_x(), a1, u2.get_x(), a2, u3.get_x(), a3);
        my_y = prototype.linear_combination(u1.get_y(), a1, u2.get_y(), a2, u3.get_y(), a3);
        my_z = prototype.linear_combination(u1.get_z(), a1, u2.get_z(), a2, u3.get_z(), a3);
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
    Field_Vector_3D(const double& a1, const Field_Vector_3D<T>& u1, const double& a2, const Field_Vector_3D<T>& u2, const double& a3, const Field_Vector_3D<T> u3) 
    {
        const T prototype = u1.get_x();
        my_x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x());
        my_y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y());
        my_z = prototype.linear_combination(a1, u1.get_z(), a2, u2.get_z(), a3, u3.get_z());
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
    Field_Vector_3D(const T& a1, const Field_Vector_3D<T>& u1, const T a2, const Field_Vector_3D<T>& u2, const T a3, const Field_Vector_3D<T>& u3, const T& a4, const Field_Vector_3D<T>& u4) 
    {
        const T prototype = a1;
        my_x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x(), a4, u4.get_x());
        my_y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y(), a4, u4.get_y());
        my_z = prototype.linear_combination(a1, u1.get_z(), a2, u2.get_z(), a3, u3.get_z(), a4, u4.get_z());
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
    Field_Vector_3D(const T& a1, const Vector_3D u1, const T a2, const Vector_3D u2, const T a3, const Vector_3D u3, const T& a4, const Vector_3D u4) 
    {
        const T prototype = a1;
        my_x = prototype.linear_combination(u1.get_x(), a1, u2.get_x(), a2, u3.get_x(), a3, u4.get_x(), a4);
        my_y = prototype.linear_combination(u1.get_y(), a1, u2.get_y(), a2, u3.get_y(), a3, u4.get_y(), a4);
        my_z = prototype.linear_combination(u1.get_z(), a1, u2.get_z(), a2, u3.get_z(), a3, u4.get_z(), a4);
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
    Field_Vector_3D(const double& a1, const Field_Vector_3D<T>& u1, const double& a2, const Field_Vector_3D<T>& u2, const double& a3, const Field_Vector_3D<T>& u3, const double& a4, const Field_Vector_3D<T>& u4) 
    {
        const T prototype = u1.get_x();
        my_x = prototype.linear_combination(a1, u1.get_x(), a2, u2.get_x(), a3, u3.get_x(), a4, u4.get_x());
        my_y = prototype.linear_combination(a1, u1.get_y(), a2, u2.get_y(), a3, u3.get_y(), a4, u4.get_y());
        my_z = prototype.linear_combination(a1, u1.get_z(), a2, u2.get_z(), a3, u3.get_z(), a4, u4.get_z());
    }

    /** Build a {@link Field_Vector_3D} from a {@link Vector_3D}.
     * @param field field for the components
     * @param v vector to convert
     */
    Field_Vector_3D(const Field<T> field, const Vector_3D v) 
    {
        my_x = field.get_zero().add(v.get_x());
        my_y = field.get_zero().add(v.get_y());
        my_z = field.get_zero().add(v.get_z());
    }

    /** Get NULL vector (coordinates: 0, 0, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_zero(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.ZERO);
    }

    /** Get first canonical vector (coordinates: 1, 0, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_plus_i(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.PLUS_I);
    }

    /** Get opposite of the first canonical vector (coordinates: -1, 0, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_minus_i(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.MINUS_I);
    }

    /** Get second canonical vector (coordinates: 0, 1, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_plus_j(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.PLUS_J);
    }

    /** Get opposite of the second canonical vector (coordinates: 0, -1, 0).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_minus_j(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.MINUS_J);
    }

    /** Get third canonical vector (coordinates: 0, 0, 1).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_plus_k(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.PLUS_K);
    }

    /** Get opposite of the third canonical vector (coordinates: 0, 0, -1).
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_minus_k(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.MINUS_K);
    }

    /** Get a vector with all coordinates set to NaN.
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_nan(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.NaN);
    }

    /** Get a vector with all coordinates set to positive infinity.
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_positive_infinity(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.POSITIVE_INFINITY);
    }

    /** Get a vector with all coordinates set to negative infinity.
     * @param field field for the components
     * @return a vector
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> get_negative_infinity(const Field<T> field) 
    {
        return Field_Vector_3D<>(field, Vector_3D.NEGATIVE_INFINITY);
    }

    /** Get the abscissa of the vector.
     * @return abscissa of the vector
     * @see #Field_Vector_3D(Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element)
     */
    T get_x() 
    {
        return x;
    }

    /** Get the ordinate of the vector.
     * @return ordinate of the vector
     * @see #Field_Vector_3D(Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element)
     */
    T get_y() 
    {
        return y;
    }

    /** Get the height of the vector.
     * @return height of the vector
     * @see #Field_Vector_3D(Calculus_Field_Element, Calculus_Field_Element, Calculus_Field_Element)
     */
    T get_z() 
    {
        return z;
    }

    /** Get the vector coordinates as a dimension 3 array.
     * @return vector coordinates
     * @see #Field_Vector_3D(Calculus_Field_Element[])
     */
    std::vector<T> to_array() 
    {
        const std::vector<T> array = Math_Arrays::build_array(x.get_field(), 3);
        array[0] = x;
        array[1] = y;
        array[2] = z;
        return array;
    }

    /** Convert to a constant vector without extra field parts.
     * @return a constant vector
     */
    Vector_3D to_vector_3d() 
    {
        return Vector_3D(x.get_real(), y.get_real(), z.get_real());
    }

    /** Get the L<sub>1</sub> norm for the vector.
     * @return L<sub>1</sub> norm for the vector
     */
    T get_norm1() 
    {
        return x.abs().add(y.abs()).add(z.abs());
    }

    /** Get the L<sub>2</sub> norm for the vector.
     * @return Euclidean norm for the vector
     */
    T get_norm() 
    {
        // there are no cancellation problems here, so we use the straightforward formula
        return x.multiply(x).add(y.multiply(y)).add(z.multiply(z)).sqrt();
    }

    /** Get the square of the norm for the vector.
     * @return square of the Euclidean norm for the vector
     */
    T get_norm_sq() 
    {
        // there are no cancellation problems here, so we use the straightforward formula
        return x.multiply(x).add(y.multiply(y)).add(z.multiply(z));
    }

    /** Get the L<sub>&infin;</sub> norm for the vector.
     * @return L<sub>&infin;</sub> norm for the vector
     */
    T get_norm_inf() 
    {
        return std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
    }

    /** Get the azimuth of the vector.
     * @return azimuth (&alpha;) of the vector, between -&pi; and +&pi;
     * @see #Field_Vector_3D(Calculus_Field_Element, Calculus_Field_Element)
     */
    T get_alpha() 
    {
        return y.atan2(x);
    }

    /** Get the elevation of the vector.
     * @return elevation (&delta;) of the vector, between -&pi;/2 and +&pi;/2
     * @see #Field_Vector_3D(Calculus_Field_Element, Calculus_Field_Element)
     */
    T get_delta() 
    {
        return z.divide(get_norm()).asin();
    }

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a vector
     */
    Field_Vector_3D<T> add(const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(x.add(v.x), y.add(v.y), z.add(v.z));
    }

    /** Add a vector to the instance.
     * @param v vector to add
     * @return a vector
     */
    Field_Vector_3D<T> add(const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.add(v.get_x()), y.add(v.get_y()), z.add(v.get_z()));
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    Field_Vector_3D<T> add(const T factor, const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(x.get_field().get_one(), this, factor, v);
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    Field_Vector_3D<T> add(const T factor, const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.add(factor.multiply(v.get_x())), y.add(factor.multiply(v.get_y())), z.add(factor.multiply(v.get_z())));
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    Field_Vector_3D<T> add(const double factor, const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(1.0, this, factor, v);
    }

    /** Add a scaled vector to the instance.
     * @param factor scale factor to apply to v before adding it
     * @param v vector to add
     * @return a vector
     */
    Field_Vector_3D<T> add(const double factor, const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.add(factor * v.get_x()), y.add(factor * v.get_y()), z.add(factor * v.get_z()));
    }

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a vector
     */
    Field_Vector_3D<T> subtract(const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(x.subtract(v.x), y.subtract(v.y), z.subtract(v.z));
    }

    /** Subtract a vector from the instance.
     * @param v vector to subtract
     * @return a vector
     */
    Field_Vector_3D<T> subtract(const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.subtract(v.get_x()), y.subtract(v.get_y()), z.subtract(v.get_z()));
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    Field_Vector_3D<T> subtract(const T factor, const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(x.get_field().get_one(), this, factor.negate(), v);
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    Field_Vector_3D<T> subtract(const T factor, const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.subtract(factor.multiply(v.get_x())), y.subtract(factor.multiply(v.get_y())), z.subtract(factor.multiply(v.get_z())));
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    Field_Vector_3D<T> subtract(const double factor, const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(1.0, this, -factor, v);
    }

    /** Subtract a scaled vector from the instance.
     * @param factor scale factor to apply to v before subtracting it
     * @param v vector to subtract
     * @return a vector
     */
    Field_Vector_3D<T> subtract(const double factor, const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.subtract(factor * v.get_x()), y.subtract(factor * v.get_y()), z.subtract(factor * v.get_z()));
    }

    /** Get a normalized vector aligned with the instance.
     * @return a normalized vector
     * @exception Math_Runtime_Exception if the norm is zero
     */
    Field_Vector_3D<T> normalize() Math_Runtime_Exception 
    {
        const T s = get_norm();
        if (s.get_real() == 0) 
        {
            throw Math_Runtime_Exception(Localized_Geometry_Formats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR);
        }
        return scalar_multiply(s.reciprocal());
    }

    /** Get a vector orthogonal to the instance.
     * <p>There are an infinite number of normalized vectors orthogonal
     * to the instance. This method picks up one of them almost
     * arbitrarily. It is useful when one needs to compute a reference
     * frame with one of the axes in a predefined direction. The
     * following example shows how to build a frame having the k axis
     * aligned with the known vector u :
     * <pre><code>
     *   Vector_3D k = u.normalize();
     *   Vector_3D i = k.orthogonal();
     *   Vector_3D j = Vector_3D.cross_product(k, i);
     * </code></pre></p>
     * @return a normalized vector orthogonal to the instance
     * @exception Math_Runtime_Exception if the norm of the instance is NULL
     */
    Field_Vector_3D<T> orthogonal() Math_Runtime_Exception 
    {

        const double threshold = 0.6 * get_norm().get_real();
        if (threshold == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        if (std::abs(x.get_real()) <= threshold) 
        {
            const T inverse  = y.multiply(y).add(z.multiply(z)).sqrt().reciprocal();
            return Field_Vector_3D<T>(inverse.get_field().get_zero(), inverse.multiply(z), inverse.multiply(y).negate());
        }
        
        if (std::abs(y.get_real()) <= threshold) 
        {
            const T inverse  = x.multiply(x).add(z.multiply(z)).sqrt().reciprocal();
            return Field_Vector_3D<T>(inverse.multiply(z).negate(), inverse.get_field().get_zero(), inverse.multiply(x));
        }

        const T inverse  = x.multiply(x).add(y.multiply(y)).sqrt().reciprocal();
        return Field_Vector_3D<T>(inverse.multiply(y), inverse.multiply(x).negate(), inverse.get_field().get_zero());
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
    static  T angle(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2)
        Math_Runtime_Exception 
        {

        const T norm_product = v1.get_norm().multiply(v2.get_norm());
        if (norm_product.get_real() == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        const T dot = dot_product(v1, v2);
        const double threshold = norm_product.get_real() * 0.9999;
        if ((dot.get_real() < -threshold) || (dot.get_real() > threshold)) 
        {
            // the vectors are almost aligned, compute using the sine
            Field_Vector_3D<T> v3 = cross_product(v1, v2);
            if (dot.get_real() >= 0) 
            {
                return v3.get_norm().divide(norm_product).asin();
            }
            return v3.get_norm().divide(norm_product).asin().subtract(dot.get_pi()).negate();
        }

        // the vectors are sufficiently separated to use the cosine
        return dot.divide(norm_product).acos();
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
    static  T angle(const Field_Vector_3D<T> v1, const Vector_3D v2)
    {

        const T norm_product = v1.get_norm().multiply(v2.get_norm());
        if (norm_product.get_real() == 0) 
        {
            throw Math_Runtime_Exception(hipparchus::exception::Localized_Core_Formats_Type::ZERO_NORM);
        }

        const T dot = dot_product(v1, v2);
        const double threshold = norm_product.get_real() * 0.9999;
        if ((dot.get_real() < -threshold) || (dot.get_real() > threshold)) 
        {
            // the vectors are almost aligned, compute using the sine
            Field_Vector_3D<T> v3 = cross_product(v1, v2);
            if (dot.get_real() >= 0) 
            {
                return v3.get_norm().divide(norm_product).asin();
            }
            return v3.get_norm().divide(norm_product).asin().subtract(dot.get_pi()).negate();
        }

        // the vectors are sufficiently separated to use the cosine
        return dot.divide(norm_product).acos();
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
    static  T angle(const Vector_3D v1, const Field_Vector_3D<T> v2)
    {
        return angle(v2, v1);
    }

    /** Get the opposite of the instance.
     * @return a vector which is opposite to the instance
     */
    Field_Vector_3D<T> negate() 
    {
        return Field_Vector_3D<T>(x.negate(), y.negate(), z.negate());
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a vector
     */
    Field_Vector_3D<T> scalar_multiply(const T a) 
    {
        return Field_Vector_3D<T>(x.multiply(a), y.multiply(a), z.multiply(a));
    }

    /** Multiply the instance by a scalar.
     * @param a scalar
     * @return a vector
     */
    Field_Vector_3D<T> scalar_multiply(const double& a) 
    {
        return Field_Vector_3D<T>(x.multiply(a), y.multiply(a), z.multiply(a));
    }

    /**
     * Returns true if any coordinate of this vector is NaN; false otherwise
     * @return  true if any coordinate of this vector is NaN; false otherwise
     */
    bool is_nan() 
    {
        return std::isnan(x.get_real()) || std::isnan(y.get_real()) || std::isnan(z.get_real());
    }

    /**
     * Returns true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     * @return  true if any coordinate of this vector is infinite and none are NaN;
     * false otherwise
     */
    bool is_infinite() 
    {
        return !is_nan() && (std::isinfinite(x.get_real()) || std::isinfinite(y.get_real()) || std::isinfinite(z.get_real()));
    }

    /**
     * Test for the equality of two 3D vectors.
     * <p>
     * If all coordinates of two 3D vectors are exactly the same, and none of their
     * {@link Calculus_Field_Element#get_real() real part} are <code>NaN</code>, the
     * two 3D vectors are considered to be equal.
     * </p>
     * <p>
     * <code>NaN</code> coordinates are considered to affect globally the vector
     * and be equals to each other - i.e, if either (or all) real part of the
     * coordinates of the 3D vector are <code>NaN</code>, the 3D vector is <code>NaN</code>.
     * </p>
     *
     * @param other Object to test for equality to this
     * @return true if two 3D vector objects are equal, false if
     *         object is NULL, not an instance of Field_Vector_3D, or
     *         not equal to this Field_Vector_3D instance
     *
     */
    //override
    bool equals(const Object& other) 
    {
        if (*this == other) 
        {
            return true;
        }

        if (dynamic_cast<const Field_Vector_3D*>(*other) != nullptr)
        {
            ////@Suppress_Warnings("unchecked")
            const Field_Vector_3D<T> rhs = (Field_Vector_3D<T>) other;
            if (rhs.is_nan()) 
            {
                return this.is_nan();
            }

            return x.equals(rhs.x) && y.equals(rhs.y) && z.equals(rhs.z);

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
    int hash_code() 
    {
        if (is_nan()) 
        {
            return 409;
        }
        return 311 * (107 * x.hash_code() + 83 * y.hash_code() +  z.hash_code());
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
    T dot_product(const Field_Vector_3D<T> v) 
    {
        return x.linear_combination(x, v.x, y, v.y, z, v.z);
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
    T dot_product(const Vector_3D v) 
    {
        return x.linear_combination(v.get_x(), x, v.get_y(), y, v.get_z(), z);
    }

    /** Compute the cross-product of the instance with another vector.
     * @param v other vector
     * @return the cross product this ^ v as a Vector_3D
     */
    Field_Vector_3D<T> cross_product(const Field_Vector_3D<T> v) 
    {
        return Field_Vector_3D<T>(x.linear_combination(y, v.z, z.negate(), v.y), y.linear_combination(z, v.x, x.negate(), v.z), z.linear_combination(x, v.y, y.negate(), v.x));
    }

    /** Compute the cross-product of the instance with another vector.
     * @param v other vector
     * @return the cross product this ^ v as a Vector_3D
     */
    Field_Vector_3D<T> cross_product(const Vector_3D v) 
    {
        return Field_Vector_3D<T>(x.linear_combination(v.get_z(), y, -v.get_y(), z), y.linear_combination(v.get_x(), z, -v.get_z(), x), z.linear_combination(v.get_y(), x, -v.get_x(), y));
    }

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    T distance1(const Field_Vector_3D<T> v) 
    {
        const T dx = v.x.subtract(x).abs();
        const T dy = v.y.subtract(y).abs();
        const T dz = v.z.subtract(z).abs();
        return dx.add(dy).add(dz);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>1</sub> norm
     */
    T distance1(const Vector_3D v) 
    {
        const T dx = x.subtract(v.get_x()).abs();
        const T dy = y.subtract(v.get_y()).abs();
        const T dz = z.subtract(v.get_z()).abs();
        return dx.add(dy).add(dz);
    }

    /** Compute the distance between the instance and another vector according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>2</sub> norm
     */
    T distance(const Field_Vector_3D<T> v) 
    {
        const T dx = v.x.subtract(x);
        const T dy = v.y.subtract(y);
        const T dz = v.z.subtract(z);
        return dx.multiply(dx).add(dy.multiply(dy)).add(dz.multiply(dz)).sqrt();
    }

    /** Compute the distance between the instance and another vector according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>2</sub> norm
     */
    T distance(const Vector_3D v) 
    {
        const T dx = x.subtract(v.get_x());
        const T dy = y.subtract(v.get_y());
        const T dz = z.subtract(v.get_z());
        return dx.multiply(dx).add(dy.multiply(dy)).add(dz.multiply(dz)).sqrt();
    }

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    T distance_inf(const Field_Vector_3D<T> v) 
    {
        const T dx = v.x.subtract(x).abs();
        const T dy = v.y.subtract(y).abs();
        const T dz = v.z.subtract(z).abs();
        if (dx.get_real() <= dy.get_real()) 
        {
            if (dy.get_real() <= dz.get_real()) 
            {
                return dz;
            }
            return dy;
        }

        if (dx.get_real() <= dz.get_real()) 
        {
            return dz;
        }
        return dx;
    }

    /** Compute the distance between the instance and another vector according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the distance between the instance and p according to the L<sub>&infin;</sub> norm
     */
    T distance_inf(const Vector_3D v) 
    {
        const T dx = x.subtract(v.get_x()).abs();
        const T dy = y.subtract(v.get_y()).abs();
        const T dz = z.subtract(v.get_z()).abs();
        if (dx.get_real() <= dy.get_real()) 
        {
            if (dy.get_real() <= dz.get_real()) 
            {
                return dz;
            }
            return dy;
        }

        if (dx.get_real() <= dz.get_real()) 
        {
            return dz;
        }
        return dx;
    }

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    T distance_sq(const Field_Vector_3D<T> v) 
    {
        const T dx = v.x.subtract(x);
        const T dy = v.y.subtract(y);
        const T dz = v.z.subtract(z);
        return dx.multiply(dx).add(dy.multiply(dy)).add(dz.multiply(dz));
    }

    /** Compute the square of the distance between the instance and another vector.
     * <p>Calling this method is equivalent to calling:
     * <code>q.subtract(p).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v second vector
     * @return the square of the distance between the instance and p
     */
    T distance_sq(const Vector_3D v) 
    {
        const T dx = x.subtract(v.get_x());
        const T dy = y.subtract(v.get_y());
        const T dz = z.subtract(v.get_z());
        return dx.multiply(dx).add(dy.multiply(dy)).add(dz.multiply(dz));
    }

    /** Compute the dot-product of two vectors.
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the dot product v1.v2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T dot_product(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2) 
    {
        return v1.dot_product(v2);
    }

    /** Compute the dot-product of two vectors.
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the dot product v1.v2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T dot_product(const Field_Vector_3D<T> v1, const Vector_3D v2) 
    {
        return v1.dot_product(v2);
    }

    /** Compute the dot-product of two vectors.
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the dot product v1.v2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T dot_product(const Vector_3D v1, const Field_Vector_3D<T> v2) 
    {
        return v2.dot_product(v1);
    }

    /** Compute the cross-product of two vectors.
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the cross product v1 ^ v2 as a Vector
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> cross_product(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2) 
    {
        return v1.cross_product(v2);
    }

    /** Compute the cross-product of two vectors.
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the cross product v1 ^ v2 as a Vector
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> cross_product(const Field_Vector_3D<T> v1, const Vector_3D v2) 
    {
        return v1.cross_product(v2);
    }

    /** Compute the cross-product of two vectors.
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the cross product v1 ^ v2 as a Vector
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> cross_product(const Vector_3D v1, const Field_Vector_3D<T> v2) 
    {
        return Field_Vector_3D<T>(v2.x.linear_combination(v1.get_y(), v2.z, -v1.get_z(), v2.y), v2.y.linear_combination(v1.get_z(), v2.x, -v1.get_x(), v2.z), v2.z.linear_combination(v1.get_x(), v2.y, -v1.get_y(), v2.x));
    }

    /** Compute the distance between two vectors according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>1</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance1(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2) 
    {
        return v1.distance1(v2);
    }

    /** Compute the distance between two vectors according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>1</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance1(const Field_Vector_3D<T> v1, const Vector_3D v2) 
    {
        return v1.distance1(v2);
    }

    /** Compute the distance between two vectors according to the L<sub>1</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm1()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>1</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance1(const Vector_3D v1, const Field_Vector_3D<T> v2) 
    {
        return v2.distance1(v1);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2) 
    {
        return v1.distance(v2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance(const Field_Vector_3D<T> v1, const Vector_3D v2) 
    {
        return v1.distance(v2);
    }

    /** Compute the distance between two vectors according to the L<sub>2</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>2</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance(const Vector_3D v1, const Field_Vector_3D<T> v2) 
    {
        return v2.distance(v1);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>&infin;</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance_inf(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2) 
    {
        return v1.distance_inf(v2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>&infin;</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance_inf(const Field_Vector_3D<T> v1, const Vector_3D v2) 
    {
        return v1.distance_inf(v2);
    }

    /** Compute the distance between two vectors according to the L<sub>&infin;</sub> norm.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm_inf()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the distance between v1 and v2 according to the L<sub>&infin;</sub> norm
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance_inf(const Vector_3D v1, const Field_Vector_3D<T> v2) 
    {
        return v2.distance_inf(v1);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between v1 and v2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance_sq(const Field_Vector_3D<T> v1, const Field_Vector_3D<T> v2) 
    {
        return v1.distance_sq(v2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between v1 and v2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance_sq(const Field_Vector_3D<T> v1, const Vector_3D v2) 
    {
        return v1.distance_sq(v2);
    }

    /** Compute the square of the distance between two vectors.
     * <p>Calling this method is equivalent to calling:
     * <code>v1.subtract(v2).get_norm_sq()</code> except that no intermediate
     * vector is built</p>
     * @param v1 first vector
     * @param v2 second vector
     * @param <T> the type of the field elements
     * @return the square of the distance between v1 and v2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance_sq(const Vector_3D v1, const Field_Vector_3D<T> v2) 
    {
        return v2.distance_sq(v1);
    }

    /** Get a string representation of this vector.
     * @return a string representation of this vector
     */
    //override
    std::string to_string() const 
    {
        return Vector_3DFormat.get_vector_3d_format().format(to_vector_3d());
    }

    /** Get a string representation of this vector.
     * @param format the custom format for components
     * @return a string representation of this vector
     */
    std::string to_string(const Number_Format format) 
    {
        return Vector_3DFormat(format).format(to_vector_3d());
    }
};