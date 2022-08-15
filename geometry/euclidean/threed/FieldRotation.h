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

//import org.hipparchus.Calculus_Field_Element;
//import org.hipparchus.Field;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Localized_Geometry_Formats;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Field_Sin_Cos;
//import org.hipparchus.util.Math_Arrays;
#include <vector>
#include <type_traits>
#include "../../../core/CalculusFieldElement.hpp"

/**
 * This class is a re-implementation of {@link Rotation} using {@link Calculus_Field_Element}.
 * <p>Instance of this class are guaranteed to be immutable.</p>
 *
 * @param <T> the type of the field elements
 * @see Field_Vector_3D
 * @see Rotation_Order
 */
template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
class Field_Rotation  
{

private:

    /** Scalar coordinate of the quaternion. */
    const T my_q_0;

    /** First coordinate of the vectorial part of the quaternion. */
    const T my_q_1;

    /** Second coordinate of the vectorial part of the quaternion. */
    const T my_q_2;

    /** Third coordinate of the vectorial part of the quaternion. */
    const T my_q_3;

public:
    /** Build a rotation from the quaternion coordinates.
     * <p>A rotation can be built from a <em>normalized</em> quaternion, * i.e. a quaternion for which q<sub>0</sub><sup>2</sup> +
     * q<sub>1</sub><sup>2</sup> + q<sub>2</sub><sup>2</sup> +
     * q<sub>3</sub><sup>2</sup> = 1. If the quaternion is not normalized, * the constructor can normalize it in a preprocessing step.</p>
     * <p>Note that some conventions put the scalar part of the quaternion
     * as the 4<sup>th</sup> component and the vector part as the first three
     * components. This is <em>not</em> our convention. We put the scalar part
     * as the first component.</p>
     * @param my_q_0 scalar part of the quaternion
     * @param my_q_1 first coordinate of the vectorial part of the quaternion
     * @param my_q_2 second coordinate of the vectorial part of the quaternion
     * @param my_q_3 third coordinate of the vectorial part of the quaternion
     * @param needs_normalization if true, the coordinates are considered
     * not to be normalized, a normalization preprocessing step is performed
     * before using them
     */
    Field_Rotation(const T& q_0, const T& q_1, const T& q_2, const T& q_3, const bool needs_normalization) 
    {

        if (needs_normalization) 
        {
            // normalization preprocessing
            const T inv = q_0.multiply(q_0).add(q_1.multiply(q_1)).add(q_2.multiply(q_2)).add(q_3.multiply(q_3)).sqrt().reciprocal();
            my_q_0 = inv.multiply(q_0);
            my_q_1 = inv.multiply(q_1);
            my_q_2 = inv.multiply(q_2);
            my_q_3 = inv.multiply(q_3);
        }
        else 
        {
            my_q_0 = q_0;
            my_q_1 = q_1;
            my_q_2 = q_2;
            my_q_3 = q_3;
        }
    }

    /** Build a rotation from an axis and an angle.
     * <p>We use the convention that angles are oriented according to
     * the effect of the rotation on vectors around the axis. That means
     * that if (i, j, k) is a direct frame and if we first provide +k as
     * the axis and &pi;/2 as the angle to this constructor, and then
     * {@link #apply_to(Field_Vector_3D) apply} the instance to +i, we will get
     * +j.</p>
     * <p>Another way to represent our convention is to say that a rotation
     * of angle &theta; about the unit vector (x, y, z) is the same as the
     * rotation build from quaternion components { cos(-&theta;/2), * x * sin(-&theta;/2), y * sin(-&theta;/2), z * sin(-&theta;/2) }.
     * Note the minus sign on the angle!</p>
     * <p>On the one hand this convention is consistent with a vectorial
     * perspective (moving vectors in fixed frames), on the other hand it
     * is different from conventions with a frame perspective (fixed vectors
     * viewed from different frames) like the ones used for example in spacecraft
     * attitude community or in the graphics community.</p>
     * @param axis axis around which to rotate
     * @param angle rotation angle.
     * @param convention convention to use for the semantics of the angle
     * @exception  if the axis norm is zero
     */
    Field_Rotation(const Field_Vector_3D<T>& axis, const T& angle, const Rotation_Convention& convention)
    {
        const T norm = axis.get_norm();
        if (norm.get_real() == 0) 
        {
            throw (Localized_Geometry_Formats.ZERO_NORM_FOR_ROTATION_AXIS);
        }

        const T half_angle = angle.multiply(convention == Rotation_Convention.VECTOR_OPERATOR ? -0.5 : 0.5);
        const Field_Sin_Cos<T> sin_cos = Sin_Cos(half_angle);
        const T coeff = sin_cos.sin().divide(norm);

        my_q_0 = sin_cos.cos();
        my_q_1 = coeff.multiply(axis.get_x());
        my_q_2 = coeff.multiply(axis.get_y());
        my_q_3 = coeff.multiply(axis.get_z());

    }

    /** Build a {@link Field_Rotation} from a {@link Rotation}.
     * @param field field for the components
     * @param r rotation to convert
     */
    Field_Rotation(const Field<T>& field, const Rotation& r) 
    {
        my_q_0 = field.get_zero().add(r.get_my_q_0());
        my_q_1 = field.get_zero().add(r.get_my_q_1());
        my_q_2 = field.get_zero().add(r.get_my_q_2());
        my_q_3 = field.get_zero().add(r.get_my_q_3());
    }

    /** Build a rotation from a 3X3 matrix.

     * <p>Rotation matrices are orthogonal matrices, i.e. unit matrices
     * (which are matrices for which m.m<sup>T</sup> = I) with real
     * coefficients. The module of the determinant of unit matrices is
     * 1, among the orthogonal 3X3 matrices, only the ones having a
     * positive determinant (+1) are rotation matrices.</p>

     * <p>When a rotation is defined by a matrix with truncated values
     * (typically when it is extracted from a technical sheet where only
     * four to five significant digits are available), the matrix is not
     * orthogonal anymore. This constructor handles this case
     * transparently by using a copy of the given matrix and applying a
     * correction to the copy in order to perfect its orthogonality. If
     * the Frobenius norm of the correction needed is above the given
     * threshold, then the matrix is considered to be too far from a
     * true rotation matrix and an exception is thrown.<p>

     * @param m rotation matrix
     * @param threshold convergence threshold for the iterative
     * orthogonality correction (convergence is reached when the
     * difference between two steps of the Frobenius norm of the
     * correction is below this threshold)

     * @exception  if the matrix is not a 3X3
     * matrix, or if it cannot be transformed into an orthogonal matrix
     * with the given threshold, or if the determinant of the resulting
     * orthogonal matrix is negative

     */
    Field_Rotation(const std::vector<std::vector<T>>& m, const double& threshold)
    {
        // dimension check
        if ((m.size() != 3) || (m[0].size() != 3) || (m[1].size() != 3) || (m[2].size() != 3)) 
        {
            throw (Localized_Geometry_Formats.ROTATION_MATRIX_DIMENSIONS, m.size(), m[0].size());
        }

        // compute a "close" orthogonal matrix
        const std::vector<std::vector<T>> ort = orthogonalize_matrix(m, threshold);

        // check the sign of the determinant
        const T d0 = ort[1][1].multiply(ort[2][2]).subtract(ort[2][1].multiply(ort[1][2]));
        const T d1 = ort[0][1].multiply(ort[2][2]).subtract(ort[2][1].multiply(ort[0][2]));
        const T d2 = ort[0][1].multiply(ort[1][2]).subtract(ort[1][1].multiply(ort[0][2]));
        const T det =
                ort[0][0].multiply(d0).subtract(ort[1][0].multiply(d1)).add(ort[2][0].multiply(d2));
        if (det.get_real() < 0.0) 
        {
            throw (Localized_Geometry_Formats.CLOSEST_ORTHOGONAL_MATRIX_HAS_NEGATIVE_DETERMINANT, det);
        }

        const std::vector<T> quat = mat2quat(ort);
        my_q_0 = quat[0];
        my_q_1 = quat[1];
        my_q_2 = quat[2];
        my_q_3 = quat[3];

    }

    /** Build the rotation that transforms a pair of vectors into another pair.

     * <p>Except for possible scale factors, if the instance were applied to
     * the pair (u<sub>1</sub>, u<sub>2</sub>) it will produce the pair
     * (v<sub>1</sub>, v<sub>2</sub>).</p>

     * <p>If the angular separation between u<sub>1</sub> and u<sub>2</sub> is
     * not the same as the angular separation between v<sub>1</sub> and
     * v<sub>2</sub>, then a corrected v'<sub>2</sub> will be used rather than
     * v<sub>2</sub>, the corrected vector will be in the (&pm;v<sub>1</sub>, * +v<sub>2</sub>) half-plane.</p>

     * @param u1 first vector of the origin pair
     * @param u2 second vector of the origin pair
     * @param v1 desired image of u1 by the rotation
     * @param v2 desired image of u2 by the rotation
     * @exception Math_Runtime_Exception if the norm of one of the vectors is zero, * or if one of the pair is degenerated (i.e. the vectors of the pair are collinear)
     */
    Field_Rotation(Field_Vector_3D<T> u1, Field_Vector_3D<T> u2, Field_Vector_3D<T> v1, Field_Vector_3D<T> v2)
    {
        // build orthonormalized base from u1, u2
        // this fails when vectors are NULL or collinear, which is forbidden to define a rotation
        const Field_Vector_3D<T> u3 = Field_Vector_3D.cross_product(u1, u2).normalize();
        u2 = Field_Vector_3D.cross_product(u3, u1).normalize();
        u1 = u1.normalize();

        // build an orthonormalized base from v1, v2
        // this fails when vectors are NULL or collinear, which is forbidden to define a rotation
        const Field_Vector_3D<T> v3 = Field_Vector_3D.cross_product(v1, v2).normalize();
        v2 = Field_Vector_3D.cross_product(v3, v1).normalize();
        v1 = v1.normalize();

        // buid a matrix transforming the first base into the second one
        const std::vector<std::vector<T>> array = Math_Arrays::build_array(u1.get_x().get_field(), 3, 3);
        array[0][0] = u1.get_x().multiply(v1.get_x()).add(u2.get_x().multiply(v2.get_x())).add(u3.get_x().multiply(v3.get_x()));
        array[0][1] = u1.get_y().multiply(v1.get_x()).add(u2.get_y().multiply(v2.get_x())).add(u3.get_y().multiply(v3.get_x()));
        array[0][2] = u1.get_z().multiply(v1.get_x()).add(u2.get_z().multiply(v2.get_x())).add(u3.get_z().multiply(v3.get_x()));
        array[1][0] = u1.get_x().multiply(v1.get_y()).add(u2.get_x().multiply(v2.get_y())).add(u3.get_x().multiply(v3.get_y()));
        array[1][1] = u1.get_y().multiply(v1.get_y()).add(u2.get_y().multiply(v2.get_y())).add(u3.get_y().multiply(v3.get_y()));
        array[1][2] = u1.get_z().multiply(v1.get_y()).add(u2.get_z().multiply(v2.get_y())).add(u3.get_z().multiply(v3.get_y()));
        array[2][0] = u1.get_x().multiply(v1.get_z()).add(u2.get_x().multiply(v2.get_z())).add(u3.get_x().multiply(v3.get_z()));
        array[2][1] = u1.get_y().multiply(v1.get_z()).add(u2.get_y().multiply(v2.get_z())).add(u3.get_y().multiply(v3.get_z()));
        array[2][2] = u1.get_z().multiply(v1.get_z()).add(u2.get_z().multiply(v2.get_z())).add(u3.get_z().multiply(v3.get_z()));

        std::vector<T> quat = mat2quat(array);
        my_q_0 = quat[0];
        my_q_1 = quat[1];
        my_q_2 = quat[2];
        my_q_3 = quat[3];
    }

    /** Build one of the rotations that transform one vector into another one.

     * <p>Except for a possible scale factor, if the instance were
     * applied to the vector u it will produce the vector v. There is an
     * infinite number of such rotations, this constructor choose the
     * one with the smallest associated angle (i.e. the one whose axis
     * is orthogonal to the (u, v) plane). If u and v are collinear, an
     * arbitrary rotation axis is chosen.</p>

     * @param u origin vector
     * @param v desired image of u by the rotation
     * @exception Math_Runtime_Exception if the norm of one of the vectors is zero
     */
    public Field_Rotation(const Field_Vector_3D<T>& u, const Field_Vector_3D<T>& v) 
    {
        const T norm_product = u.get_norm().multiply(v.get_norm());
        if (norm_product.get_real() == 0) 
        {
            throw Math_Runtime_Exception(Localized_Geometry_Formats.ZERO_NORM_FOR_ROTATION_DEFINING_VECTOR);
        }

        const T dot = Field_Vector_3D.dot_product(u, v);

        if (dot.get_real() < ((2.0e-15 - 1.0) * norm_product.get_real())) 
        {
            // special case u = -v: we select a PI angle rotation around
            // an arbitrary vector orthogonal to u
            const Field_Vector_3D<T> w = u.orthogonal();
            my_q_0 = norm_product.get_field().get_zero();
            my_q_1 = w.get_x().negate();
            my_q_2 = w.get_y().negate();
            my_q_3 = w.get_z().negate();
        }
        else 
        {
            // general case: (u, v) defines a plane, we select
            // the shortest possible rotation: axis orthogonal to this plane
            my_q_0 = dot.divide(norm_product).add(1.0).multiply(0.5).sqrt();
            const T coeff = my_q_0.multiply(norm_product).multiply(2.0).reciprocal();
            const Field_Vector_3D<T> q = Field_Vector_3D.cross_product(v, u);
            my_q_1 = coeff.multiply(q.get_x());
            my_q_2 = coeff.multiply(q.get_y());
            my_q_3 = coeff.multiply(q.get_z());
        }

    }

    /** Build a rotation from three Cardan or Euler elementary rotations.

     * <p>Cardan rotations are three successive rotations around the
     * canonical axes X, Y and Z, each axis being used once. There are
     * 6 such sets of rotations (XYZ, XZY, YXZ, YZX, ZXY and ZYX). Euler
     * rotations are three successive rotations around the canonical
     * axes X, Y and Z, the first and last rotations being around the
     * same axis. There are 6 such sets of rotations (XYX, XZX, YXY, * YZY, ZXZ and ZYZ), the most popular one being ZXZ.</p>
     * <p>Beware that many people routinely use the term Euler angles even
     * for what really are Cardan angles (this confusion is especially
     * widespread in the aerospace business where Roll, Pitch and Yaw angles
     * are often wrongly tagged as Euler angles).</p>

     * @param order order of rotations to compose, from left to right
     * (i.e. we will use {@code r1.compose(r2.compose(r3, convention), convention)})
     * @param convention convention to use for the semantics of the angle
     * @param alpha1 angle of the first elementary rotation
     * @param alpha2 angle of the second elementary rotation
     * @param alpha3 angle of the third elementary rotation
     */
    Field_Rotation(const Rotation_Order& order, const Rotation_Convention& convention, const T& alpha1, const T& alpha2, const T& alpha3) 
    {
        const Field<T> field = alpha1.get_field();
        const Field_Rotation<T> r1 = Field_Rotation<>(new Field_Vector_3D<>(field, order.get_a1()), alpha1, convention);
        const Field_Rotation<T> r2 = Field_Rotation<>(new Field_Vector_3D<>(field, order.get_a2()), alpha2, convention);
        const Field_Rotation<T> r3 = Field_Rotation<>(new Field_Vector_3D<>(field, order.get_a3()), alpha3, convention);
        const Field_Rotation<T> composed = r1.compose(r2.compose(r3, convention), convention);
        my_q_0 = composed.my_q_0;
        my_q_1 = composed.my_q_1;
        my_q_2 = composed.my_q_2;
        my_q_3 = composed.my_q_3;
    }

    /** Get identity rotation.
     * @param field field for the components
     * @return a rotation
     * @param <T> the type of the field elements
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Rotation<T> get_identity(const Field<T> field) 
    {
        return Field_Rotation<>(field, Rotation.IDENTITY);
    }

    /** Apply the instance to another rotation.
 * <p>
 * Calling this method is equivalent to call
 * {@link #compose(Rotation, Rotation_Convention)
 * compose(r, Rotation_Convention.VECTOR_OPERATOR)}.
 * </p>
 * @param r rotation to apply the rotation to
 * @return a rotation which is the composition of r by the instance
 */
    Field_Rotation<T> apply_to(const Rotation& r)
    {
        return compose(r, Rotation_Convention.VECTOR_OPERATOR);
    }

    /** Compose the instance with another rotation.
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#VECTOR_OPERATOR vector operator} convention, * applying the instance to a rotation is computing the composition
     * in an order compliant with the following rule : let {@code u} be any
     * vector and {@code v} its image by {@code r1} (i.e.
     * {@code r1.apply_to(u) = v}). Let {@code w} be the image of {@code v} by
     * rotation {@code r2} (i.e. {@code r2.apply_to(v) = w}). Then
     * {@code w = comp.apply_to(u)}, where
     * {@code comp = r2.compose(r1, Rotation_Convention.VECTOR_OPERATOR)}.
     * </p>
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#FRAME_TRANSFORM frame transform} convention, * the application order will be reversed. So keeping the exact same
     * meaning of all {@code r1}, {@code r2}, {@code u}, {@code v}, {@code w}
     * and  {@code comp} as above, {@code comp} could also be computed as
     * {@code comp = r1.compose(r2, Rotation_Convention.FRAME_TRANSFORM)}.
     * </p>
     * @param r rotation to apply the rotation to
     * @param convention convention to use for the semantics of the angle
     * @return a rotation which is the composition of r by the instance
     */
    Field_Rotation<T> compose(const Rotation& r, const Rotation_Convention& convention)
    {
        return convention == Rotation_Convention.VECTOR_OPERATOR
            ? compose_internal(r)
            : apply_to(r, this);
    }

    /** Apply a rotation to another rotation.
 * Applying a rotation to another rotation is computing the composition
 * in an order compliant with the following rule : let u be any
 * vector and v its image by r_inner (i.e. r_inner.apply_to(u) = v), let w be the image
 * of v by r_outer (i.e. r_outer.apply_to(v) = w), then w = comp.apply_to(u), * where comp = apply_to(r_outer, r_inner).
 * @param r1 rotation to apply
 * @param r_inner rotation to apply the rotation to
 * @param <T> the type of the field elements
 * @return a rotation which is the composition of r by the instance
 */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Rotation<T> apply_to(const Rotation r1, const Field_Rotation<T> r_inner)
    {
        return Field_Rotation<T>(r_inner.my_q_0.multiply(r1.get_my_q_0()).subtract(r_inner.my_q_1.multiply(r1.get_my_q_1()).add(r_inner.my_q_2.multiply(r1.get_my_q_2())).add(r_inner.my_q_3.multiply(r1.get_my_q_3()))), r_inner.my_q_1.multiply(r1.get_my_q_0()).add(r_inner.my_q_0.multiply(r1.get_my_q_1())).add(r_inner.my_q_2.multiply(r1.get_my_q_3()).subtract(r_inner.my_q_3.multiply(r1.get_my_q_2()))), r_inner.my_q_2.multiply(r1.get_my_q_0()).add(r_inner.my_q_0.multiply(r1.get_my_q_2())).add(r_inner.my_q_3.multiply(r1.get_my_q_1()).subtract(r_inner.my_q_1.multiply(r1.get_my_q_3()))), r_inner.my_q_3.multiply(r1.get_my_q_0()).add(r_inner.my_q_0.multiply(r1.get_my_q_3())).add(r_inner.my_q_1.multiply(r1.get_my_q_2()).subtract(r_inner.my_q_2.multiply(r1.get_my_q_1()))), false);
    }

    /** Apply the inverse of the instance to another rotation.
     * <p>
     * Calling this method is equivalent to call
     * {@link #compose_inverse(Field_Rotation, Rotation_Convention)
     * compose_inverse(r, Rotation_Convention.VECTOR_OPERATOR)}.
     * </p>
     * @param r rotation to apply the rotation to
     * @return a rotation which is the composition of r by the inverse
     * of the instance
     */
    Field_Rotation<T> apply_inverse_to(const Field_Rotation<T>& r)
    {
        return compose_inverse(r, Rotation_Convention.VECTOR_OPERATOR);
    }

    /** Compose the inverse of the instance with another rotation.
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#VECTOR_OPERATOR vector operator} convention, * applying the inverse of the instance to a rotation is computing
     * the composition in an order compliant with the following rule :
     * let {@code u} be any vector and {@code v} its image by {@code r1}
     * (i.e. {@code r1.apply_to(u) = v}). Let {@code w} be the inverse image
     * of {@code v} by {@code r2} (i.e. {@code r2.apply_inverse_to(v) = w}).
     * Then {@code w = comp.apply_to(u)}, where
     * {@code comp = r2.compose_inverse(r1)}.
     * </p>
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#FRAME_TRANSFORM frame transform} convention, * the application order will be reversed, which means it is the
     * <em>innermost</em> rotation that will be reversed. So keeping the exact same
     * meaning of all {@code r1}, {@code r2}, {@code u}, {@code v}, {@code w}
     * and  {@code comp} as above, {@code comp} could also be computed as
     * {@code comp = r1.revert().compose_inverse(r2.revert(), Rotation_Convention.FRAME_TRANSFORM)}.
     * </p>
     * @param r rotation to apply the rotation to
     * @param convention convention to use for the semantics of the angle
     * @return a rotation which is the composition of r by the inverse
     * of the instance
     */
    Field_Rotation<T> compose_inverse(const Field_Rotation<T>& r, const Rotation_Convention& convention)
    {
        return convention == Rotation_Convention.VECTOR_OPERATOR
            ? compose_inverse_internal(r)
            : r.compose_internal(revert());
    }

    /** Apply the inverse of the instance to another rotation.
 * <p>
 * Calling this method is equivalent to call
 * {@link #compose_inverse(Rotation, Rotation_Convention)
 * compose_inverse(r, Rotation_Convention.VECTOR_OPERATOR)}.
 * </p>
 * @param r rotation to apply the rotation to
 * @return a rotation which is the composition of r by the inverse
 * of the instance
 */
    Field_Rotation<T> apply_inverse_to(const Rotation& r)
    {
        return compose_inverse(r, Rotation_Convention.VECTOR_OPERATOR);
    }

    /** Compose the inverse of the instance with another rotation.
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#VECTOR_OPERATOR vector operator} convention, * applying the inverse of the instance to a rotation is computing
     * the composition in an order compliant with the following rule :
     * let {@code u} be any vector and {@code v} its image by {@code r1}
     * (i.e. {@code r1.apply_to(u) = v}). Let {@code w} be the inverse image
     * of {@code v} by {@code r2} (i.e. {@code r2.apply_inverse_to(v) = w}).
     * Then {@code w = comp.apply_to(u)}, where
     * {@code comp = r2.compose_inverse(r1)}.
     * </p>
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#FRAME_TRANSFORM frame transform} convention, * the application order will be reversed, which means it is the
     * <em>innermost</em> rotation that will be reversed. So keeping the exact same
     * meaning of all {@code r1}, {@code r2}, {@code u}, {@code v}, {@code w}
     * and  {@code comp} as above, {@code comp} could also be computed as
     * {@code comp = r1.revert().compose_inverse(r2.revert(), Rotation_Convention.FRAME_TRANSFORM)}.
     * </p>
     * @param r rotation to apply the rotation to
     * @param convention convention to use for the semantics of the angle
     * @return a rotation which is the composition of r by the inverse
     * of the instance
     */
    Field_Rotation<T> compose_inverse(const Rotation& r, const Rotation_Convention& convention)
    {
        return convention == Rotation_Convention.VECTOR_OPERATOR
            ? compose_inverse_internal(r)
            : apply_to(r, revert());
    }

    /** Apply the inverse of a rotation to another rotation.
     * Applying the inverse of a rotation to another rotation is computing
     * the composition in an order compliant with the following rule :
     * let u be any vector and v its image by r_inner (i.e. r_inner.apply_to(u) = v), * let w be the inverse image of v by r_outer
     * (i.e. r_outer.apply_inverse_to(v) = w), then w = comp.apply_to(u), where
     * comp = apply_inverse_to(r_outer, r_inner).
     * @param r_outer rotation to apply the rotation to
     * @param r_inner rotation to apply the rotation to
     * @param <T> the type of the field elements
     * @return a rotation which is the composition of r by the inverse
     * of the instance
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Rotation<T> apply_inverse_to(const Rotation& r_outer, const Field_Rotation<T>& r_inner)
    {
        return Field_Rotation<T>(r_inner.my_q_0.multiply(r_outer.get_my_q_0()).add(r_inner.my_q_1.multiply(r_outer.get_my_q_1()).add(r_inner.my_q_2.multiply(r_outer.get_my_q_2())).add(r_inner.my_q_3.multiply(r_outer.get_my_q_3()))).negate(), r_inner.my_q_0.multiply(r_outer.get_my_q_1()).add(r_inner.my_q_2.multiply(r_outer.get_my_q_3()).subtract(r_inner.my_q_3.multiply(r_outer.get_my_q_2()))).subtract(r_inner.my_q_1.multiply(r_outer.get_my_q_0())), r_inner.my_q_0.multiply(r_outer.get_my_q_2()).add(r_inner.my_q_3.multiply(r_outer.get_my_q_1()).subtract(r_inner.my_q_1.multiply(r_outer.get_my_q_3()))).subtract(r_inner.my_q_2.multiply(r_outer.get_my_q_0())), r_inner.my_q_0.multiply(r_outer.get_my_q_3()).add(r_inner.my_q_1.multiply(r_outer.get_my_q_2()).subtract(r_inner.my_q_2.multiply(r_outer.get_my_q_1()))).subtract(r_inner.my_q_3.multiply(r_outer.get_my_q_0())), false);
    }

    /** Compute the <i>distance</i> between two rotations.
     * <p>The <i>distance</i> is intended here as a way to check if two
     * rotations are almost similar (i.e. they transform vectors the same way)
     * or very different. It is mathematically defined as the angle of
     * the rotation r that prepended to one of the rotations gives the other
     * one:</p>
     * <pre>
     *        r<sub>1</sub>(r) = r<sub>2</sub>
     * </pre>
     * <p>This distance is an angle between 0 and &pi;. Its value is the smallest
     * possible upper bound of the angle in radians between r<sub>1</sub>(v)
     * and r<sub>2</sub>(v) for all possible vectors v. This upper bound is
     * reached for some v. The distance is equal to 0 if and only if the two
     * rotations are identical.</p>
     * <p>Comparing two rotations should always be done using this value rather
     * than for example comparing the components of the quaternions. It is much
     * more stable, and has a geometric meaning. Also comparing quaternions
     * components is error prone since for example quaternions (0.36, 0.48, -0.48, -0.64)
     * and (-0.36, -0.48, 0.48, 0.64) represent exactly the same rotation despite
     * their components are different (they are exact opposites).</p>
     * @param r1 first rotation
     * @param r2 second rotation
     * @param <T> the type of the field elements
     * @return <i>distance</i> between r1 and r2
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  T distance(const Field_Rotation<T> r1, const Field_Rotation<T> r2)
    {
        return r1.compose_inverse_internal(r2).get_angle();
    }

    /** Revert a rotation.
     * Build a rotation which reverse the effect of another
     * rotation. This means that if r(u) = v, then r.revert(v) = u. The
     * instance is not changed.
     * @return a rotation whose effect is the reverse of the effect
     * of the instance
     */
    Field_Rotation<T> revert() 
    {
        return Field_Rotation<T>(my_q_0.negate(), my_q_1, my_q_2, my_q_3, false);
    }

    /** Get the scalar coordinate of the quaternion.
     * @return scalar coordinate of the quaternion
     */
    T get_my_q_0() const
    {
        return my_q_0;
    }

    /** Get the first coordinate of the vectorial part of the quaternion.
     * @return first coordinate of the vectorial part of the quaternion
     */
    T get_my_q_1() const
    {
        return my_q_1;
    }

    /** Get the second coordinate of the vectorial part of the quaternion.
     * @return second coordinate of the vectorial part of the quaternion
     */
    T get_my_q_2() const
    {
        return my_q_2;
    }

    /** Get the third coordinate of the vectorial part of the quaternion.
     * @return third coordinate of the vectorial part of the quaternion
     */
    T get_my_q_3() const
    {
        return my_q_3;
    }

    /** Get the normalized axis of the rotation.
     * <p>
     * Note that as {@link #get_angle()} always returns an angle
     * between 0 and &pi;, changing the convention changes the
     * direction of the axis, not the sign of the angle.
     * </p>
     * @param convention convention to use for the semantics of the angle
     * @return normalized axis of the rotation
     * @see #Field_Rotation(Field_Vector_3D, Calculus_Field_Element, Rotation_Convention)
     */
    Field_Vector_3D<T> get_axis(const Rotation_Convention convention) 
    {
        const T squared_sine = my_q_1.multiply(my_q_1).add(my_q_2.multiply(my_q_2)).add(my_q_3.multiply(my_q_3));
        if (squared_sine.get_real() == 0) 
        {
            const Field<T> field = squared_sine.get_field();
            return Field_Vector_3D<T>(convention == Rotation_Convention.VECTOR_OPERATOR ? field.get_one(): field.get_one().negate(), field.get_zero(), field.get_zero());
        }
        else 
        {
            const double sgn = convention == Rotation_Convention.VECTOR_OPERATOR ? +1 : -1;
            if (my_q_0.get_real() < 0) 
            {
                T inverse = squared_sine.sqrt().reciprocal().multiply(sgn);
                return Field_Vector_3D<T>(my_q_1.multiply(inverse), my_q_2.multiply(inverse), my_q_3.multiply(inverse));
            }
            const T inverse = squared_sine.sqrt().reciprocal().negate().multiply(sgn);
            return Field_Vector_3D<T>(my_q_1.multiply(inverse), my_q_2.multiply(inverse), my_q_3.multiply(inverse));
        }
    }

    /** Get the angle of the rotation.
     * @return angle of the rotation (between 0 and &pi;)
     * @see #Field_Rotation(Field_Vector_3D, Calculus_Field_Element, Rotation_Convention)
     */
    T get_angle() 
    {
        if ((my_q_0.get_real() < -0.1) || (my_q_0.get_real() > 0.1)) 
        {
            return my_q_1.multiply(my_q_1).add(my_q_2.multiply(my_q_2)).add(my_q_3.multiply(my_q_3)).sqrt().asin().multiply(2);
        }
        else if (my_q_0.get_real() < 0) 
        {
            return my_q_0.negate().acos().multiply(2);
        }
        return my_q_0.acos().multiply(2);
    }

    /** Get the Cardan or Euler angles corresponding to the instance.

     * <p>The equations show that each rotation can be defined by two
     * different values of the Cardan or Euler angles set. For example
     * if Cardan angles are used, the rotation defined by the angles
     * a<sub>1</sub>, a<sub>2</sub> and a<sub>3</sub> is the same as
     * the rotation defined by the angles &pi; + a<sub>1</sub>, &pi;
     * - a<sub>2</sub> and &pi; + a<sub>3</sub>. This method implements
     * the following arbitrary choices:</p>
     * <ul>
     *   <li>for Cardan angles, the chosen set is the one for which the
     *   second angle is between -&pi;/2 and &pi;/2 (i.e its cosine is
     *   positive),</li>
     *   <li>for Euler angles, the chosen set is the one for which the
     *   second angle is between 0 and &pi; (i.e its sine is positive).</li>
     * </ul>

     * <p>Cardan and Euler angle have a very disappointing drawback: all
     * of them have singularities. This means that if the instance is
     * too close to the singularities corresponding to the given
     * rotation order, it will be impossible to retrieve the angles. For
     * Cardan angles, this is often called gimbal lock. There is
     * <em>nothing</em> to do to prevent this, it is an intrinsic problem
     * with Cardan and Euler representation (but not a problem with the
     * rotation itself, which is perfectly well defined). For Cardan
     * angles, singularities occur when the second angle is close to
     * -&pi;/2 or +&pi;/2, for Euler angle singularities occur when the
     * second angle is close to 0 or &pi;, this implies that the identity
     * rotation is always singular for Euler angles!</p>

     * @param order rotation order to use
     * @param convention convention to use for the semantics of the angle
     * @return an array of three angles, in the order specified by the set
     * @exception Math_Illegal_State_Exception if the rotation is
     * singular with respect to the angles set specified
     */
    std::vector<T> get_angles(const Rotation_Order order, Rotation_Convention convention)
    {
        if (convention == Rotation_Convention.VECTOR_OPERATOR) 
        {
            if (order == Rotation_Order.XYZ) 
            {

                // r (+K) coordinates are :
                //  sin (theta), -cos (theta) sin (phi), cos (theta) cos (phi)
                // (-r) (+I) coordinates are :
                // cos (psi) cos (theta), -sin (psi) cos (theta), sin (theta)
                const // and we can choose to have theta in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(vector(0, 0, 1));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(1, 0, 0));
                if  ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_y().negate().atan2(v1.get_z()), v2.get_z().asin(), v2.get_y().negate().atan2(v2.get_x()));

            }
            else if (order == Rotation_Order.XZY) 
            {

                // r (+J) coordinates are :
                // -sin (psi), cos (psi) cos (phi), cos (psi) sin (phi)
                // (-r) (+I) coordinates are :
                // cos (theta) cos (psi), -sin (psi), sin (theta) cos (psi)
                // and we can choose to have psi in the interval [-PI/2 ; +PI/2]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 1, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(1, 0, 0));
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_z().atan2(v1.get_y()), v2.get_y().asin().negate(), v2.get_z().atan2(v2.get_x()));

            }
            else if (order == Rotation_Order.YXZ) 
            {

                // r (+K) coordinates are :
                //  cos (phi) sin (theta), -sin (phi), cos (phi) cos (theta)
                // (-r) (+J) coordinates are :
                // sin (psi) cos (phi), cos (psi) cos (phi), -sin (phi)
                // and we can choose to have phi in the interval [-PI/2 ; +PI/2]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 0, 1));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 1, 0));
                if ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_x().atan2(v1.get_z()), v2.get_z().asin().negate(), v2.get_x().atan2(v2.get_y()));

            }
            else if (order == Rotation_Order.YZX) 
            {

                // r (+I) coordinates are :
                // cos (psi) cos (theta), sin (psi), -cos (psi) sin (theta)
                // (-r) (+J) coordinates are :
                // sin (psi), cos (phi) cos (psi), -sin (phi) cos (psi)
                // and we can choose to have psi in the interval [-PI/2 ; +PI/2]
                const Field_Vector_3D<T> v1 = apply_to(vector(1, 0, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 1, 0));
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_z().negate().atan2(v1.get_x()), v2.get_x().asin(), v2.get_z().negate().atan2(v2.get_y()));

            }
            else if (order == Rotation_Order.ZXY) 
            {

                // r (+J) coordinates are :
                // -cos (phi) sin (psi), cos (phi) cos (psi), sin (phi)
                // (-r) (+K) coordinates are :
                // -sin (theta) cos (phi), sin (phi), cos (theta) cos (phi)
                // and we can choose to have phi in the interval [-PI/2 ; +PI/2]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 1, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 0, 1));
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_x().negate().atan2(v1.get_y()), v2.get_y().asin(), v2.get_x().negate().atan2(v2.get_z()));

            }
            else if (order == Rotation_Order.ZYX) 
            {

                // r (+I) coordinates are :
                //  cos (theta) cos (psi), cos (theta) sin (psi), -sin (theta)
                // (-r) (+K) coordinates are :
                // -sin (theta), sin (phi) cos (theta), cos (phi) cos (theta)
                // and we can choose to have theta in the interval [-PI/2 ; +PI/2]
                const Field_Vector_3D<T> v1 = apply_to(vector(1, 0, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 0, 1));
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_y().atan2(v1.get_x()), v2.get_x().asin().negate(), v2.get_y().atan2(v2.get_z()));

            }
            else if (order == Rotation_Order.XYX) 
            {

                // r (+I) coordinates are :
                //  cos (theta), sin (phi1) sin (theta), -cos (phi1) sin (theta)
                // (-r) (+I) coordinates are :
                // cos (theta), sin (theta) sin (phi2), sin (theta) cos (phi2)
                // and we can choose to have theta in the interval [0 ; PI]
                const Field_Vector_3D<T> v1 = apply_to(vector(1, 0, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(1, 0, 0));
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_y().atan2(v1.get_z().negate()), v2.get_x().acos(), v2.get_y().atan2(v2.get_z()));

            }
            else if (order == Rotation_Order.XZX) 
            {

                // r (+I) coordinates are :
                //  cos (psi), cos (phi1) sin (psi), sin (phi1) sin (psi)
                // (-r) (+I) coordinates are :
                // cos (psi), -sin (psi) cos (phi2), sin (psi) sin (phi2)
                // and we can choose to have psi in the interval [0 ; PI]
                const Field_Vector_3D<T> v1 = apply_to(vector(1, 0, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(1, 0, 0));
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_z().atan2(v1.get_y()), v2.get_x().acos(), v2.get_z().atan2(v2.get_y().negate()));

            }
            else if (order == Rotation_Order.YXY) 
            {

                // r (+J) coordinates are :
                //  sin (theta1) sin (phi), cos (phi), cos (theta1) sin (phi)
                // (-r) (+J) coordinates are :
                // sin (phi) sin (theta2), cos (phi), -sin (phi) cos (theta2)
                // and we can choose to have phi in the interval [0 ; PI]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 1, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 1, 0));
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_x().atan2(v1.get_z()), v2.get_y().acos(), v2.get_x().atan2(v2.get_z().negate()));

            }
            else if (order == Rotation_Order.YZY) 
            {

                // r (+J) coordinates are :
                //  -cos (theta1) sin (psi), cos (psi), sin (theta1) sin (psi)
                // (-r) (+J) coordinates are :
                // sin (psi) cos (theta2), cos (psi), sin (psi) sin (theta2)
                // and we can choose to have psi in the interval [0 ; PI]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 1, 0));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 1, 0));
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_z().atan2(v1.get_x().negate()), v2.get_y().acos(), v2.get_z().atan2(v2.get_x()));

            }
            else if (order == Rotation_Order.ZXZ) 
            {

                // r (+K) coordinates are :
                //  sin (psi1) sin (phi), -cos (psi1) sin (phi), cos (phi)
                // (-r) (+K) coordinates are :
                // sin (phi) sin (psi2), sin (phi) cos (psi2), cos (phi)
                // and we can choose to have phi in the interval [0 ; PI]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 0, 1));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 0, 1));
                if ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_x().atan2(v1.get_y().negate()), v2.get_z().acos(), v2.get_x().atan2(v2.get_y()));

            }
            else 
            { // last possibility is ZYZ

                // r (+K) coordinates are :
                //  cos (psi1) sin (theta), sin (psi1) sin (theta), cos (theta)
                // (-r) (+K) coordinates are :
                // -sin (theta) cos (psi2), sin (theta) sin (psi2), cos (theta)
                // and we can choose to have theta in the interval [0 ; PI]
                const Field_Vector_3D<T> v1 = apply_to(vector(0, 0, 1));
                const Field_Vector_3D<T> v2 = apply_inverse_to(vector(0, 0, 1));
                if ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v1.get_y().atan2(v1.get_x()), v2.get_z().acos(), v2.get_y().atan2(v2.get_x().negate()));

            }
        }
        else 
        {
            if (order == Rotation_Order.XYZ) 
            {

                // r (Vector_3D.plus_i) coordinates are :
                //  cos (theta) cos (psi), -cos (theta) sin (psi), sin (theta)
                // (-r) (Vector_3D.plus_k) coordinates are :
                // sin (theta), -sin (phi) cos (theta), cos (phi) cos (theta)
                // and we can choose to have theta in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_I);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_K);
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_y().negate().atan2(v2.get_z()), v2.get_x().asin(), v1.get_y().negate().atan2(v1.get_x()));

            }
            else if (order == Rotation_Order.XZY) 
            {

                // r (Vector_3D.plus_i) coordinates are :
                // cos (psi) cos (theta), -sin (psi), cos (psi) sin (theta)
                // (-r) (Vector_3D.plus_j) coordinates are :
                // -sin (psi), cos (phi) cos (psi), sin (phi) cos (psi)
                // and we can choose to have psi in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_I);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_J);
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_z().atan2(v2.get_y()), v2.get_x().asin().negate(), v1.get_z().atan2(v1.get_x()));

            }
            else if (order == Rotation_Order.YXZ) 
            {

                // r (Vector_3D.plus_j) coordinates are :
                // cos (phi) sin (psi), cos (phi) cos (psi), -sin (phi)
                // (-r) (Vector_3D.plus_k) coordinates are :
                // sin (theta) cos (phi), -sin (phi), cos (theta) cos (phi)
                // and we can choose to have phi in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_J);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_K);
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_x().atan2(v2.get_z()), v2.get_y().asin().negate(), v1.get_x().atan2(v1.get_y()));

            }
            else if (order == Rotation_Order.YZX) 
            {

                // r (Vector_3D.plus_j) coordinates are :
                // sin (psi), cos (psi) cos (phi), -cos (psi) sin (phi)
                // (-r) (Vector_3D.plus_i) coordinates are :
                // cos (theta) cos (psi), sin (psi), -sin (theta) cos (psi)
                // and we can choose to have psi in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_J);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_I);
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_z().negate().atan2(v2.get_x()), v2.get_y().asin(), v1.get_z().negate().atan2(v1.get_y()));

            }
            else if (order == Rotation_Order.ZXY) 
            {

                // r (Vector_3D.plus_k) coordinates are :
                //  -cos (phi) sin (theta), sin (phi), cos (phi) cos (theta)
                // (-r) (Vector_3D.plus_j) coordinates are :
                // -sin (psi) cos (phi), cos (psi) cos (phi), sin (phi)
                // and we can choose to have phi in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_K);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_J);
                if ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_x().negate().atan2(v2.get_y()), v2.get_z().asin(), v1.get_x().negate().atan2(v1.get_z()));

            }
            else if (order == Rotation_Order.ZYX) 
            {

                // r (Vector_3D.plus_k) coordinates are :
                //  -sin (theta), cos (theta) sin (phi), cos (theta) cos (phi)
                // (-r) (Vector_3D.plus_i) coordinates are :
                // cos (psi) cos (theta), sin (psi) cos (theta), -sin (theta)
                // and we can choose to have theta in the interval [-PI/2 ; +PI/2]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_K);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_I);
                if  ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_y().atan2(v2.get_x()), v2.get_z().asin().negate(), v1.get_y().atan2(v1.get_z()));

            }
            else if (order == Rotation_Order.XYX) 
            {

                // r (Vector_3D.plus_i) coordinates are :
                //  cos (theta), sin (phi2) sin (theta), cos (phi2) sin (theta)
                // (-r) (Vector_3D.plus_i) coordinates are :
                // cos (theta), sin (theta) sin (phi1), -sin (theta) cos (phi1)
                // and we can choose to have theta in the interval [0 ; PI]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_I);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_I);
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_y().atan2(v2.get_z().negate()), v2.get_x().acos(), v1.get_y().atan2(v1.get_z()));

            }
            else if (order == Rotation_Order.XZX) 
            {

                // r (Vector_3D.plus_i) coordinates are :
                //  cos (psi), -cos (phi2) sin (psi), sin (phi2) sin (psi)
                // (-r) (Vector_3D.plus_i) coordinates are :
                // cos (psi), sin (psi) cos (phi1), sin (psi) sin (phi1)
                // and we can choose to have psi in the interval [0 ; PI]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_I);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_I);
                if ((v2.get_x().get_real() < -0.9999999999) || (v2.get_x().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_z().atan2(v2.get_y()), v2.get_x().acos(), v1.get_z().atan2(v1.get_y().negate()));

            }
            else if (order == Rotation_Order.YXY) 
            {

                // r (Vector_3D.plus_j) coordinates are :
                // sin (phi) sin (theta2), cos (phi), -sin (phi) cos (theta2)
                // (-r) (Vector_3D.plus_j) coordinates are :
                //  sin (theta1) sin (phi), cos (phi), cos (theta1) sin (phi)
                // and we can choose to have phi in the interval [0 ; PI]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_J);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_J);
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_x().atan2(v2.get_z()), v2.get_y().acos(), v1.get_x().atan2(v1.get_z().negate()));

            }
            else if (order == Rotation_Order.YZY) 
            {

                // r (Vector_3D.plus_j) coordinates are :
                // sin (psi) cos (theta2), cos (psi), sin (psi) sin (theta2)
                // (-r) (Vector_3D.plus_j) coordinates are :
                //  -cos (theta1) sin (psi), cos (psi), sin (theta1) sin (psi)
                // and we can choose to have psi in the interval [0 ; PI]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_J);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_J);
                if ((v2.get_y().get_real() < -0.9999999999) || (v2.get_y().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_z().atan2(v2.get_x().negate()), v2.get_y().acos(), v1.get_z().atan2(v1.get_x()));

            }
            else if (order == Rotation_Order.ZXZ) 
            {

                // r (Vector_3D.plus_k) coordinates are :
                // sin (phi) sin (psi2), sin (phi) cos (psi2), cos (phi)
                // (-r) (Vector_3D.plus_k) coordinates are :
                //  sin (psi1) sin (phi), -cos (psi1) sin (phi), cos (phi)
                // and we can choose to have phi in the interval [0 ; PI]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_K);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_K);
                if ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_x().atan2(v2.get_y().negate()), v2.get_z().acos(), v1.get_x().atan2(v1.get_y()));

            }
            else
            { // last possibility is ZYZ

                // r (Vector_3D.plus_k) coordinates are :
                // -sin (theta) cos (psi2), sin (theta) sin (psi2), cos (theta)
                // (-r) (Vector_3D.plus_k) coordinates are :
                //  cos (psi1) sin (theta), sin (psi1) sin (theta), cos (theta)
                // and we can choose to have theta in the interval [0 ; PI]
                Field_Vector_3D<T> v1 = apply_to(Vector_3D.PLUS_K);
                Field_Vector_3D<T> v2 = apply_inverse_to(Vector_3D.PLUS_K);
                if ((v2.get_z().get_real() < -0.9999999999) || (v2.get_z().get_real() > 0.9999999999)) 
                {
                    throw Math_Illegal_State_Exception(Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
                }
                return build_array(v2.get_y().atan2(v2.get_x()), v2.get_z().acos(), v1.get_y().atan2(v1.get_x().negate()));

            }
        }
    }

    /** Get the 3X3 matrix corresponding to the instance
     * @return the matrix corresponding to the instance
     */
    std::vector<std::vector<T>> get_matrix()
    {
        // products
        const T my_q_0my_q_0 = my_q_0.multiply(my_q_0);
        const T my_q_0my_q_1 = my_q_0.multiply(my_q_1);
        const T my_q_0my_q_2 = my_q_0.multiply(my_q_2);
        const T my_q_0my_q_3 = my_q_0.multiply(my_q_3);
        const T my_q_1my_q_1 = my_q_1.multiply(my_q_1);
        const T my_q_1my_q_2 = my_q_1.multiply(my_q_2);
        const T my_q_1my_q_3 = my_q_1.multiply(my_q_3);
        const T my_q_2my_q_2 = my_q_2.multiply(my_q_2);
        const T my_q_2my_q_3 = my_q_2.multiply(my_q_3);
        const T my_q_3my_q_3 = my_q_3.multiply(my_q_3);

        // create the matrix
        const std::vector<std::vector<T>> m = Math_Arrays::build_array(my_q_0.get_field(), 3, 3);

        m[0][0] = my_q_0my_q_0.add(my_q_1my_q_1).multiply(2).subtract(1);
        m[1][0] = my_q_1my_q_2.subtract(my_q_0my_q_3).multiply(2);
        m[2][0] = my_q_1my_q_3.add(my_q_0my_q_2).multiply(2);

        m[0][1] = my_q_1my_q_2.add(my_q_0my_q_3).multiply(2);
        m[1][1] = my_q_0my_q_0.add(my_q_2my_q_2).multiply(2).subtract(1);
        m[2][1] = my_q_2my_q_3.subtract(my_q_0my_q_1).multiply(2);

        m[0][2] = my_q_1my_q_3.subtract(my_q_0my_q_2).multiply(2);
        m[1][2] = my_q_2my_q_3.add(my_q_0my_q_1).multiply(2);
        m[2][2] = my_q_0my_q_0.add(my_q_3my_q_3).multiply(2).subtract(1);

        return m;
    }

    /** Convert to a constant vector without derivatives.
     * @return a constant vector
     */
    Rotation to_rotation()
    {
        return Rotation(my_q_0.get_real(), my_q_1.get_real(), my_q_2.get_real(), my_q_3.get_real(), false);
    }

    /** Apply the rotation to a vector.
     * @param u vector to apply the rotation to
     * @return a vector which is the image of u by the rotation
     */
    Field_Vector_3D<T> apply_to(const Field_Vector_3D<T> u)
    {

        const T x = u.get_x();
        const T y = u.get_y();
        const T z = u.get_z();

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));

        return Field_Vector_3D<T>(my_q_0.multiply(x.multiply(my_q_0).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x), my_q_0.multiply(y.multiply(my_q_0).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y), my_q_0.multiply(z.multiply(my_q_0).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z));

    }

    /** Apply the rotation to a vector.
     * @param u vector to apply the rotation to
     * @return a vector which is the image of u by the rotation
     */
    Field_Vector_3D<T> apply_to(const Vector_3D& u)
    {
        const double x = u.get_x();
        const double y = u.get_y();
        const double z = u.get_z();

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));

        return Field_Vector_3D<T>(my_q_0.multiply(my_q_0.multiply(x).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x), my_q_0.multiply(my_q_0.multiply(y).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y), my_q_0.multiply(my_q_0.multiply(z).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z));
    }

    /** Apply the rotation to a vector stored in an array.
     * @param in an array with three items which stores vector to rotate
     * @param out an array with three items to put result to (it can be the same
     * array as in)
     */
    void apply_to(const std::vector<T> in, const std::vector<T> out)
    {

        const T x = in[0];
        const T y = in[1];
        const T z = in[2];

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));

        out[0] = my_q_0.multiply(x.multiply(my_q_0).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x);
        out[1] = my_q_0.multiply(y.multiply(my_q_0).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y);
        out[2] = my_q_0.multiply(z.multiply(my_q_0).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z);

    }

    /** Apply the rotation to a vector stored in an array.
     * @param in an array with three items which stores vector to rotate
     * @param out an array with three items to put result to
     */
    void apply_to(const std::vector<double> in, const std::vector<T> out)
    {

        const double x = in[0];
        const double y = in[1];
        const double z = in[2];

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));

        out[0] = my_q_0.multiply(my_q_0.multiply(x).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x);
        out[1] = my_q_0.multiply(my_q_0.multiply(y).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y);
        out[2] = my_q_0.multiply(my_q_0.multiply(z).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z);

    }

    /** Apply a rotation to a vector.
     * @param r rotation to apply
     * @param u vector to apply the rotation to
     * @param <T> the type of the field elements
     * @return a vector which is the image of u by the rotation
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> apply_to(const Rotation r, const Field_Vector_3D<T> u)
    {

        const T x = u.get_x();
        const T y = u.get_y();
        const T z = u.get_z();

        const T s = x.multiply(r.get_my_q_1()).add(y.multiply(r.get_my_q_2())).add(z.multiply(r.get_my_q_3()));

        return Field_Vector_3D<T>(x.multiply(r.get_my_q_0()).subtract(z.multiply(r.get_my_q_2()).subtract(y.multiply(r.get_my_q_3()))).multiply(r.get_my_q_0()).add(s.multiply(r.get_my_q_1())).multiply(2).subtract(x), y.multiply(r.get_my_q_0()).subtract(x.multiply(r.get_my_q_3()).subtract(z.multiply(r.get_my_q_1()))).multiply(r.get_my_q_0()).add(s.multiply(r.get_my_q_2())).multiply(2).subtract(y), z.multiply(r.get_my_q_0()).subtract(y.multiply(r.get_my_q_1()).subtract(x.multiply(r.get_my_q_2()))).multiply(r.get_my_q_0()).add(s.multiply(r.get_my_q_3())).multiply(2).subtract(z));

    }

    /** Apply the inverse of the rotation to a vector.
     * @param u vector to apply the inverse of the rotation to
     * @return a vector which such that u is its image by the rotation
     */
    Field_Vector_3D<T> apply_inverse_to(const Field_Vector_3D<T> u)
    {

        const T x = u.get_x();
        const T y = u.get_y();
        const T z = u.get_z();

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));
        const T m0 = my_q_0.negate();

        return Field_Vector_3D<T>(m0.multiply(x.multiply(m0).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x), m0.multiply(y.multiply(m0).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y), m0.multiply(z.multiply(m0).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z));

    }

    /** Apply the inverse of the rotation to a vector.
     * @param u vector to apply the inverse of the rotation to
     * @return a vector which such that u is its image by the rotation
     */
    Field_Vector_3D<T> apply_inverse_to(const Vector_3D& u)
    {

        const double x = u.get_x();
        const double y = u.get_y();
        const double z = u.get_z();

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));
        const T m0 = my_q_0.negate();

        return Field_Vector_3D<T>(m0.multiply(m0.multiply(x).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x), m0.multiply(m0.multiply(y).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y), m0.multiply(m0.multiply(z).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z));

    }

    /** Apply the inverse of the rotation to a vector stored in an array.
     * @param in an array with three items which stores vector to rotate
     * @param out an array with three items to put result to (it can be the same
     * array as in)
     */
    void apply_inverse_to(const std::vector<T> in, const std::vector<T> out)
    {

        const T x = in[0];
        const T y = in[1];
        const T z = in[2];

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));
        const T m0 = my_q_0.negate();

        out[0] = m0.multiply(x.multiply(m0).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x);
        out[1] = m0.multiply(y.multiply(m0).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y);
        out[2] = m0.multiply(z.multiply(m0).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z);

    }

    /** Apply the inverse of the rotation to a vector stored in an array.
     * @param in an array with three items which stores vector to rotate
     * @param out an array with three items to put result to
     */
    void apply_inverse_to(const std::vector<double> in, const std::vector<T> out)
    {

        const double x = in[0];
        const double y = in[1];
        const double z = in[2];

        const T s = my_q_1.multiply(x).add(my_q_2.multiply(y)).add(my_q_3.multiply(z));
        const T m0 = my_q_0.negate();

        out[0] = m0.multiply(m0.multiply(x).subtract(my_q_2.multiply(z).subtract(my_q_3.multiply(y)))).add(s.multiply(my_q_1)).multiply(2).subtract(x);
        out[1] = m0.multiply(m0.multiply(y).subtract(my_q_3.multiply(x).subtract(my_q_1.multiply(z)))).add(s.multiply(my_q_2)).multiply(2).subtract(y);
        out[2] = m0.multiply(m0.multiply(z).subtract(my_q_1.multiply(y).subtract(my_q_2.multiply(x)))).add(s.multiply(my_q_3)).multiply(2).subtract(z);

    }

    /** Apply the inverse of a rotation to a vector.
     * @param r rotation to apply
     * @param u vector to apply the inverse of the rotation to
     * @param <T> the type of the field elements
     * @return a vector which such that u is its image by the rotation
     */
    template<typename T, typename std::enable_if<std::is_base_of<Calculus_Field_Element<T>, T>::value>::type* = nullptr>
    static  Field_Vector_3D<T> apply_inverse_to(const Rotation r, const Field_Vector_3D<T> u)
    {

        const T x = u.get_x();
        const T y = u.get_y();
        const T z = u.get_z();

        const T s = x.multiply(r.get_my_q_1()).add(y.multiply(r.get_my_q_2())).add(z.multiply(r.get_my_q_3()));
        const double m0 = -r.get_my_q_0();

        return Field_Vector_3D<T>(x.multiply(m0).subtract(z.multiply(r.get_my_q_2()).subtract(y.multiply(r.get_my_q_3()))).multiply(m0).add(s.multiply(r.get_my_q_1())).multiply(2).subtract(x), y.multiply(m0).subtract(x.multiply(r.get_my_q_3()).subtract(z.multiply(r.get_my_q_1()))).multiply(m0).add(s.multiply(r.get_my_q_2())).multiply(2).subtract(y), z.multiply(m0).subtract(y.multiply(r.get_my_q_1()).subtract(x.multiply(r.get_my_q_2()))).multiply(m0).add(s.multiply(r.get_my_q_3())).multiply(2).subtract(z));
    }

    /** Apply the instance to another rotation.
     * <p>
     * Calling this method is equivalent to call
     * {@link #compose(Field_Rotation, Rotation_Convention)
     * compose(r, Rotation_Convention.VECTOR_OPERATOR)}.
     * </p>
     * @param r rotation to apply the rotation to
     * @return a rotation which is the composition of r by the instance
     */
    Field_Rotation<T> apply_to(const Field_Rotation<T> r)
    {
        return compose(r, Rotation_Convention.VECTOR_OPERATOR);
    }

    /** Compose the instance with another rotation.
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#VECTOR_OPERATOR vector operator} convention, * applying the instance to a rotation is computing the composition
     * in an order compliant with the following rule : let {@code u} be any
     * vector and {@code v} its image by {@code r1} (i.e.
     * {@code r1.apply_to(u) = v}). Let {@code w} be the image of {@code v} by
     * rotation {@code r2} (i.e. {@code r2.apply_to(v) = w}). Then
     * {@code w = comp.apply_to(u)}, where
     * {@code comp = r2.compose(r1, Rotation_Convention.VECTOR_OPERATOR)}.
     * </p>
     * <p>
     * If the semantics of the rotations composition corresponds to a
     * {@link Rotation_Convention#FRAME_TRANSFORM frame transform} convention, * the application order will be reversed. So keeping the exact same
     * meaning of all {@code r1}, {@code r2}, {@code u}, {@code v}, {@code w}
     * and  {@code comp} as above, {@code comp} could also be computed as
     * {@code comp = r1.compose(r2, Rotation_Convention.FRAME_TRANSFORM)}.
     * </p>
     * @param r rotation to apply the rotation to
     * @param convention convention to use for the semantics of the angle
     * @return a rotation which is the composition of r by the instance
     */
    Field_Rotation<T> compose(const Field_Rotation<T> r, const Rotation_Convention convention)
    {
        return convention == Rotation_Convention.VECTOR_OPERATOR
            ? compose_internal(r)
            : r.compose_internal(this);
    }

private:

    /** Convert an orthogonal rotation matrix to a quaternion.
     * @param ort orthogonal rotation matrix
     * @return quaternion corresponding to the matrix
     */
    std::vector<T> mat2quat(const std::vector<std::vector<T>> ort)
    {

        const std::vector<T> quat = Math_Arrays::build_array(ort[0][0].get_field(), 4);

        // There are different ways to compute the quaternions elements
        // from the matrix. They all involve computing one element from
        // the diagonal of the matrix, and computing the three other ones
        // using a formula involving a division by the first element, // which unfortunately can be zero. sin_ce the norm of the
        // quaternion is 1, we know at least one element has an absolute
        // value greater or equal to 0.5, so it is always possible to
        // select the right formula and avoid division by zero and even
        // numerical inaccuracy. Checking the elements in turn and using
        // the first one greater than 0.45 is safe (this leads to a simple
        // test since qi = 0.45 implies 4 qi^2 - 1 = -0.19)
        T s = ort[0][0].add(ort[1][1]).add(ort[2][2]);
        if (s.get_real() > -0.19)
        {
            // compute my_q_0 and deduce my_q_1, my_q_2 and my_q_3
            quat[0] = s.add(1.0).sqrt().multiply(0.5);
            T inv = quat[0].reciprocal().multiply(0.25);
            quat[1] = inv.multiply(ort[1][2].subtract(ort[2][1]));
            quat[2] = inv.multiply(ort[2][0].subtract(ort[0][2]));
            quat[3] = inv.multiply(ort[0][1].subtract(ort[1][0]));
        }
        else
        {
            s = ort[0][0].subtract(ort[1][1]).subtract(ort[2][2]);
            if (s.get_real() > -0.19)
            {
                // compute my_q_1 and deduce my_q_0, my_q_2 and my_q_3
                quat[1] = s.add(1.0).sqrt().multiply(0.5);
                T inv = quat[1].reciprocal().multiply(0.25);
                quat[0] = inv.multiply(ort[1][2].subtract(ort[2][1]));
                quat[2] = inv.multiply(ort[0][1].add(ort[1][0]));
                quat[3] = inv.multiply(ort[0][2].add(ort[2][0]));
            }
            else
            {
                s = ort[1][1].subtract(ort[0][0]).subtract(ort[2][2]);
                if (s.get_real() > -0.19)
                {
                    // compute my_q_2 and deduce my_q_0, my_q_1 and my_q_3
                    quat[2] = s.add(1.0).sqrt().multiply(0.5);
                    T inv = quat[2].reciprocal().multiply(0.25);
                    quat[0] = inv.multiply(ort[2][0].subtract(ort[0][2]));
                    quat[1] = inv.multiply(ort[0][1].add(ort[1][0]));
                    quat[3] = inv.multiply(ort[2][1].add(ort[1][2]));
                }
                else
                {
                    // compute my_q_3 and deduce my_q_0, my_q_1 and my_q_2
                    s = ort[2][2].subtract(ort[0][0]).subtract(ort[1][1]);
                    quat[3] = s.add(1.0).sqrt().multiply(0.5);
                    T inv = quat[3].reciprocal().multiply(0.25);
                    quat[0] = inv.multiply(ort[0][1].subtract(ort[1][0]));
                    quat[1] = inv.multiply(ort[0][2].add(ort[2][0]));
                    quat[2] = inv.multiply(ort[2][1].add(ort[1][2]));
                }
            }
        }

        return quat;
    }

    /** Create a dimension 3 array.
     * @param a0 first array element
     * @param a1 second array element
     * @param a2 third array element
     * @return array
     */
    std::vector<T> build_array(const T a0, const T a1, const T a2) 
    {
        const std::vector<T> array = Math_Arrays::build_array(a0.get_field(), 3);
        array[0] = a0;
        array[1] = a1;
        array[2] = a2;
        return array;
    }

    /** Create a constant vector.
     * @param x abscissa
     * @param y ordinate
     * @param z height
     * @return a constant vector
     */
    Field_Vector_3D<T> vector(const double& x, const double& y, const double& z) 
    {
        const T zero = my_q_0.get_field().get_zero();
        return Field_Vector_3D<T>(zero.add(x), zero.add(y), zero.add(z));
    }

    /** Compose the instance with another rotation using vector operator convention.
     * @param r rotation to apply the rotation to
     * @return a rotation which is the composition of r by the instance
     * using vector operator convention
     */
    Field_Rotation<T> compose_internal(const Field_Rotation<T> r) 
    {
        return Field_Rotation<T>(r.my_q_0.multiply(my_q_0).subtract(r.my_q_1.multiply(my_q_1).add(r.my_q_2.multiply(my_q_2)).add(r.my_q_3.multiply(my_q_3))), r.my_q_1.multiply(my_q_0).add(r.my_q_0.multiply(my_q_1)).add(r.my_q_2.multiply(my_q_3).subtract(r.my_q_3.multiply(my_q_2))), r.my_q_2.multiply(my_q_0).add(r.my_q_0.multiply(my_q_2)).add(r.my_q_3.multiply(my_q_1).subtract(r.my_q_1.multiply(my_q_3))), r.my_q_3.multiply(my_q_0).add(r.my_q_0.multiply(my_q_3)).add(r.my_q_1.multiply(my_q_2).subtract(r.my_q_2.multiply(my_q_1))), false);
    }

    /** Compose the instance with another rotation using vector operator convention.
     * @param r rotation to apply the rotation to
     * @return a rotation which is the composition of r by the instance
     * using vector operator convention
     */
    Field_Rotation<T> compose_internal(const Rotation r) 
    {
        return Field_Rotation<T>(my_q_0.multiply(r.get_my_q_0()).subtract(my_q_1.multiply(r.get_my_q_1()).add(my_q_2.multiply(r.get_my_q_2())).add(my_q_3.multiply(r.get_my_q_3()))), my_q_0.multiply(r.get_my_q_1()).add(my_q_1.multiply(r.get_my_q_0())).add(my_q_3.multiply(r.get_my_q_2()).subtract(my_q_2.multiply(r.get_my_q_3()))), my_q_0.multiply(r.get_my_q_2()).add(my_q_2.multiply(r.get_my_q_0())).add(my_q_1.multiply(r.get_my_q_3()).subtract(my_q_3.multiply(r.get_my_q_1()))), my_q_0.multiply(r.get_my_q_3()).add(my_q_3.multiply(r.get_my_q_0())).add(my_q_2.multiply(r.get_my_q_1()).subtract(my_q_1.multiply(r.get_my_q_2()))), false);
    }

    /** Compose the inverse of the instance with another rotation
     * using vector operator convention.
     * @param r rotation to apply the rotation to
     * @return a rotation which is the composition of r by the inverse
     * of the instance using vector operator convention
     */
    Field_Rotation<T> compose_inverse_internal(Field_Rotation<T> r) 
    {
        return Field_Rotation<T>(r.my_q_0.multiply(my_q_0).add(r.my_q_1.multiply(my_q_1).add(r.my_q_2.multiply(my_q_2)).add(r.my_q_3.multiply(my_q_3))).negate(), r.my_q_0.multiply(my_q_1).add(r.my_q_2.multiply(my_q_3).subtract(r.my_q_3.multiply(my_q_2))).subtract(r.my_q_1.multiply(my_q_0)), r.my_q_0.multiply(my_q_2).add(r.my_q_3.multiply(my_q_1).subtract(r.my_q_1.multiply(my_q_3))).subtract(r.my_q_2.multiply(my_q_0)), r.my_q_0.multiply(my_q_3).add(r.my_q_1.multiply(my_q_2).subtract(r.my_q_2.multiply(my_q_1))).subtract(r.my_q_3.multiply(my_q_0)), false);
    }

    /** Compose the inverse of the instance with another rotation
     * using vector operator convention.
     * @param r rotation to apply the rotation to
     * @return a rotation which is the composition of r by the inverse
     * of the instance using vector operator convention
     */
    Field_Rotation<T> compose_inverse_internal(Rotation& r) 
    {
        return Field_Rotation<T>(my_q_0.multiply(r.get_my_q_0()).add(my_q_1.multiply(r.get_my_q_1()).add(my_q_2.multiply(r.get_my_q_2())).add(my_q_3.multiply(r.get_my_q_3()))).negate(), my_q_1.multiply(r.get_my_q_0()).add(my_q_3.multiply(r.get_my_q_2()).subtract(my_q_2.multiply(r.get_my_q_3()))).subtract(my_q_0.multiply(r.get_my_q_1())), my_q_2.multiply(r.get_my_q_0()).add(my_q_1.multiply(r.get_my_q_3()).subtract(my_q_3.multiply(r.get_my_q_1()))).subtract(my_q_0.multiply(r.get_my_q_2())), my_q_3.multiply(r.get_my_q_0()).add(my_q_2.multiply(r.get_my_q_1()).subtract(my_q_1.multiply(r.get_my_q_2()))).subtract(my_q_0.multiply(r.get_my_q_3())), false);
    }

    /** Perfect orthogonality on a 3X3 matrix.
     * @param m initial matrix (not exactly orthogonal)
     * @param threshold convergence threshold for the iterative
     * orthogonality correction (convergence is reached when the
     * difference between two steps of the Frobenius norm of the
     * correction is below this threshold)
     * @return an orthogonal matrix close to m
     * @exception  if the matrix cannot be
     * orthogonalized with the given threshold after 10 iterations
     */
    std::vector<std::vector<T>> orthogonalize_matrix(const std::vector<std::vector<T>>& m, const double& threshold)
    {
        T x00 = m[0][0];
        T x01 = m[0][1];
        T x02 = m[0][2];
        T x10 = m[1][0];
        T x11 = m[1][1];
        T x12 = m[1][2];
        T x20 = m[2][0];
        T x21 = m[2][1];
        T x22 = m[2][2];
        double fn{};
        double fn1;

        const std::vector<std::vector<T>> o = Math_Arrays::build_array(m[0][0].get_field(), 3, 3);

        // iterative correction: Xn+1 = Xn - 0.5 * (Xn.Mt.Xn - M)
        int i;
        for (i = 0; i < 11; ++i) 
        {

            // Mt.Xn
            const T mx00 = m[0][0].multiply(x00).add(m[1][0].multiply(x10)).add(m[2][0].multiply(x20));
            const T mx10 = m[0][1].multiply(x00).add(m[1][1].multiply(x10)).add(m[2][1].multiply(x20));
            const T mx20 = m[0][2].multiply(x00).add(m[1][2].multiply(x10)).add(m[2][2].multiply(x20));
            const T mx01 = m[0][0].multiply(x01).add(m[1][0].multiply(x11)).add(m[2][0].multiply(x21));
            const T mx11 = m[0][1].multiply(x01).add(m[1][1].multiply(x11)).add(m[2][1].multiply(x21));
            const T mx21 = m[0][2].multiply(x01).add(m[1][2].multiply(x11)).add(m[2][2].multiply(x21));
            const T mx02 = m[0][0].multiply(x02).add(m[1][0].multiply(x12)).add(m[2][0].multiply(x22));
            const T mx12 = m[0][1].multiply(x02).add(m[1][1].multiply(x12)).add(m[2][1].multiply(x22));
            const T mx22 = m[0][2].multiply(x02).add(m[1][2].multiply(x12)).add(m[2][2].multiply(x22));

            // Xn+1
            o[0][0] = x00.subtract(x00.multiply(mx00).add(x01.multiply(mx10)).add(x02.multiply(mx20)).subtract(m[0][0]).multiply(0.5));
            o[0][1] = x01.subtract(x00.multiply(mx01).add(x01.multiply(mx11)).add(x02.multiply(mx21)).subtract(m[0][1]).multiply(0.5));
            o[0][2] = x02.subtract(x00.multiply(mx02).add(x01.multiply(mx12)).add(x02.multiply(mx22)).subtract(m[0][2]).multiply(0.5));
            o[1][0] = x10.subtract(x10.multiply(mx00).add(x11.multiply(mx10)).add(x12.multiply(mx20)).subtract(m[1][0]).multiply(0.5));
            o[1][1] = x11.subtract(x10.multiply(mx01).add(x11.multiply(mx11)).add(x12.multiply(mx21)).subtract(m[1][1]).multiply(0.5));
            o[1][2] = x12.subtract(x10.multiply(mx02).add(x11.multiply(mx12)).add(x12.multiply(mx22)).subtract(m[1][2]).multiply(0.5));
            o[2][0] = x20.subtract(x20.multiply(mx00).add(x21.multiply(mx10)).add(x22.multiply(mx20)).subtract(m[2][0]).multiply(0.5));
            o[2][1] = x21.subtract(x20.multiply(mx01).add(x21.multiply(mx11)).add(x22.multiply(mx21)).subtract(m[2][1]).multiply(0.5));
            o[2][2] = x22.subtract(x20.multiply(mx02).add(x21.multiply(mx12)).add(x22.multiply(mx22)).subtract(m[2][2]).multiply(0.5));

            // correction on each elements
            const double corr00 = o[0][0].get_real() - m[0][0].get_real();
            const double corr01 = o[0][1].get_real() - m[0][1].get_real();
            const double corr02 = o[0][2].get_real() - m[0][2].get_real();
            const double corr10 = o[1][0].get_real() - m[1][0].get_real();
            const double corr11 = o[1][1].get_real() - m[1][1].get_real();
            const double corr12 = o[1][2].get_real() - m[1][2].get_real();
            const double corr20 = o[2][0].get_real() - m[2][0].get_real();
            const double corr21 = o[2][1].get_real() - m[2][1].get_real();
            const double corr22 = o[2][2].get_real() - m[2][2].get_real();

            // Frobenius norm of the correction
            fn1 = corr00 * corr00 + corr01 * corr01 + corr02 * corr02 +
                  corr10 * corr10 + corr11 * corr11 + corr12 * corr12 +
                  corr20 * corr20 + corr21 * corr21 + corr22 * corr22;

            // convergence test
            if (std::abs(fn1 - fn) <= threshold) 
            {
                return o;
            }

            // prepare next iteration
            x00 = o[0][0];
            x01 = o[0][1];
            x02 = o[0][2];
            x10 = o[1][0];
            x11 = o[1][1];
            x12 = o[1][2];
            x20 = o[2][0];
            x21 = o[2][1];
            x22 = o[2][2];
            fn  = fn1;
        }
        // the algorithm did not converge after 10 iterations
        throw (Localized_Geometry_Formats.UNABLE_TO_ORTHOGONOLIZE_MATRIX, i - 1);
    }
};