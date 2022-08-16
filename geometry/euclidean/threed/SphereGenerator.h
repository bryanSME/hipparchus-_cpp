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

  //import java.util.Arrays;
  //import java.util.List;

  //import org.hipparchus.fraction.Big_Fraction;
  //import org.hipparchus.geometry.enclosing.Enclosing_Ball;
  //import org.hipparchus.geometry.enclosing.Support_Ball_Generator;
  //import org.hipparchus.geometry.euclidean.twod.Disk_Generator;
  //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
  //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
  //import org.hipparchus.util.FastMath;

  /** Class generating an enclosing ball from its support points.
   */
class Sphere_Generator : Support_Ball_Generator<Euclidean_3D, Vector_3D>
{
	/** {@inherit_doc} */
	//override
	public Enclosing_Ball<Euclidean_3D, Vector_3D> ball_on_support(const List<Vector_3D> support)
	{
		if (support.is_empty())
		{
			return Enclosing_Ball<>(Vector_3D.ZERO, -INFINITY);
		}
		else
		{
			const Vector_3D v_a = support.get(0);
			if (support.size() < 2)
			{
				return Enclosing_Ball<>(v_a, 0, v_a);
			}
			else
			{
				const Vector_3D v_b = support.get(1);
				if (support.size() < 3)
				{
					const Vector_3D center = Vector_3D(0.5, v_a, 0.5, v_b);

					// we could have computed r directly from the v_a and v_b
					// (it was done this way up to Hipparchus 1.0), but as center
					// is approximated in the computation above, it is better to
					// take the const value of center and compute r from the distances
					// to center of all support points, using a max to ensure all support
					// points belong to the ball
					// see <https://github.com/Hipparchus-Math/hipparchus/issues/20>
					const double r = std::max(Vector_3D.distance(v_a, center), Vector_3D.distance(v_b, center));
					return Enclosing_Ball<>(center, r, v_a, v_b);
				}
				else
				{
					const Vector_3D v_c = support.get(2);
					if (support.size() < 4)
					{
						// delegate to 2D disk generator
						const Plane p = Plane(v_a, v_b, v_c, 1.0e-10 * (v_a.get_norm1() + v_b.get_norm1() + v_c.get_norm1()));
						const Enclosing_Ball<Euclidean_2D, Vector_2D> disk =
							Disk_Generator().ball_on_support(Arrays.as_list(p.to_sub_space(v_a), p.to_sub_space(v_b), p.to_sub_space(v_c)));

						// convert back to 3D
						const Vector_3D center = p.to_space(disk.get_center());

						// we could have computed r directly from the v_a and v_b
						// (it was done this way up to Hipparchus 1.0), but as center
						// is approximated in the computation above, it is better to
						// take the const value of center and compute r from the distances
						// to center of all support points, using a max to ensure all support
						// points belong to the ball
						// see <https://github.com/Hipparchus-Math/hipparchus/issues/20>
						const double r = std::max(Vector_3D.distance(v_a, center), std::max(Vector_3D.distance(v_b, center), Vector_3D.distance(v_c, center)));
						return Enclosing_Ball<>(center, r, v_a, v_b, v_c);
					}
					else
					{
						const Vector_3D vD = support.get(3);
						// a sphere is 3D can be defined as:
						// (1)   (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 = r^2
						// which can be written:
						// (2)   (x^2 + y^2 + z^2) - 2 x_0 x - 2 y_0 y - 2 z_0 z + (x_0^2 + y_0^2 + z_0^2 - r^2) = 0
						// or simply:
						// (3)   (x^2 + y^2 + z^2) + a x + b y + c z + d = 0
						// with sphere center coordinates -a/2, -b/2, -c/2
						// If the sphere exists, a b, c and d are a non zero solution to
						// [ (x^2  + y^2  + z^2)    x    y   z    1 ]   [ 1 ]   [ 0 ]
						// [ (x_a^2 + y_a^2 + zA^2)   x_a   y_a  zA   1 ]   [ a ]   [ 0 ]
						// [ (x_b^2 + yB^2 + zB^2)   x_b   yB  zB   1 ] * [ b ] = [ 0 ]
						// [ (x_c^2 + yC^2 + zC^2)   x_c   yC  zC   1 ]   [ c ]   [ 0 ]
						// [ (xD^2 + yD^2 + zD^2)   xD   yD  zD   1 ]   [ d ]   [ 0 ]
						// So the determinant of the matrix is zero. Computing this determinant
						// by expanding it using the minors m_ij of first row leads to
						// (4)   m_11 (x^2 + y^2 + z^2) - m_12 x + m_13 y - m_14 z + m_15 = 0
						// So by identifying equations (2) and (4) we get the coordinates
						// of center as:
						//      x_0 = +m_12 / (2 m_11)
						//      y_0 = -m_13 / (2 m_11)
						//      z_0 = +m_14 / (2 m_11)
						// Note that the minors m_11, m_12, m_13 and m_14 all have the last column
						// filled with 1.0, hence simplifying the computation
						const std::vector<Big_Fraction>c2 =
						{
							Big_Fraction(v_a.get_x()), Big_Fraction(v_b.get_x()), Big_Fraction(v_c.get_x()), Big_Fraction(vD.get_x())
						};
						const std::vector<Big_Fraction>c3 =
						{
							Big_Fraction(v_a.get_y()), Big_Fraction(v_b.get_y()), Big_Fraction(v_c.get_y()), Big_Fraction(vD.get_y())
						};
						const std::vector<Big_Fraction>c4 =
						{
							Big_Fraction(v_a.get_z()), Big_Fraction(v_b.get_z()), Big_Fraction(v_c.get_z()), Big_Fraction(vD.get_z())
						};
						const std::vector<Big_Fraction>c1 =
						{
							c2[0].multiply(c2[0]).add(c3[0].multiply(c3[0])).add(c4[0].multiply(c4[0])), c2[1].multiply(c2[1]).add(c3[1].multiply(c3[1])).add(c4[1].multiply(c4[1])), c2[2].multiply(c2[2]).add(c3[2].multiply(c3[2])).add(c4[2].multiply(c4[2])), c2[3].multiply(c2[3]).add(c3[3].multiply(c3[3])).add(c4[3].multiply(c4[3]))
						};
						const Big_Fraction two_m11 = minor(c2, c3, c4).multiply(2);
						const Big_Fraction m12 = minor(c1, c3, c4);
						const Big_Fraction m13 = minor(c1, c2, c4);
						const Big_Fraction m14 = minor(c1, c2, c3);
						const Vector_3D center = Vector_3D(m12.divide(two_m11).double_value(), -m13.divide(two_m11).double_value(), m14.divide(two_m11).double_value());

						// we could have computed r directly from the minors above
						// (it was done this way up to Hipparchus 1.0), but as center
						// is approximated in the computation above, it is better to
						// take the const value of center and compute r from the distances
						// to center of all support points, using a max to ensure all support
						// points belong to the ball
						// see <https://github.com/Hipparchus-Math/hipparchus/issues/20>
						const double r = std::max(Vector_3D.distance(v_a, center), std::max(Vector_3D.distance(v_b, center), std::max(Vector_3D.distance(v_c, center), Vector_3D.distance(vD, center))));
						return Enclosing_Ball<>(center, r, v_a, v_b, v_c, vD);
					}
				}
			}
		}
	}

	/** Compute a dimension 4 minor, when 4<sup>th</sup> column is known to be filled with 1.0.
	 * @param c1 first column
	 * @param c2 second column
	 * @param c3 third column
	 * @return value of the minor computed has an exact fraction
	 */
	private Big_Fraction minor(const std::vector<Big_Fraction>c1, const std::vector<Big_Fraction>c2, const std::vector<Big_Fraction>c3)
	{
		return      c2[0].multiply(c3[1]).multiply(c1[2].subtract(c1[3])).
			add(c2[0].multiply(c3[2]).multiply(c1[3].subtract(c1[1]))).
			add(c2[0].multiply(c3[3]).multiply(c1[1].subtract(c1[2]))).
			add(c2[1].multiply(c3[0]).multiply(c1[3].subtract(c1[2]))).
			add(c2[1].multiply(c3[2]).multiply(c1[0].subtract(c1[3]))).
			add(c2[1].multiply(c3[3]).multiply(c1[2].subtract(c1[0]))).
			add(c2[2].multiply(c3[0]).multiply(c1[1].subtract(c1[3]))).
			add(c2[2].multiply(c3[1]).multiply(c1[3].subtract(c1[0]))).
			add(c2[2].multiply(c3[3]).multiply(c1[0].subtract(c1[1]))).
			add(c2[3].multiply(c3[0]).multiply(c1[2].subtract(c1[1]))).
			add(c2[3].multiply(c3[1]).multiply(c1[0].subtract(c1[2]))).
			add(c2[3].multiply(c3[2]).multiply(c1[1].subtract(c1[0])));
	}
}
