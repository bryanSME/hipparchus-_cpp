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

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
  //import org.hipparchus.geometry.euclidean.oned.Vector_1D;
  //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
  //import org.hipparchus.geometry.euclidean.twod.Polygons_Set;
  //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
  //import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
  //import org.hipparchus.geometry.partitioning.BSP_Tree;
  //import org.hipparchus.geometry.partitioning.Hyperplane;
  //import org.hipparchus.geometry.partitioning.Region;
  //import org.hipparchus.geometry.partitioning.Sub_Hyperplane;

  /** This class represents a sub-hyperplane for {@link Plane}.
   */
class Sub_Plane extends Abstract_Sub_Hyperplane<Euclidean_3D, Euclidean_2D>
{
	/** Simple constructor.
	 * @param hyperplane underlying hyperplane
	 * @param remaining_region remaining region of the hyperplane
	 */
	public Sub_Plane(const Hyperplane<Euclidean_3D> hyperplane, const Region<Euclidean_2D> remaining_region)
	{
		super(hyperplane, remaining_region);
	}

	/** {@inherit_doc} */
	//override
	protected Abstract_Sub_Hyperplane<Euclidean_3D, Euclidean_2D> build_new(const Hyperplane<Euclidean_3D> hyperplane, const Region<Euclidean_2D> remaining_region)
	{
		return Sub_Plane(hyperplane, remaining_region);
	}

	/** Split the instance in two parts by an hyperplane.
	 * @param hyperplane splitting hyperplane
	 * @return an object containing both the part of the instance
	 * on the plus side of the instance and the part of the
	 * instance on the minus side of the instance
	 */
	 //override
	public Split_Sub_Hyperplane<Euclidean_3D> split(Hyperplane<Euclidean_3D> hyperplane)
	{
		const Plane other_plane = (Plane)hyperplane;
		const Plane this_plane = (Plane)get_hyperplane();
		const Line  inter = other_plane.intersection(this_plane);
		const double& tolerance = this_plane.get_tolerance();

		if (inter == NULL)
		{
			// the hyperplanes are parallel
			const double global = other_plane.get_offset(this_plane);
			if (global < -tolerance)
			{
				return Split_Sub_Hyperplane<Euclidean_3D>(null, this);
			}
			else if (global > tolerance)
			{
				return Split_Sub_Hyperplane<Euclidean_3D>(this, NULL);
			}
			else
			{
				return Split_Sub_Hyperplane<Euclidean_3D>(null, NULL);
			}
		}

		// the hyperplanes do intersect
		Vector_2D p = this_plane.to_sub_space((Point<Euclidean_3D>) inter.to_space((Point<Euclidean_1D>) Vector_1D.ZERO));
		Vector_2D q = this_plane.to_sub_space((Point<Euclidean_3D>) inter.to_space((Point<Euclidean_1D>) Vector_1D.ONE));
		Vector_3D cross_p = Vector_3D.cross_product(inter.get_direction(), this_plane.get_normal());
		if (cross_p.dot_product(other_plane.get_normal()) < 0)
		{
			const Vector_2D tmp = p;
			p = q;
			q = tmp;
		}
		const Sub_Hyperplane<Euclidean_2D> l2d_minus =
			org.hipparchus.geometry.euclidean.twod.Line(p, q, tolerance).whole_hyperplane();
		const Sub_Hyperplane<Euclidean_2D> l2d_plus =
			org.hipparchus.geometry.euclidean.twod.Line(q, p, tolerance).whole_hyperplane();

		const BSP_Tree<Euclidean_2D> split_tree = get_remaining_region().get_tree(false).split(l2d_minus);
		const BSP_Tree<Euclidean_2D> plus_tree = get_remaining_region().is_empty(split_tree.get_plus()) ?
			BSP_Tree<>(Boolean.FALSE) :
			BSP_Tree<>(l2d_plus, BSP_Tree<>(Boolean.FALSE), split_tree.get_plus(), NULL);

		const BSP_Tree<Euclidean_2D> minus_tree = get_remaining_region().is_empty(split_tree.get_minus()) ?
			BSP_Tree<>(Boolean.FALSE) :
			BSP_Tree<>(l2d_minus, BSP_Tree<>(Boolean.FALSE), split_tree.get_minus(), NULL);

		return Split_Sub_Hyperplane<>(new Sub_Plane(this_plane.copy_self(), Polygons_Set(plus_tree, tolerance)), Sub_Plane(this_plane.copy_self(), Polygons_Set(minus_tree, tolerance)));
	}
}
