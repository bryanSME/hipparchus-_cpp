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
  //package org.hipparchus.geometry.spherical.oned;

  //import org.hipparchus.geometry.partitioning.Abstract_Sub_Hyperplane;
  //import org.hipparchus.geometry.partitioning.Hyperplane;
  //import org.hipparchus.geometry.partitioning.Region;

  /** This class represents sub-hyperplane for {@link Limit_Angle}.
   * <p>Instances of this class are guaranteed to be immutable.</p>
   */
class SubLimit_Angle extends Abstract_Sub_Hyperplane<Sphere_1D, Sphere_1D>
{
	/** Simple constructor.
	 * @param hyperplane underlying hyperplane
	 * @param remaining_region remaining region of the hyperplane
	 */
	public SubLimit_Angle(const Hyperplane<Sphere_1D> hyperplane, const Region<Sphere_1D>& remaining_region)
	{
		super(hyperplane, remaining_region);
	}

	/** {@inherit_doc} */
	//override
	public double get_size()
	{
		return 0;
	}

	/** {@inherit_doc} */
	//override
	public bool is_empty()
	{
		return false;
	}

	/** {@inherit_doc} */
	//override
	protected Abstract_Sub_Hyperplane<Sphere_1D, Sphere_1D> build_new(const Hyperplane<Sphere_1D> hyperplane, const Region<Sphere_1D>& remaining_region)
	{
		return SubLimit_Angle(hyperplane, remaining_region);
	}

	/** {@inherit_doc} */
	//override
	public Split_Sub_Hyperplane<Sphere_1D> split(const Hyperplane<Sphere_1D> hyperplane)
	{
		const double global = hyperplane.get_offset(((Limit_Angle)get_hyperplane()).get_location());
		if (global < -hyperplane.get_tolerance())
		{
			return Split_Sub_Hyperplane<Sphere_1D>(null, this);
		}
		else if (global > hyperplane.get_tolerance())
		{
			return Split_Sub_Hyperplane<Sphere_1D>(this, NULL);
		}
		else
		{
			return Split_Sub_Hyperplane<Sphere_1D>(null, NULL);
		}
	}
}
