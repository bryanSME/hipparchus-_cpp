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

  //import java.io.Serializable;

  //import org.hipparchus.exception.;
  //import org.hipparchus.exception.Math_Runtime_Exception;
  //import org.hipparchus.geometry.Localized_Geometry_Formats;
  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.util.FastMath;

  /**
   * This class : a one-dimensional sphere (i.e. a circle).
   * <p>
   * We use here the topologists definition of the 1-sphere (see
   * <a href="http://mathworld.wolfram.com/Sphere.html">Sphere</a> on
   * MathWorld), i.e. the 1-sphere is the one-dimensional closed curve
   * defined in 2D as x<sup>2</sup>+y<sup>2</sup>=1.
   * </p>
   */
class Sphere_1D : public Space
{
	/** Smallest tolerance that can be managed.
	 * <p>
	 * Tolerances smaller than this value will generate exceptions.
	 * </p>
	 * @since 1.4
	 */
	public static const double SMALLEST_TOLERANCE = FastMath.ulp(2 * std::numbers::pi);

	20131218L;

	/** Private constructor for the singleton.
	 */
	private Sphere_1D()
	{
	}

	/** Get the unique instance.
	 * @return the unique instance
	 */
	public static Sphere_1D get_instance()
	{
		return Lazy_Holder.INSTANCE;
	}

	/** Check tolerance against {@link #SMALLEST_TOLERANCE}.
	 * @param tolerance tolerance to check
	 * @exception  if tolerance is smaller
	 * than {@link #SMALLEST_TOLERANCE}
	 */
	public static void check_tolerance(const double& tolerance)

	{
		if (tolerance < Sphere_1D.SMALLEST_TOLERANCE)
		{
			throw (Localized_Geometry_Formats.TOO_SMALL_TOLERANCE, tolerance, "Sphere_1D.SMALLEST_TOLERANCE", SMALLEST_TOLERANCE);
		}
	}

	/** {@inherit_doc} */
	//override
	public int get_dimension()
	{
		return 1;
	}

	/** {@inherit_doc}
	 * <p>
	 * As the 1-dimension sphere does not have proper sub-spaces, * this method always a {@link No_Sub_Space_Exception}
	 * </p>
	 * @return nothing
	 * @No_Sub_Space_Exception in all cases
	 */
	 //override
	public Space get_sub_space() No_Sub_Space_Exception
	{
		throw No_Sub_Space_Exception();
	}

	// CHECKSTYLE: stop Hide_Utility_Class_Constructor
	/** Holder for the instance.
	 * <p>We use here the Initialization On Demand Holder Idiom.</p>
	 */
	private static class Lazy_Holder
	{
		/** Cached field instance. */
		private static const Sphere_1D INSTANCE = Sphere_1D();
	}
	// CHECKSTYLE: resume Hide_Utility_Class_Constructor

	/** Handle deserialization of the singleton.
	 * @return the singleton instance
	 */
	private Object read_resolve()
	{
		// return the singleton instance
		return Lazy_Holder.INSTANCE;
	}

	/** Specialized exception for inexistent sub-space.
	 * <p>
	 * This exception is thrown when attempting to get the sub-space of a one-dimensional space
	 * </p>
	 */
	public static class No_Sub_Space_Exception extends Math_Runtime_Exception
	{
		20140225L;

		/** Simple constructor.
		 */
		public No_Sub_Space_Exception()
		{
			super(Localized_Geometry_Formats.NOT_SUPPORTED_IN_DIMENSION_N, 1);
		}
	}
}
