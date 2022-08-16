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

  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;

  /**
   * This class : a three-dimensional space.
   */
class Euclidean_3D, Space
{
	6249091865814886817L;

	/** Private constructor for the singleton.
	 */
	private Euclidean_3D()
	{
	}

	/** Get the unique instance.
	 * @return the unique instance
	 */
	public static Euclidean_3D get_instance()
	{
		return Lazy_Holder.INSTANCE;
	}

	/** {@inherit_doc} */
	//override
	public int get_dimension()
	{
		return 3;
	}

	/** {@inherit_doc} */
	//override
	public Euclidean_2D get_sub_space()
	{
		return Euclidean_2D.get_instance();
	}

	// CHECKSTYLE: stop Hide_Utility_Class_Constructor
	/** Holder for the instance.
	 * <p>We use here the Initialization On Demand Holder Idiom.</p>
	 */
	private static class Lazy_Holder
	{
		/** Cached field instance. */
		private static const Euclidean_3D INSTANCE = Euclidean_3D();
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
}
