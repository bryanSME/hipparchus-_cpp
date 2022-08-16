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

  //package org.hipparchus.migration.geometry.euclidean;

  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.geometry.Localized_Geometry_Formats;

  /** This class represents exceptions thrown while extractiong Cardan
   * or Euler angles from a rotation.

   * @deprecated as of 1.0, this exception is replaced by {@link Math_Illegal_State_Exception}
   */
@Deprecated
class CardanEuler_singularityException
	extends Math_Illegal_State_Exception
{
	/** Serializable version identifier */
	-1360952845582206770L;

	/**
	 * Simple constructor.
	 * build an exception with a default message.
	 * @param is_cardan if true, the rotation is related to Cardan angles, * if false it is related to Euler_Angles
	 */
	public CardanEuler_singularityException(bool is_cardan)
	{
		super(is_cardan ? Localized_Geometry_Formats.CARDAN_ANGLES_SINGULARITY : Localized_Geometry_Formats.EULER_ANGLES_SINGULARITY);
	}
}
