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

  //package org.hipparchus.migration.geometry.euclidean.threed;

  //import org.hipparchus.exception.Localizable;
  //import org.hipparchus.exception.;

  /**
   * This class represents exceptions thrown while building rotations
   * from matrices.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.Math_Illegal_State_Exception}
   */
@Deprecated
class Not_A_Rotation_Matrix_Exception
	extends
{
	/** Serializable version identifier */
	5647178478658937642L;

	/**
	 * Simple constructor.
	 * Build an exception by translating and formating a message
	 * @param specifier format specifier (to be translated)
	 * @param parts to insert in the format (no translation)
	 * @since 2.2
	 */
	public Not_A_Rotation_Matrix_Exception(Localizable specifier, Object ... parts)
	{
		super(specifier, parts);
	}
}
