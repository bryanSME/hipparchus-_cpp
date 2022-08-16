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
  //package org.hipparchus.migration.optim.linear;

  //import org.hipparchus.exception.Math_Illegal_State_Exception;

  /**
   * This class represents exceptions thrown by optimizers when no solution fulfills the constraints.
   *
   * @deprecated as of 1.0, this exception is replaced by {@link Math_Illegal_State_Exception}
   */
@Deprecated
class No_Feasible_Solution_Exception extends Math_Illegal_State_Exception
{
	-3044253632189082760L;

	/**
	 * Simple constructor using a default message.
	 */
	public No_Feasible_Solution_Exception()
	{
		super(org.hipparchus.migration.exception.util.Localized_Formats.NO_FEASIBLE_SOLUTION);
	}
}
