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
  //package org.hipparchus.analysis.interpolation;

  //import org.hipparchus.analysis.Trivariate_Function;
  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.util.Math_Arrays;
  //import org.hipparchus.util.Math_Utils;
#include <vector>
#include "Tricubic_Function.h"

  /**
   * Function that : the
   * <a href="http://en.wikipedia.org/wiki/Tricubic_interpolation">
   * tricubic spline interpolation</a>, as proposed in
   * <blockquote>
   *  Tricubic interpolation in three dimensions<br>
   *  F. Lekien and J. Marsden<br>
   *  <em>Int. J. Numer. Meth. Eng</em> 2005; <b>63</b>:455-471<br>
   * </blockquote>
   *
   */
class Tricubic_Interpolating_Function : Trivariate_Function
{
private:
	/**
	 * Matrix to compute the spline coefficients from the function values
	 * and function derivatives values
	 */
	static inline const std::vector<std::vector<double>> AINV
	{
		{ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ -3,3,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 2,-2,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ -3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,-3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 9,-9,-9,9,0,0,0,0,6,3,-6,-3,0,0,0,0,6,-6,3,-3,0,0,0,0,0,0,0,0,0,0,0,0,4,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ -6,6,6,-6,0,0,0,0,-3,-3,3,3,0,0,0,0,-4,4,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,-2,-2,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 2,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ -6,6,6,-6,0,0,0,0,-4,-2,4,2,0,0,0,0,-3,3,-3,3,0,0,0,0,0,0,0,0,0,0,0,0,-2,-1,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 4,-4,-4,4,0,0,0,0,2,2,-2,-2,0,0,0,0,2,-2,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,0,0,0,0,-2,-1,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,0,0,0,0,1,1,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,-1,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-9,-9,9,0,0,0,0,0,0,0,0,0,0,0,0,6,3,-6,-3,0,0,0,0,6,-6,3,-3,0,0,0,0,4,2,2,1,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6,6,6,-6,0,0,0,0,0,0,0,0,0,0,0,0,-3,-3,3,3,0,0,0,0,-4,4,-2,2,0,0,0,0,-2,-2,-1,-1,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6,6,6,-6,0,0,0,0,0,0,0,0,0,0,0,0,-4,-2,4,2,0,0,0,0,-3,3,-3,3,0,0,0,0,-2,-1,-2,-1,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,-4,-4,4,0,0,0,0,0,0,0,0,0,0,0,0,2,2,-2,-2,0,0,0,0,2,-2,2,-2,0,0,0,0,1,1,1,1,0,0,0,0 },
		{-3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,-3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 9,-9,0,0,-9,9,0,0,6,3,0,0,-6,-3,0,0,0,0,0,0,0,0,0,0,6,-6,0,0,3,-3,0,0,0,0,0,0,0,0,0,0,4,2,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ -6,6,0,0,6,-6,0,0,-3,-3,0,0,3,3,0,0,0,0,0,0,0,0,0,0,-4,4,0,0,-2,2,0,0,0,0,0,0,0,0,0,0,-2,-2,0,0,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0,0,-1,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-9,0,0,-9,9,0,0,0,0,0,0,0,0,0,0,6,3,0,0,-6,-3,0,0,0,0,0,0,0,0,0,0,6,-6,0,0,3,-3,0,0,4,2,0,0,2,1,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6,6,0,0,6,-6,0,0,0,0,0,0,0,0,0,0,-3,-3,0,0,3,3,0,0,0,0,0,0,0,0,0,0,-4,4,0,0,-2,2,0,0,-2,-2,0,0,-1,-1,0,0 },
		{ 9,0,-9,0,-9,0,9,0,0,0,0,0,0,0,0,0,6,0,3,0,-6,0,-3,0,6,0,-6,0,3,0,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,2,0,2,0,1,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,9,0,-9,0,-9,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,3,0,-6,0,-3,0,6,0,-6,0,3,0,-3,0,0,0,0,0,0,0,0,0,4,0,2,0,2,0,1,0 },
		{ -27,27,27,-27,27,-27,-27,27,-18,-9,18,9,18,9,-18,-9,-18,18,-9,9,18,-18,9,-9,-18,18,18,-18,-9,9,9,-9,-12,-6,-6,-3,12,6,6,3,-12,-6,12,6,-6,-3,6,3,-12,12,-6,6,-6,6,-3,3,-8,-4,-4,-2,-4,-2,-2,-1 },
		{ 18,-18,-18,18,-18,18,18,-18,9,9,-9,-9,-9,-9,9,9,12,-12,6,-6,-12,12,-6,6,12,-12,-12,12,6,-6,-6,6,6,6,3,3,-6,-6,-3,-3,6,6,-6,-6,3,3,-3,-3,8,-8,4,-4,4,-4,2,-2,4,4,2,2,2,2,1,1 },
		{ -6,0,6,0,6,0,-6,0,0,0,0,0,0,0,0,0,-3,0,-3,0,3,0,3,0,-4,0,4,0,-2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,-2,0,-1,0,-1,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,-6,0,6,0,6,0,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,-3,0,3,0,3,0,-4,0,4,0,-2,0,2,0,0,0,0,0,0,0,0,0,-2,0,-2,0,-1,0,-1,0 },
		{ 18,-18,-18,18,-18,18,18,-18,12,6,-12,-6,-12,-6,12,6,9,-9,9,-9,-9,9,-9,9,12,-12,-12,12,6,-6,-6,6,6,3,6,3,-6,-3,-6,-3,8,4,-8,-4,4,2,-4,-2,6,-6,6,-6,3,-3,3,-3,4,2,4,2,2,1,2,1 },
		{ -12,12,12,-12,12,-12,-12,12,-6,-6,6,6,6,6,-6,-6,-6,6,-6,6,6,-6,6,-6,-8,8,8,-8,-4,4,4,-4,-3,-3,-3,-3,3,3,3,3,-4,-4,4,4,-2,-2,2,2,-4,4,-4,4,-2,2,-2,2,-2,-2,-2,-2,-1,-1,-1,-1 }, 
		{ 2,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,2,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ -6,6,0,0,6,-6,0,0,-4,-2,0,0,4,2,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,-3,3,0,0,0,0,0,0,0,0,0,0,-2,-1,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 4,-4,0,0,-4,4,0,0,2,2,0,0,-2,-2,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,2,-2,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0 }, 
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6,6,0,0,6,-6,0,0,0,0,0,0,0,0,0,0,-4,-2,0,0,4,2,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,-3,3,0,0,-2,-1,0,0,-2,-1,0,0 },
		{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,-4,0,0,-4,4,0,0,0,0,0,0,0,0,0,0,2,2,0,0,-2,-2,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,2,-2,0,0,1,1,0,0,1,1,0,0 }, 
		{ -6,0,6,0,6,0,-6,0,0,0,0,0,0,0,0,0,-4,0,-2,0,4,0,2,0,-3,0,3,0,-3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,-1,0,-2,0,-1,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,-6,0,6,0,6,0,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4,0,-2,0,4,0,2,0,-3,0,3,0,-3,0,3,0,0,0,0,0,0,0,0,0,-2,0,-1,0,-2,0,-1,0 },
		{ 18,-18,-18,18,-18,18,18,-18,12,6,-12,-6,-12,-6,12,6,12,-12,6,-6,-12,12,-6,6,9,-9,-9,9,9,-9,-9,9,8,4,4,2,-8,-4,-4,-2,6,3,-6,-3,6,3,-6,-3,6,-6,3,-3,6,-6,3,-3,4,2,2,1,4,2,2,1 },
		{ -12,12,12,-12,12,-12,-12,12,-6,-6,6,6,6,6,-6,-6,-8,8,-4,4,8,-8,4,-4,-6,6,6,-6,-6,6,6,-6,-4,-4,-2,-2,4,4,2,2,-3,-3,3,3,-3,-3,3,3,-4,4,-2,2,-4,4,-2,2,-2,-2,-1,-1,-2,-2,-1,-1 }, 
		{ 4,0,-4,0,-4,0,4,0,0,0,0,0,0,0,0,0,2,0,2,0,-2,0,-2,0,2,0,-2,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0 },
		{ 0,0,0,0,0,0,0,0,4,0,-4,0,-4,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,0,-2,0,-2,0,2,0,-2,0,2,0,-2,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0 },
		{ -12,12,12,-12,12,-12,-12,12,-8,-4,8,4,8,4,-8,-4,-6,6,-6,6,6,-6,6,-6,-6,6,6,-6,-6,6,6,-6,-4,-2,-4,-2,4,2,4,2,-4,-2,4,2,-4,-2,4,2,-3,3,-3,3,-3,3,-3,3,-2,-1,-2,-1,-2,-1,-2,-1 },
		{ 8,-8,-8,8,-8,8,8,-8,4,4,-4,-4,-4,-4,4,4,4,-4,4,-4,-4,4,-4,4,4,-4,-4,4,4,-4,-4,4,2,2,2,2,-2,-2,-2,-2,2,2,-2,-2,2,2,-2,-2,2,-2,2,-2,2,-2,2,-2,1,1,1,1,1,1,1,1 }
	};

	/** Samples x-coordinates */
	const std::vector<double> my_xval;
	/** Samples y-coordinates */
	const std::vector<double> my_yval;
	/** Samples z-coordinates */
	const std::vector<double> my_zval;
	/** Set of cubic splines pacthing the whole data grid */
	std::vector<std::vector<std::vector<Tricubic_Function>>> my_splines;

public:
	/**
	 * @param x Sample values of the x-coordinate, in increasing order.
	 * @param y Sample values of the y-coordinate, in increasing order.
	 * @param z Sample values of the y-coordinate, in increasing order.
	 * @param f Values of the function on every grid point.
	 * @param dFdX Values of the partial derivative of function with respect to x on every grid point.
	 * @param d_fd_y Values of the partial derivative of function with respect to y on every grid point.
	 * @param d_fd_z Values of the partial derivative of function with respect to z on every grid point.
	 * @param d2FdXdY Values of the cross partial derivative of function on every grid point.
	 * @param d2_fd_xd_z Values of the cross partial derivative of function on every grid point.
	 * @param d2_fd_yd_z Values of the cross partial derivative of function on every grid point.
	 * @param d3_fd_xd_yd_z Values of the cross partial derivative of function on every grid point.
	 * @ if any of the arrays has zero length.
	 * @ if the various arrays do not contain the expected number of elements.
	 * @ if {@code x}, {@code y} or {@code z} are not strictly increasing.
	 */
	Tricubic_Interpolating_Function(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<std::vector<std::vector<double>>>& f, std::vector<std::vector<std::vector<double>>>& dFdX, std::vector < std::vector<std::vector<double>>>& d_fd_y, std::vector < std::vector<std::vector<double>>>& d_fd_z, std::vector<std::vector<std::vector<double>>>& d2FdXdY, std::vector<std::vector<std::vector<double>>>& d2_fd_xd_z, std::vector<std::vector<std::vector<double>>>& d2_fd_yd_z, std::vector<std::vector<std::vector<double>>>& d3_fd_xd_yd_z)
	:
		my_xval{ x },
		my_yval{ y },
		my_zval{ z }
	
	{
		const int x_len = x.size();
		const int y_len = y.size();
		const int z_len = z.size();

		if (x_len == 0 || y_len == 0 || z.size() == 0 || f.size() == 0 || f[0].size() == 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::NO_DATA);
		}
		if (x_len != f.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, f.size());
		}
		if (x_len != dFdX.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, dFdX.size());
		}
		if (x_len != d_fd_y.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, d_fd_y.size());
		}
		if (x_len != d_fd_z.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, d_fd_z.size());
		}
		if (x_len != d2FdXdY.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, d2FdXdY.size());
		}
		if (x_len != d2_fd_xd_z.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, d2_fd_xd_z.size());
		}
		if (x_len != d2_fd_yd_z.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, d2_fd_yd_z.size());
		}
		if (x_len != d3_fd_xd_yd_z.size())
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, x_len, d3_fd_xd_yd_z.size());
		}

		//Math_Arrays::check_order(x);
		//Math_Arrays::check_order(y);
		//Math_Arrays::check_order(z);

		const int last_i = x_len - 1;
		const int last_j = y_len - 1;
		const int last_k = z_len - 1;
		my_splines = std::vector<std::vector<std::vector<Tricubic_Function>>>(last_i, std::vector<std::vector<Tricubic_Function>>(last_j, std::vector<Tricubic_Function>(last_k)));

		for (int i{}; i < last_i; i++)
		{
			if (f[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, f[i].size(), y_len);
			}
			if (dFdX[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, dFdX[i].size(), y_len);
			}
			if (d_fd_y[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d_fd_y[i].size(), y_len);
			}
			if (d_fd_z[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d_fd_z[i].size(), y_len);
			}
			if (d2FdXdY[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d2FdXdY[i].size(), y_len);
			}
			if (d2_fd_xd_z[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d2_fd_xd_z[i].size(), y_len);
			}
			if (d2_fd_yd_z[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d2_fd_yd_z[i].size(), y_len);
			}
			if (d3_fd_xd_yd_z[i].size() != y_len)
			{
				throw std::exception("not implemented");
				//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d3_fd_xd_yd_z[i].size(), y_len);
			}

			const int ip1 = i + 1;
			const double x_r = my_xval[ip1] - my_xval[i];
			for (int j{}; j < last_j; j++)
			{
				if (f[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, f[i][j].size(), z_len);
				}
				if (dFdX[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, dFdX[i][j].size(), z_len);
				}
				if (d_fd_y[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d_fd_y[i][j].size(), z_len);
				}
				if (d_fd_z[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d_fd_z[i][j].size(), z_len);
				}
				if (d2FdXdY[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d2FdXdY[i][j].size(), z_len);
				}
				if (d2_fd_xd_z[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d2_fd_xd_z[i][j].size(), z_len);
				}
				if (d2_fd_yd_z[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d2_fd_yd_z[i][j].size(), z_len);
				}
				if (d3_fd_xd_yd_z[i][j].size() != z_len)
				{
					throw std::exception("not implemented");
					//throw (hipparchus::exception::Localized_Core_Formats_Type::DIMENSIONS_MISMATCH, d3_fd_xd_yd_z[i][j].size(), z_len);
				}

				const int jp1 = j + 1;
				const double y_r = my_yval[jp1] - my_yval[j];
				const double x_ry_r = x_r * y_r;
				for (int k{}; k < last_k; k++)
				{
					const int& kp1 = k + 1;
					const double zR = my_zval[kp1] - my_zval[k];
					const double x_rzR = x_r * zR;
					const double y_rz_r = y_r * zR;
					const double x_ry_rz_r = x_r * y_rz_r;

					const auto beta = std::vector<double>
					{
						f[i][j][k], f[ip1][j][k], f[i][jp1][k], f[ip1][jp1][k], f[i][j][kp1], f[ip1][j][kp1], f[i][jp1][kp1], f[ip1][jp1][kp1],
						dFdX[i][j][k] * x_r, dFdX[ip1][j][k] * x_r, dFdX[i][jp1][k] * x_r, dFdX[ip1][jp1][k] * x_r, dFdX[i][j][kp1] * x_r, dFdX[ip1][j][kp1] * x_r, dFdX[i][jp1][kp1] * x_r, dFdX[ip1][jp1][kp1] * x_r,
						d_fd_y[i][j][k] * y_r, d_fd_y[ip1][j][k] * y_r, d_fd_y[i][jp1][k] * y_r, d_fd_y[ip1][jp1][k] * y_r, d_fd_y[i][j][kp1] * y_r, d_fd_y[ip1][j][kp1] * y_r, d_fd_y[i][jp1][kp1] * y_r, d_fd_y[ip1][jp1][kp1] * y_r,
						d_fd_z[i][j][k] * zR, d_fd_z[ip1][j][k] * zR, d_fd_z[i][jp1][k] * zR, d_fd_z[ip1][jp1][k] * zR, d_fd_z[i][j][kp1] * zR, d_fd_z[ip1][j][kp1] * zR, d_fd_z[i][jp1][kp1] * zR, d_fd_z[ip1][jp1][kp1] * zR,
						d2FdXdY[i][j][k] * x_ry_r, d2FdXdY[ip1][j][k] * x_ry_r, d2FdXdY[i][jp1][k] * x_ry_r, d2FdXdY[ip1][jp1][k] * x_ry_r, d2FdXdY[i][j][kp1] * x_ry_r, d2FdXdY[ip1][j][kp1] * x_ry_r, d2FdXdY[i][jp1][kp1] * x_ry_r, d2FdXdY[ip1][jp1][kp1] * x_ry_r,
						d2_fd_xd_z[i][j][k] * x_rzR, d2_fd_xd_z[ip1][j][k] * x_rzR, d2_fd_xd_z[i][jp1][k] * x_rzR, d2_fd_xd_z[ip1][jp1][k] * x_rzR, d2_fd_xd_z[i][j][kp1] * x_rzR, d2_fd_xd_z[ip1][j][kp1] * x_rzR, d2_fd_xd_z[i][jp1][kp1] * x_rzR, d2_fd_xd_z[ip1][jp1][kp1] * x_rzR,
						d2_fd_yd_z[i][j][k] * y_rz_r, d2_fd_yd_z[ip1][j][k] * y_rz_r, d2_fd_yd_z[i][jp1][k] * y_rz_r, d2_fd_yd_z[ip1][jp1][k] * y_rz_r, d2_fd_yd_z[i][j][kp1] * y_rz_r, d2_fd_yd_z[ip1][j][kp1] * y_rz_r, d2_fd_yd_z[i][jp1][kp1] * y_rz_r, d2_fd_yd_z[ip1][jp1][kp1] * y_rz_r,
						d3_fd_xd_yd_z[i][j][k] * x_ry_rz_r, d3_fd_xd_yd_z[ip1][j][k] * x_ry_rz_r, d3_fd_xd_yd_z[i][jp1][k] * x_ry_rz_r, d3_fd_xd_yd_z[ip1][jp1][k] * x_ry_rz_r, d3_fd_xd_yd_z[i][j][kp1] * x_ry_rz_r, d3_fd_xd_yd_z[ip1][j][kp1] * x_ry_rz_r, d3_fd_xd_yd_z[i][jp1][kp1] * x_ry_rz_r, d3_fd_xd_yd_z[ip1][jp1][kp1] * x_ry_rz_r, };

					throw std::exception("not implemented");
					//my_splines[i][j][k] = Tricubic_Function(compute_coefficients(beta));
				}
			}
		}
	}

	/**
	 * {@inherit_doc}
	 *
	 * @ if any of the variables is outside its interpolation range.
	 */
	 //override
	double value(const double& x, const double& y, const double& z)
	{
		const int i = search_index(x, my_xval);
		if (i == -1)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, x, xval[0], xval[xval.size() - 1]);
		}
		const int j = search_index(y, my_yval);
		if (j == -1)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, y, yval[0], yval[yval.size() - 1]);
		}
		const int& k = search_index(z, my_zval);
		if (k == -1)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::OUT_OF_RANGE_SIMPLE, z, zval[0], zval[zval.size() - 1]);
		}

		const double xN = (x - my_xval[i]) / (my_xval[i + 1] - my_xval[i]);
		const double yN = (y - my_yval[j]) / (my_yval[j + 1] - my_yval[j]);
		const double zN = (z - my_zval[k]) / (my_zval[k + 1] - my_zval[k]);

		return my_splines[i][j][k].value(xN, yN, zN);
	}

	/**
	 * Indicates whether a point is within the interpolation range.
	 *
	 * @param x First coordinate.
	 * @param y Second coordinate.
	 * @param z Third coordinate.
	 * @return {@code true} if (x, y, z) is a valid point.
	 */
	bool is_valid_point(const double& x, const double& y, const double& z)
	{
		return !(x < my_xval[0] ||
			x > my_xval[my_xval.size() - 1] ||
			y < my_yval[0] ||
			y > my_yval[my_yval.size() - 1] ||
			z < my_zval[0] ||
			z > my_zval[my_zval.size() - 1]);
	}

private:
	/**
	 * @param c Coordinate.
	 * @param val Coordinate samples.
	 * @return the index in {@code val} corresponding to the interval containing {@code c}, or {@code -1}
	 *   if {@code c} is out of the range defined by the end values of {@code val}.
	 */
	int search_index(const double& c, const std::vector<double>& val)
	{
		if (c < val[0])
		{
			return -1;
		}

		const int max = val.size();
		for (int i{ 1 }; i < max; i++)
		{
			if (c <= val[i])
			{
				return i - 1;
			}
		}

		return -1;
	}

	/**
	 * Compute the spline coefficients from the list of function values and
	 * function partial derivatives values at the four corners of a grid
	 * element. They must be specified in the following order:
	 * <ul>
	 *  <li>f(0,0,0)</li>
	 *  <li>f(1,0,0)</li>
	 *  <li>f(0,1,0)</li>
	 *  <li>f(1,1,0)</li>
	 *  <li>f(0,0,1)</li>
	 *  <li>f(1,0,1)</li>
	 *  <li>f(0,1,1)</li>
	 *  <li>f(1,1,1)</li>
	 *
	 *  <li>f<sub>x</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>x</sub>(1,1,1)</li>
	 *
	 *  <li>f<sub>y</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>y</sub>(1,1,1)</li>
	 *
	 *  <li>f<sub>z</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>z</sub>(1,1,1)</li>
	 *
	 *  <li>f<sub>xy</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>xy</sub>(1,1,1)</li>
	 *
	 *  <li>f<sub>xz</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>xz</sub>(1,1,1)</li>
	 *
	 *  <li>f<sub>yz</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>yz</sub>(1,1,1)</li>
	 *
	 *  <li>f<sub>xyz</sub>(0,0,0)</li>
	 *  <li>... <em>(same order as above)</em></li>
	 *  <li>f<sub>xyz</sub>(1,1,1)</li>
	 * </ul>
	 * where the subscripts indicate the partial derivative with respect to
	 * the corresponding variable(s).
	 *
	 * @param beta List of function values and function partial derivatives values.
	 * @return the spline coefficients.
	 */
	std::vector<double> compute_coefficients(const std::vector<double>& beta)
	{
		const int sz{ 64 };
		auto a = std::vector<double>(sz);

		for (int i{}; i < sz; i++)
		{
			double result{};
			const std::vector<double> row = AINV[i];
			for (int j{}; j < sz; j++)
			{
				result += row[j] * beta[j];
			}
			a[i] = result;
		}

		return a;
	}
};