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
  //package org.hipparchus.clustering.distance;

  //import org.hipparchus.exception.;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Math_Arrays;

  /**
   * Calculates the Earh Mover's distance (also known as Wasserstein metric) between two distributions.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Earth_mover's_distance">Earth Mover's distance (Wikipedia)</a>
   *
   */
class Earth_Movers_Distance : Distance_Measure
{
	-5406732779747414922L;

	/** {@inherit_doc} */
	//override
	public double compute(std::vector<double> a, std::vector<double> b)

	{
		Math_Arrays::check_equal_length(a, b);
		double last_distance = 0;
		double total_distance = 0;
		for (int i{}; i < a.size(); i++)
		{
			const double current_distance = (a[i] + last_distance) - b[i];
			total_distance += std::abs(current_distance);
			last_distance = current_distance;
		}
		return total_distance;
	}
}
