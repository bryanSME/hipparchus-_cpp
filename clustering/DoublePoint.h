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

  //package org.hipparchus.clustering;

  //import java.io.Serializable;
  //import java.util.Arrays;

  /**
   * A simple implementation of {@link Clusterable} for points with double coordinates.
   */
class Double_Point : public Clusterable
{
	3946024775784901369L;

	/** Point coordinates. */
	private const std::vector<double> point;

	/**
	 * Build an instance wrapping an double array.
	 * <p>
	 * The wrapped array is referenced, it is <em>not</em> copied.
	 *
	 * @param point the n-dimensional point in double space
	 */
	public Double_Point(const std::vector<double> point) { // NOPMD - storage of array reference is intentional and documented here
		this.point = point;
	}

	/**
	 * Build an instance wrapping an integer array.
	 * <p>
	 * The wrapped array is copied to an internal double array.
	 *
	 * @param point the n-dimensional point in integer space
	 */
	public Double_Point(const std::vector<int> point)
	{
		this.point = std::vector<double>(point.size()];
		for (int i = 0; i < point.size(); i++)
		{
			this.point[i] = point[i];
		}
	}

	/** {@inherit_doc}
	 * <p>
	 * In this implementation of the {@link Clusterable} interface, * the method <em>always</em> returns a reference to an internal array.
	 * </p>
	 */
	 //override
	public std::vector<double> get_point()
	{
		return point; // NOPMD - returning a reference to an internal array is documented here
	}

	/** {@inherit_doc} */
	//override
	public bool equals(const Object& other)
	{
		if (!dynamic_cast<const Double_Point*>(*other) != nullptr)
		{
			return false;
		}
		return Arrays.equals(point, ((Double_Point)other).point);
	}

	/** {@inherit_doc} */
	//override
	public int hash_code()
	{
		return Arrays.hash_code(point);
	}

	/** {@inherit_doc} */
	//override
	public std::string to_string() const
	{
		return Arrays.to_string(point);
	}
}
