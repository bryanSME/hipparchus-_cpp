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

  //import org.hipparchus.geometry.Point;
  //import org.hipparchus.geometry.Space;

#include <vector>

/** This class represents a ball enclosing some points.
 * @param <S> Space type.
 * @param <P> Point type.
 * @see Space
 * @see Point
 * @see Encloser
 */
class Enclosing_Ball<S extends Space, P extends Point<S>>
{
private:
	/** Center of the ball. */
	const P my_center;

	/** Radius of the ball. */
	const double my_radius;

	/** Support points used to define the ball. */
	const std::vector<P> my_support;

public:
	/** Simple constructor.
	 * @param center center of the ball
	 * @param radius radius of the ball
	 * @param support support points used to define the ball
	 */
	 //@Safe_Varargs
	Enclosing_Ball(const P& center, const double& radius, P ... support)
		:
		my_center{ center },
		my_radius{ radius },
		my_support{ support }
	{}

	/** Get the center of the ball.
	 * @return center of the ball
	 */
	P get_center() const
	{
		return my_center;
	}

	/** Get the radius of the ball.
	 * @return radius of the ball (can be negative if the ball is empty)
	 */
	double get_radius() const
	{
		return my_radius;
	}

	/** Get the support points used to define the ball.
	 * @return support points used to define the ball
	 */
	std::vector<P> get_support() const
	{
		return my_support;
	}

	/** Get the number of support points used to define the ball.
	 * @return number of support points used to define the ball
	 */
	int get_support_size() const
	{
		return my_support.size();
	}

	/** Check if a point is within the ball or at boundary.
	 * @param point point to test
	 * @return true if the point is within the ball or at boundary
	 */
	bool contains(const P& point) const
	{
		return point.distance(my_center) <= my_radius;
	}

	/** Check if a point is within an enlarged ball or at boundary.
	 * @param point point to test
	 * @param margin margin to consider
	 * @return true if the point is within the ball enlarged
	 * by the margin or at boundary
	 */
	bool contains(const P& point, const double& margin) const
	{
		return point.distance(my_center) <= my_radius + margin;
	}
}
