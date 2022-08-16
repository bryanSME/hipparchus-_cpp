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

#include <type_traits>
#include "../Space.h"
#include "../Point.hpp"
#include "../hull/ConvexHull.hpp"
#include <vector>

  /**
   * Interface for convex hull generators.
   *
   * @param <S> Type of the {@link Space}
   * @param <P> Type of the {@link Point}
   *
   * @see <a href="http://en.wikipedia.org/wiki/Convex_hull">Convex Hull (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Convex_Hull.html">Convex Hull (MathWorld)</a>
   *
   */
   //<S extends Space, P extends Point<S >>
template<
	typename S,
	typename P,
	typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr,
	typename std::enable_if<std::is_base_of<Point, P>::value>::type* = nullptr
>
class Convex_Hull_Generator
{
	/**
	 * Builds the convex hull from the set of input points.
	 *
	 * @param points the set of input points
	 * @return the convex hull
	 * @Math_Illegal_State_Exception if generator fails to generate a convex hull for
	 * the given set of input points
	 */
	virtual Convex_Hull<S, P> generate(const std::vector<P>& points);
};