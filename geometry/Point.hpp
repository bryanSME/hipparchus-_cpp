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
 
#include<type_traits>
#include "Space.h"

  /** This interface represents a generic geometrical point.
   * @param <S> Type of the space.
   * @see Space
   * @see Vector
   */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Point
{
	/** Get the space to which the point belongs.
	 * @return containing space
	 */
	virtual Space get_space() = 0;

	/**
	 * Returns true if any coordinate of this point is NaN; false otherwise
	 * @return  true if any coordinate of this point is NaN; false otherwise
	 */
	virtual bool is_nan() = 0;

	/** Compute the distance between the instance and another point.
	 * @param p second point
	 * @return the distance between the instance and p
	 */
	virtual double distance(Point<S> p) = 0;
};