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
  //package org.hipparchus.geometry.partitioning;
#include "../Space.h"
#include <type_traits>
#include "../../core/CalculusFieldElement.hpp"
//import org.hipparchus.geometry.Space;

/** This interface represents the remaining parts of an hyperplane after
 * other parts have been chopped off.

 * <p>sub-hyperplanes are obtained when parts of an {@link
 * Hyperplane hyperplane} are chopped off by other hyperplanes that
 * intersect it. The remaining part is a convex region. Such objects
 * appear in {@link BSP_Tree BSP trees} as the intersection of a cut
 * hyperplane with the convex region which it splits, the chopping
 * hyperplanes are the cut hyperplanes closer to the tree root.</p>

 * <p>
 * Note that this interface is <em>not</em> intended to be implemented
 * by Hipparchus users, it is only intended to be implemented
 * within the library itself. New methods may be added even for minor
 * versions, which breaks compatibility for external implementations.
 * </p>

 * @param <S> Type of the embedding space.

 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Sub_Hyperplane
{
	/** Copy the instance.
	 * <p>The instance created is completely independent of the original
	 * one. A deep copy is used, none of the underlying objects are
	 * shared (except for the nodes attributes and immutable
	 * objects).</p>
	 * @return a sub-hyperplane, copy of the instance
	 */
	virtual Sub_Hyperplane<S> copy_self() = 0;

	/** Get the underlying hyperplane.
	 * @return underlying hyperplane
	 */
	virtual Hyperplane<S> get_hyperplane() = 0;

	/** Check if the instance is empty.
	 * @return true if the instance is empty
	 */
	virtual bool is_empty() = 0;

	/** Get the size of the instance.
	 * @return the size of the instance (this is a length in 1D, an area
	 * in 2D, a volume in 3D ...)
	 */
	virtual double get_size() = 0;

	/** Split the instance in two parts by an hyperplane.
	 * @param hyperplane splitting hyperplane
	 * @return an object containing both the part of the instance
	 * on the plus side of the hyperplane and the part of the
	 * instance on the minus side of the hyperplane
	 */
	virtual Split_Sub_Hyperplane<S> split(Hyperplane<S> hyperplane) = 0;

	/** Compute the union of the instance and another sub-hyperplane.
	 * @param other other sub-hyperplane to union (<em>must</em> be in the
	 * same hyperplane as the instance)
	 * @return a sub-hyperplane, union of the instance and other
	 */
	virtual Sub_Hyperplane<S> reunite(Sub_Hyperplane<S> other) = 0;

	/** Class holding the results of the {@link #split split} method.
	 * @param <U> Type of the embedding space.
	 */
	template<typename U, typename std::enable_if<std::is_base_of<Space, U>::value>::type* = nullptr>
	class Split_Sub_Hyperplane
	{
	private:
		/** Part of the sub-hyperplane on the plus side of the splitting hyperplane. */
		const Sub_Hyperplane<U> my_plus;

		/** Part of the sub-hyperplane on the minus side of the splitting hyperplane. */
		const Sub_Hyperplane<U> my_minus;

	public:
		/** Build a Split_Sub_Hyperplane from its parts.
		 * @param plus part of the sub-hyperplane on the plus side of the
		 * splitting hyperplane
		 * @param minus part of the sub-hyperplane on the minus side of the
		 * splitting hyperplane
		 */
		Split_Sub_Hyperplane(const Sub_Hyperplane<U> plus, const Sub_Hyperplane<U> minus)
			: my_plus{ plus }, my_minus{ minus } {};

		/** Get the part of the sub-hyperplane on the plus side of the splitting hyperplane.
		 * @return part of the sub-hyperplane on the plus side of the splitting hyperplane
		 */
		Sub_Hyperplane<U> get_plus() const
		{
			return my_plus;
		}

		/** Get the part of the sub-hyperplane on the minus side of the splitting hyperplane.
		 * @return part of the sub-hyperplane on the minus side of the splitting hyperplane
		 */
		Sub_Hyperplane<U> get_minus() const
		{
			return my_minus;
		}

		/** Get the side of the split sub-hyperplane with respect to its splitter.
		 * @return {@link Side#PLUS} if only {@link #get_plus()} is neither NULL nor empty, * {@link Side#MINUS} if only {@link #get_minus()} is neither NULL nor empty, * {@link Side#BOTH} if both {@link #get_plus()} and {@link #get_minus()}
		 * are neither NULL nor empty or {@link Side#HYPER} if both {@link #get_plus()} and
		 * {@link #get_minus()} are either NULL or empty
		 */
		Side get_side()
		{
			if (plus != NULL && !plus.is_empty())
			{
				if (minus != NULL && !minus.is_empty())
				{
					return Side.BOTH;
				}
				return Side.PLUS;
			}
			if (minus != NULL && !minus.is_empty())
			{
				return Side.MINUS;
			}
			return Side.HYPER;
		}
	}
};