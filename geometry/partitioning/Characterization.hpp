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

//import org.hipparchus.exception.Math_Runtime_Exception;
//import org.hipparchus.geometry.Space;

/** Cut sub-hyperplanes characterization with respect to inside/outside cells.
 * @see Boundary_Builder
 * @param <S> Type of the space.
 */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Characterization
{
private:
	/** Part of the cut sub-hyperplane that touch outside cells. */
	Sub_Hyperplane<S> my_outside_touching;

	/** Part of the cut sub-hyperplane that touch inside cells. */
	Sub_Hyperplane<S> my_inside_touching;

	/** Nodes that were used to split the outside touching part. */
	const Nodes_Set<S> my_outside_splitters;

	/** Nodes that were used to split the inside touching part. */
	const Nodes_Set<S> my_inside_splitters;

	/** Filter the parts of an hyperplane belonging to the boundary.
	 * <p>The filtering consist in splitting the specified
	 * sub-hyperplane into several parts lying in inside and outside
	 * cells of the tree. The principle is to call this method twice for
	 * each cut sub-hyperplane in the tree, once on the plus node and
	 * once on the minus node. The parts that have the same flag
	 * (inside/inside or outside/outside) do not belong to the boundary
	 * while parts that have different flags (inside/outside or
	 * outside/inside) do belong to the boundary.</p>
	 * @param node current BSP tree node
	 * @param sub sub-hyperplane to characterize
	 * @param splitters nodes that did split the current one
	 */
	void characterize(const BSP_Tree<S>& node, const Sub_Hyperplane<S>& sub, const List<BSP_Tree<S>>& splitters)
	{
		if (node.get_cut() == NULL)
		{
			// we have reached a leaf node
			if (const auto inside = static_cast<bool>(node.get_attribute()); inside)
			{
				add_inside_touching(sub, splitters);
			}
			else
			{
				add_outside_touching(sub, splitters);
			}
			return;
		}

		const Hyperplane<S> hyperplane = node.get_cut().get_hyperplane();
		const Sub_Hyperplane.Split_Sub_Hyperplane<S> split = sub.split(hyperplane);
		switch (split.get_side())
		{
		case PLUS:
			characterize(node.get_plus(), sub, splitters);
			break;
		case MINUS:
			characterize(node.get_minus(), sub, splitters);
			break;
		case BOTH:
			splitters.add(node);
			characterize(node.get_plus(), split.get_plus(), splitters);
			characterize(node.get_minus(), split.get_minus(), splitters);
			splitters.remove(splitters.size() - 1);
			break;
		default:
			// this should not happen
			throw Math_Runtime_Exception.create_internal_error();
		}
	}

	/** Add a part of the cut sub-hyperplane known to touch an outside cell.
	 * @param sub part of the cut sub-hyperplane known to touch an outside cell
	 * @param splitters sub-hyperplanes that did split the current one
	 */
	void add_outside_touching(const Sub_Hyperplane<S>& sub, const List<BSP_Tree<S>>& splitters)
	{
		my_outside_touching = my_outside_touching == NULL
			? sub
			: my_outside_touching.reunite(sub);
		my_outside_splitters.add_all(splitters);
	}

	/** Add a part of the cut sub-hyperplane known to touch an inside cell.
	 * @param sub part of the cut sub-hyperplane known to touch an inside cell
	 * @param splitters sub-hyperplanes that did split the current one
	 */
	void add_inside_touching(const Sub_Hyperplane<S>& sub, const List<BSP_Tree<S>>& splitters)
	{
		my_inside_touching = my_inside_touching == NULL
			? sub
			: my_inside_touching.reunite(sub);

		my_inside_splitters.add_all(splitters);
	}

public:
	/** Simple constructor.
	 * <p>Characterization consists in splitting the specified
	 * sub-hyperplane into several parts lying in inside and outside
	 * cells of the tree. The principle is to compute characterization
	 * twice for each cut sub-hyperplane in the tree, once on the plus
	 * node and once on the minus node. The parts that have the same flag
	 * (inside/inside or outside/outside) do not belong to the boundary
	 * while parts that have different flags (inside/outside or
	 * outside/inside) do belong to the boundary.</p>
	 * @param node current BSP tree node
	 * @param sub sub-hyperplane to characterize
	 */
	Characterization(const BSP_Tree<S>& node, const Sub_Hyperplane<S>& sub)
	{
		my_outside_touching = NULL;
		my_inside_touching = NULL;
		my_outside_splitters = Nodes_Set<>();
		my_inside_splitters = Nodes_Set<>();
		characterize(node, sub, Array_list<>());
	}

	/** Check if the cut sub-hyperplane touches outside cells.
	 * @return true if the cut sub-hyperplane touches outside cells
	 */
	bool touch_outside()
	{
		return my_outside_touching != NULL && !my_outside_touching.empty();
	}

	/** Get all the parts of the cut sub-hyperplane known to touch outside cells.
	 * @return parts of the cut sub-hyperplane known to touch outside cells
	 * (may be NULL or empty)
	 */
	Sub_Hyperplane<S> outside_touching()
	{
		return my_outside_touching;
	}

	/** Get the nodes that were used to split the outside touching part.
	 * <p>
	 * Splitting nodes are internal nodes (i.e. they have a non-null
	 * cut sub-hyperplane).
	 * </p>
	 * @return nodes that were used to split the outside touching part
	 */
	Nodes_Set<S> get_outside_splitters() const
	{
		return my_outside_splitters;
	}

	/** Check if the cut sub-hyperplane touches inside cells.
	 * @return true if the cut sub-hyperplane touches inside cells
	 */
	bool touch_inside() const
	{
		return my_inside_touching != NULL && !my_inside_touching.empty();
	}

	/** Get all the parts of the cut sub-hyperplane known to touch inside cells.
	 * @return parts of the cut sub-hyperplane known to touch inside cells
	 * (may be NULL or empty)
	 */
	Sub_Hyperplane<S> inside_touching() const
	{
		return my_inside_touching;
	}

	/** Get the nodes that were used to split the inside touching part.
	 * <p>
	 * Splitting nodes are internal nodes (i.e. they have a non-null
	 * cut sub-hyperplane).
	 * </p>
	 * @return nodes that were used to split the inside touching part
	 */
	Nodes_Set<S> get_inside_splitters() const
	{
		return my_inside_splitters;
	}
};