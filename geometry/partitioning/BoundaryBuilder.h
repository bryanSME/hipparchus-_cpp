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

  //import org.hipparchus.geometry.Space;

  /** Visitor building boundary shell tree.
   * <p>
   * The boundary shell is represented as {@link Boundary_Attribute boundary attributes}
   * at each internal node.
   * </p>
   * @param <S> Type of the space.
   */
class Boundary_Builder<S extends Space> : BSP_Tree_Visitor<S>
{
	/** {@inherit_doc} */
	//override
	public Order visit_order(BSP_Tree<S> node)
	{
		return Order.PLUS_MINUS_SUB;
	}

	/** {@inherit_doc} */
	//override
	public void visit_internal_node(BSP_Tree<S> node)
	{
		Sub_Hyperplane<S> plus_outside = NULL;
		Sub_Hyperplane<S> plus_inside = NULL;
		Nodes_Set<S>      splitters = NULL;

		// characterize the cut sub-hyperplane, // first with respect to the plus sub-tree
		const Characterization<S> plus_char = Characterization<>(node.get_plus(), node.get_cut().copy_self());

		if (plus_char.touch_outside())
		{
			// plus_char.outside_touching() corresponds to a subset of the cut sub-hyperplane
			// known to have outside cells on its plus side, we want to check if parts
			// of this subset do have inside cells on their minus side
			const Characterization<S> minus_char = Characterization<>(node.get_minus(), plus_char.outside_touching());
			if (minus_char.touch_inside())
			{
				// this part belongs to the boundary, // it has the outside on its plus side and the inside on its minus side
				plus_outside = minus_char.inside_touching();
				splitters = Nodes_Set<>();
				splitters.add_all(minus_char.get_inside_splitters());
				splitters.add_all(plus_char.get_outside_splitters());
			}
		}

		if (plus_char.touch_inside())
		{
			// plus_char.inside_touching() corresponds to a subset of the cut sub-hyperplane
			// known to have inside cells on its plus side, we want to check if parts
			// of this subset do have outside cells on their minus side
			const Characterization<S> minus_char = Characterization<>(node.get_minus(), plus_char.inside_touching());
			if (minus_char.touch_outside())
			{
				// this part belongs to the boundary, // it has the inside on its plus side and the outside on its minus side
				plus_inside = minus_char.outside_touching();
				if (splitters == NULL)
				{
					splitters = Nodes_Set<>();
				}
				splitters.add_all(minus_char.get_outside_splitters());
				splitters.add_all(plus_char.get_inside_splitters());
			}
		}

		if (splitters != NULL)
		{
			// the parent nodes are natural splitters for boundary sub-hyperplanes
			for (BSP_Tree<S> up = node.get_parent(); up != NULL; up = up.get_parent())
			{
				splitters.add(up);
			}
		}

		// set the boundary attribute at non-leaf nodes
		node.set_attribute(new Boundary_Attribute<>(plus_outside, plus_inside, splitters));
	}

	/** {@inherit_doc} */
	//override
	public void visit_leaf_node(BSP_Tree<S> node)
	{
	}
}
