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

#include "../Space.h"
#include <type_traits>

  /** Visitor computing the boundary size.
   * @param <S> Type of the space.
   */
template<typename S, typename std::enable_if<std::is_base_of<Space, S>::value>::type* = nullptr>
class Boundary_size_visitor<S> : public BSP_Tree_Visitor<S>
{
	/** Size of the boundary. */
	private double boundary_size;

	/** Simple constructor.
	 */
	Boundary_size_visitor()
	{
		boundary_size = 0;
	}

	/** {@inherit_doc}*/
	//override
	public Order visit_order(const BSP_Tree<S> node)
	{
		return Order.MINUS_SUB_PLUS;
	}

	/** {@inherit_doc}*/
	//override
	public void visit_internal_node(const BSP_Tree<S> node)
	{
		//@Suppress_Warnings("unchecked")
		const Boundary_Attribute<S> attribute =
			(Boundary_Attribute<S>) node.get_attribute();
		if (attribute.get_plus_outside() != NULL)
		{
			boundary_size += attribute.get_plus_outside().get_size();
		}
		if (attribute.get_plus_inside() != NULL)
		{
			boundary_size += attribute.get_plus_inside().get_size();
		}
	}

	/** {@inherit_doc}*/
	//override
	public void visit_leaf_node(const BSP_Tree<S> node)
	{
	}

	/** Get the size of the boundary.
	 * @return size of the boundary
	 */
	public double get_size()
	{
		return boundary_size;
	}
}
