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

  //import java.util.Array_list;
  //import java.util.Iterator;
  //import java.util.List;

  //import org.hipparchus.geometry.Space;

  /** Set of {@link BSP_Tree BSP tree} nodes.
   * @see Boundary_Attribute
   * @param <S> Type of the space.
   */
class Nodes_Set<S extends Space> : Iterable<BSP_Tree<S>>
{
	/** List of sub-hyperplanes. */
	private const List<BSP_Tree<S>> list;

	/** Simple constructor.
	 */
	public Nodes_Set()
	{
		list = Array_list<>();
	}

	/** Add a node if not already known.
	 * @param node node to add
	 */
	public void add(const BSP_Tree<S> node)
	{
		for (const BSP_Tree<S> existing : list)
		{
			if (node == existing)
			{
				// the node is already known, don't add it
				return;
			}
		}

		// the node was not known, add it
		list.add(node);
	}

	/** Add nodes if they are not already known.
	 * @param iterator nodes iterator
	 */
	public void add_all(const Iterable<BSP_Tree<S>> iterator)
	{
		for (const BSP_Tree<S> node : iterator)
		{
			add(node);
		}
	}

	/** {@inherit_doc} */
	//override
	public Iterator<BSP_Tree<S>> iterator()
	{
		return list.iterator();
	}
}
