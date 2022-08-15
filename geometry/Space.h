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

  /** This interface represents a generic space, with affine and vectorial counterparts.
   * @see Vector
   */
class Space
{
	/** Get the dimension of the space.
	 * @return dimension of the space
	 */
	virtual int get_dimension() = 0;

	/** Get the n-1 dimension subspace of this space.
	 * @return n-1 dimension sub-space of this space
	 * @see #get_dimension()
	 * @exception Math_Runtime_Exception for dimension-1 spaces
	 * which do not have sub-spaces
	 */
	virtual Space get_sub_space() = 0;
};