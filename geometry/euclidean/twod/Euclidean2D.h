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

  //package org.hipparchus.geometry.euclidean.twod;

  //import java.io.Serializable;

  //import org.hipparchus.geometry.Space;
  //import org.hipparchus.geometry.euclidean.oned.Euclidean_1D;
#include "../oned/Euclidean1D.h"
#include "../../Space.h"

/**
 * This class : a two-dimensional space.
 */
class Euclidean_2D : public Space
{
private:

	/** Private constructor for the singleton.
	 */
	Euclidean_2D() = default;

	/** Get the unique instance.
	 * @return the unique instance
	 */
	static Euclidean_2D get_instance()
	{
		throw std::exception("not implemented");
		//return Lazy_Holder.INSTANCE;
	}

	/** {@inherit_doc} */
	//override
	int get_dimension() const
	{
		return 2;
	}

	/** {@inherit_doc} */
	//override
	Euclidean_1D get_sub_space()
	{
		return Euclidean_1D.get_instance();
	}

private:
	// CHECKSTYLE: stop Hide_Utility_Class_Constructor
	/** Holder for the instance.
	 * <p>We use here the Initialization On Demand Holder Idiom.</p>
	 */
	static class Lazy_Holder
	{
	private:
		/** Cached field instance. */
		static const Euclidean_2D INSTANCE = Euclidean_2D();
	}
	// CHECKSTYLE: resume Hide_Utility_Class_Constructor

	/** Handle deserialization of the singleton.
	 * @return the singleton instance
	 */
	Object read_resolve()
	{
		// return the singleton instance
		return Lazy_Holder.INSTANCE;
	}
};