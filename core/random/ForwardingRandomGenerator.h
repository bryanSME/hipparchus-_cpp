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

#include <vector>
#include "RandomGenerator.h"

  /**
   * A random generator which forwards all its method calls to another random generator.
   * <p>
   * Subclasses should //override one or more methods to modify the behavior of the backing
   * random generator as desired per the decorator pattern.
   */
class Forwarding_Random_Generator : public Random_Generator
{
protected:
	/**
	 * Returns the backing delegate instance that methods are forwarded to.
	 * <p>
	 * Concrete subclasses //override this method to supply the instance being decorated.
	 *
	 * @return the delegate instance
	 */
	virtual Random_Generator delegate() = 0;

public:
	/** {@inherit_doc} */
	//override
	void set_seed(const int& seed)
	{
		delegate().set_seed(seed);
	}

	/** {@inherit_doc} */
	//override
	void set_seed(const std::vector<int>& seed)
	{
		delegate().set_seed(seed);
	}

	/** {@inherit_doc} */
	//override
	void set_seed(long seed)
	{
		delegate().set_seed(seed);
	}

	/** {@inherit_doc} */
	//override
	void next_bytes(std::vector<std::byte>bytes)
	{
		delegate().next_bytes(bytes);
	}

	/** {@inherit_doc} */
	//override
	void next_bytes(const std::vector<std::byte>& bytes, const int& offset, const int& length)
	{
		delegate().next_bytes(bytes, offset, length);
	}

	/** {@inherit_doc} */
	//override
	int next_int()
	{
		return delegate().next_int();
	}

	/** {@inherit_doc} */
	//override
	int next_int(const int& n)
	{
		return delegate().next_int(n);
	}

	/** {@inherit_doc} */
	//override
	long next_long()
	{
		return delegate().next_long();
	}

	/** {@inherit_doc} */
	//override
	long next_long(const long& n)
	{
		return delegate().next_long(n);
	}

	/** {@inherit_doc} */
	//override
	bool next_boolean()
	{
		return delegate().next_boolean();
	}

	/** {@inherit_doc} */
	//override
	float next_float()
	{
		return delegate().next_float();
	}

	/** {@inherit_doc} */
	//override
	double next_double()
	{
		return delegate().next_double();
	}

	/** {@inherit_doc} */
	//override
	double next_gaussian()
	{
		return delegate().next_gaussian();
	}
};