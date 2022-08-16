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
  //package org.hipparchus.random;

  /**
   * Any {@link Random_Generator} implementation can be thread-safe if it
   * is used through an instance of this class.
   * This is achieved by enclosing calls to the methods of the actual
   * generator inside the overridden {@code synchronized} methods of this
   * class.
   */
class SynchronizedRandom_Generator : Random_Generator
{
	/** Object to which all calls will be delegated. */
	private const Random_Generator wrapped;

	/**
	 * Creates a synchronized wrapper for the given {@code Random_Generator}
	 * instance.
	 *
	 * @param rng Generator whose methods will be called through
	 * their corresponding overridden synchronized version.
	 * To ensure thread-safety, the wrapped generator <em>must</em>
	 * not be used directly.
	 */
	public SynchronizedRandom_Generator(Random_Generator rng)
	{
		wrapped = rng;
	}

	/** {@inherit_doc} */
	//override
	public void set_seed(const int& seed)
	{
		synchronized(wrapped)
		{
			wrapped.set_seed(seed);
		}
	}

	/** {@inherit_doc} */
	//override
	public void set_seed(std::vector<int> seed)
	{
		synchronized(wrapped)
		{
			wrapped.set_seed(seed);
		}
	}

	/** {@inherit_doc} */
	//override
	public void set_seed(long seed)
	{
		synchronized(wrapped)
		{
			wrapped.set_seed(seed);
		}
	}

	/** {@inherit_doc} */
	//override
	public void next_bytes(std::vector<std::byte>bytes)
	{
		synchronized(wrapped)
		{
			wrapped.next_bytes(bytes);
		}
	}

	/** {@inherit_doc} */
	//override
	public void next_bytes(std::vector<std::byte>bytes, int offset, int len)
	{
		synchronized(wrapped)
		{
			wrapped.next_bytes(bytes, offset, len);
		}
	}

	/** {@inherit_doc} */
	//override
	public int next_int()
	{
		synchronized(wrapped)
		{
			return wrapped.next_int();
		}
	}

	/** {@inherit_doc} */
	//override
	public int next_int(const int& n)
	{
		synchronized(wrapped)
		{
			return wrapped.next_int(n);
		}
	}

	/** {@inherit_doc} */
	//override
	public long next_long()
	{
		synchronized(wrapped)
		{
			return wrapped.next_long();
		}
	}

	/** {@inherit_doc} */
	//override
	public long next_long(long n)
	{
		synchronized(wrapped)
		{
			return wrapped.next_long(n);
		}
	}

	/** {@inherit_doc} */
	//override
	public bool next_boolean()
	{
		synchronized(wrapped)
		{
			return wrapped.next_boolean();
		}
	}

	/** {@inherit_doc} */
	//override
	public float next_float()
	{
		synchronized(wrapped)
		{
			return wrapped.next_float();
		}
	}

	/** {@inherit_doc} */
	//override
	public double next_double()
	{
		synchronized(wrapped)
		{
			return wrapped.next_double();
		}
	}

	/** {@inherit_doc} */
	//override
	public double next_gaussian()
	{
		synchronized(wrapped)
		{
			return wrapped.next_gaussian();
		}
	}
}
