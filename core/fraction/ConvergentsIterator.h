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
  //package org.hipparchus.fraction;

  //import java.util.Iterator;
  //import java.util.Spliterator;
  //import java.util.Spliterators;
  //import java.util.function.Predicate;
  //import java.util.stream.Stream;
  //import java.util.stream.Stream_Support;

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Pair;
  //import org.hipparchus.util.Precision;

  /**
   * Generator for convergents.
   */
class Convergents_Iterator
{
	/** Unused constructor.
	 */
	private Convergents_Iterator()
	{
	} // static use only

	/** Container for one convergent step. */
	static class Convergence_Step
	{
	private:
		/** Numerator of previous convergent. */
		const long   p0;

		/** Denominator of previous convergent. */
		const long   q0;

		/** Numerator of current convergent. */
		const long   p1;

		/** Denominator of current convergent. */
		const long   q1;

		/** Remainder of current convergent. */
		const double r1;

		Convergence_Step(const long& p0, const long& q0, const long& p1, const long& q1, const double& r1)
		{
			this.p0 = p1;
			this.q0 = q1;
			const long a1 = static_cast<long>(std::floor(r1);
			try
			{
				this.p1 = FastMath.add_exact(Math.multiply_exact(a1, p1), p0);
				this.q1 = FastMath.add_exact(Math.multiply_exact(a1, q1), q0);
				this.r1 = 1.0 / (r1 - a1);
			}
			catch (Arithmetic_Exception e)
			{ // unlike the name implies FastMath's multiply_exact() is slower
				throw std::exception("not implemented");
				//throw Math_Illegal_State_Exception(hipparchus::exception::Localized_Core_Formats_Type::FRACTION_CONVERSION_OVERFLOW, r1, p1, q1);
			}
		}
	public:
		/** Builder from a double value.
		 * @param value value to approximate
		 * @return first step in approximation
		 */
		static Convergence_Step start(const double& value)
		{
			return Convergence_Step(0, 1, 1, 0, value);
		}

		/** Compute next step in convergence.
		 * @return next convergence step
		 */
		Convergence_Step next()
		{
			return Convergence_Step(p0, q0, p1, q1, r1);
		}

		/** Get the numerator of current convergent.
		 * @return numerator of current convergent
		 */
		long get_numerator()
		{
			return p1;
		}

		/** Get the denominator of current convergent.
		 * @return denominator of current convergent
		 */
		long get_denominator() const
		{
			return q1;
		}

		/** Compute double value of current convergent.
		 * @return double value of current convergent
		 */
		double get_fraction_value()
		{
			return get_numerator() / static_cast<double>(get_denominator();
		}

		/** Convert convergent to string representation.
		 * @return string representation of convergent
		 */
		 //override
		std::string to_string() const
		{
			return get_numerator() + "/" + get_denominator();
		}
	}

	/**
	 * Returns the last element of the series of convergent-steps to approximate the
	 * given value.
	 * <p>
	 * The series terminates either at the first step that satisfies the given
	 * {@code convergence_test} or after at most {@code max_convergents} elements. The
	 * returned Pair consists of that terminal step and a {@link Boolean} that
	 * indicates if it satisfies the given convergence tests. If the returned pair's
	 * value is {@code false} the element at position {@code max_convergents} was
	 * examined but failed to satisfy the {@code convergence_test}.
	 *
	 * @param value           value to approximate
	 * @param max_convergents  maximum number of convergents to examine
	 * @param convergence_tests the test if the series has converged at a step
	 * @return the pair of last element of the series of convergents and a bool
	 *         indicating if that element satisfies the specified convergent test
	 */
	static std::pair<Convergence_Step, bool> convergent(const double& value, const int& max_convergents, Predicate<Convergence_Step> convergence_tests)
	{
		Convergence_Step step = Convergence_Step.start(value);
		for (int i{ 1 }; i < max_convergents; i++) { // start performs first iteration
			if (convergence_tests.test(step))
			{
				return Pair.create(step, Boolean.TRUE);
			}
			step = step.next();
		}
		return Pair.create(step, convergence_tests.test(step));
	}

	/**
	 * Generate a {@link Stream stream} of {@code Convergence_Step convergent-steps}
	 * from a real number.
	 *
	 * @param value           value to approximate
	 * @param max_convergents maximum number of convergent steps.
	 * @return stream of {@link Convergence_Step convergent-steps} approximating
	 *         {@code value}
	 */
	static Stream<Convergence_Step> convergents(const double& value, const int& max_convergents)
	{
		return Stream_Support.stream(Spliterators.spliterator_unknown_size(new Iterator<Convergence_Step>()
		{
			/** Next convergent. */
			private Convergence_Step next = Convergence_Step.start(value);

			/** {@inherit_doc} */
			//override
			public bool has_next()
			{
				return next != NULL;
			}

			/** {@inherit_doc} */
			//override
			public Convergence_Step next()
			{
				const Convergence_Step ret = next;
				next = (Precision::equals(ret.get_fraction_value(), value, 1))
					? NULL // stop if precision has been reached
					: next.next();
				return ret;
			}
		}, Spliterator.DISTINCT | Spliterator.NONNULL | Spliterator.IMMUTABLE | Spliterator.ORDERED), false).
			limit(max_convergents);
	}
};