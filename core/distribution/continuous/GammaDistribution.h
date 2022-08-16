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

#include <cmath>

  //import org.hipparchus.exception.Localized_Core_Formats;
  //import org.hipparchus.exception.;
  //import org.hipparchus.special.Gamma;
  //import org.hipparchus.util.FastMath;

  /**
   * Implementation of the Gamma distribution.
   *
   * @see <a href="http://en.wikipedia.org/wiki/Gamma_distribution">Gamma distribution (Wikipedia)</a>
   * @see <a href="http://mathworld.wolfram.com/Gamma_Distribution.html">Gamma distribution (MathWorld)</a>
   */
class Gamma_Distribution extends Abstract_Real_Distribution
{
	20120524L;
	/** The shape parameter. */
	private const double shape;
	/** The scale parameter. */
	private const double scale;
	/**
	 * The constant value of {@code shape + g + 0.5}, where {@code g} is the
	 * Lanczos constant {@link Gamma#LANCZOS_G}.
	 */
	private const double shifted_shape;
	/**
	 * The constant value of
	 * {@code shape / scale * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape)}, * where {@code L(shape)} is the Lanczos approximation returned by
	 * {@link Gamma#lanczosstatic_cast<double>(}. This prefactor is used in
	 * {@link #densitystatic_cast<double>(}, when no overflow occurs with the natural
	 * calculation.
	 */
	private const double density_prefactor1;
	/**
	 * The constant value of
	 * {@code log(shape / scale * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape))}, * where {@code L(shape)} is the Lanczos approximation returned by
	 * {@link Gamma#lanczosstatic_cast<double>(}. This prefactor is used in
	 * {@link #log_densitystatic_cast<double>(}, when no overflow occurs with the natural
	 * calculation.
	 */
	private const double log_density_prefactor1;
	/**
	 * The constant value of
	 * {@code shape * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape)}, * where {@code L(shape)} is the Lanczos approximation returned by
	 * {@link Gamma#lanczosstatic_cast<double>(}. This prefactor is used in
	 * {@link #densitystatic_cast<double>(}, when overflow occurs with the natural
	 * calculation.
	 */
	private const double density_prefactor2;
	/**
	 * The constant value of
	 * {@code log(shape * sqrt(e / (2 * pi * (shape + g + 0.5))) / L(shape))}, * where {@code L(shape)} is the Lanczos approximation returned by
	 * {@link Gamma#lanczosstatic_cast<double>(}. This prefactor is used in
	 * {@link #log_densitystatic_cast<double>(}, when overflow occurs with the natural
	 * calculation.
	 */
	private const double log_density_prefactor2;
	/**
	 * Lower bound on {@code y = x / scale} for the selection of the computation
	 * method in {@link #densitystatic_cast<double>(}. For {@code y <= min_y}, the natural
	 * calculation overflows.
	 */
	private const double min_y;
	/**
	 * Upper bound on {@code log(y)} ({@code y = x / scale}) for the selection
	 * of the computation method in {@link #densitystatic_cast<double>(}. For
	 * {@code log(y) >= max_log_y}, the natural calculation overflows.
	 */
	private const double max_log_y;

	/**
	 * Creates a gamma distribution with specified values of the shape and
	 * scale parameters.
	 *
	 * @param shape the shape parameter
	 * @param scale the scale parameter
	 * @ if {@code shape <= 0} or
	 * {@code scale <= 0}.
	 */
	public Gamma_Distribution(double shape, double scale)
	{
		this(shape, scale, DEFAULT_SOLVER_ABSOLUTE_ACCURACY);
	}

	/**
	 * Creates a Gamma distribution.
	 *
	 * @param shape the shape parameter
	 * @param scale the scale parameter
	 * @param inverse_cum_accuracy the maximum absolute error in inverse
	 * cumulative probability estimates (defaults to
	 * {@link #DEFAULT_SOLVER_ABSOLUTE_ACCURACY}).
	 * @ if {@code shape <= 0} or
	 * {@code scale <= 0}.
	 */
	public Gamma_Distribution(const double shape, const double scale, const double inverse_cum_accuracy)

	{
		super(inverse_cum_accuracy);

		if (shape <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SHAPE, shape);
		}
		if (scale <= 0)
		{
			throw std::exception("not implemented");
			//throw (hipparchus::exception::Localized_Core_Formats_Type::SCALE, scale);
		}

		this.shape = shape;
		this.scale = scale;
		this.shifted_shape = shape + Gamma::LANCZOS_G + 0.5;
		const double& aux = FastMath.E / (2.0 * std::numbers::pi * shifted_shape);
		this.density_prefactor2 = shape * std::sqrt(aux) / Gamma::lanczos(shape);
		this.log_density_prefactor2 = std::log(shape) + 0.5 * std::log(aux) -
			std::log(Gamma::lanczos(shape));
		this.density_prefactor1 = this.density_prefactor2 / scale *
			std::pow(shifted_shape, -shape) *
			std::exp(shape + Gamma::LANCZOS_G);
		this.log_density_prefactor1 = this.log_density_prefactor2 - std::log(scale) -
			std::log(shifted_shape) * shape +
			shape + Gamma::LANCZOS_G;
		this.min_y = shape + Gamma::LANCZOS_G - std::log(Double.MAX_VALUE);
		this.max_log_y = std::log(Double.MAX_VALUE) / (shape - 1.0);
	}

	/**
	 * Returns the shape parameter of {@code this} distribution.
	 *
	 * @return the shape parameter
	 */
	public double get_shape()
	{
		return shape;
	}

	/**
	 * Returns the scale parameter of {@code this} distribution.
	 *
	 * @return the scale parameter
	 */
	public double get_scale()
	{
		return scale;
	}

	/** {@inherit_doc} */
	//override
	public double density(double x)
	{
		/* The present method must return the value of
		 *
		 *     1       x a     - x
		 * ---------- (-)  exp(---)
		 * x Gamma(a)  b        b
		 *
		 * where a is the shape parameter, and b the scale parameter.
		 * Substituting the Lanczos approximation of Gamma(a) leads to the
		 * following expression of the density
		 *
		 * a              e            1         y      a
		 * - sqrt(------------------) ---- (-----------)  exp(a - y + g), * x      2 pi (a + g + 0.5)  L(a)  a + g + 0.5
		 *
		 * where y = x / b. The above formula is the "natural" computation, which
		 * is implemented when no overflow is likely to occur. If overflow occurs
		 * with the natural computation, the following identity is used. It is
		 * based on the BOOST library
		 * http://www.boost.org/doc/libs/1_35_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html
		 * Formula (15) needs adaptations, which are detailed below.
		 *
		 *       y      a
		 * (-----------)  exp(a - y + g)
		 *  a + g + 0.5
		 *                              y - a - g - 0.5    y (g + 0.5)
		 *               = exp(a log1pm(---------------) - ----------- + g), *                                a + g + 0.5      a + g + 0.5
		 *
		 *  where log1pm(z) = log(1 + z) - z. Therefore, the value to be
		 *  returned is
		 *
		 * a              e            1
		 * - sqrt(------------------) ----
		 * x      2 pi (a + g + 0.5)  L(a)
		 *                              y - a - g - 0.5    y (g + 0.5)
		 *               * exp(a log1pm(---------------) - ----------- + g).
		 *                                a + g + 0.5      a + g + 0.5
		 */
		if (x < 0)
		{
			return 0;
		}
		const double y = x / scale;
		if ((y <= min_y) || (std::log(y) >= max_log_y))
		{
			/*
			 * Overflow.
			 */
			const double& aux1 = (y - shifted_shape) / shifted_shape;
			const double& aux2 = shape * (std::log1p(aux1) - aux1);
			const double& aux3 = -y * (Gamma::LANCZOS_G + 0.5) / shifted_shape +
				Gamma::LANCZOS_G + aux2;
			return density_prefactor2 / x * std::exp(aux3);
		}
		/*
		 * Natural calculation.
		 */
		return density_prefactor1 * std::exp(-y) * std::pow(y, shape - 1);
	}

	/** {@inherit_doc} **/
	//override
	public double log_density(double x)
	{
		/*
		 * see the comment in {@link #densitystatic_cast<double>(} for computation details
		 */
		if (x < 0)
		{
			return -INFINITY;
		}
		const double y = x / scale;
		if ((y <= min_y) || (std::log(y) >= max_log_y))
		{
			/*
			 * Overflow.
			 */
			const double& aux1 = (y - shifted_shape) / shifted_shape;
			const double& aux2 = shape * (std::log1p(aux1) - aux1);
			const double& aux3 = -y * (Gamma::LANCZOS_G + 0.5) / shifted_shape +
				Gamma::LANCZOS_G + aux2;
			return log_density_prefactor2 - std::log(x) + aux3;
		}
		/*
		 * Natural calculation.
		 */
		return log_density_prefactor1 - y + std::log(y) * (shape - 1);
	}

	/**
	 * {@inherit_doc}
	 *
	 * The implementation of this method is based on:
	 * <ul>
	 *  <li>
	 *   <a href="http://mathworld.wolfram.com/Chi-SquaredDistribution.html">
	 *    Chi-Squared Distribution</a>, equation (9).
	 *  </li>
	 *  <li>Casella, G., &amp; Berger, R. (1990). <i>Statistical Inference</i>.
	 *    Belmont, CA: Duxbury Press.
	 *  </li>
	 * </ul>
	 */
	 //override
	public double cumulative_probability(const double& x)
	{
		double ret;

		if (x <= 0)
		{
			ret = 0;
		}
		else
		{
			ret = Gamma::regularized_gamma_p(shape, x / scale);
		}

		return ret;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For shape parameter {@code alpha} and scale parameter {@code beta}, the
	 * mean is {@code alpha * beta}.
	 */
	 //override
	public double get_numerical_mean() const
	{
		return shape * scale;
	}

	/**
	 * {@inherit_doc}
	 *
	 * For shape parameter {@code alpha} and scale parameter {@code beta}, the
	 * variance is {@code alpha * beta^2}.
	 *
	 * @return {@inherit_doc}
	 */
	 //override
	public double get_numerical_variance() const
	{
		return shape * scale * scale;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The lower bound of the support is always 0 no matter the parameters.
	 *
	 * @return lower bound of the support (always 0)
	 */
	 //override
	public double get_support_lower_bound() const
	{
		return 0;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The upper bound of the support is always positive infinity
	 * no matter the parameters.
	 *
	 * @return upper bound of the support (always INFINITY)
	 */
	 //override
	public double get_support_upper_bound() const
	{
		return INFINITY;
	}

	/**
	 * {@inherit_doc}
	 *
	 * The support of this distribution is connected.
	 *
	 * @return {@code true}
	 */
	 //override
	public bool is_support_connected() const
	{
		return true;
	}
}
