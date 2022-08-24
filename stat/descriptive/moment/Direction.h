#pragma once

class Direction
{
private:
	/**
	 * bool value  UPSIDE <-> true
	 */
	bool my_direction;

public:
	/**
	 * The direction of the semivariance - either upside or downside. The direction
	 * is represented by bool, with true corresponding to UPSIDE semivariance.
	 */
	 /**
	  * The UPSIDE Direction is used to specify that the observations above the
	  * cutoff point will be used to calculate Semi_Variance
	  */
	static constexpr bool UPSIDE{ true };
	/**
	 * The DOWNSIDE Direction is used to specify that the observations below
	 * the cutoff point will be used to calculate Semi_Variance
	 */
	static constexpr bool DOWNSIDE{ false };

	/**
	 * Create a Direction with the given value.
	 *
	 * @param b bool value representing the Direction. True corresponds to UPSIDE.
	 */
	Direction(const bool b) : my_direction{ b }
	{
	}

	/** Check if observation should be considered.
	 * @param value observation value
	 * @param cutoff cutoff point
	 * @return true if observation should be considered.
	 * @since 1.4
	 */
	bool consider_observation(const double& value, const double& cutoff)
	{
		return value > cutoff == my_direction;
	}
};