/** Class holding the results of the {@link #split split} method.
	 */
class Split
{
private:
	/** Part of the arcs set on the plus side of the splitting arc. */
	const Arcs_Set my_plus;

	/** Part of the arcs set on the minus side of the splitting arc. */
	const Arcs_Set my_minus;

	/** Build a Split from its parts.
	 * @param plus part of the arcs set on the plus side of the
	 * splitting arc
	 * @param minus part of the arcs set on the minus side of the
	 * splitting arc
	 */
	Split(const Arcs_Set& plus, const Arcs_Set& minus) : my_plus{ plus }, my_minus{ minus }
	{
	}

public:
	/** Get the part of the arcs set on the plus side of the splitting arc.
	 * @return part of the arcs set on the plus side of the splitting arc
	 */
	Arcs_Set get_plus() const
	{
		return my_plus;
	}

	/** Get the part of the arcs set on the minus side of the splitting arc.
	 * @return part of the arcs set on the minus side of the splitting arc
	 */
	Arcs_Set get_minus() const
	{
		return my_minus;
	}

	/** Get the side of the split arc with respect to its splitter.
	 * @return {@link Side#PLUS} if only {@link #get_plus()} returns non-null, * {@link Side#MINUS} if only {@link #get_minus()} returns non-null, * {@link Side#BOTH} if both {@link #get_plus()} and {@link #get_minus()}
	 * return non-null or {@link Side#HYPER} if both {@link #get_plus()} and
	 * {@link #get_minus()} return NULL
	 */
	Side get_side()
	{
		if (my_plus != NULL)
		{
			return my_minus != NULL
				? Side::BOTH
				: Side::PLUS;
		}
		return minus != NULL
			? Side::MINUS
			: Side::HYPER;
	}
};