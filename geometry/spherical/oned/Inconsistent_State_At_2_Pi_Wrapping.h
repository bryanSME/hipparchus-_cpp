/** Specialized exception for inconsistent BSP tree state inconsistency.
 * <p>
 * This exception is thrown at {@link Arcs_Set} construction time when the
 * {@link org.hipparchus.geometry.partitioning.Region.Location inside/outside}
 * state is not consistent at the 0, \(2 \pi \) crossing.
 * </p>
 */
class Inconsistent_State_At_2_Pi_Wrapping
{
public:
	/** Simple constructor.
	 */
	Inconsistent_State_At_2_Pi_Wrapping()
	{
		super(Localized_Geometry_Formats.INCONSISTENT_STATE_AT_2_PI_WRAPPING);
	}
}