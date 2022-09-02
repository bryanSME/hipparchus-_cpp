#pragma once

#include <vector>

/** Local iterator for sub-arcs. */
class Sub_Arcs_Iterator : std::vector<double>::Iterator
{
private:
	/** Start of the first arc. */
	const BSP_Tree<Sphere_1D> first_start;

	/** Current node. */
	BSP_Tree<Sphere_1D> current;

	/** Sub-arc no yet returned. */
	std::vector<double> pending;

	/** Walk the tree to select the pending sub-arc.
	 */
	void select_pending()
	{
		// look for the start of the arc
		BSP_Tree<Sphere_1D> start = current;
		while (start != NULL && !is_arc_start(start))
		{
			start = next_internal_node(start);
		}

		if (start == NULL)
		{
			// we have exhausted the iterator
			current = NULL;
			pending = NULL;
			return;
		}

		// look for the end of the arc
		BSP_Tree<Sphere_1D> end = start;
		while (end != NULL && !is_arc_end(end))
		{
			end = next_internal_node(end);
		}

		if (end != NULL)
		{
			// we have identified the arc
			pending = std::vector<double>
			{
				get_angle(start), get_angle(end)
			};

			// prepare search for next arc
			current = end;
		}
		else
		{
			// the const arc wraps around 2\pi, its end is before the first start
			end = first_start;
			while (end != NULL && !is_arc_end(end))
			{
				end = previous_internal_node(end);
			}
			if (end == NULL)
			{
				// this should never happen
				throw Math_Runtime_Exception.create_internal_error();
			}

			// we have identified the last arc
			pending = std::vector<double>
			{
				get_angle(start), get_angle(end) + Math_Utils::TWO_PI
			};

			// there won't be any other arcs
			current = NULL;
		}
	}

public:
	/** Simple constructor.
	 */
	Sub_Arcs_Iterator()
	{
		first_start = get_first_arc_start();
		current = first_start;

		if (first_start == NULL)
		{
			// all the leaf tree nodes share the same inside/outside status
			if ((Boolean)get_first_leaf(get_tree(false)).get_attribute())
			{
				// it is an inside node, it represents the full circle
				pending = std::vector<double>
				{
					0, Math_Utils::TWO_PI
				};
			}
			else
			{
				pending = NULL;
			}
		}
		else
		{
			select_pending();
		}
	}



	/** {@inherit_doc} */
	//override
	bool has_next()
	{
		return pending != NULL;
	}

	/** {@inherit_doc} */
	//override
	std::vector<double> next()
	{
		if (pending == NULL)
		{
			throw No_Such_Element_Exception();
		}
		const std::vector<double> next = pending;
		select_pending();
		return next;
	}

	/** {@inherit_doc} */
	//override
	void remove() const
	{
		throw Unsupported_Operation_Exception();
	}
};