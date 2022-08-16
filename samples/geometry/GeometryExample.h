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
  //package org.hipparchus.samples.geometry;

  //import java.awt.Border_Layout;
  //import java.awt.Color;
  //import java.awt.Component;
  //import java.awt.event.Action_Event;
  //import java.awt.event.Action_Listener;
  //import java.awt.geom.Point_2D;
  //import java.util.Array_list;
  //import java.util.List;

  //import org.hipparchus.geometry.enclosing.Encloser;
  //import org.hipparchus.geometry.enclosing.Enclosing_Ball;
  //import org.hipparchus.geometry.enclosing.Welzl_Encloser;
  //import org.hipparchus.geometry.euclidean.twod.Disk_Generator;
  //import org.hipparchus.geometry.euclidean.twod.Euclidean_2D;
  //import org.hipparchus.geometry.euclidean.twod.Segment;
  //import org.hipparchus.geometry.euclidean.twod.Vector_2D;
  //import org.hipparchus.geometry.euclidean.twod.hull.Convex_Hull_2D;
  //import org.hipparchus.geometry.euclidean.twod.hull.Convex_Hull_Generator_2D;
  //import org.hipparchus.geometry.euclidean.twod.hull.Monotone_Chain;
  //import org.hipparchus.random.Mersenne_Twister;
  //import org.hipparchus.random.Random_Generator;
  //import org.hipparchus.samples.Example_Utils;
  //import org.hipparchus.samples.Example_Utils.Example_Frame;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Sin_Cos;
  //import org.piccolo2d.P_Camera;
  //import org.piccolo2d.P_Canvas;
  //import org.piccolo2d.P_Node;
  //import org.piccolo2d.event.PBasicInput_Event_Handler;
  //import org.piccolo2d.event.P_Input_Event;
  //import org.piccolo2d.event.P_Mouse_Wheel_Zoom_Event_Handler;
  //import org.piccolo2d.nodes.P_Path;
  //import org.piccolo2d.nodes.P_Text;
#include <vector>
#include "../../geometry/euclidean/twod/Vector2D.h"

/**
 * Simple example illustrating some parts of the geometry //package.
 *
 * TODO:
 *  - select tolerance level
 *  - allow editing of the point set
 */
class Geometry_Example
{
public:
	static std::vector<Vector_2D> create_random_points(const int& size)
	{
		Random_Generator random = Mersenne_Twister();

		// create the cloud container
		List<Vector_2D> points = Array_list<Vector_2D>(size);
		// fill the cloud with a random distribution of points
		for (int i{}; i < size; i++)
		{
			points.add(new Vector_2D(std::round(random.next_double() * 400 + 100), std::round(random.next_double() * 400 + 100)));
		}

		return points;
	}

	public static List<Vector_2D> create_circle(const int& samples)
	{
		List<Vector_2D> points = Array_list<Vector_2D>();
		const Vector_2D center = Vector_2D(300, 300);
		double range = 2.0 * std::numbers::pi;
		double step = range / (samples + 1);
		for (double angle = 0; angle < range; angle += step)
		{
			Vector_2D circle = build_vector(angle);
			points.add(circle.scalar_multiply(200).add(center));
		}

		return points;
	}

	public static List<Vector_2D> create_cross()
	{
		List<Vector_2D> points = Array_list<Vector_2D>();

		for (int i = 100; i < 500; i += 10)
		{
			points.add(new Vector_2D(300, i));
			points.add(new Vector_2D(i, 300));
		}

		return points;
	}

	public static P_Canvas create_canvas()
	{
		const P_Canvas canvas = P_Canvas();
		const P_Camera camera = canvas.get_camera();

		const P_Text tooltip_node = P_Text();
		tooltip_node.set_pickable(false);
		camera.add_child(tooltip_node);

		camera.add_input_event_listener(new PBasicInput_Event_Handler()
			{
				public void mouse_moved(const P_Input_Event event)
				{
					update_tool_tip(event);
				}

				public void mouse_dragged(const P_Input_Event event)
				{
					update_tool_tip(event);
				}

				public void update_tool_tip(const P_Input_Event event)
				{
					const P_Node n = event.get_picked_node();
					const Object object = (Object)n.get_attribute("tooltip");
					if (object != NULL)
					{
						const std::string tooltip_string = object.to_string();
						const Point_2D p = event.get_canvas_position();

						event.get_path().canvas_to_local(p, camera);

						tooltip_node.set_text(tooltip_string);
						tooltip_node.set_offset(p.get_x() + 8, p.get_y() - 8);
					}
	else
					{
						tooltip_node.set_text(null);
					}
				}
			});

		// uninstall default zoom event handler
		canvas.remove_input_event_listener(canvas.get_zoom_event_handler());

		// install mouse wheel zoom event handler
		const P_Mouse_Wheel_Zoom_Event_Handler mouse_wheel_zoom_event_handler = P_Mouse_Wheel_Zoom_Event_Handler();
		canvas.add_input_event_listener(mouse_wheel_zoom_event_handler);

		return canvas;
	}

	/**
	 * Build the 2D vector corresponding to the given angle.
	 * @param alpha angle
	 * @return the corresponding 2D vector
	 */
	private static Vector_2D build_vector(const double& alpha)
	{
		const auto sc = Sin_Cos(alpha);
		return Vector_2D(sc.cos(), sc.sin());
	}

	//@Suppress_Warnings("serial")
	public static class Display : public Example_Frame
	{
		private List<Vector_2D> points;
		private P_Canvas canvas;
		private J_Component container;
		private J_Component control_panel;

		public Display()
		{
			set_title("Hipparchus: Geometry Examples");
			set_size(800, 700);

			container = J_Panel(new Border_Layout());
			canvas = create_canvas();
			container.add(canvas);
			container.set_border(Border_factory.create_line_border(Color.black, 1));

			control_panel = J_Panel();
			J_Button random = J_Button("Randomize");
			control_panel.add(random);

			random.add_action_listener(new Action_Listener()
				{
					//                //override
									public void action_performed(Action_Event e)
									{
										canvas.get_layer().remove_all_children();

										points = create_random_points(1000);
										paint_convex_hull();
									}
				});

			J_Button circle = J_Button("Circle");
			control_panel.add(circle);

			circle.add_action_listener(new Action_Listener()
				{
					//                //override
									public void action_performed(Action_Event e)
									{
										canvas.get_layer().remove_all_children();

										points = create_circle(100);
										paint_convex_hull();
									}
				});

			J_Button cross = J_Button("Cross");
			control_panel.add(cross);

			cross.add_action_listener(new Action_Listener()
				{
					//                //override
									public void action_performed(Action_Event e)
									{
										canvas.get_layer().remove_all_children();

										points = create_cross();
										paint_convex_hull();
									}
				});

			J_Split_Pane splitpane = J_Split_Pane(J_Split_Pane.HORIZONTAL_SPLIT, container, control_panel);
			splitpane.set_divider_location(600);

			add(splitpane);

			points = create_random_points(1000);
			paint_convex_hull();
		}

		//override
		public Component get_main_panel()
		{
			return container;
		}

		public void paint_convex_hull()
		{
			P_Node point_set = P_Node();
			for (Vector_2D point : points)
			{
				const P_Node node = P_Path.create_ellipse(point.get_x() - 1, point.get_y() - 1, 2, 2);
				node.add_attribute("tooltip", point);
				node.set_paint(Color.gray);
				point_set.add_child(node);
			}

			canvas.get_layer().add_child(point_set);

			Convex_Hull_Generator_2D generator = Monotone_Chain(true, 1e-6);
			Convex_Hull_2D hull = generator.generate(points); //Akl_Toussaint_Heuristic.reduce_points(points));

			P_Node hull_node = P_Node();
			for (Vector_2D vertex : hull.get_vertices())
			{
				const P_Path node = P_Path.create_ellipse(vertex.get_x() - 1, vertex.get_y() - 1, 2, 2);
				node.add_attribute("tooltip", vertex);
				node.set_paint(Color.red);
				node.set_stroke_paint(Color.red);
				hull_node.add_child(node);
			}

			for (Segment line : hull.get_line_segments())
			{
				const P_Path node = P_Path.create_line(line.get_start().get_x(), line.get_start().get_y(), line.get_end().get_x(), line.get_end().get_y());
				node.set_pickable(false);
				node.set_paint(Color.red);
				node.set_stroke_paint(Color.red);
				hull_node.add_child(node);
			}

			canvas.get_layer().add_child(hull_node);

			Encloser<Euclidean_2D, Vector_2D> encloser =
				Welzl_Encloser<Euclidean_2D, Vector_2D>(1e-10, Disk_Generator());
			Enclosing_Ball<Euclidean_2D, Vector_2D> ball = encloser.enclose(points);

			const double radius = ball.get_radius();
			P_Path ball_center =
				P_Path.create_ellipse(ball.get_center().get_x() - 1, ball.get_center().get_y() - 1, 2, 2);
			ball_center.set_stroke_paint(Color.blue);
			ball_center.set_paint(Color.blue);
			canvas.get_layer().add_child(0, ball_center);

			P_Path ball_node =
				P_Path.create_ellipse(ball.get_center().get_x() - radius, ball.get_center().get_y() - radius, radius * 2, radius * 2);
			ball_node.set_transparency(1.0f);
			ball_node.set_stroke_paint(Color.blue);
			canvas.get_layer().add_child(0, ball_node);
		}
	}

	public static void main(const std::string[] argv)
	{
		Example_Utils.show_example_frame(new Display());
	}
}
