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
  //package org.hipparchus.samples;

  //import java.awt.Basic_Stroke;
  //import java.awt.Color;
  //import java.awt.Component;
  //import java.awt.Font;
  //import java.awt.Grid_Bag_Constraints;
  //import java.awt.Grid_Bag_Layout;
  //import java.awt.Insets;
  //import java.util.Array_list;
  //import java.util.List;

  //import javax.swing.Border_factory;
  //import javax.swing.Box_layout;
  //import javax.swing.J_Component;
  //import javax.swing.J_Label;
  //import javax.swing.J_Panel;
  //import javax.swing.J_Scroll_Pane;

  //import org.hipparchus.distribution.Real_Distribution;
  //import org.hipparchus.distribution.continuous.Beta_Distribution;
  //import org.hipparchus.distribution.continuous.Cauchy_Distribution;
  //import org.hipparchus.distribution.continuous.Chi_Squared_Distribution;
  //import org.hipparchus.distribution.continuous.Exponential_Distribution;
  //import org.hipparchus.distribution.continuous.F_Distribution;
  //import org.hipparchus.distribution.continuous.Gamma_Distribution;
  //import org.hipparchus.distribution.continuous.Levy_Distribution;
  //import org.hipparchus.distribution.continuous.Log_Normal_Distribution;
  //import org.hipparchus.distribution.continuous.Normal_Distribution;
  //import org.hipparchus.distribution.continuous.Pareto_Distribution;
  //import org.hipparchus.distribution.continuous.T_Distribution;
  //import org.hipparchus.distribution.continuous.Weibull_Distribution;
  //import org.hipparchus.samples.Example_Utils.Example_Frame;
  //import org.hipparchus.util.FastMath;

  //import com.xeiam.xchart.Chart;
  //import com.xeiam.xchart.Chart_Builder;
  //import com.xeiam.xchart.Series;
  //import com.xeiam.xchart.Series_marker;
  //import com.xeiam.xchart.Style_Manager.Chart_Type;
  //import com.xeiam.xchart.Style_Manager.Legend_Position;
  //import com.xeiam.xchart.X_Chart_Panel;

#include <string>
#include <vector>

/**
 * Displays pdf/cdf for real distributions.
 */
class Real_Distribution_comparison
{
public:
	static void add_p_d_f_series(Chart chart, Real_Distribution distribution, std::string desc, int lower_bound, int upper_bound)
	{
		// generates Log data
		List<Number> x_data = Array_list<Number>();
		List<Number> y_data = Array_list<Number>();
		int samples = 100;
		double step_size = (upper_bound - lower_bound) / static_cast<double>(samples;
		for (double x = lower_bound; x <= upper_bound; x += step_size)
		{
			try
			{
				double density = distribution.density(x);
				if (!std::isinf(density) && !std::isnan(density))
				{
					x_data.add(x);
					y_data.add(density);
				}
			}
			catch (Exception e)
			{
				// ignore
				// some distributions may reject certain values depending on the parameter settings
			}
		}

		Series series = chart.add_series(desc, x_data, y_data);
		series.set_marker(Series_marker.NONE);
		series.set_line_style(new Basic_Stroke(1.2f));
	}

	static void add_c_d_f_series(Chart chart, Real_Distribution distribution, std::string desc, int lower_bound, int upper_bound)
	{
		// generates Log data
		List<Number> x_data = Array_list<Number>();
		List<Number> y_data = Array_list<Number>();
		int samples = 100;
		double step_size = (upper_bound - lower_bound) / static_cast<double>(samples;
		for (double x = lower_bound; x <= upper_bound; x += step_size)
		{
			double density = distribution.cumulative_probability(x);
			if (!std::isinf(density) && !std::isnan(density))
			{
				x_data.add(x);
				y_data.add(density);
			}
		}

		Series series = chart.add_series(desc, x_data, y_data);
		series.set_marker(Series_marker.NONE);
		series.set_line_style(new Basic_Stroke(1.2f));
	}

	static Chart create_chart(std::string title, int min_x, int max_x, Legend_Position position)
	{
		Chart chart = Chart_Builder().width(235).height(200).build();

		// Customize Chart
		chart.set_chart_title(title);
		chart.get_style_manager().set_chart_title_visible(true);
		chart.get_style_manager().set_chart_title_font(new Font("Arial", Font.PLAIN, 10));
		chart.get_style_manager().set_legend_position(position);
		chart.get_style_manager().set_legend_visible(true);
		chart.get_style_manager().set_legend_font(new Font("Arial", Font.PLAIN, 10));
		chart.get_style_manager().set_legend_padding(6);
		chart.get_style_manager().set_legend_series_line_length(6);
		chart.get_style_manager().set_axis_tick_labels_font(new Font("Arial", Font.PLAIN, 9));

		chart.get_style_manager().set_x_axis_min(min_x);
		chart.get_style_manager().set_x_axis_max(max_x);
		chart.get_style_manager().set_chart_background_color(Color.white);
		chart.get_style_manager().set_chart_padding(4);

		chart.get_style_manager().set_chart_type(Chart_Type.Line);
		return chart;
	}

	static J_Component create_component(std::string distribution_name, int min_x, int max_x, std::string[] series_text, Real_Distribution... series)
	{
		J_Component container = J_Panel();
		container.set_layout(new Box_layout(container, Box_layout.PAGE_AXIS));

		container.add(new J_Label(distribution_name));

		Chart chart = create_chart("PDF", min_x, max_x, Legend_Position.Inside_N_E);
		int i = 0;
		for (Real_Distribution d : series)
		{
			add_p_d_f_series(chart, d, series_text[i++], min_x, max_x);
		}
		container.add(new X_Chart_Panel(chart));

		chart = create_chart("CDF", min_x, max_x, Legend_Position.Inside_S_E);
		i = 0;
		for (Real_Distribution d : series)
		{
			add_c_d_f_series(chart, d, series_text[i++], min_x, max_x);
		}
		container.add(new X_Chart_Panel(chart));

		container.set_border(Border_factory.create_line_border(Color.black, 1));
		return container;
	}

	////@Suppress_Warnings("serial")
	static class Display extends Example_Frame
	{
	private:
		J_Component my_container;

	public:
		Display()
		{
			set_title("Hipparchus: Real distributions overview");
			set_size(1320, 920);

			container = J_Panel();
			container.set_layout(new Grid_Bag_Layout());

			Grid_Bag_Constraints c = Grid_Bag_Constraints();
			c.fill = Grid_Bag_Constraints.VERTICAL;
			c.gridx = 0;
			c.gridy = 0;
			c.insets = Insets(2, 2, 2, 2);

			J_Component comp = NULL;

			comp = create_component("Normal", -5, 5, std::string[]{ "μ=0,σ\u00B2=0.2", "μ=0,σ\u00B2=1", "μ=0,σ\u00B2=5", "μ=-2,σ\u00B2=0.5" }, Normal_Distribution(0, std::sqrt(0.2)), Normal_Distribution(), Normal_Distribution(0, std::sqrt(5)), Normal_Distribution(-2, std::sqrt(0.5)));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Beta", 0, 1, std::string[]{ "α=β=0.5", "α=5,β=1", "α=1,β=3", "α=2,β=2", "α=2,β=5" }, Beta_Distribution(0.5, 0.5), Beta_Distribution(5, 1), Beta_Distribution(1, 3), Beta_Distribution(2, 2), Beta_Distribution(2, 5));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Cauchy", -5, 5, std::string[]{ "x=0,γ=0.5", "x=0,γ=1", "x=0,γ=2", "x=-2,γ=1" }, Cauchy_Distribution(0, 0.5), Cauchy_Distribution(0, 1), Cauchy_Distribution(0, 2), Cauchy_Distribution(-2, 1));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Chi_Squared", 0, 5, std::string[]{ "k=1", "k=2", "k=3", "k=4", "k=6" }, Chi_Squared_Distribution(1), Chi_Squared_Distribution(2), Chi_Squared_Distribution(3), Chi_Squared_Distribution(4), Chi_Squared_Distribution(6));
			container.add(comp, c);

			c.gridy++;
			c.gridx = 0;
			comp = create_component("Exponential", 0, 5, std::string[]{ "λ=0.5", "λ=1", "λ=1.5", "λ=2.5" }, Exponential_Distribution(0.5), Exponential_Distribution(1), Exponential_Distribution(1.5), Exponential_Distribution(2.5));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Fisher-Snedecor", 0, 5, std::string[]{ "d1=1,d2=1", "d1=2,d2=1", "d1=5,d2=2", "d1=100,d2=1", "d1=100,d2=100" }, F_Distribution(1, 1), F_Distribution(2, 1), F_Distribution(5, 2), F_Distribution(100, 1), F_Distribution(100, 100));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Gamma", 0, 20, std::string[]{ "k=1,θ=2", "k=2,θ=2", "k=3,θ=2", "k=5,θ=1", "k=9,θ=0.5" }, Gamma_Distribution(1, 2), Gamma_Distribution(2, 2), Gamma_Distribution(3, 2), Gamma_Distribution(5, 1), Gamma_Distribution(9, 0.5));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Levy", 0, 3, std::string[]{ "c=0.5", "c=1", "c=2", "c=4", "c=8" }, Levy_Distribution(0, 0.5), Levy_Distribution(0, 1), Levy_Distribution(0, 2), Levy_Distribution(0, 4), Levy_Distribution(0, 8));
			container.add(comp, c);

			c.gridy++;
			c.gridx = 0;
			comp = create_component("Log-Normal", 0, 3, std::string[]{ "μ=0,σ\u00B2=10", "μ=0,σ\u00B2=1.5", "μ=0,σ\u00B2=1", "μ=0,σ\u00B2=0.5", "μ=0,σ\u00B2=0.25", "μ=0,σ\u00B2=0.125" }, Log_Normal_Distribution(0, 10), Log_Normal_Distribution(0, 1.5), Log_Normal_Distribution(0, 1), Log_Normal_Distribution(0, 0.5), Log_Normal_Distribution(0, 0.25), Log_Normal_Distribution(0, 0.125));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Pareto", 0, 5, std::string[]{ "x=1,α=1", "x=1,α=2", "x=1,α=3", "x=1,α=10" }, Pareto_Distribution(1, 1), Pareto_Distribution(1, 2), Pareto_Distribution(1, 3), Pareto_Distribution(1, 10));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Student-T", -5, 5, std::string[]{ "df=1", "df=2", "df=5", "df=10000" }, T_Distribution(1), T_Distribution(2), T_Distribution(5), T_Distribution(10000));
			container.add(comp, c);

			c.gridx++;
			comp = create_component("Weibull", 0, 3, std::string[]{ "λ=0.5,k=1", "λ=1,k=1", "λ=1.5,k=1", "λ=5,k=1" }, Weibull_Distribution(0.5, 1), Weibull_Distribution(1, 1), Weibull_Distribution(1.5, 1), Weibull_Distribution(5, 1));
			container.add(comp, c);

			J_Scroll_Pane scroll_pane = J_Scroll_Pane(container);
			add(scroll_pane);
		}

		//override
		Component get_main_panel() const
		{
			return my_container;
		}
	}

	static void main(std::string[] args)
	{
		Example_Utils.show_example_frame(new Display());
	}
}
