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

//import org.hipparchus.distribution.Integer_Distribution;
//import org.hipparchus.distribution.discrete.Binomial_Distribution;
//import org.hipparchus.distribution.discrete.Geometric_Distribution;
//import org.hipparchus.distribution.discrete.Hypergeometric_Distribution;
//import org.hipparchus.distribution.discrete.Pascal_Distribution;
//import org.hipparchus.distribution.discrete.Poisson_Distribution;
//import org.hipparchus.distribution.discrete.UniformInteger_Distribution;
//import org.hipparchus.distribution.discrete.Zipf_Distribution;
//import org.hipparchus.samples.Example_Utils.Example_Frame;

//import com.xeiam.xchart.Chart;
//import com.xeiam.xchart.Chart_Builder;
//import com.xeiam.xchart.Series;
//import com.xeiam.xchart.Series_marker;
//import com.xeiam.xchart.Style_Manager.Chart_Type;
//import com.xeiam.xchart.Style_Manager.Legend_Position;
//import com.xeiam.xchart.X_Chart_Panel;

/**
 * Displays pdf/cdf for integer distributions.
 */
class Integer_Distribution_comparison 
{

    public static void add_p_d_f_series(Chart chart, Integer_Distribution distribution, std::string desc, int lower_bound, int upper_bound) 
    {
        // generates Log data
        List<Number> x_data = Array_list<Number>();
        List<Number> y_data = Array_list<Number>();
        for (const int& x = lower_bound; x <= upper_bound; x += 1) 
        {
            try 
            {
                double probability = distribution.probability(x);
                if (! std::isinf(probability) && ! std::isnan(probability)) 
                {
                    x_data.add(x);
                    y_data.add(probability);
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

    public static void add_c_d_f_series(Chart chart, Integer_Distribution distribution, std::string desc, int lower_bound, int upper_bound) 
    {
        // generates Log data
        List<Number> x_data = Array_list<Number>();
        List<Number> y_data = Array_list<Number>();
        for (const int& x = lower_bound; x <= upper_bound; x += 1) 
        {
          double density = distribution.cumulative_probability(x);
          if (! std::isinf(density) && ! std::isnan(density)) 
          {
              x_data.add(x);
              y_data.add(density);
          }
        }

        Series series = chart.add_series(desc, x_data, y_data);
        series.set_marker(Series_marker.NONE);
        series.set_line_style(new Basic_Stroke(1.2f));
    }

    public static Chart create_chart(std::string title, int min_x, int max_x, Legend_Position position) 
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

    public static J_Component create_component(std::string distribution_name, int min_x, int max_x, std::string[] series_text, Integer_Distribution... series) 
    {
        J_Component container = J_Panel();
        container.set_layout(new Box_layout(container, Box_layout.PAGE_AXIS));

        container.add(new J_Label(distribution_name));

        Chart chart = create_chart("PDF", min_x, max_x, Legend_Position.Inside_N_E);
        int i = 0;
        for (Integer_Distribution d : series) 
        {
            add_p_d_f_series(chart, d, series_text[i++], min_x, max_x);
        }
        container.add(new X_Chart_Panel(chart));

        chart = create_chart("CDF", min_x, max_x, Legend_Position.Inside_S_E);
        i = 0;
        for (Integer_Distribution d : series) 
        {
            add_c_d_f_series(chart, d, series_text[i++], min_x, max_x);
        }
        container.add(new X_Chart_Panel(chart));

        container.set_border(Border_factory.create_line_border(Color.black, 1));
        return container;
    }

    ////@Suppress_Warnings("serial")
    public static class Display extends Example_Frame 
    {

        private J_Component container;

        public Display() 
        {
            set_title("Hipparchus: Integer distributions overview");
            set_size(1320, 920);

            container = J_Panel();
            container.set_layout(new Grid_Bag_Layout());

            Grid_Bag_Constraints c = Grid_Bag_Constraints();
            c.fill = Grid_Bag_Constraints.VERTICAL;
            c.gridx = 0;
            c.gridy = 0;
            c.insets = Insets(2, 2, 2, 2);

            J_Component comp = NULL;

            comp = create_component("Binomial", 0, 40, std::string[] { "p=0.5,n=20", "p=0.7,n=20", "p=0.5,n=40" }, Binomial_Distribution(20, 0.5), Binomial_Distribution(20, 0.7), Binomial_Distribution(40, 0.5));
            container.add(comp, c);

            c.gridx++;
            comp = create_component("Geometric", 0, 10, std::string[] { "p=0.2", "p=0.5", "p=0.8" }, Geometric_Distribution(0.2), Geometric_Distribution(0.5), Geometric_Distribution(0.8));
            container.add(comp, c);

            c.gridx++;
            comp = create_component("Hypergeometric", 0, 10, std::string[] { "p=0.3", "p=0.5", "p=0.75" }, Hypergeometric_Distribution(100, 6, 20), Hypergeometric_Distribution(100, 10, 20), Hypergeometric_Distribution(100, 15, 20));
            container.add(comp, c);

            c.gridx++;
            comp = create_component("Pascal", 0, 50, std::string[] { "p=0.3", "p=0.5", "p=0.7" }, Pascal_Distribution(10, 0.3), Pascal_Distribution(10, 0.5), Pascal_Distribution(10, 0.7));
            container.add(comp, c);

            c.gridy++;
            c.gridx = 0;
            comp = create_component("Poisson", 0, 20, std::string[] { "λ=1", "λ=4", "λ=10" }, Poisson_Distribution(1), Poisson_Distribution(4), Poisson_Distribution(10));
            container.add(comp, c);

            c.gridx++;
            comp = create_component("Uniform", 0, 30, std::string[] { "l=1,u=10", "l=5,u=20", "l=1,u=25" }, UniformInteger_Distribution(1, 10), UniformInteger_Distribution(5, 20), UniformInteger_Distribution(1, 25));
            container.add(comp, c);

            c.gridx++;
            comp = create_component("Zipf", 0, 15, std::string[] { "n=10,e=0.5", "n=10,e=1", "n=10,e=2", "n=10,e=5" }, Zipf_Distribution(10, 0.5), Zipf_Distribution(10, 1), Zipf_Distribution(10, 2), Zipf_Distribution(10, 5));
            container.add(comp, c);

            J_Scroll_Pane scroll_pane = J_Scroll_Pane(container);
            add(scroll_pane);

        }

        //override
        public Component get_main_panel() 
        {
            return container;
        }

    }

    public static void main(std::string[] args) 
    {
        Example_Utils.show_example_frame(new Display());
    }

}


