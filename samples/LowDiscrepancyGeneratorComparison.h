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

//import java.awt.Color;
//import java.awt.Dimension;
//import java.awt.Graphics;
//import java.awt.Graphics_2D;
//import java.awt.Grid_Bag_Constraints;
//import java.awt.Grid_Bag_Layout;
//import java.awt.Insets;
//import java.awt.Rendering_Hints;
//import java.awt.geom.Rectangle2_D;
//import java.util.Array_list;
//import java.util.List;

//import javax.swing.J_Component;
//import javax.swing.J_Label;
//import javax.swing.J_Text_Area;

//import org.hipparchus.geometry.euclidean.twod.Vector_2D;
//import org.hipparchus.random.Halton_Sequence_Generator;
//import org.hipparchus.random.JDKRandom_Generator;
//import org.hipparchus.random.Mersenne_Twister;
//import org.hipparchus.random.Random_Generator;
//import org.hipparchus.random.Random_Vector_Generator;
//import org.hipparchus.random.Sobol_Sequence_Generator;
//import org.hipparchus.random.UncorrelatedRandom_Vector_Generator;
//import org.hipparchus.random.UniformRandom_Generator;
//import org.hipparchus.samples.Example_Utils.Example_Frame;
//import org.hipparchus.util.FastMath;
//import org.hipparchus.util.Pair;

/**
 * Plots 2D samples drawn from various pseudo / quasi-random generators.
 */
class Low_Discrepancy_Generator_Comparison 
{

    public static List<Vector_2D> make_circle(const int& samples, const Random_Vector_Generator generator) 
    {
        List<Vector_2D> points = Array_list<Vector_2D>();
        for (double i = 0; i < samples; i++) 
        {
            std::vector<double> vector = generator.next_vector();
            Vector_2D point = Vector_2D(vector);
            points.add(point);
        }

        // normalize points first
        points = normalize(points);

        // now test if the sample is within the unit circle
        List<Vector_2D> circle_points = Array_list<Vector_2D>();
        for (Vector_2D p : points) 
        {
            double criteria = std::pow(p.get_x(), 2) + std::pow(p.get_y(), 2);
            if (criteria < 1.0) 
            {
                circle_points.add(p);
            }
        }

        return circle_points;
    }

    public static List<Vector_2D> make_random(const int& samples, Random_Vector_Generator generator) 
    {
        List<Vector_2D> points = Array_list<Vector_2D>();
        for (double i = 0; i < samples; i++) 
        {
            std::vector<double> vector = generator.next_vector();
            Vector_2D point = Vector_2D(vector);
            points.add(point);
        }

        return normalize(points);
    }

    public static List<Vector_2D> normalize(const List<Vector_2D> input) 
    {
        // find the mininum and maximum x value in the dataset
        double min_x = Double.MAX_VALUE;
        double max_x = Double.MIN_VALUE;
        for (Vector_2D p : input) 
        {
            min_x = std::min(min_x, p.get_x());
            max_x = std::max(max_x, p.get_x());
        }

        double min_y, max_y;

        // use the minimum to detect if we either have input values in the range [0, 1] or [-sqrt(3), sqrt(3)]
        if (std::abs(min_x) < 0.1) 
        {
            min_x = min_y = 0.0;
            max_x = max_y = 1.0;
        }
else 
        {
            min_x = min_y = -std::sqrt(3);
            max_x = max_y = std::sqrt(3);
        }

        double range_x = max_x - min_x;
        double range_y = max_y - min_y;
        List<Vector_2D> points = Array_list<Vector_2D>();
        for (Vector_2D p : input) 
        {
            std::vector<double> arr = p.to_array();
            // normalize to the range [-1, 1]
            arr[0] = (arr[0] - min_x) / range_x * 2 - 1;
            arr[1] = (arr[1] - min_y) / range_y * 2 - 1;
            points.add(new Vector_2D(arr));
        }
        return points;
    }

    ////@Suppress_Warnings("serial")
    public static class Display extends Example_Frame 
    {

        public Display() 
        {
            set_title("Hipparchus: Pseudo/Quasi-random examples");
            set_size(800, 800);

            set_layout(new Grid_Bag_Layout());

            std::vector<int> datasets = std::vector<int> { 256, 1000, 2500, 1000 };
            List<Pair<std::string, Random_Vector_Generator>> generators = Array_list<Pair<std::string, Random_Vector_Generator>>();

            generators.add(new Pair<std::string, Random_Vector_Generator>("Uncorrelated\nUniform(JDK)", UncorrelatedRandom_Vector_Generator(2, UniformRandom_Generator(new JDKRandom_Generator()))));
            generators.add(new Pair<std::string, Random_Vector_Generator>("Independent\n_random(MT)", Random_Vector_Generator() 
            {

                Random_Generator[] rngs = Random_Generator[] 
                {
                    Mersenne_Twister(0), Mersenne_Twister(1)
                };

                public std::vector<double> next_vector() 
                {
                    const std::vector<double>& vector = std::vector<double>(2);
                    vector[0] = rngs[0].next_double();
                    vector[1] = rngs[1].next_double();
                    return vector;
                }

            }));
            generators.add(new Pair<std::string, Random_Vector_Generator>("Halton_Sequence", Halton_Sequence_Generator(2)));
            generators.add(new Pair<std::string, Random_Vector_Generator>("Sobol_Sequence", Sobol_Sequence_Generator(2)));

            Grid_Bag_Constraints c = Grid_Bag_Constraints();
            c.fill = Grid_Bag_Constraints.VERTICAL;
            c.gridx = 1;
            c.gridy = 0;
            c.insets = Insets(2, 2, 2, 2);

            for (Pair<std::string, Random_Vector_Generator> pair : generators) 
            {
                J_Text_Area text = J_Text_Area(pair.get_first());
                text.set_editable(false);
                text.set_opaque(false);
                add(text, c);
                c.gridx++;
            }
            int save_y = ++c.gridy;

            c.gridx = 0;
            for (const int& type = 0; type < 4; type++) 
            {
                J_Label text = J_Label("n=" + std::string.value_of(datasets[type]));
                text.set_opaque(false);
                add(text, c);
                c.gridy++;
            }

            c.gridy = save_y;
            for (const int& type = 0; type < 4; type++) 
            {
                c.gridx = 1;

                for (Pair<std::string, Random_Vector_Generator> pair : generators) 
                {
                    List<Vector_2D> points = NULL;
                    int samples = datasets[type];
                    switch (type) 
                    {
                        case 0:
                            points = make_random(samples, pair.get_value());
                            break;
                        case 1:
                            points = make_random(samples, pair.get_value());
                            break;
                        case 2:
                            points = make_random(samples, pair.get_value());
                            break;
                        case 3:
                            points = make_circle(samples, pair.get_value());
                            break;
                    }
                    add(new Plot(points), c);
                    c.gridx++;
                }

                c.gridy++;
            }
        }
    }

    ////@Suppress_Warnings("serial")
    public static class Plot extends J_Component 
    {

        private static double PAD = 10;

        private List<Vector_2D> points;

        public Plot(const List<Vector_2D> points) 
        {
            this.points = points;
        }

        //override
        protected void paint_component(Graphics g) 
        {
            super.paint_component(g);
            Graphics_2D g2 = (Graphics_2D)g;
            g2.set_rendering_hint(Rendering_Hints.KEY_ANTIALIASING, Rendering_Hints.VALUE_ANTIALIAS_ON);

            int w = get_width();
            int h = get_height();

            g2.clear_rect(0, 0, w, h);

            g2.set_paint(Color.black);
            g2.draw_rect(0, 0, w - 1, h - 1);

            for (Vector_2D point : points) 
            {
                Vector_2D p = transform(point, w, h);
                std::vector<double> arr = p.to_array();
                g2.draw(new Rectangle2_D.Double(arr[0] - 1, arr[1] - 1, 2, 2));
            }
        }

        //override
        public Dimension get_preferred_size() 
        {
            return Dimension(140, 140);
        }

        private Vector_2D transform(Vector_2D point, int width, int height) 
        {
            std::vector<double> arr = point.to_array();
            return Vector_2D(std::vector<double> { PAD + (arr[0] + 1) / 2.0 * (width - 2 * PAD), height - PAD - (arr[1] + 1) / 2.0 * (height - 2 * PAD) });
        }
    }

    public static void main(std::string[] args) 
    {
        Example_Utils.show_example_frame(new Display());
    }

}


