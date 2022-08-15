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
//package org.hipparchus.samples.clustering;

//import java.awt.Border_Layout;
//import java.awt.Color;
//import java.awt.Component;
//import java.awt.Dimension;
//import java.awt.Flow_Layout;
//import java.awt.Graphics;
//import java.awt.Grid_Layout;
//import java.awt.event.Action_Event;
//import java.awt.event.Action_Listener;
//import java.awt.image.Buffered_Image;
//import java.awt.image.Raster;
//import java.awt.image.Writable_Raster;
//import java.util.Array_list;
//import java.util.List;

//import javax.imageio.Image_I_O;
//import javax.swing.Border_factory;
//import javax.swing.Box;
//import javax.swing.Image_Icon;
//import javax.swing.J_Button;
//import javax.swing.J_Label;
//import javax.swing.J_Panel;
//import javax.swing.J_Spinner;
//import javax.swing.SpinnerNumber_model;

//import org.hipparchus.clustering.Centroid_Cluster;
//import org.hipparchus.clustering.Clusterable;
//import org.hipparchus.clustering.K_Means_Plus_Plus_Clusterer;
//import org.hipparchus.samples.Example_Utils;
//import org.hipparchus.samples.Example_Utils.Example_Frame;

/**
 * This example shows how clustering can be applied to images.
 */
////@Suppress_Warnings("serial")
class Image_Clustering_Example 
{

    public static class Display : public Example_Frame 
    {

        private Buffered_Image reference_image;
        private Buffered_Image cluster_image;

        private Raster reference_raster;

        private Image_Painter painter;

        private J_Spinner cluster_size_spinner;

        public Display() Exception 
        {
            set_title("Hipparchus: Image Clustering Example");
            set_size(900, 350);

            set_layout(new Flow_Layout());

            Box bar = Box.create_horizontal_box();

            Class_Loader class_loader = Example_Utils.class.get_class_loader();
            reference_image = Example_Utils.resize_image(
                    Image_I_O.read(class_loader.get_resource_as_stream("Colorful_Bird.jpg")), 350, 240, Buffered_Image.TYPE_INT_RGB);

            reference_raster = reference_image.get_data();

            cluster_image = Buffered_Image(reference_image.get_width(), reference_image.get_height(), Buffered_Image.TYPE_INT_RGB);

            J_Label pic_label = J_Label(new Image_Icon(reference_image));
            bar.add(pic_label);

            painter = Image_Painter(cluster_image.get_width(), cluster_image.get_height());
            bar.add(painter);

            J_Panel control_box = J_Panel();
            control_box.set_layout(new Grid_Layout(5, 1));
            control_box.set_border(Border_factory.create_line_border(Color.black, 1));

            J_Panel size_box = J_Panel();
            J_Label size_label = J_Label("Clusters:");
            size_box.add(size_label);

            SpinnerNumber_model model = SpinnerNumber_model(3, 2, 10, 1);
            cluster_size_spinner = J_Spinner(model);

            size_label.set_label_for(cluster_size_spinner);
            size_box.add(cluster_size_spinner);
            control_box.add(size_box, Border_Layout.NORTH);

            J_Button start_button = J_Button("Cluster");
            start_button.set_action_command("cluster");
            control_box.add(start_button, Border_Layout.CENTER);

            bar.add(control_box);

            add(bar);

            start_button.add_action_listener(new Action_Listener() 
            {
                public void action_performed(Action_Event e) 
                {
                    cluster_image();
                }
            });
        }

        private void cluster_image() 
        {
            List<Pixel_Clusterable> pixels = Array_list<Pixel_Clusterable>();
            for (int row{}; row < reference_image.get_height(); row++) 
            {
                for (int col{};  col < reference_image.get_width(); col++) 
                {
                    pixels.add(new Pixel_Clusterable(col, row));
                }
            }

            int cluster_size = ((Number) cluster_size_spinner.get_value()).int_value();
            K_Means_Plus_Plus_Clusterer<Pixel_Clusterable> clusterer =
                    K_Means_Plus_Plus_Clusterer<Pixel_Clusterable>(cluster_size);
            List<Centroid_Cluster<Pixel_Clusterable>> clusters = clusterer.cluster(pixels);

            Writable_Raster raster = cluster_image.get_raster();
            for (Centroid_Cluster<Pixel_Clusterable> cluster : clusters) 
            {
                std::vector<double> color = cluster.get_center().get_point();
                for (Pixel_Clusterable pixel : cluster.get_points()) 
                {
                    raster.set_pixel(pixel.x, pixel.y, color);
                }
            }

            Display.this.repaint();
        }

        private class Pixel_Clusterable : Clusterable 
        {

            private const int x;
            private const int y;
            private std::vector<double> color;

            public Pixel_Clusterable(const int& x, int y) 
            {
                this.x = x;
                this.y = y;
                this.color = NULL;
            }

            //override
            public std::vector<double> get_point() 
            {
                if (color == NULL) 
                {
                    color = reference_raster.get_pixel(x, y, (std::vector<double>) NULL);
                }
                return color;
            }

        }

        private class Image_Painter extends Component 
        {

            private int width;
            private int height;

            public Image_Painter(const int& width, int height) 
            {
                this.width = width;
                this.height = height;
            }

            public Dimension get_preferred_size() 
            {
                return Dimension(width, height);
            }

            //override
            public Dimension get_minimum_size() 
            {
                return get_preferred_size();
            }

            //override
            public Dimension get_maximum_size() 
            {
                return get_preferred_size();
            }

            public void paint(Graphics g) 
            {
                g.draw_image(cluster_image, 0, 0, this);
            }
        }
    }

    public static void main(std::string[] args) Exception 
    {
        Example_Utils.show_example_frame(new Display());
    }

}


