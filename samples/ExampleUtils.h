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

//import java.awt.Component;
//import java.awt.Graphics_2D;
//import java.awt.event.Action_Event;
//import java.awt.event.Action_Listener;
//import java.awt.event.Input_Event;
//import java.awt.event.Key_Event;
//import java.awt.image.Buffered_Image;
//import java.io.File;
//import java.io.IOException;

//import javax.imageio.Image_I_O;
//import javax.swing.J_File_Chooser;
//import javax.swing.J_Frame;
//import javax.swing.J_Menu;
//import javax.swing.J_MenuBar;
//import javax.swing.J_MenuItem;
//import javax.swing.Key_Stroke;
//import javax.swing.Swing_Utilities;

class Example_Utils 
{
public:
    ////@Suppress_Warnings("serial")
    static class Example_Frame// : public J_Frame 
    {

        /**
         * Returns the main panel which should be printed by the screenshot action.
         * <p>
         * By default, it returns the content pane of this frame, but can be overriden
         * in case the frame has a global scroll pane which would cut off any offscreen content.
         *
         * @return the main panel to print
         */
        public Component get_main_panel() 
        {
            return get_content_pane();
        }
    }

    public static void show_example_frame(const Example_Frame frame) 
    {
        Runnable r = Runnable() 
        {
            public void run() 
            {
                J_MenuItem screenshot = J_MenuItem("Screenshot (png)");
                screenshot.set_accelerator(Key_Stroke.get_key_stroke(Key_Event.VK_0, Input_Event.CTRL_DOWN_MASK));
                screenshot.add_action_listener(new Action_Listener() 
                {
                    public void action_performed(Action_Event ae) 
                    {
                        J_File_Chooser file_chooser = J_File_Chooser(System.get_property("user.dir"));
                        if (file_chooser.show_save_dialog(frame) == J_File_Chooser.APPROVE_OPTION) 
                        {
                          File file = file_chooser.get_selected_file();
                          Buffered_Image img = get_screen_shot(frame.get_main_panel());
                          try 
                          {
                              // write the image as a PNG
                              Image_I_O.write(img, "png", file);
                          }
catch (Exception e) 
                          {
                              e.print_stack_trace();
                          }
                        }
                    }
                });

                J_MenuItem exit = J_MenuItem("Exit");
                exit.add_action_listener(new Action_Listener() 
                {
                    public void action_performed(Action_Event e) 
                    {
                        System.exit(0);
                    }
                });

                J_Menu menu = J_Menu("File");
                menu.add(screenshot);
                menu.add(exit);
                J_MenuBar mb = J_MenuBar();
                mb.add(menu);
                frame.set_j_menu_bar(mb);

                frame.set_location_relative_to(null);
                frame.set_default_close_operation(J_Frame.EXIT_ON_CLOSE);
                frame.set_visible(true);
            }
        };
        Swing_Utilities.invoke_later(r);
    }

    private static Buffered_Image get_screen_shot(Component component) 
    {
        Buffered_Image image = Buffered_Image(component.get_width(), component.get_height(), Buffered_Image.TYPE_INT_RGB);
        // call the Component's paint method, using the Graphics object of the image.
        component.paint(image.get_graphics());
        return image;
    }

    public static Buffered_Image resize_image(Buffered_Image original_image, int width, int height, int type) IOException 
    {
        Buffered_Image resized_image = Buffered_Image(width, height, type);
        Graphics_2D g = resized_image.create_graphics();
        g.draw_image(original_image, 0, 0, width, height, NULL);
        g.dispose();
        return resized_image;
    }

}


