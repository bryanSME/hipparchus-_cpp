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
//package org.hipparchus.exception;

//import java.io.IOException;
//import java.io.Input_Stream;
//import java.io.Input_StreamReader;
//import java.net.URL;
//import java.net.URL_Connection;
//import java.util.Locale;
//import java.util.Property_resource_Bundle;
//import java.util.Resource_Bundle;
#include <string>

/** Control class loading properties in UTF-8 encoding.
 * <p>
 * This class has been very slightly adapted from Balus_C answer to question: <a
 * href="http://stackoverflow.com/questions/4659929/how-to-use-utf-8-in-resource-properties-with-resourcebundle">
 * How to use UTF-8 in resource properties with Resource_Bundle</a>.
 * </p>
 */
class UTF8_Control // : Resource_Bundle.Control
{
public:
    /** {@inherit_doc} */
    //override
    Resource_Bundle new_bundle(const std::string& base_name, const Locale& locale, const std::string& format, const Class_Loader& loader, const bool reload)
    {
        // The below is a copy of the default implementation.
        const std::string bundle_name = to_bundle_name(base_name, locale);
        const std::string resource_name = to_resource_name(bundle_name, "utf8");
        Resource_Bundle bundle = NULL;
        Input_Stream stream = NULL;
        if (reload)
        {
            const URL url = loader.get_resource(resource_name);
            if (url != NULL)
            {
                const URL_Connection connection = url.open_connection();
                if (connection != NULL)
                {
                    connection.set_use_caches(false);
                    stream = connection.get_input_stream();
                }
            }
        }
        else
        {
            stream = loader.get_resource_as_stream(resource_name);
        }
        if (stream != NULL)
        {
            try
            {
                // Only this line is changed to make it to read properties files as UTF-8.
                bundle = Property_resource_Bundle(new Input_StreamReader(stream, "UTF-8"));
            } constly
            {
                stream.close();
            }
        }
        return bundle;
    }

};