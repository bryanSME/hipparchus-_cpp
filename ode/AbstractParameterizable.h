#pragma once
/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
//package org.hipparchus.ode;

//import java.util.Array_list;
//import java.util.Collection;
//import java.util.List;

//import org.hipparchus.exception.;

/** This virtual class provides boilerplate parameters list.
 *
 */

class Abstract_Parameterizable : Parameterizable 
{

   /** List of the parameters names. */
    private const List<std::string> parameters_names;

    /** Simple constructor.
     * @param names names of the supported parameters
     */
    protected Abstract_Parameterizable(const std::string ... names) 
    {
        parameters_names = Array_list<>();
        for (const std::string name : names) 
        {
            parameters_names.add(name);
        }
    }

    /** Simple constructor.
     * @param names names of the supported parameters
     */
    protected Abstract_Parameterizable(const Collection<std::string> names) 
    {
        parameters_names = Array_list<>();
        parameters_names.add_all(names);
    }

    /** {@inherit_doc} */
    //override
    public List<std::string> get_parameters_names() 
    {
        return parameters_names;
    }

    /** {@inherit_doc} */
    //override
    public bool is_supported(const std::string name) 
    {
        for (const std::string supported_name : parameters_names) 
        {
            if (supported_name.equals(name)) 
            {
                return true;
            }
        }
        return false;
    }

    /** Check if a parameter is supported and throw an Illegal_Argument_Exception if not.
     * @param name name of the parameter to check
     * @exception  if the parameter is not supported
     * @see #is_supported(std::string)
     */
    public void complain_if_not_supported(const std::string name)
         
        {
        if (!is_supported(name)) 
        {
            throw (Localized_ODE_Formats.UNKNOWN_PARAMETER, name);
        }
    }

}


