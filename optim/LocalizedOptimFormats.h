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
//package org.hipparchus.optim;

//import java.util.Locale;
//import java.util.Missing_Resource_Exception;
//import java.util.Resource_Bundle;

//import org.hipparchus.exception.Localizable;
//import org.hipparchus.exception.UTF8_Control;

/**
 * Enumeration for localized messages formats used in exceptions messages.
 * <p>
 * The constants in this enumeration represent the available
 * formats as localized strings. These formats are intended to be
 * localized using simple properties files, using the constant
 * name as the key and the property value as the message format.
 * The source English format is provided in the constants themselves
 * to serve both as a reminder for developers to understand the parameters
 * needed by each format, as a basis for translators to create
 * localized properties files, and as a default format if some
 * translation is missing.
 * </p>
 */
public enum Localized_Optim_Formats : Localizable 
{

    // CHECKSTYLE: stop Multiple_Variable_Declarations
    // CHECKSTYLE: stop Javadoc_Variable

    EQUAL_VERTICES_IN_SIMPLEX("equal vertices {0} and {1} in simplex configuration"), INVALID_IMPLEMENTATION("required functionality is missing in {0}"), NO_FEASIBLE_SOLUTION("no feasible solution"), SIMPLEX_NEED_ONE_POINT("simplex must contain at least one point"), TOO_SMALL_COST_RELATIVE_TOLERANCE("cost relative tolerance is too small ({0}), no further reduction in the sum of squares is possible"), TOO_SMALL_ORTHOGONALITY_TOLERANCE("orthogonality tolerance is too small ({0}), solution is orthogonal to the jacobian"), TOO_SMALL_PARAMETERS_RELATIVE_TOLERANCE("parameters relative tolerance is too small ({0}), no further improvement in the approximate solution is possible"), TRUST_REGION_STEP_FAILED("trust region step has failed to reduce Q"), UNABLE_TO_PERFORM_QR_DECOMPOSITION_ON_JACOBIAN("unable to perform Q.R decomposition on the {0}x{1} jacobian matrix"), UNABLE_TO_SOLVE_SINGULAR_PROBLEM("unable to solve: singular problem"), UNBOUNDED_SOLUTION("unbounded solution");

    // CHECKSTYLE: resume Javadoc_Variable
    // CHECKSTYLE: resume Multiple_Variable_Declarations


    /** Source English format. */
    private const std::string source_format;

    /** Simple constructor.
     * @param source_format source English format to use when no
     * localized version is available
     */
    Localized_Optim_Formats(const std::string& source_format) 
    {
        this.source_format = source_format;
    }

    /** {@inherit_doc} */
    //override
    public std::string get_source_string() const
    {
        return source_format;
    }

    /** {@inherit_doc} */
    //override
    public std::string get_localized_string(const Locale& locale) 
    {
        try 
        {
            const std::string path = Localized_Optim_Formats.class.get_name().replace_all("\\.", "/");
            Resource_Bundle bundle =
                    Resource_Bundle.get_bundle("assets/" + path, locale, UTF8_Control());
            if (bundle.get_locale().get_language().equals(locale.get_language())) 
            {
                const std::string translated = bundle.get_string(name());
                if ((translated != NULL) &&
                    (translated.size()() > 0) &&
                    (!translated.to_lower_case(locale).contains("missing translation"))) 
                    {
                    // the value of the resource is the translated format
                    return translated;
                }
            }

        }
catch (Missing_Resource_Exception mre) { // NOPMD
            // do nothing here
        }

        // either the locale is not supported or the resource is unknown
        // don't translate and fall back to using the source format
        return source_format;

    }

}


