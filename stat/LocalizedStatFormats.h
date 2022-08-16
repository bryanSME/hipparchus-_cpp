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
 //package org.hipparchus.stat;

 //import java.util.Locale;
 //import java.util.Missing_Resource_Exception;
 //import java.util.Resource_Bundle;

 //import org.hipparchus.exception.Localizable;
 //import org.hipparchus.exception.UTF8_Control;
#include "../core/exception/Localizable.h"
#include <string>

/**
 * Enumeration for localized messages formats used in exceptions messages.
 * <p>
 * The constants in this enumeration represent the available formats as
 * localized strings. These formats are intended to be localized using simple
 * properties files, using the constant name as the key and the property value
 * as the message format. The source English format is provided in the constants
 * themselves to serve both as a reminder for developers to understand the
 * parameters needed by each format, as a basis for translators to create
 * localized properties files, and as a default format if some translation is
 * missing.
 * </p>
 */
class Localized_Stat_Formats : public Localizable
{
	// CHECKSTYLE: stop Multiple_Variable_Declarations
	// CHECKSTYLE: stop Javadoc_Variable
	//TODO: This needs to extend the enum @see LocalizedCoreFormats::Localized_Core_Format_Types
	enum foo
	{
		/*TIES_ARE_NOT_ALLOWED("Ties are not allowed."),
		INSUFFICIENT_DATA_FOR_T_STATISTIC("insufficient data for t statistic, needs at least 2, got {0}"),
		NOT_ENOUGH_DATA_REGRESSION("the number of observations is not sufficient to conduct regression"),
		INVALID_REGRESSION_OBSERVATION("length of regressor array = {0} does not match the number of variables = {1} in the model"),
		NOT_ENOUGH_DATA_FOR_NUMBER_OF_PREDICTORS("not enough data ({0} rows) for this many predictors ({1} predictors)"),
		NOT_SUPPORTED_NAN_STRATEGY("NaN strategy {0} not supported"),
		NO_REGRESSORS("Regression model must include at least one regressor"),
		COVARIANCE_MATRIX("covariance matrix"),
		OUT_OF_BOUNDS_QUANTILE_VALUE("out of bounds quantile value: {0}, must be in (0, 100]"),
		OUT_OF_BOUNDS_CONFIDENCE_LEVEL("out of bounds confidence level {0}, must be between {1} and {2}"),
		OUT_OF_BOUND_SIGNIFICANCE_LEVEL("out of bounds significance level {0}, must be between {1} and {2}"),
		SIGNIFICANCE_LEVEL("significance level ({0})"),
		TOO_MANY_REGRESSORS("too many regressors ({0}) specified, only {1} in the model"),
		TWO_OR_MORE_CATEGORIES_REQUIRED("two or more categories required, got {0}"),
		TWO_OR_MORE_VALUES_IN_CATEGORY_REQUIRED("two or more values required in each category, one has {0}")*/
	};

private:
	/** Source English format. */
	const std::string my_source_format;

public:

	/**
	 * Simple constructor.
	 *
	 * @param source_format source English format to use when no localized
	 *        version is available
	 */
	Localized_Stat_Formats(const std::string& source_format) : my_source_format{ source_format } {};

	/** {@inherit_doc} */
	//override
	std::string get_source_string() const
	{
		return my_source_format;
	}

	/** {@inherit_doc} */
	//override
	std::string get_localized_string(const std::locale& locale)
	{
		throw std::exception("Not Implemented");
		//try
		//{
		//const std::string path = Localized_Stat_Formats.class.get_name()
		//    .replace_all("\\.", "/");
		//Resource_Bundle bundle = Resource_Bundle
		//    .get_bundle("assets/" + path, locale, UTF8_Control());
		//if (bundle.get_locale().get_language().equals(locale.get_language()))
		//{
		//    const std::string translated = bundle.get_string(name());
		//    if ((translated != NULL) && (translated.size()() > 0) &&
		//        (!translated.to_lower_case(locale).contains("missing translation")))
		//        {
		//        // the value of the resource is the translated format
		//        return translated;
		//    }
		//}

		//}
		//catch (Missing_Resource_Exception mre) { // NOPMD
		//    // do nothing here
		//}

		// either the locale is not supported or the resource is unknown
		// don't translate and fall back to using the source format
		return my_source_format;
	}
};