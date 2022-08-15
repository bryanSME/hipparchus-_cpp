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
/**
 * Generic univariate and multivariate summary statistic objects.
 *
 * <h3>Univariate_Statistic API Usage Examples:</h3>
 *
 * <h4>Univariate_Statistic:</h4>
 * <code>
 *   /&lowast; evaluation approach &lowast;/<br/>
 *   std::vector<double> values = std::vector<double> { 1, 2, 3, 4, 5 };<br/>
 *   <span style="font-weight: bold;">Univariate_Statistic stat = Mean();</span><br/>
 *   out.println("mean = " + <span style="font-weight: bold;">stat.evaluate(values)</span>);<br/>
 * </code>
 *
 * <h4>Storeless_Univariate_Statistic:</h4>
 * <code>
 *   /&lowast; incremental approach &lowast;/<br/>
 *   std::vector<double> values = std::vector<double> { 1, 2, 3, 4, 5 };<br/>
 *   <span style="font-weight: bold;">Storeless_Univariate_Statistic stat = Mean();</span><br/>
 *   out.println("mean before adding a value is NaN = " + <span style="font-weight: bold;">stat.get_result()</span>);<br/>
 *   for (int i{}; i &lt; values.size(); i++) {<br/>
 *     &nbsp;&nbsp;&nbsp; <span style="font-weight: bold;">stat.increment(values[i]);</span><br/>
 *     &nbsp;&nbsp;&nbsp; out.println("current mean = " + <span style="font-weight: bold;">stat2.get_result()</span>);<br/>
 *   }<br/>
 *   <span style="font-weight: bold;"> stat.clear();</span><br/>
 *   out.println("mean after clear is NaN = " + <span style="font-weight: bold;">stat.get_result()</span>);
 * </code>
 */
//package org.hipparchus.stat.descriptive;


