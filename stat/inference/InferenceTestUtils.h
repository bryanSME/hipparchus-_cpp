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
//package org.hipparchus.stat.inference;

//import java.util.Collection;

//import org.hipparchus.distribution.Real_Distribution;
//import org.hipparchus.exception.;
//import org.hipparchus.exception.Math_Illegal_State_Exception;
//import org.hipparchus.exception.Null_Argument_Exception;
//import org.hipparchus.stat.descriptive.Statistical_Summary;

/**
 * A collection of static methods to create inference test instances or to
 * perform inference tests.
 */
class Inference_Test_Utils  
{

    /** Singleton T_Test instance. */
    private static const T_Test T_TEST = T_Test();

    /** Singleton Chi_Square_Test instance. */
    private static const Chi_Square_Test CHI_SQUARE_TEST = Chi_Square_Test();

    /** Singleton OneWay_anova instance. */
    private static const OneWay_anova ONE_WAY_ANANOVA = OneWay_anova();

    /** Singleton G-Test instance. */
    private static const G_Test G_TEST = G_Test();

    /** Singleton K-S test instance */
    private static const Kolmogorov_Smirnov_Test KS_TEST = Kolmogorov_Smirnov_Test();

    /**
     * Prevent instantiation.
     */
    private Inference_Test_Utils() 
    {
        super();
    }

    // CHECKSTYLE: stop Javadoc_Method_Check

    /**
     * @see T_Test#homoscedastic_t(std::vector<double>, std::vector<double>)
     */
    public static double homoscedastic_t(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception 
        {
        return T_TEST.homoscedastic_t(sample1, sample2);
    }

    /**
     * @see T_Test#homoscedastic_t(Statistical_Summary, Statistical_Summary)
     */
    public static double homoscedastic_t(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception 
        {
        return T_TEST.homoscedastic_t(sample_stats1, sample_stats2);
    }

    /**
     * @see T_Test#homoscedastic_t_test(std::vector<double>, std::vector<double>, double)
     */
    public static bool homoscedastic_t_test(const std::vector<double> sample1, const std::vector<double> sample2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.homoscedastic_t_test(sample1, sample2, alpha);
    }

    /**
     * @see T_Test#homoscedastic_t_test(std::vector<double>, std::vector<double>)
     */
    public static double homoscedastic_t_test(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.homoscedastic_t_test(sample1, sample2);
    }

    /**
     * @see T_Test#homoscedastic_t_test(Statistical_Summary, Statistical_Summary)
     */
    public static double homoscedastic_t_test(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.homoscedastic_t_test(sample_stats1, sample_stats2);
    }

    /**
     * @see T_Test#paired_t(std::vector<double>, std::vector<double>)
     */
    public static double paired_t(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception 
        {
        return T_TEST.paired_t(sample1, sample2);
    }

    /**
     * @see T_Test#paired_t_test(std::vector<double>, std::vector<double>, double)
     */
    public static bool paired_t_test(const std::vector<double> sample1, const std::vector<double> sample2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.paired_t_test(sample1, sample2, alpha);
    }

    /**
     * @see T_Test#paired_t_test(std::vector<double>, std::vector<double>)
     */
    public static double paired_t_test(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.paired_t_test(sample1, sample2);
    }

    /**
     * @see T_Test#t(double, std::vector<double>)
     */
    public static double t(const double& mu, const std::vector<double> observed)
        , Null_Argument_Exception 
        {
        return T_TEST.t(mu, observed);
    }

    /**
     * @see T_Test#t(double, Statistical_Summary)
     */
    public static double t(const double& mu, const Statistical_Summary sample_stats)
        , Null_Argument_Exception 
        {
        return T_TEST.t(mu, sample_stats);
    }

    /**
     * @see T_Test#t(std::vector<double>, std::vector<double>)
     */
    public static double t(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception 
        {
        return T_TEST.t(sample1, sample2);
    }

    /**
     * @see T_Test#t(Statistical_Summary, Statistical_Summary)
     */
    public static double t(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception 
        {
        return T_TEST.t(sample_stats1, sample_stats2);
    }

    /**
     * @see T_Test#t_test(double, std::vector<double>, double)
     */
    public static bool t_test(const double& mu, const std::vector<double> sample, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(mu, sample, alpha);
    }

    /**
     * @see T_Test#t_test(double, std::vector<double>)
     */
    public static double t_test(const double& mu, const std::vector<double> sample)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(mu, sample);
    }

    /**
     * @see T_Test#t_test(double, Statistical_Summary, double)
     */
    public static bool t_test(const double& mu, const Statistical_Summary sample_stats, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(mu, sample_stats, alpha);
    }

    /**
     * @see T_Test#t_test(double, Statistical_Summary)
     */
    public static double t_test(const double& mu, const Statistical_Summary sample_stats)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(mu, sample_stats);
    }

    /**
     * @see T_Test#t_test(std::vector<double>, std::vector<double>, double)
     */
    public static bool t_test(const std::vector<double> sample1, const std::vector<double> sample2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(sample1, sample2, alpha);
    }

    /**
     * @see T_Test#t_test(std::vector<double>, std::vector<double>)
     */
    public static double t_test(const std::vector<double> sample1, const std::vector<double> sample2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(sample1, sample2);
    }

    /**
     * @see T_Test#t_test(Statistical_Summary, Statistical_Summary, double)
     */
    public static bool t_test(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(sample_stats1, sample_stats2, alpha);
    }

    /**
     * @see T_Test#t_test(Statistical_Summary, Statistical_Summary)
     */
    public static double t_test(const Statistical_Summary sample_stats1, const Statistical_Summary sample_stats2)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return T_TEST.t_test(sample_stats1, sample_stats2);
    }

    /**
     * @see Chi_Square_Test#chi_square(std::vector<double>, long[])
     */
    public static double chi_square(const std::vector<double> expected, const std::vector<long> observed)
         
        {
        return CHI_SQUARE_TEST.chi_square(expected, observed);
    }

    /**
     * @see Chi_Square_Test#chi_square(long[][])
     */
    public static double chi_square(const std::vector<std::vector<long>> counts)
        , Null_Argument_Exception 
        {
        return CHI_SQUARE_TEST.chi_square(counts);
    }

    /**
     * @see Chi_Square_Test#chi_square_test(std::vector<double>, long[], double)
     */
    public static bool chi_square_test(const std::vector<double> expected, const std::vector<long> observed, const double& alpha)
        , Math_Illegal_State_Exception 
        {
        return CHI_SQUARE_TEST.chi_square_test(expected, observed, alpha);
    }

    /**
     * @see Chi_Square_Test#chi_square_test(std::vector<double>, long[])
     */
    public static double chi_square_test(const std::vector<double> expected, const std::vector<long> observed)
        , Math_Illegal_State_Exception 
        {
        return CHI_SQUARE_TEST.chi_square_test(expected, observed);
    }

    /**
     * @see Chi_Square_Test#chi_square_test(long[][], double)
     */
    public static bool chi_square_test(const std::vector<std::vector<long>> counts, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return CHI_SQUARE_TEST.chi_square_test(counts, alpha);
    }

    /**
     * @see Chi_Square_Test#chi_square_test(long[][])
     */
    public static double chi_square_test(const std::vector<std::vector<long>> counts)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return CHI_SQUARE_TEST.chi_square_test(counts);
    }

    /**
     * @see Chi_Square_Test#chi_square_data_sets_comparison(long[], long[])
     */
    public static double chi_square_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2)
         
        {
        return CHI_SQUARE_TEST.chi_square_data_sets_comparison(observed1, observed2);
    }

    /**
     * @see Chi_Square_Test#chi_square_test_data_sets_comparison(long[], long[])
     */
    public static double chi_square_test_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2)
        , Math_Illegal_State_Exception 
        {
        return CHI_SQUARE_TEST.chi_square_test_data_sets_comparison(observed1, observed2);
    }

    /**
     * @see Chi_Square_Test#chi_square_test_data_sets_comparison(long[], long[], double)
     */
    public static bool chi_square_test_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2, const double& alpha)
        , Math_Illegal_State_Exception 
        {
        return CHI_SQUARE_TEST.chi_square_test_data_sets_comparison(observed1, observed2, alpha);
    }

    /**
     * @see OneWay_anova#anova_f_value(Collection)
     */
    public static double one_way_anova_f_value(const Collection<std::vector<double>> category_data)
        , Null_Argument_Exception 
        {
        return ONE_WAY_ANANOVA.anova_f_value(category_data);
    }

    /**
     * @see OneWay_anova#anova_p_value(Collection)
     */
    public static double one_way_anova_p_value(const Collection<std::vector<double>> category_data)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return ONE_WAY_ANANOVA.anova_p_value(category_data);
    }

    /**
     * @see OneWay_anova#anova_test(Collection,double)
     */
    public static bool one_way_anova_test(const Collection<std::vector<double>> category_data, const double& alpha)
        , Null_Argument_Exception, Math_Illegal_State_Exception 
        {
        return ONE_WAY_ANANOVA.anova_test(category_data, alpha);
    }

     /**
     * @see G_Test#g(std::vector<double>, long[])
     */
    public static double g(const std::vector<double> expected, const std::vector<long> observed)
         
        {
        return G_TEST.g(expected, observed);
    }

    /**
     * @see G_Test#g_test( std::vector<double>,  std::vector<long> )
     */
    public static double g_test(const std::vector<double> expected, const std::vector<long> observed)
        , Math_Illegal_State_Exception 
        {
        return G_TEST.g_test(expected, observed);
    }

    /**
     * @see G_Test#g_test_intrinsic(std::vector<double>, std::vector<long> )
     */
    public static double g_test_intrinsic(const std::vector<double> expected, const std::vector<long> observed)
        , Math_Illegal_State_Exception 
        {
        return G_TEST.g_test_intrinsic(expected, observed);
    }

     /**
     * @see G_Test#g_test( std::vector<double>,long[],double)
     */
    public static bool g_test(const std::vector<double> expected, const std::vector<long> observed, const double& alpha)
        , Math_Illegal_State_Exception 
        {
        return G_TEST.g_test(expected, observed, alpha);
    }

    /**
     * @see G_Test#g_data_sets_comparison(long[], long[])
     */
    public static double g_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2)
         
        {
        return G_TEST.g_data_sets_comparison(observed1, observed2);
    }

    /**
     * @see G_Test#root_log_likelihood_ratio(long, long, long, long)
     */
    public static double root_log_likelihood_ratio(const long k11, const long k12, const long k21, const long k22)
         
        {
        return G_TEST.root_log_likelihood_ratio(k11, k12, k21, k22);
    }


    /**
     * @see G_Test#g_test_data_sets_comparison(long[], long[])
     */
    public static double g_test_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2)
        , Math_Illegal_State_Exception 
        {
        return G_TEST.g_test_data_sets_comparison(observed1, observed2);
    }

    /**
     * @see G_Test#g_test_data_sets_comparison(long[],long[],double)
     */
    public static bool g_test_data_sets_comparison(const std::vector<long> observed1, const std::vector<long> observed2, const double& alpha)
        , Math_Illegal_State_Exception 
        {
        return G_TEST.g_test_data_sets_comparison(observed1, observed2, alpha);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_statistic(Real_Distribution, std::vector<double>)
     */
    public static double kolmogorov_smirnov_statistic(Real_Distribution dist, std::vector<double> data)
            , Null_Argument_Exception 
            {
        return KS_TEST.kolmogorov_smirnov_statistic(dist, data);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_test(Real_Distribution, std::vector<double>)
     */
    public static double kolmogorov_smirnov_test(Real_Distribution dist, std::vector<double> data)
            , Null_Argument_Exception 
            {
        return KS_TEST.kolmogorov_smirnov_test(dist, data);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_test(Real_Distribution, std::vector<double>, bool)
     */
    public static double kolmogorov_smirnov_test(Real_Distribution dist, std::vector<double> data, bool strict)
            , Null_Argument_Exception 
            {
        return KS_TEST.kolmogorov_smirnov_test(dist, data, strict);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_test(Real_Distribution, std::vector<double>, double)
     */
    public static bool kolmogorov_smirnov_test(Real_Distribution dist, std::vector<double> data, double alpha)
            , Null_Argument_Exception 
            {
        return KS_TEST.kolmogorov_smirnov_test(dist, data, alpha);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_statistic(std::vector<double>, std::vector<double>)
     */
    public static double kolmogorov_smirnov_statistic(std::vector<double> x, std::vector<double> y)
            , Null_Argument_Exception 
            {
        return KS_TEST.kolmogorov_smirnov_statistic(x, y);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_test(std::vector<double>, std::vector<double>)
     */
    public static double kolmogorov_smirnov_test(std::vector<double> x, std::vector<double> y)
            , Null_Argument_Exception 
            {
        return KS_TEST.kolmogorov_smirnov_test(x, y);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#kolmogorov_smirnov_test(std::vector<double>, std::vector<double>, bool)
     */
    public static double kolmogorov_smirnov_test(std::vector<double> x, std::vector<double> y, bool strict)
            , Null_Argument_Exception  
            {
        return KS_TEST.kolmogorov_smirnov_test(x, y, strict);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#exact_p(double, int, int, bool)
     */
    public static double exact_p(double d, int m, int n, bool strict) 
    {
        return KS_TEST.exact_p(d, n, m, strict);
    }

    /**
     * @see Kolmogorov_Smirnov_Test#approximate_p(double, int, int)
     */
    public static double approximate_p(double d, int n, int m) 
    {
        return KS_TEST.approximate_p(d, n, m);
    }

    // CHECKSTYLE: resume Javadoc_Method_Check

}


