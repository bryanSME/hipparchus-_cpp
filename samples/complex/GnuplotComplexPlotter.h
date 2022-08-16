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

 /** Plotter for complex functions using gnuplot.
  */
  //package org.hipparchus.samples.complex;

  //import java.io.File;
  //import java.io.IOException;
  //import java.io.Print_Stream;
  //import java.nio.charset.Standard_Charsets;
  //import java.util.Array_list;
  //import java.util.List;
  //import java.util.Locale;

  //import org.hipparchus.analysis.integration.Iterative_Legendre_Gauss_Integrator;
  //import org.hipparchus.complex.std::complex<double>;
  //import org.hipparchus.complex.std::complex<double>_Univariate_Integrator;
  //import org.hipparchus.exception.Math_Illegal_State_Exception;
  //import org.hipparchus.special.elliptic.jacobi.Field_Jacobi_Elliptic;
  //import org.hipparchus.special.elliptic.jacobi.Jacobi_Elliptic_Builder;
  //import org.hipparchus.special.elliptic.legendre.Legendre_Elliptic_Integral;
  //import org.hipparchus.util.FastMath;
  //import org.hipparchus.util.Ryu_Double;
#include <complex>
#include <string>
#include <vector>

class Gnuplot_Complex_Plotter
{
private:
	/** Elliptic integrals characteristic. */
	std::complex<double> n;

	/** Elliptic integrals parameter. */
	std::complex<double> m;

	/** Jacobi functions computer. */
	Field_Jacobi_Elliptic<std::complex<double>> jacobi;

	/** Functions to plot. */
	List<Predefined> functions;

	/** Output directory (display to terminal if NULL). */
	File output;

	/** Domain coloring. */
	Domain_coloring coloring;

	/** Output width. */
	int width;

	/** Output height. */
	int height;

	/** Min x. */
	double x_min;

	/** Max x. */
	double x_max;

	/** Min y. */
	double y_min;

	/** Max y. */
	double y_max;

	/** Max z. */
	double z_max;

	/** View X rotation. */
	double view_x_rot;

	/** View Z rotation. */
	double view_z_rot;

	/** Indicator for 3D surfaces. */
	bool use3_d;

	/** Maximum number of integrands evaluations for each integral evaluation. */
	int max_eval;

	/** Integrator for numerically integrated elliptical integrals. */
	const std::complex<double>_Univariate_Integrator integrator;

	/** Default constructor.
	 */
	Gnuplot_Complex_Plotter()
	{
		n = std::complex<double>(3.4, 1.3);
		m = std::complex<double>(0.64);
		jacobi = Jacobi_Elliptic_Builder.build(m);
		functions = Array_list<>();
		output = NULL;
		width = 800;
		height = 800;
		coloring = SawTooth_phaseModuleValue(1.0, 0.7, 1.0, 15);
		x_min = -7;
		x_max = +7;
		y_min = -7;
		y_max = +7;
		z_max = +7;
		view_x_rot = 60;
		view_z_rot = 30;
		use3_d = false;
		max_eval = 100000;
		integrator = std::complex<double>_Univariate_Integrator(new Iterative_Legendre_Gauss_Integrator(24, 1.0e-6, 1.0e-6));
	}

	/** Display usage.
	 * @param status exit code
	 */
	static void usage(const int& status)
	{
		System.err.println("usage: java org.hipparchus.samples.complex.Gnuplot_Complex_Plotter" +
			" [--help]" +
			" [--output-dir directory]" +
			" [--3d]" +
			" [--view x_rot z_rot]" +
			" [--color {classical|enhanced-module|enhanced-phase-module}]" +
			" [--xmin x_min] [--xmax x_max] [--ymin y_min] [--ymax y_max] [--zmax z_max]" +
			" [--m mRe mIm] [--n n_re nIm] [--maxeval max_eval]" +
			" --function {id|sn|cn|dn|cs|...|sin|cos|...} [--function ...]");
		System.exit(status);
	}

	/** Interface for evaluating complex functions. */
	static interface Evaluator
	{
		/** Evaluate complex function.
		 * @param plotter associated plotter
		 * @param z free variable
		 */
		std::complex<double> value(Gnuplot_Complex_Plotter plotter, std::complex<double> z);
	}

public:
	/** Main program.
	 * @param args program arguments
	 */
	static void main(std::vector<std::string>& args)
	{
		const Gnuplot_Complex_Plotter plotter = Gnuplot_Complex_Plotter();
		try
		{
			for (int i{}; i < args.size(); ++i)
			{
				switch (args[i])
				{
				case "--help":
					usage(0);
					break;
				case "--output-dir":
					plotter.output = File(args[++i]);
					if (!(plotter.output.exists() && plotter.output.is_directory() && plotter.output.can_write()))
					{
						System.err.format(Locale.US, "cannot generate output file in %s%n", plotter.output.get_absolute_path());
						System.exit(1);
					}
					break;
				case "--width":
					plotter.width = Integer.parse_int(args[++i]);
					break;
				case "--height":
					plotter.height = Integer.parse_int(args[++i]);
					break;
				case "--color":
					switch (args[++i])
					{
					case "classical":
						plotter.coloring = Continuous_moduleValue(1.0);
						break;
					case "enhanced-module":
						plotter.coloring = Saw_Tooth_Module_Value(1.0);
						break;
					case "enhanced-phase-module":
						plotter.coloring = SawTooth_phaseModuleValue(1.0, 0.7, 1.0, 15);
						break;
					default:
						usage(1);
					}
					break;
				case "--3d":
					plotter.use3_d = true;
					break;
				case "--view":
					plotter.view_x_rot = Double.parse_double(args[++i]);
					plotter.view_z_rot = Double.parse_double(args[++i]);
					break;
				case "--xmin":
					plotter.x_min = Double.parse_double(args[++i]);
					break;
				case "--xmax":
					plotter.x_max = Double.parse_double(args[++i]);
					break;
				case "--ymin":
					plotter.y_min = Double.parse_double(args[++i]);
					break;
				case "--ymax":
					plotter.y_max = Double.parse_double(args[++i]);
					break;
				case "--zmax":
					plotter.z_max = Double.parse_double(args[++i]);
					break;
				case "--m":
				{
					plotter.m = std::complex<double>(Double.parse_double(args[++i]), Double.parse_double(args[++i]));
					plotter.jacobi = Jacobi_Elliptic_Builder.build(plotter.m);
					break;
				}
				case "--n":
				{
					plotter.n = std::complex<double>(Double.parse_double(args[++i]), Double.parse_double(args[++i]));
					break;
				}
				case "--maxeval":
				{
					plotter.max_eval = Integer.parse_int(args[++i]);
					break;
				}
				case "--function":
					try
					{
						plotter.functions.add(Predefined.value_of(args[++i]));
					}
					catch (Illegal_Argument_Exception iae)
					{
						System.err.format(Locale.US, "unknown function %s, known functions:%n", args[i]);
						for (const Predefined predefined : Predefined.values())
						{
							System.err.format(Locale.US, " %s", predefined.name());
						}
						System.err.format(Locale.US, "%n");
						System.exit(1);
					}
					break;
				default:
					usage(1);
				}
			}
			if (plotter.functions.is_empty())
			{
				usage(1);
			}

			plotter.plot();
		}
		catch (Index_Out_Of_Bounds_Exception iobe)
		{
			usage(1);
		}
		catch (IOException ioe)
		{
			System.err.println(ioe.get_localized_message());
			System.exit(1);
		}

		System.exit(0);
	}

	/** Plot the function.
	 * @IOException if gnuplot process cannot be run
	 */
	void plot() IOException
	{
		for (const Predefined predefined : functions)
		{
			const Process_Builder pb = Process_Builder("gnuplot").
				redirect_output(Process_Builder.Redirect.INHERIT).
				redirect_error(Process_Builder.Redirect.INHERIT);
			pb.environment().remove("XDG_SESSION_TYPE");
			const Process gnuplot = pb.start();
			try (Print_Stream out = Print_Stream(gnuplot.get_output_stream(), false, Standard_Charsets.UTF_8.name()))
			{
				if (output == NULL)
				{
					out.format(Locale.US, "set terminal qt size %d, %d title 'complex plotter'%n", width, height);
				}
				else
				{
					out.format(Locale.US, "set terminal pngcairo size %d, %d%n", width, height);
					out.format(Locale.US, "set output '%s'%n", File(output, predefined.name() + ".png").get_absolute_path());
				}
				out.format(Locale.US, "set xrange [%f : %f]%n", x_min, x_max);
				out.format(Locale.US, "set yrange [%f : %f]%n", y_min, y_max);
				if (use3_d)
				{
					out.format(Locale.US, "set zrange [%f : %f]%n", 0.0, z_max);
				}
				out.format(Locale.US, "set xlabel 'Re(z)'%n");
				out.format(Locale.US, "set ylabel 'Im(z)'%n");
				out.format(Locale.US, "set key off%n");
				out.format(Locale.US, "unset colorbox%n");
				out.format(Locale.US, "set title '%s'%n", predefined.title(m));
				out.format(Locale.US, "$data <<EOD%n");

				for (int i{}; i < width; ++i)
				{
					const double x = x_min + i * (x_max - x_min) / (width - 1);
					for (int j{}; j < height; ++j)
					{
						const double y = y_min + j * (y_max - y_min) / (height - 1);
						std::complex<double> z = std::complex<double>.value_of(x, y);
						std::complex<double> fz;
						try
						{
							fz = predefined.evaluator.value(this, z);
						}
						catch (Math_Illegal_State_Exception e)
						{
							fz = std::complex<double>.NaN;
						}
						out.format(Locale.US, "%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f%n", z.get_real_part(), z.get_imaginary_part(), fz.norm(), coloring.hue(fz), coloring.saturation(fz), coloring.value(fz));
					}
					out.format(Locale.US, "%n");
				}
				out.format(Locale.US, "EOD%n");
				if (use3_d)
				{
					out.format(Locale.US, "set view %f, %f%n", view_x_rot, view_z_rot);
					out.format(Locale.US, "splot $data using 1:2:3:(hsv2rgb($4,$5,$6)) with pm3d lc rgb variable%n");
				}
				else
				{
					out.format(Locale.US, "set view map scale 1%n");
					out.format(Locale.US, "splot $data using 1:2:(hsv2rgb($4,$5,$6)) with pm3d lc rgb variable%n");
				}
				if (output == NULL)
				{
					out.format(Locale.US, "pause mouse close%n");
				}
				else
				{
					System.out.format(Locale.US, "output written to %s%n", File(output, predefined.name() + ".png").get_absolute_path());
				}
			}
		}
	}

	/** Predefined complex functions for plotting. */
	private static enum Predefined
	{
		id((plotter, z)->z), sn((plotter, z)->plotter.jacobi.values_n(z).sn()), cn((plotter, z)->plotter.jacobi.values_n(z).cn()), dn((plotter, z)->plotter.jacobi.values_n(z).dn()), cs((plotter, z)->plotter.jacobi.values_s(z).cs()), ds((plotter, z)->plotter.jacobi.values_s(z).ds()), ns((plotter, z)->plotter.jacobi.values_s(z).ns()), dc((plotter, z)->plotter.jacobi.values_c(z).dc()), nc((plotter, z)->plotter.jacobi.values_c(z).nc()), sc((plotter, z)->plotter.jacobi.values_c(z).sc()), nd((plotter, z)->plotter.jacobi.values_d(z).nd()), sd((plotter, z)->plotter.jacobi.values_d(z).sd()), cd((plotter, z)->plotter.jacobi.values_d(z).cd()), arcsn((plotter, z)->plotter.jacobi.arcsn(z)), arccn((plotter, z)->plotter.jacobi.arccn(z)), arcdn((plotter, z)->plotter.jacobi.arcdn(z)), arccs((plotter, z)->plotter.jacobi.arccs(z)), arcds((plotter, z)->plotter.jacobi.arcds(z)), arcns((plotter, z)->plotter.jacobi.arcns(z)), arcdc((plotter, z)->plotter.jacobi.arcdc(z)), arcnc((plotter, z)->plotter.jacobi.arcnc(z)), arcsc((plotter, z)->plotter.jacobi.arcsc(z)), arcnd((plotter, z)->plotter.jacobi.arcnd(z)), arcsd((plotter, z)->plotter.jacobi.arcsd(z)), arccd((plotter, z)->plotter.jacobi.arccd(z)), K((plotter, z)->Legendre_Elliptic_Integral.big_k(z)), K_Prime((plotter, z)->Legendre_Elliptic_Integral.big_k_prime(z)), Fzm((plotter, z)->Legendre_Elliptic_Integral.big_f(z, plotter.m)), integrated_fzm((plotter, z)->Legendre_Elliptic_Integral.big_f(z, plotter.m, plotter.integrator, plotter.max_eval)), E((plotter, z)->Legendre_Elliptic_Integral.big_e(z)), Ezm((plotter, z)->Legendre_Elliptic_Integral.big_e(z, plotter.m)), integrated_ezm((plotter, z)->Legendre_Elliptic_Integral.big_e(z, plotter.m, plotter.integrator, plotter.max_eval)), Pi((plotter, z)->Legendre_Elliptic_Integral.big_pi(plotter.n, z)), Pizm((plotter, z)->Legendre_Elliptic_Integral.big_pi(plotter.n, z, plotter.m)), integrated_pizm((plotter, z)->Legendre_Elliptic_Integral.big_pi(plotter.n, z, plotter.m, plotter.integrator, plotter.max_eval)), sin((plotter, z)->std::sin(z)), cos((plotter, z)->std::cos(z)), tan((plotter, z)->std::tan(z)), asin((plotter, z)->std::asin(z)), acos((plotter, z)->std::acos(z)), atan((plotter, z)->std::atan(z)), sinh((plotter, z)->std::sinh(z)), cosh((plotter, z)->std::cosh(z)), tanh((plotter, z)->std::tanh(z)), asinh((plotter, z)->std::asinh(z)), acosh((plotter, z)->std::acosh(z)), atanh((plotter, z)->std::atanh(z));

	/** Function evaluator. */
	private const Evaluator evaluator;

	/** Simple constructor.
	 * @param evaluator function evaluator
	 */
	private Predefined(const Evaluator evaluator)
	{
		this.evaluator = evaluator;
	}

	/** Get plot title.
	 * @param m elliptic parameter
	 * @return plot title
	 */
	public std::string title(const std::complex<double> m)
	{
		if (name().ends_with("zm"))
		{
			return name().substring(0, name().size()() - 2) +
				"(z, m = " +
				Ryu_Double.double_to_string(m.get_real_part()) +
				(m.get_imaginary() >= 0 ? " + " : " - ") +
				Ryu_Double.double_to_string(std::abs(m.get_imaginary_part())) +
				"i)";
		}
		else
		{
			return name() + "(z)";
		}
	}
	}
};