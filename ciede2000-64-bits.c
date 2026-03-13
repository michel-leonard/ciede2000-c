#include <math.h>

// This function written in C99 is not affiliated with the CIE (International Commission on Illumination),
// and is released into the public domain. It is provided "as is" without any warranty, express or implied.

// The functional CIE ΔE2000 implementation, which operates on two L*a*b* colors, and returns their difference.
// "l" ranges from 0 to 100, while "a" and "b" are unbounded and commonly clamped to the range of -128 to 127.
static double ciede2000_functional(double l1, double a1, double b1, double l2, double a2, double b2, double kl, double kc, double kh, int canonical) {
	// Working in C99 with the CIEDE2000 color-difference formula.
	// kl, kc, kh are parametric factors to be adjusted according to
	// different viewing parameters such as textures, backgrounds...

	const double pi_1 = 3.14159265358979323846;
	const double pi_3 = 1.04719755119659774615;

	// 1. Compute chroma magnitudes ... a and b usually range from -128 to +127
	const double a1_sq = a1 * a1;
	const double b1_sq = b1 * b1;
	const double c1_orig = sqrt(a1_sq + b1_sq);

	const double a2_sq = a2 * a2;
	const double b2_sq = b2 * b2;
	const double c2_orig = sqrt(a2_sq + b2_sq);

	// 2. Compute chroma mean and apply G compensation
	const double c_avg = 0.5 * (c1_orig + c2_orig);
	const double c_avg_3 = c_avg * c_avg * c_avg;
	const double c_avg_7 = c_avg_3 * c_avg_3 * c_avg;
	const double g_denom = c_avg_7 + 6103515625.0;
	const double g_ratio = c_avg_7 / g_denom;
	const double g_sqrt = sqrt(g_ratio);
	const double g_factor = 1.0 + 0.5 * (1.0 - g_sqrt);

	// 3. Apply G correction to a components, compute corrected chroma
	const double a1_prime = a1 * g_factor;
	const double c1_prime_sq = a1_prime * a1_prime + b1 * b1;
	const double c1_prime = sqrt(c1_prime_sq);
	const double a2_prime = a2 * g_factor;
	const double c2_prime_sq = a2_prime * a2_prime + b2 * b2;
	const double c2_prime = sqrt(c2_prime_sq);

	// 4. Compute hue angles in radians, adjust for negatives and wrap
	const double h1_raw = atan2(b1, a1_prime);
	const double h2_raw = atan2(b2, a2_prime);
	const double h1_adj = h1_raw + (double) (h1_raw < 0.0) * 2.0 * pi_1;
	const double h2_adj = h2_raw + (double) (h2_raw < 0.0) * 2.0 * pi_1;
	const double delta_h = fabs(h1_adj - h2_adj);
	const double h_mean_raw = 0.5 * (h1_adj + h2_adj);
	const double h_diff_raw = 0.5 * (h2_adj - h1_adj);

	// Check if hue mean wraps around pi (180 deg)
	const double hue_wrap = (double) (pi_1 + 1E-14 < delta_h);

	// The part where most programmers get it wrong
	double h_mean;
	if (canonical && pi_1 + 1E-14 < h_mean_raw)
		// Gaurav Sharma, OpenJDK, ...
		h_mean = h_mean_raw - hue_wrap * pi_1;
	else
		// Bruce Lindbloom, Netflix’s VMAF, ...
		h_mean = h_mean_raw + hue_wrap * pi_1;

	// Michel Leonard 2026 - When mean wraps, difference wraps too
	const double h_diff = h_diff_raw + hue_wrap * pi_1;

	// 5. Compute hue rotation correction factor R_T
	const double c_bar = 0.5 * (c1_prime + c2_prime);
	const double c_bar_3 = c_bar * c_bar * c_bar;
	const double c_bar_7 = c_bar_3 * c_bar_3 * c_bar;
	const double rc_denom = c_bar_7 + 6103515625.0;
	const double R_C = sqrt(c_bar_7 / rc_denom);

	const double theta = 36.0 * h_mean - 55.0 * pi_1;
	const double theta_denom = -25.0 * pi_1 * pi_1;
	const double exp_argument = theta * theta / theta_denom;
	const double exp_term = exp(exp_argument);
	const double delta_theta = pi_3 * exp_term;
	const double sin_term = sin(delta_theta);

	// Rotation factor ... cross-effect between chroma and hue
	const double R_T = -2.0 * R_C * sin_term;

	// 6. Compute lightness term ... L nominally ranges from 0 to 100
	const double l_avg = 0.5 * (l1 + l2);
	const double l_delta_sq = (l_avg - 50.0) * (l_avg - 50.0);
	const double l_delta = l2 - l1;

	// Adaptation to the non-linearity of light perception ... S_L
	const double sl_num = 0.015 * l_delta_sq;
	const double sl_denom = sqrt(20.0 + l_delta_sq);
	const double S_L = 1.0 + sl_num / sl_denom;
	const double L_term = l_delta / (kl * S_L);

	// 7. Compute chroma-related trig terms and factor T
	const double trig_1 = 0.17 * sin(h_mean + pi_3);
	const double trig_2 = 0.24 * sin(2.0 * h_mean + 0.5 * pi_1);
	const double trig_3 = 0.32 * sin(3.0 * h_mean + 1.6 * pi_3);
	const double trig_4 = 0.2 * sin(4.0 * h_mean + 0.15 * pi_1);
	const double T = 1.0 - trig_1 + trig_2 + trig_3 - trig_4;

	const double c_sum = c1_prime + c2_prime;
	const double c_product = c1_prime * c2_prime;
	const double c_geo_mean = sqrt(c_product);

	// 8. Compute hue difference and scaling factor S_H
	const double sin_h_diff = sin(h_diff);
	const double S_H = 1.0 + 0.0075 * c_sum * T;
	const double H_term = 2.0 * c_geo_mean * sin_h_diff / (kh * S_H);

	// 9. Compute chroma difference and scaling factor S_C
	const double c_delta = c2_prime - c1_prime;
	const double S_C = 1.0 + 0.0225 * c_sum;
	const double C_term = c_delta / (kc * S_C);

	// 10. Combine lightness, chroma, hue, and interaction terms
	const double L_part = L_term * L_term;
	const double C_part = C_term * C_term;
	const double H_part = H_term * H_term;
	const double interaction = C_term * H_term * R_T;
	const double delta_e_squared = L_part + C_part + H_part + interaction;
	const double delta_e = sqrt(delta_e_squared);

	return delta_e; // Given a tolerance of 3.6e-13
}

// If you remove the constant 1E-14, the code will continue to work, but CIEDE2000
// interoperability between all programming languages will no longer be guaranteed.
