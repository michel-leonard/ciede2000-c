#include <quadmath.h>

// This function written in C99 is not affiliated with the CIE (International Commission on Illumination),
// and is released into the public domain. It is provided "as is" without any warranty, express or implied.

// The functional CIE ΔE2000 implementation, which operates on two L*a*b* colors, and returns their difference.
// "l" ranges from 0 to 100, while "a" and "b" are unbounded and commonly clamped to the range of -128 to 127.
static __float128 ciede2000_functional(__float128 l1, __float128 a1, __float128 b1, __float128 l2, __float128 a2, __float128 b2, __float128 kl, __float128 kc, __float128 kh, int canonical) {
	// Working in C99 with the CIEDE2000 color-difference formula.
	// kl, kc, kh are parametric factors to be adjusted according to
	// different viewing parameters such as textures, backgrounds...

	const __float128 pi_1 = 3.1415926535897932384626433832795028841972q;
	const __float128 pi_3 = 1.0471975511965977461542144610931676280657q;

	// 1. Compute chroma magnitudes ... a and b usually range from -128 to +127
	const __float128 a1_sq = a1 * a1;
	const __float128 b1_sq = b1 * b1;
	const __float128 c1_orig = sqrtq(a1_sq + b1_sq);

	const __float128 a2_sq = a2 * a2;
	const __float128 b2_sq = b2 * b2;
	const __float128 c2_orig = sqrtq(a2_sq + b2_sq);

	// 2. Compute chroma mean and apply G compensation
	const __float128 c_avg = 0.5q * (c1_orig + c2_orig);
	const __float128 c_avg_3 = c_avg * c_avg * c_avg;
	const __float128 c_avg_7 = c_avg_3 * c_avg_3 * c_avg;
	const __float128 g_denom = c_avg_7 + 6103515625.0q;
	const __float128 g_ratio = c_avg_7 / g_denom;
	const __float128 g_sqrt = sqrtq(g_ratio);
	const __float128 g_factor = 1.0q + 0.5q * (1.0q - g_sqrt);

	// 3. Apply G correction to a components, compute corrected chroma
	const __float128 a1_prime = a1 * g_factor;
	const __float128 c1_prime_sq = a1_prime * a1_prime + b1 * b1;
	const __float128 c1_prime = sqrtq(c1_prime_sq);
	const __float128 a2_prime = a2 * g_factor;
	const __float128 c2_prime_sq = a2_prime * a2_prime + b2 * b2;
	const __float128 c2_prime = sqrtq(c2_prime_sq);

	// 4. Compute hue angles in radians, adjust for negatives and wrap
	const __float128 h1_raw = atan2q(b1, a1_prime);
	const __float128 h2_raw = atan2q(b2, a2_prime);
	const __float128 h1_adj = h1_raw + (__float128) (h1_raw < 0.0q) * 2.0q * pi_1;
	const __float128 h2_adj = h2_raw + (__float128) (h2_raw < 0.0q) * 2.0q * pi_1;
	const __float128 delta_h = fabsq(h1_adj - h2_adj);
	const __float128 h_mean_raw = 0.5q * (h1_adj + h2_adj);
	const __float128 h_diff_raw = 0.5q * (h2_adj - h1_adj);

	// Check if hue mean wraps around pi (180 deg)
	const __float128 hue_wrap = (__float128) (pi_1 + 1E-14q < delta_h);

	// The part where most programmers get it wrong
	__float128 h_mean;
	if (canonical && pi_1 + 1E-14q < h_mean_raw)
		// Gaurav Sharma, OpenJDK, ...
		h_mean = h_mean_raw - hue_wrap * pi_1;
	else
		// Bruce Lindbloom, Netflix’s VMAF, ...
		h_mean = h_mean_raw + hue_wrap * pi_1;

	// Michel Leonard 2026 - When mean wraps, difference wraps too
	const __float128 h_diff = h_diff_raw + hue_wrap * pi_1;

	// 5. Compute hue rotation correction factor R_T
	const __float128 c_bar = 0.5q * (c1_prime + c2_prime);
	const __float128 c_bar_3 = c_bar * c_bar * c_bar;
	const __float128 c_bar_7 = c_bar_3 * c_bar_3 * c_bar;
	const __float128 rc_denom = c_bar_7 + 6103515625.0q;
	const __float128 R_C = sqrtq(c_bar_7 / rc_denom);

	const __float128 theta = 36.0q * h_mean - 55.0q * pi_1;
	const __float128 theta_denom = -25.0q * pi_1 * pi_1;
	const __float128 exp_argument = theta * theta / theta_denom;
	const __float128 exp_term = expq(exp_argument);
	const __float128 delta_theta = pi_3 * exp_term;
	const __float128 sin_term = sinq(delta_theta);

	// Rotation factor ... cross-effect between chroma and hue
	const __float128 R_T = -2.0q * R_C * sin_term;

	// 6. Compute lightness term ... L nominally ranges from 0 to 100
	const __float128 l_avg = 0.5q * (l1 + l2);
	const __float128 l_delta_sq = (l_avg - 50.0q) * (l_avg - 50.0q);
	const __float128 l_delta = l2 - l1;

	// Adaptation to the non-linearity of light perception ... S_L
	const __float128 sl_num = 0.015q * l_delta_sq;
	const __float128 sl_denom = sqrtq(20.0q + l_delta_sq);
	const __float128 S_L = 1.0q + sl_num / sl_denom;
	const __float128 L_term = l_delta / (kl * S_L);

	// 7. Compute chroma-related trig terms and factor T
	const __float128 trig_1 = 0.17q * sinq(h_mean + pi_3);
	const __float128 trig_2 = 0.24q * sinq(2.0q * h_mean + 0.5q * pi_1);
	const __float128 trig_3 = 0.32q * sinq(3.0q * h_mean + 1.6q * pi_3);
	const __float128 trig_4 = 0.2q * sinq(4.0q * h_mean + 0.15q * pi_1);
	const __float128 T = 1.0q - trig_1 + trig_2 + trig_3 - trig_4;

	const __float128 c_sum = c1_prime + c2_prime;
	const __float128 c_product = c1_prime * c2_prime;
	const __float128 c_geo_mean = sqrtq(c_product);

	// 8. Compute hue difference and scaling factor S_H
	const __float128 sin_h_diff = sinq(h_diff);
	const __float128 S_H = 1.0q + 0.0075q * c_sum * T;
	const __float128 H_term = 2.0q * c_geo_mean * sin_h_diff / (kh * S_H);

	// 9. Compute chroma difference and scaling factor S_C
	const __float128 c_delta = c2_prime - c1_prime;
	const __float128 S_C = 1.0q + 0.0225q * c_sum;
	const __float128 C_term = c_delta / (kc * S_C);

	// 10. Combine lightness, chroma, hue, and interaction terms
	const __float128 L_part = L_term * L_term;
	const __float128 C_part = C_term * C_term;
	const __float128 H_part = H_term * H_term;
	const __float128 interaction = C_term * H_term * R_T;
	const __float128 delta_e_squared = L_part + C_part + H_part + interaction;
	const __float128 delta_e = sqrtq(delta_e_squared);

	return delta_e; // Given a tolerance of 3.6e-31
}

// If you remove the constant 1E-14, the code will continue to work, but CIEDE2000
// interoperability between all programming languages will no longer be guaranteed.
