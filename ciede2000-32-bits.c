#include <math.h>

// This function written in C99 is not affiliated with the CIE (International Commission on Illumination),
// and is released into the public domain. It is provided "as is" without any warranty, express or implied.

// The functional CIE ΔE2000 implementation, which operates on two L*a*b* colors, and returns their difference.
// "l" ranges from 0 to 100, while "a" and "b" are unbounded and commonly clamped to the range of -128 to 127.
static float ciede2000_functional(float l1, float a1, float b1, float l2, float a2, float b2, float kl, float kc, float kh, int canonical) {
	// Working in C99 with the CIEDE2000 color-difference formula.
	// kl, kc, kh are parametric factors to be adjusted according to
	// different viewing parameters such as textures, backgrounds...

	const float pi_1 = 3.14159265358979323846f;
	const float pi_3 = 1.04719755119659774615f;

	// 1. Compute chroma magnitudes ... a and b usually range from -128 to +127
	const float a1_sq = a1 * a1;
	const float b1_sq = b1 * b1;
	const float c1_orig = sqrtf(a1_sq + b1_sq);

	const float a2_sq = a2 * a2;
	const float b2_sq = b2 * b2;
	const float c2_orig = sqrtf(a2_sq + b2_sq);

	// 2. Compute chroma mean and apply G compensation
	const float c_avg = 0.5f * (c1_orig + c2_orig);
	const float c_avg_3 = c_avg * c_avg * c_avg;
	const float c_avg_7 = c_avg_3 * c_avg_3 * c_avg;
	const float g_denom = c_avg_7 + 6103515625.0f;
	const float g_ratio = c_avg_7 / g_denom;
	const float g_sqrt = sqrtf(g_ratio);
	const float g_factor = 1.0f + 0.5f * (1.0f - g_sqrt);

	// 3. Apply G correction to a components, compute corrected chroma
	const float a1_prime = a1 * g_factor;
	const float c1_prime_sq = a1_prime * a1_prime + b1 * b1;
	const float c1_prime = sqrtf(c1_prime_sq);
	const float a2_prime = a2 * g_factor;
	const float c2_prime_sq = a2_prime * a2_prime + b2 * b2;
	const float c2_prime = sqrtf(c2_prime_sq);

	// 4. Compute hue angles in radians, adjust for negatives and wrap
	const float h1_raw = atan2f(b1, a1_prime);
	const float h2_raw = atan2f(b2, a2_prime);
	const float h1_adj = h1_raw + (float) (h1_raw < 0.0f) * 2.0f * pi_1;
	const float h2_adj = h2_raw + (float) (h2_raw < 0.0f) * 2.0f * pi_1;
	const float delta_h = fabsf(h1_adj - h2_adj);
	const float h_mean_raw = 0.5f * (h1_adj + h2_adj);
	const float h_diff_raw = 0.5f * (h2_adj - h1_adj);

	// Check if hue mean wraps around pi (180 deg)
	const float hue_wrap = (float) (pi_1 + 1E-6f < delta_h);

	// The part where most programmers get it wrong
	float h_mean;
	if (canonical && pi_1 + 1E-6f < h_mean_raw)
		// Gaurav Sharma, OpenJDK, ...
		h_mean = h_mean_raw - hue_wrap * pi_1;
	else
		// Bruce Lindbloom, Netflix’s VMAF, ...
		h_mean = h_mean_raw + hue_wrap * pi_1;

	// Michel Leonard 2026 - When mean wraps, difference wraps too
	const float h_diff = h_diff_raw + hue_wrap * pi_1;

	// 5. Compute hue rotation correction factor R_T
	const float c_bar = 0.5f * (c1_prime + c2_prime);
	const float c_bar_3 = c_bar * c_bar * c_bar;
	const float c_bar_7 = c_bar_3 * c_bar_3 * c_bar;
	const float rc_denom = c_bar_7 + 6103515625.0f;
	const float R_C = sqrtf(c_bar_7 / rc_denom);

	const float theta = 36.0f * h_mean - 55.0f * pi_1;
	const float theta_denom = -25.0f * pi_1 * pi_1;
	const float exp_argument = theta * theta / theta_denom;
	const float exp_term = expf(exp_argument);
	const float delta_theta = pi_3 * exp_term;
	const float sin_term = sinf(delta_theta);

	// Rotation factor ... cross-effect between chroma and hue
	const float R_T = -2.0f * R_C * sin_term;

	// 6. Compute lightness term ... L nominally ranges from 0 to 100
	const float l_avg = 0.5f * (l1 + l2);
	const float l_delta_sq = (l_avg - 50.0f) * (l_avg - 50.0f);
	const float l_delta = l2 - l1;

	// Adaptation to the non-linearity of light perception ... S_L
	const float sl_num = 0.015f * l_delta_sq;
	const float sl_denom = sqrtf(20.0f + l_delta_sq);
	const float S_L = 1.0f + sl_num / sl_denom;
	const float L_term = l_delta / (kl * S_L);

	// 7. Compute chroma-related trig terms and factor T
	const float trig_1 = 0.17f * sinf(h_mean + pi_3);
	const float trig_2 = 0.24f * sinf(2.0f * h_mean + 0.5f * pi_1);
	const float trig_3 = 0.32f * sinf(3.0f * h_mean + 1.6f * pi_3);
	const float trig_4 = 0.2f * sinf(4.0f * h_mean + 0.15f * pi_1);
	const float T = 1.0f - trig_1 + trig_2 + trig_3 - trig_4;

	const float c_sum = c1_prime + c2_prime;
	const float c_product = c1_prime * c2_prime;
	const float c_geo_mean = sqrtf(c_product);

	// 8. Compute hue difference and scaling factor S_H
	const float sin_h_diff = sinf(h_diff);
	const float S_H = 1.0f + 0.0075f * c_sum * T;
	const float H_term = 2.0f * c_geo_mean * sin_h_diff / (kh * S_H);

	// 9. Compute chroma difference and scaling factor S_C
	const float c_delta = c2_prime - c1_prime;
	const float S_C = 1.0f + 0.0225f * c_sum;
	const float C_term = c_delta / (kc * S_C);

	// 10. Combine lightness, chroma, hue, and interaction terms
	const float L_part = L_term * L_term;
	const float C_part = C_term * C_term;
	const float H_part = H_term * H_term;
	const float interaction = C_term * H_term * R_T;
	const float delta_e_squared = L_part + C_part + H_part + interaction;
	const float delta_e = sqrtf(delta_e_squared);

	return delta_e; // Given a tolerance of 1.4e-4
}

// If you remove the constant 1E-6, the code will continue to work, but CIEDE2000
// interoperability between all programming languages will no longer be guaranteed.
