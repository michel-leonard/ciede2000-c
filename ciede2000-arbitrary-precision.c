#include <mpfr.h>

// This function written in C99 is not affiliated with the CIE (International Commission on Illumination),
// and is released into the public domain. It is provided "as is" without any warranty, express or implied.

struct delta_e {
	mpfr_t sys[80];
	struct {
		mpfr_t l;
		mpfr_t a;
		mpfr_t b;
	} color[2];
	mpfr_t kl;
	mpfr_t kc;
	mpfr_t kh;
	mpfr_t result;
	int canonical;
};

// The arbitrary precision CIE ΔE2000 implementation, which operates on two L*a*b* colors, and returns their difference
static void ciede2000_arbitrary_precision(struct delta_e *manager) {
	// Developed in year 2026 with :
	// - GMP version: 6.3.0
	// - MPFR version: 4.2.1

	static const mpfr_rnd_t rnd = MPFR_RNDN; // round to nearest (precision_required = 35.5 + 3.35 * n_correct_decimals)

	// You can copy this line to extract intermediate values after the calculation
	mpfr_t *restrict c = manager->sys, *my_zero = c, *my_25_pow_7 = ++c, *my_pi = ++c, *my_2_pi = ++c,
			*my_interoperability_pi = ++c, *my_55_pi = ++c, *my_minus_25_pi_pi = ++c, *my_pi_div_3 = ++c, *my_0015 = ++c,
			*my_017 = ++c, *my_024 = ++c, *my_pi_div_2 = ++c, *my_3_pi_div_2 = ++c, *my_032 = ++c, *my_8_pi_div_15 = ++c,
			*my_02 = ++c, *my_3_pi_div_20 = ++c, *my_00075 = ++c, *my_00225 = ++c, *X = ++c, *Y = ++c, *a1_sq = ++c,
			*b1_sq = ++c, *c1_orig = ++c, *a2_sq = ++c, *b2_sq = ++c, *c2_orig = ++c, *c_avg = ++c, *c_avg_7 = ++c,
			*g_denom = ++c, *g_ratio = ++c, *g_sqrt = ++c, *g_factor = ++c, *a1_prime = ++c, *c1_prime_sq = ++c,
			*c1_prime = ++c, *a2_prime = ++c, *c2_prime_sq = ++c, *c2_prime = ++c, *h1_adj = ++c, *h2_adj = ++c,
			*delta_h = ++c, *h_mean_raw = ++c, *h_diff_raw = ++c, *h_diff = ++c, *h_mean = ++c, *c_sum = ++c, *c_bar = ++c,
			*c_bar_7 = ++c, *rc_denom = ++c, *R_C = ++c, *theta = ++c, *exp_term = ++c, *sin_term = ++c, *R_T = ++c,
			*l_avg = ++c, *l_delta = ++c, *l_delta_sq = ++c, *sl_num = ++c, *sl_denom = ++c, *S_L = ++c, *L_term = ++c,
			*trig_1 = ++c, *trig_2 = ++c, *trig_3 = ++c, *trig_4 = ++c, *T = ++c, *c_product = ++c, *c_geo_mean = ++c,
			*sin_h_diff = ++c, *S_H = ++c, *H_term = ++c, *c_delta = ++c, *S_C = ++c, *C_term = ++c, *L_part = ++c,
			*C_part = ++c, *H_part = ++c, *interaction = ++c, *delta_e_squared = ++c;

	// This algorithm does not perform degree/radian conversions

	if (mpfr_zero_p(*my_zero)) {
		// We define the constants once and for all
		mpfr_set_ui(*my_zero, 0, rnd);
		mpfr_set_ui(*my_25_pow_7, 6103515625, rnd);
		mpfr_const_pi(*my_pi, rnd);
		mpfr_mul_2ui(*my_2_pi, *my_pi, 1, rnd);
		// This approach guarantees a stable algorithm, avoiding divergences around π
		mpfr_add_d(*my_interoperability_pi, *my_pi, 1E-14, rnd);
		mpfr_mul_ui(*my_55_pi, *my_pi, 55, rnd);
		mpfr_sqr(*my_pi_div_3, *my_pi, rnd); //tmp
		mpfr_mul_si(*my_minus_25_pi_pi, *my_pi_div_3, -25, rnd);
		mpfr_div_ui(*my_pi_div_3, *my_pi, 3, rnd);
		mpfr_set_str(*my_0015, "0.015", 10, rnd);
		mpfr_set_str(*my_017, "0.17", 10, rnd);
		mpfr_set_str(*my_024, "0.24", 10, rnd);
		mpfr_div_2ui(*my_pi_div_2, *my_pi, 1, rnd);
		mpfr_mul_ui(*my_3_pi_div_2, *my_pi_div_2, 3, rnd);
		mpfr_set_str(*my_032, "0.32", 10, rnd);
		mpfr_mul_2ui(*my_02, *my_pi, 3, rnd); // tmp
		mpfr_div_ui(*my_8_pi_div_15, *my_02, 15, rnd);
		mpfr_set_str(*my_02, "0.2", 10, rnd);
		mpfr_mul_ui(*my_00075, *my_pi, 3, rnd); // tmp
		mpfr_div_ui(*my_3_pi_div_20, *my_00075, 20, rnd);
		mpfr_set_str(*my_00075, "0.0075", 10, rnd);
		mpfr_set_str(*my_00225, "0.0225", 10, rnd);
	}

	// 1. Compute chroma magnitudes ... a and b usually range from -128 to +127
	mpfr_sqr(*a1_sq, manager->color[0].a, rnd);
	mpfr_sqr(*b1_sq, manager->color[0].b, rnd);
	mpfr_add(*X, *a1_sq, *b1_sq, rnd);
	mpfr_sqrt(*c1_orig, *X, rnd);
	//
	mpfr_sqr(*a2_sq, manager->color[1].a, rnd);
	mpfr_sqr(*b2_sq, manager->color[1].b, rnd);
	mpfr_add(*X, *a2_sq, *b2_sq, rnd);
	mpfr_sqrt(*c2_orig, *X, rnd);

	// 2. Compute chroma mean and apply G compensation
	mpfr_add(*X, *c1_orig, *c2_orig, rnd);
	mpfr_div_2ui(*c_avg, *X, 1, rnd);
	mpfr_sqr(*X, *c_avg, rnd);
	mpfr_mul(*Y, *X, *c_avg, rnd);
	mpfr_sqr(*X, *Y, rnd);
	mpfr_mul(*c_avg_7, *X, *c_avg, rnd);
	mpfr_add(*g_denom, *c_avg_7, *my_25_pow_7, rnd);
	mpfr_div(*g_ratio, *c_avg_7, *g_denom, rnd);
	mpfr_sqrt(*g_sqrt, *g_ratio, rnd);
	mpfr_si_sub(*X, 1, *g_sqrt, rnd);
	mpfr_div_2ui(*Y, *X, 1, rnd);
	mpfr_add_ui(*g_factor, *Y, 1, rnd);

	// 3. Apply G correction to a components, compute corrected chroma
	mpfr_mul(*a1_prime, manager->color[0].a, *g_factor, rnd);
	mpfr_sqr(*X, *a1_prime, rnd);
	mpfr_add(*c1_prime_sq, *X, *b1_sq, rnd);
	mpfr_sqrt(*c1_prime, *c1_prime_sq, rnd);
	mpfr_mul(*a2_prime, manager->color[1].a, *g_factor, rnd);
	mpfr_sqr(*X, *a2_prime, rnd);
	mpfr_add(*c2_prime_sq, *X, *b2_sq, rnd);
	mpfr_sqrt(*c2_prime, *c2_prime_sq, rnd);
	// The use of atan (not atan2) is clearer across programming languages
	int sa = mpfr_sgn(manager->color[0].a), sb = mpfr_sgn(manager->color[0].b);
	if (sa != 0) {
		mpfr_div(*X, manager->color[0].b, *a1_prime, rnd);
		mpfr_atan(*Y, *X, rnd);
		mpfr_add(*h1_adj, *Y, 0 < sa ? sb < 0 ? *my_2_pi : *my_zero : *my_pi, rnd);
	} else
		mpfr_set(*h1_adj, sb < 0 ? *my_3_pi_div_2 : 0 < sb ? *my_pi_div_2 : *my_pi, rnd);
	//
	sa = mpfr_sgn(manager->color[1].a), sb = mpfr_sgn(manager->color[1].b);
	if (sa != 0) {
		mpfr_div(*X, manager->color[1].b, *a2_prime, rnd);
		mpfr_atan(*Y, *X, rnd);
		mpfr_add(*h2_adj, *Y, 0 < sa ? sb < 0 ? *my_2_pi : *my_zero : *my_pi, rnd);
	} else
		mpfr_set(*h2_adj, sb < 0 ? *my_3_pi_div_2 : 0 < sb ? *my_pi_div_2 : *my_pi, rnd);
	// The atan2 polyfill (customized for ΔE2000) is complete

	// 4. Check if hue mean wraps around pi (180 deg)
	mpfr_sub(*X, *h2_adj, *h1_adj, rnd);
	mpfr_add(*Y, *h1_adj, *h2_adj, rnd);
	mpfr_div_2ui(*h_diff_raw, *X, 1, rnd);
	mpfr_div_2ui(*h_mean_raw, *Y, 1, rnd);
	mpfr_abs(*delta_h, *X, rnd);
	//  Michel Leonard 2026 - The part where most programmers get it wrong
	if (mpfr_cmp(*my_interoperability_pi, *delta_h) < 0) {
		// When mean wraps, difference wraps too
		mpfr_add(*h_diff, *h_diff_raw, *my_pi, rnd);
		if (manager->canonical && mpfr_cmp(*my_interoperability_pi, *h_mean_raw) < 0)
			// Sharma’s implementation, OpenJDK, ...
			mpfr_sub(*h_mean, *h_mean_raw, *my_pi, rnd);
		else
			// Lindbloom’s implementation, Netflix’s VMAF, ...
			mpfr_add(*h_mean, *h_mean_raw, *my_pi, rnd);
	} else {
		mpfr_set(*h_mean, *h_mean_raw, rnd);
		mpfr_set(*h_diff, *h_diff_raw, rnd);
	}

	// 5. Compute hue rotation correction factor R_T
	mpfr_add(*c_sum, *c1_prime, *c2_prime, rnd);
	mpfr_div_2ui(*c_bar, *c_sum, 1, rnd);
	mpfr_sqr(*X, *c_bar, rnd);
	mpfr_mul(*Y, *X, *c_bar, rnd);
	mpfr_sqr(*X, *Y, rnd);
	mpfr_mul(*c_bar_7, *X, *c_bar, rnd);
	mpfr_add(*rc_denom, *c_bar_7, *my_25_pow_7, rnd);
	mpfr_div(*X, *c_bar_7, *rc_denom, rnd);
	mpfr_sqrt(*R_C, *X, rnd);
	//
	mpfr_mul_ui(*X, *h_mean, 36, rnd);
	mpfr_sub(*theta, *X, *my_55_pi, rnd);
	mpfr_sqr(*X, *theta, rnd);
	mpfr_div(*Y, *X, *my_minus_25_pi_pi, rnd);
	mpfr_exp(*exp_term, *Y, rnd);
	mpfr_mul(*X, *exp_term, *my_pi_div_3, rnd);
	mpfr_sin(*sin_term, *X, rnd);
	// Rotation factor ... cross-effect between chroma and hue
	mpfr_mul(*X, *R_C, *sin_term, rnd);
	mpfr_mul_si(*R_T, *X, -2, rnd);

	// 6. Compute lightness term ... L nominally ranges from 0 to 100
	mpfr_add(*X, manager->color[0].l, manager->color[1].l, rnd);
	mpfr_div_2ui(*l_avg, *X, 1, rnd);
	mpfr_sub_ui(*X, *l_avg, 50, rnd);
	mpfr_sqr(*l_delta_sq, *X, rnd);
	mpfr_sub(*l_delta, manager->color[1].l, manager->color[0].l, rnd);
	// Adaptation to the non-linearity of light perception ... S_L
	mpfr_mul(*sl_num, *my_0015, *l_delta_sq, rnd);
	mpfr_add_ui(*X, *l_delta_sq, 20, rnd);
	mpfr_sqrt(*sl_denom, *X, rnd);
	mpfr_div(*X, *sl_num, *sl_denom, rnd);
	mpfr_add_ui(*S_L, *X, 1, rnd);
	mpfr_mul(*X, *S_L, manager->kl, rnd);
	mpfr_div(*L_term, *l_delta, *X, rnd);

	// 7. Compute chroma-related trig terms and factor T
	mpfr_set(*X, *h_mean, rnd);
	mpfr_add(*trig_1, *X, *my_pi_div_3, rnd);
	mpfr_sin(*Y, *trig_1, rnd);
	mpfr_mul(*trig_1, *my_017, *Y, rnd);
	// Stage 2
	mpfr_add(*Y, *X, *h_mean, rnd);
	mpfr_add(*trig_2, *Y, *my_pi_div_2, rnd);
	mpfr_sin(*X, *trig_2, rnd);
	mpfr_mul(*trig_2, *my_024, *X, rnd);
	// Stage 3
	mpfr_add(*X, *Y, *h_mean, rnd);
	mpfr_add(*trig_3, *X, *my_8_pi_div_15, rnd);
	mpfr_sin(*Y, *trig_3, rnd);
	mpfr_mul(*trig_3, *my_032, *Y, rnd);
	// Stage 4
	mpfr_add(*Y, *X, *h_mean, rnd);
	mpfr_add(*trig_4, *Y, *my_3_pi_div_20, rnd);
	mpfr_sin(*X, *trig_4, rnd);
	mpfr_mul(*trig_4, *my_02, *X, rnd);
	// Get T using what has been calculated in each stage
	mpfr_ui_sub(*X, 1, *trig_1, rnd);
	mpfr_add(*Y, *X, *trig_2, rnd);
	mpfr_add(*X, *Y, *trig_3, rnd);
	mpfr_sub(*T, *X, *trig_4, rnd);
	//
	mpfr_mul(*c_product, *c1_prime, *c2_prime, rnd);
	mpfr_sqrt(*c_geo_mean, *c_product, rnd);

	// 8. Compute hue difference and scaling factor S_H
	mpfr_sin(*sin_h_diff, *h_diff, rnd);
	mpfr_mul(*X, *c_sum, *T, rnd);
	mpfr_mul(*Y, *X, *my_00075, rnd);
	mpfr_add_ui(*S_H, *Y, 1, rnd);
	mpfr_mul(*X, manager->kh, *S_H, rnd);
	mpfr_div(*Y, *sin_h_diff, *X, rnd);
	mpfr_mul(*X, *Y, *c_geo_mean, rnd);
	mpfr_mul_2ui(*H_term, *X, 1, rnd);

	// 9. Compute chroma difference and scaling factor S_C
	mpfr_sub(*c_delta, *c2_prime, *c1_prime, rnd);
	mpfr_mul(*X, *c_sum, *my_00225, rnd);
	mpfr_add_ui(*S_C, *X, 1, rnd);
	mpfr_mul(*X, *S_C, manager->kc, rnd);
	mpfr_div(*C_term, *c_delta, *X, rnd);

	// 10. Combine lightness, chroma, hue, and interaction terms
	mpfr_sqr(*L_part, *L_term, rnd);
	mpfr_sqr(*C_part, *C_term, rnd);
	mpfr_sqr(*H_part, *H_term, rnd);
	mpfr_mul(*X, *H_term, *R_T, rnd);
	mpfr_mul(*interaction, *C_term, *X, rnd);
	mpfr_add(*X, *L_part, *C_part, rnd);
	mpfr_add(*Y, *X, *H_part, rnd);
	mpfr_add(*delta_e_squared, *Y, *interaction, rnd);

	// Store the CIEDE2000 metric as a result
	mpfr_sqrt(manager->result, *delta_e_squared, rnd);
}

// If you remove the constant 1E-14, the code will continue to work, but CIEDE2000
// interoperability between all programming languages will no longer be guaranteed.

// Source code tested by Michel LEONARD
// Website: ciede2000.pages-perso.free.fr

////////////////////////////////////////////////////////////////////////////////////
////////////                                                            ////////////
////////////            ALLOCATES/CLEARS THE MEMORY FOR THE             ////////////
////////////           CIEDE2000 ARBITRARY PRECISION NUMBERS            ////////////
////////////                                                            ////////////
////////////   precision_required = 35.5 + 3.35 * n_correct_decimals    ////////////
////////////                                                            ////////////
////////////////////////////////////////////////////////////////////////////////////

static void ciede2000_arbitrary_precision_memory(struct delta_e *c, mpfr_prec_t precision) {
	if (precision) {
		// Pass a non-zero precision to perform a boot
		for (size_t i = 0; i < sizeof(c->sys) / sizeof(*c->sys); ++i) mpfr_init2(c->sys[i], precision);
		mpfr_inits2(precision, c->color[0].l, c->color[0].a, c->color[0].b,
					c->color[1].l, c->color[1].a, c->color[1].b,
					c->kl, c->kc, c->kh, c->result, NULL);
		mpfr_set_ui(*c->sys, 0, MPFR_RNDN);
	} else {
		// Pass a zero precision to release the memory
		for (size_t i = 0; i < sizeof(c->sys) / sizeof(*c->sys); ++i) mpfr_clear(c->sys[i]);
		mpfr_clears(c->color[0].l, c->color[0].a, c->color[0].b,
					c->color[1].l, c->color[1].a, c->color[1].b,
					c->kl, c->kc, c->kh, c->result, NULL);
		(*c->sys)->_mpfr_d = 0;
	}
	mpfr_free_cache();
}
