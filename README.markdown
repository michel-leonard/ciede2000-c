# CIEDE2000 color difference formula in C

This page presents the CIEDE2000 color difference, implemented in the C99 programming language.

![Logo](https://raw.githubusercontent.com/michel-leonard/ciede2000-color-matching/refs/heads/main/docs/assets/images/logo.jpg)

## Our CIEDE2000 offer

These 4 production-ready files, released in 2026, contain the CIEDE2000 algorithm.

Source File|Type|Bits|Purpose|Advantage|
|:--:|:--:|:--:|:--:|:--:|
[ciede2000-32-bits.c](./ciede2000-32-bits.c)|`float`|32|General|Lightness, Speed|
[ciede2000-64-bits.c](./ciede2000-64-bits.c)|`double`|64|Scientific|Interoperability|
[ciede2000-128-bits.c](./ciede2000-128-bits.c)|`__float128`|128|Metrology|Higher Precision|
[ciede2000-arbitrary-precision.c](./ciede2000-arbitrary-precision.c)|`mpfr_t`|Unlimited|Metrology|–|

### Software Versions

- C99, C11, C17, C23
- GCC 13.3
- Clang 16
- libquadmath 14.2 (only for 128-bit precision)
- GMP 6.3 (only for arbitrary precision)
- MPFR 4.2.1 (only for arbitrary precision)

### Example Usage

We calculate the CIEDE2000 distance between two colors, first without and then with parametric factors.

```c
// You must name this file "demo.c" and include the contents of "ciede2000-64-bits.c" here

#include <stdio.h>

int main(void) {
	// Example of two L*a*b* colors
	const double l1 = 88.3, a1 = 126.1, b1 = -1.7;
	const double l2 = 89.3, a2 = 109.1, b2 = 4.6;

	double delta_e = ciede2000(l1, a1, b1, l2, a2, b2);
	printf("CIEDE2000 = %.14f\n", delta_e); // ΔE2000 = 3.39372310844506

	// Example of parametric factors used in the textile industry
	const double kl = 2.0, kc = 1.0, kh = 1.0;

	// Perform a CIEDE2000 calculation compliant with that of Gaurav Sharma
	const int canonical = 1;

	delta_e = ciede2000_with_parameters(l1, a1, b1, l2, a2, b2, kl, kc, kh, canonical);
	printf("CIEDE2000 = %.14f\n", delta_e); // ΔE2000 = 3.34905104605513
}
```

- Compile with `clang` or `gcc -std=c99 -Wall -Wextra -Wpedantic -Ofast -o ciede2000-demo demo.c -lm`
- Execute `./ciede2000-demo` to display the calculated color differences

## Public Domain Licence

You are free to use these files, even for commercial purposes.
