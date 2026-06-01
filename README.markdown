Browse : [Wolfram Language](https://github.com/michel-leonard/ciede2000-wolfram-language) · [AWK](https://github.com/michel-leonard/ciede2000-awk) · [BC](https://github.com/michel-leonard/ciede2000-basic-calculator) · [C#](https://github.com/michel-leonard/ciede2000-csharp) · [C++](https://github.com/michel-leonard/ciede2000-cpp) · **C99** · [Dart](https://github.com/michel-leonard/ciede2000-dart) · [Go](https://github.com/michel-leonard/ciede2000-go) · [JavaScript](https://github.com/michel-leonard/ciede2000-javascript) · [Java](https://github.com/michel-leonard/ciede2000-java) · [Julia](https://github.com/michel-leonard/ciede2000-julia)

# CIEDE2000 color difference formula in C

This page presents the CIEDE2000 color difference, implemented in the C99 programming language.

![Logo](https://raw.githubusercontent.com/michel-leonard/ciede2000-color-matching/refs/heads/main/docs/assets/images/logo.jpg)

## About

Here you’ll find the first rigorously correct implementation of CIEDE2000 that doesn’t use any conversion between degrees and radians. Set parameter `canonical` to obtain results in line with your existing pipeline.

`canonical`|The algorithm operates...|
|:--:|-|
`0`|in accordance with the CIEDE2000 values currently used by many industry players|
`1`|in accordance with the CIEDE2000 values provided by [this](https://hajim.rochester.edu/ece/sites/gsharma/ciede2000/) academic MATLAB function|

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

### Test Results

LEONARD’s tests are based on well-chosen L\*a\*b\* colors, with various parametric factors `kL`, `kC` and `kH`.

<details>
<summary>Display test results for 3 correct decimal places in 32-bits</summary>

```
CIEDE2000 Verification Summary :
          Compliance : [ ] CANONICAL [X] SIMPLIFIED
  First Checked Line : 40.0,0.5,-128.0,49.91,0.0,24.0,1.0,1.0,1.0,51.0186691
           Precision : 3 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 29.61 seconds
     Average Delta E : 63.58
   Average Deviation : 7.5e-06
   Maximum Deviation : 0.00015
```

```
CIEDE2000 Verification Summary :
          Compliance : [X] CANONICAL [ ] SIMPLIFIED
  First Checked Line : 40.0,0.5,-128.0,49.91,0.0,24.0,1.0,1.0,1.0,51.0184669
           Precision : 3 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 29.57 seconds
     Average Delta E : 63.58
   Average Deviation : 7.1e-06
   Maximum Deviation : 0.00019
```

</details>

<details>
<summary>Display test results for 12 correct decimal places in 64-bits</summary>

```
CIEDE2000 Verification Summary :
          Compliance : [ ] CANONICAL [X] SIMPLIFIED
  First Checked Line : -0.0,8.0,32.0,0.0,32.00003,-128.0,1.0,1.0,1.0,52.957588884756206
           Precision : 12 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 31.40 seconds
     Average Delta E : 67.12
   Average Deviation : 1.5e-14
   Maximum Deviation : 3.4e-13
```

```
CIEDE2000 Verification Summary :
          Compliance : [X] CANONICAL [ ] SIMPLIFIED
  First Checked Line : -0.0,8.0,32.0,0.0,32.00003,-128.0,1.0,1.0,1.0,52.957396265103696
           Precision : 12 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 31.25 seconds
     Average Delta E : 67.12
   Average Deviation : 1.4e-14
   Maximum Deviation : 3.4e-13
```

</details>

<details>
<summary>Display test results for 30 correct decimal places in 128-bits</summary>

```
CIEDE2000 Verification Summary :
          Compliance : [ ] CANONICAL [X] SIMPLIFIED
  First Checked Line : 40.0,0.0,127.9999992,40.0,0.0005,-7.0,1.0,1.0,1.0,39.4274155669737183315...
           Precision : 30 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 247.17 seconds
     Average Delta E : 67.15
   Average Deviation : 1.2e-32
   Maximum Deviation : 2.5e-31
```

```
CIEDE2000 Verification Summary :
          Compliance : [X] CANONICAL [ ] SIMPLIFIED
  First Checked Line : 40.0,0.0,127.9999992,40.0,0.0005,-7.0,1.0,1.0,1.0,39.4276099196716983684...
           Precision : 30 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 239.35 seconds
     Average Delta E : 67.15
   Average Deviation : 1.1e-32
   Maximum Deviation : 2.5e-31
```

</details>

<details>
<summary>Display test results for 250 correct decimal places in arbitrary precision</summary>

```
CIEDE2000 Verification Summary :
          Compliance : [X] CANONICAL [ ] SIMPLIFIED
  First Checked Line : 40.0,0.5,-128.0,49.91,0.0,24.0,1.0,1.0,1.0,51.01846301969812634838534673...
           Precision : 250 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 2162.72 seconds
     Average Delta E : 67.13
   Average Deviation : 2.5e-252
   Maximum Deviation : 5e-252
```

```
CIEDE2000 Verification Summary :
          Compliance : [ ] CANONICAL [X] SIMPLIFIED
  First Checked Line : 40.0,0.5,-128.0,49.91,0.0,24.0,1.0,1.0,1.0,51.01866090771253205230602031...
           Precision : 250 decimal digits
           Successes : 10000000
               Error : 0
            Duration : 2176.73 seconds
     Average Delta E : 67.13
   Average Deviation : 2.5e-252
   Maximum Deviation : 5e-252
```

</details>

## Public Domain Licence

You are free to use these files, even for commercial purposes.
