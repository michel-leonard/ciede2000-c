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

- GCC 13.3.0
- GMP 6.3.0
- MPFR 4.2.1

## Public Domain Licence

You are free to use these files, even for commercial purposes.
