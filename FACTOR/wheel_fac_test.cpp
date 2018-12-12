/**
 * wheel_fac_test.cpp - A wheel factorization routine and its test code.
 * WHEELFAC() is intended to match the non-user-interface part of the
 * FACTOR routine written for the Casio fx-5800P calculator, and WFSUB()
 * is intended to match the WFSUB routine.
 *
 * Written in 2018 by Ben Tesch.
 *
 * To the extent possible under law, the author has dedicated all copyright
 * and related and neighboring rights to this software to the public domain
 * worldwide. This software is distributed without any warranty.
 * The text of the CC0 Public Domain Dedication should be reproduced at the
 * end of this file. If not, see http://creativecommons.org/publicdomain/zero/1.0/
 */
#include <cinttypes>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <stdexcept>

/**
 * Calculates the integer square root of a 64 bit unsigned number.
 * From http://www.codecodex.com/wiki/Calculate_an_integer_square_root
 * with stylistic modifications and extension to 64 bits.
 */
uint64_t isqrt64(const uint64_t &number)
{
  uint64_t root = 0ull;
  uint64_t remainder = number;
  uint64_t place = 0x4000000000000000ull;

  while (place > remainder) place >>= 2u;

  while (place != 0ull)
  {
    if (remainder >= (root + place))
    {
      remainder -= root + place;
      root += place << 1u;
    }
    root >>= 1u;
    place >>= 2u;
  }

  return root;
}

/**
 * Test function.
 * Returns true if number is prime and false otherwise.
 */
bool isPrime(const uint64_t &number)
{
  if (number == 0ull) return false;
  if ((number == 1ull) | (number == 2ull)) return true;
  if ((number & 1ull) == 0ull) return false;

  uint64_t max_prime = isqrt64(number) + 1ull;
  for (uint64_t test = 3ull; test <= max_prime; test += 2ull)
  {
    if (number % test == 0ull) return false;
  }

  return true;
}

/**
 * Simulating the variables in the Casio fx-5800P with
 * globals. To maintain similarity, Z_[0] is not used.
 */
double A;
double B;
double C;
double D;
double E;
double F;
double G;
double H;
double I;
double J;
double K;
double L;
double M;
double N;
double O;
double P;
double Q;
double R;
double S;
double T;
double U;
double V;
double W;
double X;
double Y;
double Z;
double Z_[23];

/**
 * Implements the Casio fx-5800P Frac( function using
 * the C++ standard library.
 */
inline double Frac(const double &number)
{
  double scratch;
  return std::modf(number, &scratch);
}

/**
 * Subroutine of WHEELFAC()
 * Call with B a known factor of A.
 * Pulls all powers of B out of A and stores
 * B as a new factor as specified below:
 *
 * 1. E holds number of factors (ignoring powers)
 *    and must lie on the range [1,11].
 *
 *    a. E is limited to 11 because 
 *       1 * 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29
 *       is the largest set of unique factors whose
 *       product is less than or equal to 10000000000, the
 *       largest number in the range represented as exact
 *       integers in the fx-5800P.
 *
 * 2. Z_[E] holds B.
 * 
 * 3. Z_[E + 11] holds the highest power of B
 *    that is a factor of A.
 *
 * WFSUB modifies A, C, D, E, and Z_[1] through Z_[22].
 * It uses but does not modify B.
 */
void WFSUB(void)
{
  E = E + 1.0;
  if (E > 11.0) { std::printf("\nERROR: ran out of factor storage space\n\n"); A = 1.0; return; }
  Z_[static_cast<uint16_t>(E)] = B;

  do
  {
    A = D;
    Z_[static_cast<uint16_t>(E + 11.0)] = Z_[static_cast<uint16_t>(E + 11.0)] + 1.0;
    D = A / B;
  } while (Frac(D) == 0.0);
  
  C = std::trunc(std::sqrt(A));
  return;
}

/**
 * Wheel factorization with a basis of 2, 3, 5, and 7.
 * Stores factors of A in Z_[1] through Z_[22] per the
 * comment above WFSUB.
 * Modifies, A, B, C, D, E, and Z_[1] through Z_[22].
 * Calls WFSUB.
 */
void WHEELFAC(void)
{
  for (E = 1.0; E <= 22.0; E = E + 1.0) Z_[static_cast<uint16_t>(E)] = 0.0;
  Z_[1] = 1.0;
  Z_[12] = 1.0;
  E = 1.0;
  C = std::trunc(std::sqrt(A));
  
  B =  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
  B =  3.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
  B =  5.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
  B =  7.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
  B = 11.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
  
  while (true)
  {
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  8.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  8.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  6.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  4.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B + 10.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B +  2.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
    B = B + 10.0; D = A / B; if (Frac(D) == 0.0) WFSUB(); if (B > C) goto one;
  }
  
  one: 
  if (A > 1.0)
  {
    E = E + 1.0;
    Z_[static_cast<uint16_t>(E)] = A;
    A = 1.0;
    Z_[static_cast<uint16_t>(E + 11.0)] = 1.0;
  }
  
  return;
}

/**
 * Test the WHEELFAC() routine.
 * Take and validate user input to determine the range of numbers
 * used during testing.
 * The stride argument facilitates testing a large number range
 * by running multiple instances of this program in parallel.
 */
int main(int argc, char *argv[])
{
  std::printf("cmdline arguments: [start] [end] [stride]\n");
  std::printf("all arguments required. must be positive integers.\n\n");

  if (argc != 4)
  {
    std::printf("ERROR: must run with 3 arguments.\n\n");
    return 0;
  }

  int64_t istart;
  try { istart = std::stoll(argv[1]); }
  catch (const std::out_of_range &e)
  {
    std::printf("ERROR: start argument exceeds range of int64\n");
    std::printf("       pick start on range [1,9999999999]\n\n");
    return 0;
  }
  catch (const std::invalid_argument &e)
  {
    std::printf("ERROR: could not convert start argument to integer\n\n");
    return 0;
  }
  if ((istart < 1ll) | (istart > 9999999999ll))
  {
    std::printf("ERROR: start argument out of allowed range [1,9999999999]\n\n");
    return 0;
  }

  int64_t iend;
  try { iend = std::stoll(argv[2]); }
  catch (const std::out_of_range &e)
  {
    std::printf("ERROR: end argument exceeds range of int64\n");
    std::printf("       pick end on range [1,9999999999]\n\n");
    return 0;
  }
  catch (const std::invalid_argument &e)
  {
    std::printf("ERROR: could not convert end argument to integer\n\n");
    return 0;
  }
  if ((iend < 1ll) | (iend > 9999999999ll))
  {
    std::printf("ERROR: end argument out of allowed range [1,9999999999]\n\n");
    return 0;
  }
  if (iend < istart)
  {
    std::printf("ERROR: end argument %" PRIi64 " less than start argument %" PRIi64 "\n\n", iend, istart);
    return 0;
  }

  int64_t istride;
  try { istride = std::stoll(argv[3]); }
  catch (const std::out_of_range &e)
  {
    std::printf("ERROR: stride argument exceeds range of int64\n\n");
    return 0;
  }
  catch (const std::invalid_argument &e)
  {
    std::printf("ERROR: could not convert stride argument to integer\n\n");
    return 0;
  }
  if (istride < 1ll)
  {
    std::printf("ERROR: stride argument must be positive\n\n");
    return 0;
  }
  if (istride > (iend - istart))
  {
    std::printf("WARNING: stride = %" PRIi64 " exceeds (end - start) = %" PRIi64 "\n", istride, iend - istart);
    std::printf("         only start = %" PRIi64 " will be evaluated.\n\n", istart);
  }

  uint64_t total_evals = 1ull + (iend - istart) / istride;
  uint64_t n_evals = 0ull;
  uint64_t print_count = 0ull;
  for (int64_t i = istart; i <= iend; i += istride)
  {
    double i_dbl = static_cast<double>(i);
    A = i_dbl;
    WHEELFAC();
    
    double product = 1.0;
    for (uint32_t j = 1u; j <= 11u; j += 1u)
    {
      if (Z_[j] == 0.0) break;
      if (Z_[j] < 0.0) std::printf("\nERROR: negative factor %.0f for input %.0f\n\n", Z_[j], i_dbl);
      if (Z_[j] != std::trunc(Z_[j])) std::printf("\nERROR: non-integer factor %.16f for input %.0f\n\n", Z_[j], i_dbl);
      if (!isPrime(static_cast<uint64_t>(Z_[j]))) std::printf("\nERROR: composite factor %.0f for input %.0f\n\n", Z_[j], i_dbl);

      for (double k = Z_[j + 11u]; k > 0.0; k = k - 1.0) product = product * Z_[j];
    }
    
    if (product != i_dbl)
    {
      std::printf("\nERROR: incorrect factorization for input %.0f\n\n", i_dbl);
      for (uint32_t j = 1u; j <= 11u; j += 1u)
      {
        if (Z_[j] == 0.0) break;
        std::printf("%.0f^%.0f\n", Z_[j], Z_[j + 11u]);
      }
    }
    
    bool sorted = true;
    for (uint32_t j = 1u; j <= 10u; j += 1u)
    {
      if (Z_[j] == 0.0) break;
      if ((Z_[j + 1u] != 0.0) & (Z_[j + 1u] <= Z_[j]))
      {
        sorted = false;
        break;
      }
    }
    
    if (!sorted)
    {
      std::printf("\nERROR: factors not sorted for input %.0f\n\n", i_dbl);
      for (uint32_t j = 1u; j <= 11u; j += 1u)
      {
        if (Z_[j] == 0.0) break;
        std::printf("%.0f^%.0f\n", Z_[j], Z_[j + 11u]);
      }
    }

    n_evals++;
    print_count++;
    if (print_count >= 10000000ull)
    {
      print_count = 0ull;
      std::printf("progress: %" PRIu64 " / %" PRIu64 "\n", n_evals, total_evals);
    }
  }

  std::printf("\nIf no errors above, it worked.\n");
  std::printf("start = %" PRIi64 ", end = %" PRIi64 ", stride = %" PRIi64 "\n\n", istart, iend, istride);
  return 0;
}

/*
Creative Commons Legal Code

CC0 1.0 Universal

    CREATIVE COMMONS CORPORATION IS NOT A LAW FIRM AND DOES NOT PROVIDE
    LEGAL SERVICES. DISTRIBUTION OF THIS DOCUMENT DOES NOT CREATE AN
    ATTORNEY-CLIENT RELATIONSHIP. CREATIVE COMMONS PROVIDES THIS
    INFORMATION ON AN "AS-IS" BASIS. CREATIVE COMMONS MAKES NO WARRANTIES
    REGARDING THE USE OF THIS DOCUMENT OR THE INFORMATION OR WORKS
    PROVIDED HEREUNDER, AND DISCLAIMS LIABILITY FOR DAMAGES RESULTING FROM
    THE USE OF THIS DOCUMENT OR THE INFORMATION OR WORKS PROVIDED
    HEREUNDER.

Statement of Purpose

The laws of most jurisdictions throughout the world automatically confer
exclusive Copyright and Related Rights (defined below) upon the creator
and subsequent owner(s) (each and all, an "owner") of an original work of
authorship and/or a database (each, a "Work").

Certain owners wish to permanently relinquish those rights to a Work for
the purpose of contributing to a commons of creative, cultural and
scientific works ("Commons") that the public can reliably and without fear
of later claims of infringement build upon, modify, incorporate in other
works, reuse and redistribute as freely as possible in any form whatsoever
and for any purposes, including without limitation commercial purposes.
These owners may contribute to the Commons to promote the ideal of a free
culture and the further production of creative, cultural and scientific
works, or to gain reputation or greater distribution for their Work in
part through the use and efforts of others.

For these and/or other purposes and motivations, and without any
expectation of additional consideration or compensation, the person
associating CC0 with a Work (the "Affirmer"), to the extent that he or she
is an owner of Copyright and Related Rights in the Work, voluntarily
elects to apply CC0 to the Work and publicly distribute the Work under its
terms, with knowledge of his or her Copyright and Related Rights in the
Work and the meaning and intended legal effect of CC0 on those rights.

1. Copyright and Related Rights. A Work made available under CC0 may be
protected by copyright and related or neighboring rights ("Copyright and
Related Rights"). Copyright and Related Rights include, but are not
limited to, the following:

  i. the right to reproduce, adapt, distribute, perform, display,
     communicate, and translate a Work;
 ii. moral rights retained by the original author(s) and/or performer(s);
iii. publicity and privacy rights pertaining to a person's image or
     likeness depicted in a Work;
 iv. rights protecting against unfair competition in regards to a Work,
     subject to the limitations in paragraph 4(a), below;
  v. rights protecting the extraction, dissemination, use and reuse of data
     in a Work;
 vi. database rights (such as those arising under Directive 96/9/EC of the
     European Parliament and of the Council of 11 March 1996 on the legal
     protection of databases, and under any national implementation
     thereof, including any amended or successor version of such
     directive); and
vii. other similar, equivalent or corresponding rights throughout the
     world based on applicable law or treaty, and any national
     implementations thereof.

2. Waiver. To the greatest extent permitted by, but not in contravention
of, applicable law, Affirmer hereby overtly, fully, permanently,
irrevocably and unconditionally waives, abandons, and surrenders all of
Affirmer's Copyright and Related Rights and associated claims and causes
of action, whether now known or unknown (including existing as well as
future claims and causes of action), in the Work (i) in all territories
worldwide, (ii) for the maximum duration provided by applicable law or
treaty (including future time extensions), (iii) in any current or future
medium and for any number of copies, and (iv) for any purpose whatsoever,
including without limitation commercial, advertising or promotional
purposes (the "Waiver"). Affirmer makes the Waiver for the benefit of each
member of the public at large and to the detriment of Affirmer's heirs and
successors, fully intending that such Waiver shall not be subject to
revocation, rescission, cancellation, termination, or any other legal or
equitable action to disrupt the quiet enjoyment of the Work by the public
as contemplated by Affirmer's express Statement of Purpose.

3. Public License Fallback. Should any part of the Waiver for any reason
be judged legally invalid or ineffective under applicable law, then the
Waiver shall be preserved to the maximum extent permitted taking into
account Affirmer's express Statement of Purpose. In addition, to the
extent the Waiver is so judged Affirmer hereby grants to each affected
person a royalty-free, non transferable, non sublicensable, non exclusive,
irrevocable and unconditional license to exercise Affirmer's Copyright and
Related Rights in the Work (i) in all territories worldwide, (ii) for the
maximum duration provided by applicable law or treaty (including future
time extensions), (iii) in any current or future medium and for any number
of copies, and (iv) for any purpose whatsoever, including without
limitation commercial, advertising or promotional purposes (the
"License"). The License shall be deemed effective as of the date CC0 was
applied by Affirmer to the Work. Should any part of the License for any
reason be judged legally invalid or ineffective under applicable law, such
partial invalidity or ineffectiveness shall not invalidate the remainder
of the License, and in such case Affirmer hereby affirms that he or she
will not (i) exercise any of his or her remaining Copyright and Related
Rights in the Work or (ii) assert any associated claims and causes of
action with respect to the Work, in either case contrary to Affirmer's
express Statement of Purpose.

4. Limitations and Disclaimers.

 a. No trademark or patent rights held by Affirmer are waived, abandoned,
    surrendered, licensed or otherwise affected by this document.
 b. Affirmer offers the Work as-is and makes no representations or
    warranties of any kind concerning the Work, express, implied,
    statutory or otherwise, including without limitation warranties of
    title, merchantability, fitness for a particular purpose, non
    infringement, or the absence of latent or other defects, accuracy, or
    the present or absence of errors, whether or not discoverable, all to
    the greatest extent permissible under applicable law.
 c. Affirmer disclaims responsibility for clearing rights of other persons
    that may apply to the Work or any use thereof, including without
    limitation any person's Copyright and Related Rights in the Work.
    Further, Affirmer disclaims responsibility for obtaining any necessary
    consents, permissions or other rights required for any use of the
    Work.
 d. Affirmer understands and acknowledges that Creative Commons is not a
    party to this document and has no duty or obligation with respect to
    this CC0 or use of the Work.
*/