#include "Candia-v2/SplittingFn.hpp"
#include "Candia-v2/SpecialFuncs.hpp"
#include <limits>
using namespace Candia2;

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <iomanip>

namespace orig{
namespace P3GSOFT {
    double A4gluon = 0.0;  // Mimics the Fortran COMMON block
}

// -------------------------------------------------------------
// P3GGA
// -------------------------------------------------------------
double P3GGA(double Y, int nf, int IMOD)
{
    using namespace P3GSOFT;

    double YM  = 1.0 / Y;
    double Y1  = 1.0 - Y;
    double DL  = std::log(Y);
    double DL1 = std::log1p(-Y);

    double nf2 = nf * nf;
    double nf3 = nf * nf2;

    // Large-x coefficients
    A4gluon =  40880.330
             - 11714.246  * nf
             +   440.04876 * nf2
             +     7.3627750 * nf3;

    double Ccoeff = 8.5814120e4 - 1.3880515e4 * nf + 1.3511111e2 * nf2;
    double Dcoeff = 5.4482808e4 - 4.3411337e3 * nf - 2.1333333e1 * nf2;

    // 1/x * ln^a(x) terms
    double bfkl0 = -8.308617314e3;
    double bfkl1 = -1.069119905e5 - 9.963830436e2 * nf;

    // Base part
    double P3gg01 =
          bfkl0 * std::pow(DL, 3) * YM
        + bfkl1 * std::pow(DL, 2) * YM
        + Ccoeff * DL1
        + Dcoeff - A4gluon;

    double P3ggApp1 = 0.0;
    double P3ggApp2 = 0.0;

    if (nf == 3) {
        P3ggApp1 = P3gg01
            + 3.4 * bfkl1 * DL * YM
            - 345063.0 * Y1 * YM
            + 86650.0 * (1.0 + Y * Y) * Y1
            + 158160.0 * DL
            - 15741.0 * Y1 * std::pow(DL1, 2)
            - 9417.0 * Y1 * std::pow(DL1, 3);

        P3ggApp2 = P3gg01
            + 5.4 * bfkl1 * DL * YM
            - 1265632.0 * Y1 * YM
            - 656644.0 * (1.0 + Y * Y) * Y1
            - 1352233.0 * DL
            + 203298.0 * Y1 * std::pow(DL1, 2)
            + 39112.0 * Y1 * std::pow(DL1, 3);
    }
    else if (nf == 4) {
        P3ggApp1 = P3gg01
            + 3.4 * bfkl1 * DL * YM
            - 342625.0 * Y1 * YM
            + 100372.0 * (1.0 + Y * Y) * Y1
            + 189167.0 * DL
            - 29762.0 * Y1 * std::pow(DL1, 2)
            - 12102.0 * Y1 * std::pow(DL1, 3);

        P3ggApp2 = P3gg01
            + 5.4 * bfkl1 * DL * YM
            - 1271540.0 * Y1 * YM
            - 649661.0 * (1.0 + Y * Y) * Y1
            - 1334919.0 * DL
            + 191263.0 * Y1 * std::pow(DL1, 2)
            + 36867.0 * Y1 * std::pow(DL1, 3);
    }
    else if (nf == 5) {
        P3ggApp1 = P3gg01
            + 3.4 * bfkl1 * DL * YM
            - 337540.0 * Y1 * YM
            + 119366.0 * (1.0 + Y * Y) * Y1
            + 223769.0 * DL
            - 45129.0 * Y1 * std::pow(DL1, 2)
            - 15046.0 * Y1 * std::pow(DL1, 3);

        P3ggApp2 = P3gg01
            + 5.4 * bfkl1 * DL * YM
            - 1274800.0 * Y1 * YM
            - 637406.0 * (1.0 + Y * Y) * Y1
            - 1314010.0 * DL
            + 177882.0 * Y1 * std::pow(DL1, 2)
            + 34362.0 * Y1 * std::pow(DL1, 3);
    }
    else {
        throw std::invalid_argument("Error in P3GGA: invalid nf (must be 3,4,5)");
    }

	double res{};
    if (IMOD == 1)
        res = P3ggApp1;
    else if (IMOD == 2)
        res = P3ggApp2;
    else
        res = 0.5 * (P3ggApp1 + P3ggApp2);
	return res/16.0;
}

double P3GGA_2410(double y, int nf, int imod)
{
    double ym = 1.0 / y;
    double y1 = 1.0 - y;
    double dl = std::log(y);
    double dl1 = std::log(1.0 - y);

    double nf_d = static_cast<double>(nf);
    double nf2 = nf_d * nf_d;
    double nf3 = nf_d * nf2;

    // Large-x coefficients
    double a4gluon = 40880.330 - 11714.246 * nf_d + 440.04876 * nf2 + 7.3627750 * nf3;
    double ccoeff  = 8.5814120e4 - 1.3880515e4 * nf_d + 1.3511111e2 * nf2;
    double dcoeff  = 5.4482808e4 - 4.3411337e3 * nf_d - 2.1333333e1 * nf2;

    double x1l4cff = 5.6460905e1 * nf_d - 3.6213992 * nf2;
    double x1l3cff = 2.4755054e2 * nf_d - 4.0559671e1 * nf2 + 1.5802469 * nf3;

    // Small-x coefficients
    double bfkl0 = -8.3086173e3;
    double bfkl1 = -1.0691199e5 - 9.9638304e2 * nf_d;

    double x0l6cff =  1.44e2 - 2.7786008e1 * nf_d + 7.9012346e-1 * nf2;
    double x0l5cff = -1.44e2 - 1.6208066e2 * nf_d + 1.4380247e1 * nf2;
    double x0l4cff =  2.6165784e4 - 3.3447551e3 * nf_d + 9.1522635e1 * nf2 - 1.9753086e-1 * nf3;

    // Resulting part
    double p3gg01 = bfkl0 * std::pow(dl, 3) * ym
                  + bfkl1 * std::pow(dl, 2) * ym
                  + x0l6cff * std::pow(dl, 6)
                  + x0l5cff * std::pow(dl, 5)
                  + x0l4cff * std::pow(dl, 4)
                  + ccoeff * dl1
                  + dcoeff - a4gluon
                  + x1l4cff * y1 * std::pow(dl1, 4)
                  + x1l3cff * y1 * std::pow(dl1, 3);

    double p3ggapp1 = 0.0;
    double p3ggapp2 = 0.0;

    if (nf == 3) {
        p3ggapp1 = p3gg01
                 - 421311.0  * y1 * dl * ym
                 - 325557.0  * y1 * ym
                 + 1679790.0 * y1
                 - 1456863.0 * y1 * y
                 + 3246307.0 * y1 * dl
                 + 2026324.0 * std::pow(dl, 2)
                 + 549188.0  * std::pow(dl, 3)
                 + 8337.0    * y1 * dl1
                 + 26718.0   * y1 * std::pow(dl1, 2)
                 - 27049.0   * (y1 * y1) * std::pow(dl1, 3);

        p3ggapp2 = p3gg01
                 - 700113.0  * y1 * dl * ym
                 - 2300581.0 * y1 * ym
                 + 896407.0  * y1 * (1.0 + 2.0 * y)
                 - 162733.0  * y1 * (y * y)
                 - 2661862.0 * y1 * dl
                 + 196759.0  * std::pow(dl, 2)
                 - 260607.0  * std::pow(dl, 3)
                 + 84068.0   * y1 * dl1
                 + 346318.0  * y1 * std::pow(dl1, 2)
                 + 315725.0  * dl * std::pow(dl1, 2);

    } else if (nf == 4) {
        p3ggapp1 = p3gg01
                 - 437084.0  * y1 * dl * ym
                 - 361570.0  * y1 * ym
                 + 1696070.0 * y1
                 - 1457385.0 * y1 * y
                 + 3195104.0 * y1 * dl
                 + 2009021.0 * std::pow(dl, 2)
                 + 544380.0  * std::pow(dl, 3)
                 + 9938.0    * y1 * dl1
                 + 24376.0   * y1 * std::pow(dl1, 2)
                 - 22143.0   * (y1 * y1) * std::pow(dl1, 3);

        p3ggapp2 = p3gg01
                 - 706649.0  * y1 * dl * ym
                 - 2274637.0 * y1 * ym
                 + 836544.0  * y1 * (1.0 + 2.0 * y)
                 - 199929.0  * y1 * (y * y)
                 - 2683760.0 * y1 * dl
                 + 168802.0  * std::pow(dl, 2)
                 - 250799.0  * std::pow(dl, 3)
                 + 36967.0   * y1 * dl1
                 + 24530.0   * y1 * std::pow(dl1, 2)
                 - 71470.0   * (y1 * y1) * std::pow(dl1, 2);

    } else if (nf == 5) {
        p3ggapp1 = p3gg01
                 - 439426.0  * y1 * dl * ym
                 - 293679.0  * y1 * ym
                 + 1916281.0 * y1
                 - 1615883.0 * y1 * y
                 + 3648786.0 * y1 * dl
                 + 2166231.0 * std::pow(dl, 2)
                 + 594588.0  * std::pow(dl, 3)
                 + 50406.0   * y1 * dl1
                 + 24692.0   * y1 * std::pow(dl1, 2)
                 + 174067.0  * (y1 * y1) * dl1;

        p3ggapp2 = p3gg01
                 - 705978.0  * y1 * dl * ym
                 - 2192234.0 * y1 * ym
                 + 1730508.0 * y1 * y
                 + 353143.0  * y1 * (2.0 - y * y)
                 - 2602682.0 * y1 * dl
                 + 178960.0  * std::pow(dl, 2)
                 - 218133.0  * std::pow(dl, 3)
                 + 2285.0    * y1 * dl1
                 + 19295.0   * y1 * std::pow(dl1, 2)
                 - 13719.0   * (y1 * y1) * std::pow(dl1, 2);

    } else if (nf == 6) {
        p3ggapp1 = p3gg01
                 - 476018.0  * y1 * dl * ym
                 - 469289.0  * y1 * ym
                 + 2049351.0 * y1
                 - 1589000.0 * y1 * y
                 + 3185549.0 * y1 * dl
                 + 1994521.0 * std::pow(dl, 2)
                 + 527723.0  * std::pow(dl, 3)
                 - 340674.0  * y1 * dl1
                 + 22460.0   * y1 * std::pow(dl1, 2)
                 - 394556.0  * dl * dl1;

        p3ggapp2 = p3gg01
                 - 709863.0  * y1 * dl * ym
                 - 2134347.0 * y1 * ym
                 + 1605315.0 * y1 * y
                 + 360743.0  * y1 * (2.0 - y * y)
                 - 2426250.0 * y1 * dl
                 + 230631.0  * std::pow(dl, 2)
                 - 185804.0  * std::pow(dl, 3)
                 - 7992.9    * y1 * dl1
                 + 15918.0   * y1 * std::pow(dl1, 2)
                 - 32771.0   * (y1 * y1) * dl1;
    } else {
        std::cerr << " Error in function p3gga: choice of nf " << std::endl;
        std::abort();
    }

    if (imod == 1) return p3ggapp1;
    if (imod == 2) return p3ggapp2;
    return (0.5 * (p3ggapp1 + p3ggapp2))/16.0;
}

// -------------------------------------------------------------
// P3GGB
// -------------------------------------------------------------
double P3GGB(double Y, int nf, int IMOD)
{
    using namespace P3GSOFT;

    double nf2 = nf * nf;
    double nf3 = nf * nf2;

    A4gluon =  40880.330
             - 11714.246  * nf
             +   440.04876 * nf2
             +     7.3627750 * nf3;

    return A4gluon/16.0;
}

// -------------------------------------------------------------
// P3GGC
// -------------------------------------------------------------
double P3GGC(double Y, int nf, int IMOD)
{
    using namespace P3GSOFT;

    double nf2 = nf * nf;
    double nf3 = nf * nf2;

    double B4gluon =
          68587.64
        - 18143.983 * nf
        +   423.81135 * nf2
        +     0.90672154 * nf3;

    return B4gluon/16.0;
}


// -------------------------------------------------------------
// P3GQA_2512
// -------------------------------------------------------------
double P3GQA(double Y, int nf, int IMOD)
{
    double YM  = 1.0 / Y;
    double Y1  = 1.0 - Y;
    double DL  = std::log(Y);
    double DL1 = std::log(Y1);

    double nf2 = nf * nf;
    double nf3 = nf * nf2;

    // Large-x coefficients
    double x1L5cff = 1.3443073e1 - 5.4869684e-1 * nf;
    double x1L4cff = 3.7539831e2 - 3.4494742e1 * nf + 8.7791495e-1 * nf2;

    double y1L5cff = 2.2222222e1 - 5.4869684e-1 * nf;
    double y1L4cff = 6.6242163e2 - 4.7992684e1 * nf + 8.7791495e-1 * nf2;

    // Small-x x^-1 coefficients
    double bfkl0 = -8.3086173e3 / 2.25;
    double bfkl1 = (-1.0691199e5 - nf * 9.9638304e2) / 2.25;

    // Small-x double logs (x^0)
    double x0L6cff = 5.2235940e1 - 7.3744856e0 * nf;
    double x0L5cff = -2.9221399e2 + 1.8436214e0 * nf;
    double x0L4cff = 7.3106077e3 - 3.7887135e2 * nf - 3.2438957e1 * nf2;

    // Base contribution
    double P3gq01 =
          bfkl0 * YM * std::pow(DL, 3)
        + bfkl1 * YM * std::pow(DL, 2)
        + x0L6cff * std::pow(DL, 6)
        + x0L5cff * std::pow(DL, 5)
        + x0L4cff * std::pow(DL, 4)
        + x1L4cff * std::pow(DL1, 4)
        + x1L5cff * std::pow(DL1, 5)
        + y1L4cff * Y1 * std::pow(DL1, 4)
        + y1L5cff * Y1 * std::pow(DL1, 5);

    double P3gqApp1 = 0.0;
    double P3gqApp2 = 0.0;

    if (nf == 3) {
        P3gqApp1 = P3gq01
            + 3.5 * bfkl1 * YM * DL
            - 27891.0 * YM * Y1
            - 309124.0
            + 1056866.0 * Y * (2.0 - Y)
            - 124735.0 * DL
            - 16246.0 * std::pow(DL, 2)
            + 131175.0 * std::pow(DL, 3)
            + 4970.1 * std::pow(DL1, 3)
            + 60041.0 * std::pow(DL1, 2)
            + 343181.0 * DL1
            - 958330.0 * DL * DL1;

        P3gqApp2 = P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1139334.0 * YM * Y1
            + 143008.0
            - 290390.0 * Y * (2.0 - Y)
            - 659492.0 * DL
            + 303685.0 * std::pow(DL, 2)
            - 81867.0 * std::pow(DL, 3)
            + 1811.8 * std::pow(DL1, 3)
            - 465.9 * std::pow(DL1, 2)
            - 51206.0 * DL1
            + 274249.0 * DL * DL1;
    }
    else if (nf == 4) {
        P3gqApp1 = P3gq01
            + 3.5 * bfkl1 * YM * DL
            - 8302.8 * YM * Y1
            - 347706.0
            + 1105306.0 * Y * (2.0 - Y)
            - 127650.0 * DL
            - 29728.0 * std::pow(DL, 2)
            + 137537.0 * std::pow(DL, 3)
            + 4658.1 * std::pow(DL1, 3)
            + 59205.0 * std::pow(DL1, 2)
            + 345513.0 * DL1
            - 995120.0 * DL * DL1;

        P3gqApp2 = P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1129822.0 * YM * Y1
            + 108527.0
            - 254166.0 * Y * (2.0 - Y)
            - 667254.0 * DL
            + 293099.0 * std::pow(DL, 2)
            - 77437.0 * std::pow(DL, 3)
            + 1471.3 * std::pow(DL1, 3)
            - 1850.3 * std::pow(DL1, 2)
            - 52451.0 * DL1
            + 248634.0 * DL * DL1;
    }
    else if (nf == 5) {
        P3gqApp1 = P3gq01
            + 3.5 * bfkl1 * YM * DL
            + 14035.0 * YM * Y1
            - 384003.0
            + 1152711.0 * Y * (2.0 - Y)
            - 126346.0 * DL
            - 42967.0 * std::pow(DL, 2)
            + 144270.0 * std::pow(DL, 3)
            + 4385.5 * std::pow(DL1, 3)
            + 58688.0 * std::pow(DL1, 2)
            + 348988.0 * DL1
            - 1031165.0 * DL * DL1;

        P3gqApp2 = P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1117561.0 * YM * Y1
            + 76329.0
            - 218973.0 * Y * (2.0 - Y)
            - 670799.0 * DL
            + 282763.0 * std::pow(DL, 2)
            - 72633.0 * std::pow(DL, 3)
            + 1170.0 * std::pow(DL1, 3)
            - 2915.5 * std::pow(DL1, 2)
            - 52548.0 * DL1
            + 223771.0 * DL * DL1;
    }
    else if (nf == 6) {
        P3gqApp1 = P3gq01
            + 3.5 * bfkl1 * YM * DL
            + 39203.0 * YM * Y1
            - 417914.0
            + 1199042.0 * Y * (2.0 - Y)
            - 120750.0 * DL
            - 55941.0 * std::pow(DL, 2)
            + 151383.0 * std::pow(DL, 3)
            + 4149.2 * std::pow(DL1, 3)
            + 58466.0 * std::pow(DL1, 2)
            + 353589.0 * DL1
            - 1066510.0 * DL * DL1;

        P3gqApp2 = P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1102470.0 * YM * Y1
            + 46517.0
            - 184858.0 * Y * (2.0 - Y)
            - 670056.0 * DL
            + 272689.0 * std::pow(DL, 2)
            - 67453.0 * std::pow(DL, 3)
            + 905.0 * std::pow(DL1, 3)
            - 3686.2 * std::pow(DL1, 2)
            - 51523.0 * DL1
            + 199594.0 * DL * DL1;
    }
    else {
        throw std::invalid_argument("Error in P3GQA_2512: invalid nf (must be 3–6)");
    }

	double res{};
    if (IMOD == 1)
        res = P3gqApp1;
    else if (IMOD == 2)
        res = P3gqApp2;
    else
        res = 0.5 * (P3gqApp1 + P3gqApp2);
	return res/16.0;
}

// -------------------------------------------------------------
// P3QGA
// -------------------------------------------------------------
double P3QGA(double Y, int nf, int IMOD)
{
    double YM  = 1.0 / Y;
    double Y1  = 1.0 - Y;
    double DL  = std::log(Y);
    double DL1 = std::log(1.0 - Y);

    double nf2 = nf * nf;
    double nf3 = nf * nf2;

    // Large-x coefficients
    double x1L5cff = 1.8518519e0 * nf - 4.1152263e-1 * nf2;
    double x1L4cff = 3.5687794e1 * nf - 3.5116598e0 * nf2
                   - 8.2304527e-2 * nf3;

    double y1L5cff = 2.8806584e0 * nf + 8.2304527e-1 * nf2;
    double y1L4cff = -4.0511391e1 * nf + 5.5418381e0 * nf2
                   + 1.6460905e-1 * nf3;

    // Small-x coefficients
    double bfkl1   = 3.9357613e3 * nf;

    double x0L6cff = -1.9588477e1 * nf + 2.7654321e0 * nf2;
    double x0L5cff =  2.1573663e1 * nf + 1.7244444e1 * nf2;
    double x0L4cff = -2.8667643e3 * nf + 3.0122403e2 * nf2
                   + 4.1316872e0 * nf3;

    // Base contribution
    double P3QG01 =
          bfkl1 * YM * std::pow(DL, 2)
        + x0L6cff * std::pow(DL, 6)
        + x0L5cff * std::pow(DL, 5)
        + x0L4cff * std::pow(DL, 4)
        + x1L4cff * std::pow(DL1, 4)
        + x1L5cff * std::pow(DL1, 5)
        + y1L4cff * Y1 * std::pow(DL1, 4)
        + y1L5cff * Y1 * std::pow(DL1, 5);

    double P3qgApp1 = 0.0;
    double P3qgApp2 = 0.0;

    if (nf == 3) {
        P3qgApp1 = P3QG01
            + 187500.0 * YM * DL
            + 826060.0 * YM * Y1
            - 150474.0
            + 226254.0 * Y * (2.0 - Y)
            + 577733.0 * DL
            - 180747.0 * std::pow(DL, 2)
            + 95411.0  * std::pow(DL, 3)
            + 119.8    * std::pow(DL1, 3)
            + 7156.3   * std::pow(DL1, 2)
            + 45790.0  * DL1
            - 95682.0  * DL * DL1;

        P3qgApp2 = P3QG01
            + 135000.0 * YM * DL
            + 484742.0 * YM * Y1
            - 11627.0
            - 187478.0 * Y * (2.0 - Y)
            + 413512.0 * DL
            - 82500.0  * std::pow(DL, 2)
            + 29987.0  * std::pow(DL, 3)
            - 850.1    * std::pow(DL1, 3)
            - 11425.0  * std::pow(DL1, 2)
            - 75323.0  * DL1
            + 282836.0 * DL * DL1;
    }
    else if (nf == 4) {
        P3qgApp1 = P3QG01
            + 250000.0 * YM * DL
            + 1089180.0 * YM * Y1
            - 241088.0
            + 342902.0 * Y * (2.0 - Y)
            + 720081.0 * DL
            - 247071.0 * std::pow(DL, 2)
            + 126405.0 * std::pow(DL, 3)
            + 272.4    * std::pow(DL1, 3)
            + 10911.0  * std::pow(DL1, 2)
            + 60563.0  * DL1
            - 161448.0 * DL * DL1;

        P3qgApp2 = P3QG01
            + 180000.0 * YM * DL
            + 634090.0 * YM * Y1
            - 55958.0
            - 208744.0 * Y * (2.0 - Y)
            + 501120.0 * DL
            - 116073.0 * std::pow(DL, 2)
            + 39173.0  * std::pow(DL, 3)
            - 1020.8   * std::pow(DL1, 3)
            - 13864.0  * std::pow(DL1, 2)
            - 100922.0 * DL1
            + 343243.0 * DL * DL1;
    }
    else if (nf == 5) {
        P3qgApp1 = P3QG01
            + 312500.0 * YM * DL
            + 1345700.0 * YM * Y1
            - 350466.0
            + 480028.0 * Y * (2.0 - Y)
            + 837903.0 * DL
            - 315928.0 * std::pow(DL, 2)
            + 157086.0 * std::pow(DL, 3)
            + 472.7    * std::pow(DL1, 3)
            + 15415.0  * std::pow(DL1, 2)
            + 75644.0  * DL1
            - 244869.0 * DL * DL1;

        P3qgApp2 = P3QG01
            + 225000.0 * YM * DL
            + 776837.0 * YM * Y1
            - 119054.0
            - 209530.0 * Y * (2.0 - Y)
            + 564202.0 * DL
            - 152181.0 * std::pow(DL, 2)
            + 48046.0  * std::pow(DL, 3)
            - 1143.8   * std::pow(DL1, 3)
            - 15553.0  * std::pow(DL1, 2)
            - 126212.0 * DL1
            + 385995.0 * DL * DL1;
    }
    else if (nf == 6) {
        P3qgApp1 = P3QG01
            + 375000.0 * YM * DL
            + 1595330.0 * YM * Y1
            - 477729.0
            + 637552.0 * Y * (2.0 - Y)
            + 931556.0 * DL
            - 387017.0 * std::pow(DL, 2)
            + 187509.0 * std::pow(DL, 3)
            + 715.5    * std::pow(DL1, 3)
            + 20710.0  * std::pow(DL1, 2)
            + 91373.0  * DL1
            - 346374.0 * DL * DL1;

        P3qgApp2 = P3QG01
            + 270000.0 * YM * DL
            + 912695.0 * YM * Y1
            - 200034.0
            - 189918.0 * Y * (2.0 - Y)
            + 603114.0 * DL
            - 190521.0 * std::pow(DL, 2)
            + 56661.0  * std::pow(DL, 3)
            - 1224.3   * std::pow(DL1, 3)
            - 16453.0  * std::pow(DL1, 2)
            - 150856.0 * DL1
            + 410661.0 * DL * DL1;
    }
    else {
        throw std::invalid_argument("Error in P3QGA: invalid nf (must be 3–6)");
    }

	double res{};
    if (IMOD == 1)
        res = P3qgApp1;
    else if (IMOD == 2)
        res = P3qgApp2;
    else
        res = 0.5 * (P3qgApp1 + P3qgApp2);
	return res/16.0;
}

// ------------------------------------------------------------------
// P3GQA_2512
// ------------------------------------------------------------------
double P3GQA_2512(double Y, int nf, int IMOD)
{
    double YM  = 1.0 / Y;
    double Y1  = 1.0 - Y;
    double DL  = std::log(Y);
    double DL1 = std::log(Y1);

    int nf2 = nf * nf;
    int nf3 = nf * nf2;   // (kept for strict faithfulness)

    // --------------------------------------------------------------
    // Known large-x coefficients
    // --------------------------------------------------------------

    double x1L5cff = 1.3443073e1 - 5.4869684e-1 * nf;
    double x1L4cff = 3.7539831e2 - 3.4494742e1 * nf
                   + 8.7791495e-1 * nf2;

    double y1L5cff = 2.2222222e1 - 5.4869684e-1 * nf;
    double y1L4cff = 6.6242163e2 - 4.7992684e1 * nf
                   + 8.7791495e-1 * nf2;

    // --------------------------------------------------------------
    // Small-x x^-1 coefficients
    // --------------------------------------------------------------

    double bfkl0 = -8.3086173e3 / 2.25;
    double bfkl1 = (-1.0691199e5 - nf * 9.9638304e2) / 2.25;

    // --------------------------------------------------------------
    // Small-x double logs
    // --------------------------------------------------------------

    double x0L6cff =  5.2235940e1 - 7.3744856e0 * nf;
    double x0L5cff = -2.9221399e2 + 1.8436214e0 * nf;
    double x0L4cff =  7.3106077e3 - 3.7887135e2 * nf
                    - 3.2438957e1 * nf2;

    // --------------------------------------------------------------
    // Base function: P3gq01
    // --------------------------------------------------------------

    double P3gq01 =
          bfkl0 * YM * std::pow(DL, 3)
        + bfkl1 * YM * std::pow(DL, 2)
        + x0L6cff * std::pow(DL, 6)
        + x0L5cff * std::pow(DL, 5)
        + x0L4cff * std::pow(DL, 4)
        + x1L4cff * std::pow(DL1, 4)
        + x1L5cff * std::pow(DL1, 5)
        + y1L4cff * Y1 * std::pow(DL1, 4)
        + y1L5cff * Y1 * std::pow(DL1, 5);

    double P3gqApp1 = 0.0;
    double P3gqApp2 = 0.0;

    // --------------------------------------------------------------
    // nf-specific approximations
    // --------------------------------------------------------------

    if (nf == 3)
    {
        P3gqApp1 =
              P3gq01
            + 3.5 * bfkl1 * YM * DL
            - 27891.  * YM * Y1
            - 309124.
            + 1056866. * Y * (2.0 - Y)
            - 124735. * DL
            - 16246.  * std::pow(DL,2)
            + 131175. * std::pow(DL,3)
            + 4970.1  * std::pow(DL1,3)
            + 60041.  * std::pow(DL1,2)
            + 343181. * DL1
            - 958330. * DL * DL1;

        P3gqApp2 =
              P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1139334. * YM * Y1
            + 143008.
            - 290390. * Y * (2.0 - Y)
            - 659492. * DL
            + 303685. * std::pow(DL,2)
            - 81867.  * std::pow(DL,3)
            + 1811.8  * std::pow(DL1,3)
            - 465.9   * std::pow(DL1,2)
            - 51206.  * DL1
            + 274249. * DL * DL1;
    }
    else if (nf == 4)
    {
        P3gqApp1 =
              P3gq01
            + 3.5 * bfkl1 * YM * DL
            - 8302.8 * YM * Y1
            - 347706.
            + 1105306. * Y * (2.0 - Y)
            - 127650. * DL
            - 29728.  * std::pow(DL,2)
            + 137537. * std::pow(DL,3)
            + 4658.1  * std::pow(DL1,3)
            + 59205.  * std::pow(DL1,2)
            + 345513. * DL1
            - 995120. * DL * DL1;

        P3gqApp2 =
              P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1129822. * YM * Y1
            + 108527.
            - 254166. * Y * (2.0 - Y)
            - 667254. * DL
            + 293099. * std::pow(DL,2)
            - 77437.  * std::pow(DL,3)
            + 1471.3  * std::pow(DL1,3)
            - 1850.3  * std::pow(DL1,2)
            - 52451.  * DL1
            + 248634. * DL * DL1;
    }
    else if (nf == 5)
    {
        P3gqApp1 =
              P3gq01
            + 3.5 * bfkl1 * YM * DL
            + 14035. * YM * Y1
            - 384003.
            + 1152711. * Y * (2.0 - Y)
            - 126346. * DL
            - 42967.  * std::pow(DL,2)
            + 144270. * std::pow(DL,3)
            + 4385.5  * std::pow(DL1,3)
            + 58688.  * std::pow(DL1,2)
            + 348988. * DL1
            - 1031165. * DL * DL1;

        P3gqApp2 =
              P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1117561. * YM * Y1
            + 76329.
            - 218973. * Y * (2.0 - Y)
            - 670799. * DL
            + 282763. * std::pow(DL,2)
            - 72633.  * std::pow(DL,3)
            + 1170.0  * std::pow(DL1,3)
            - 2915.5  * std::pow(DL1,2)
            - 52548.  * DL1
            + 223771. * DL * DL1;
    }
    else if (nf == 6)
    {
        P3gqApp1 =
              P3gq01
            + 3.5 * bfkl1 * YM * DL
            + 39203. * YM * Y1
            - 417914.
            + 1199042. * Y * (2.0 - Y)
            - 120750. * DL
            - 55941.  * std::pow(DL,2)
            + 151383. * std::pow(DL,3)
            + 4149.2  * std::pow(DL1,3)
            + 58466.  * std::pow(DL1,2)
            + 353589. * DL1
            - 1066510. * DL * DL1;

        P3gqApp2 =
              P3gq01
            + 7.0 * bfkl1 * YM * DL
            - 1102470. * YM * Y1
            + 46517.
            - 184858. * Y * (2.0 - Y)
            - 670056. * DL
            + 272689. * std::pow(DL,2)
            - 67453.  * std::pow(DL,3)
            + 905.0   * std::pow(DL1,3)
            - 3686.2  * std::pow(DL1,2)
            - 51523.  * DL1
            + 199594. * DL * DL1;
    }
    else
    {
        throw std::runtime_error("Error in P3GQA_2512: nf must be 3..6");
    }

    // --------------------------------------------------------------
    // Return according to IMOD
    // --------------------------------------------------------------

	double res{};
    if (IMOD == 1)
        res = P3gqApp1;
    else if (IMOD == 2)
        res = P3gqApp2;
    else
        res = 0.5 * (P3gqApp1 + P3gqApp2);
	return res/16.0;
}


// -------------------------------------------------------------
// P3PSA
// -------------------------------------------------------------
double P3PSA(double Y, int nf, int IMOD)
{
    double YM  = 1.0 / Y;
    double Y1  = 1.0 - Y;
    double DL  = std::log(Y);
    double DL1 = std::log(1.0 - Y);

    double nf2 = nf * nf;
    double nf3 = nf * nf2;

    // Large-x coefficients
    double x1L4cff = -5.6460905e1 * nf + 3.6213992e0 * nf2;
    double x1L3cff = -2.4755054e2 * nf + 4.0559671e1 * nf2
                   - 1.5802469e0 * nf3;

    double y1L4cff = -1.3168724e1 * nf;
    double y1L3cff = -1.9911111e2 * nf + 1.3695473e1 * nf2;

    // Small-x coefficients
    double bfkl1   = 1.7492273e3 * nf;

    double x0L6cff = -7.5061728e0 * nf + 7.9012346e-1 * nf2;
    double x0L5cff =  2.8549794e1 * nf + 3.7925926e0 * nf2;
    double x0L4cff = -8.5480010e2 * nf + 7.7366255e1 * nf2
                   - 1.9753086e-1 * nf3;

    // Base contribution
    double P3ps01 =
          bfkl1 * std::pow(DL, 2) * YM
        + x0L6cff * std::pow(DL, 6)
        + x0L5cff * std::pow(DL, 5)
        + x0L4cff * std::pow(DL, 4)
        + x1L3cff * Y1 * std::pow(DL1, 3)
        + x1L4cff * Y1 * std::pow(DL1, 4)
        + y1L3cff * Y1 * Y1 * std::pow(DL1, 3)
        + y1L4cff * Y1 * Y1 * std::pow(DL1, 4);

    double P3psApp1 = 0.0;
    double P3psApp2 = 0.0;

    if (nf == 3) {
        P3psApp1 = P3ps01
            + 67731.0  * Y1 * DL * YM
            + 274100.0 * Y1 * YM
            - 104493.0 * Y1 * (1.0 + 2.0 * Y)
            + 34403.0  * Y1 * Y * Y
            + 353656.0 * Y1 * DL
            + 10620.0  * std::pow(DL, 2)
            + 40006.0  * std::pow(DL, 3)
            - 7412.1   * Y1 * DL1
            - 2365.1   * Y1 * std::pow(DL1, 2)
            + 1533.0   * Y1 * Y1 * std::pow(DL1, 2);

        P3psApp2 = P3ps01
            + 54593.0  * Y1 * DL * YM
            + 179748.0 * Y1 * YM
            - 195263.0 * Y1
            + 12789.0  * Y1 * Y * (1.0 + Y)
            + 4700.0   * Y1 * DL
            - 103604.0 * std::pow(DL, 2)
            - 2758.3   * std::pow(DL, 3)
            - 2801.2   * Y1 * DL1
            - 1986.9   * Y1 * std::pow(DL1, 2)
            - 6005.9   * Y1 * Y1 * std::pow(DL1, 2);
    }
    else if (nf == 4) {
        P3psApp1 = P3ps01
            + 90154.0  * Y1 * DL * YM
            + 359084.0 * Y1 * YM
            - 136319.0 * Y1 * (1.0 + 2.0 * Y)
            + 45379.0  * Y1 * Y * Y
            + 461167.0 * Y1 * DL
            + 13869.0  * std::pow(DL, 2)
            + 52525.0  * std::pow(DL, 3)
            - 7498.2   * Y1 * DL1
            - 2491.5   * Y1 * std::pow(DL1, 2)
            + 1727.2   * Y1 * Y1 * std::pow(DL1, 2);

        P3psApp2 = P3ps01
            + 72987.0  * Y1 * DL * YM
            + 235802.0 * Y1 * YM
            - 254921.0 * Y1
            + 17138.0  * Y1 * Y * (1.0 + Y)
            + 5212.9   * Y1 * DL
            - 135378.0 * std::pow(DL, 2)
            - 3350.9   * std::pow(DL, 3)
            - 1472.7   * Y1 * DL1
            - 1997.2   * Y1 * std::pow(DL1, 2)
            - 8123.3   * Y1 * Y1 * std::pow(DL1, 2);
    }
    else if (nf == 5) {
        P3psApp1 = P3ps01
            + 112481.0 * Y1 * DL * YM
            + 440555.0 * Y1 * YM
            - 166581.0 * Y1 * (1.0 + 2.0 * Y)
            + 56087.0  * Y1 * Y * Y
            + 562992.0 * Y1 * DL
            + 16882.0  * std::pow(DL, 2)
            + 64577.0  * std::pow(DL, 3)
            - 6570.1   * Y1 * DL1
            - 2365.7   * Y1 * std::pow(DL1, 2)
            + 1761.7   * Y1 * Y1 * std::pow(DL1, 2);

        P3psApp2 = P3ps01
            + 91468.0  * Y1 * DL * YM
            + 289658.0 * Y1 * YM
            - 311749.0 * Y1
            + 21521.0  * Y1 * Y * (1.0 + Y)
            + 4908.9   * Y1 * DL
            - 165795.0 * std::pow(DL, 2)
            - 3814.9   * std::pow(DL, 3)
            + 804.5    * Y1 * DL1
            - 1760.8   * Y1 * std::pow(DL1, 2)
            - 10295.0  * Y1 * Y1 * std::pow(DL1, 2);
    }
    else if (nf == 6) {
        P3psApp1 = P3ps01
            + 134701.0 * Y1 * DL * YM
            + 518318.0 * Y1 * YM
            - 195241.0 * Y1 * (1.0 + 2.0 * Y)
            + 66517.0  * Y1 * Y * Y
            + 658832.0 * Y1 * DL
            + 19605.0  * std::pow(DL, 2)
            + 76125.0  * std::pow(DL, 3)
            - 4734.5   * Y1 * DL1
            - 2035.2   * Y1 * std::pow(DL1, 2)
            + 1633.1   * Y1 * Y1 * std::pow(DL1, 2);

        P3psApp2 = P3ps01
            + 110032.0 * Y1 * DL * YM
            + 341158.0 * Y1 * YM
            - 365676.0 * Y1
            + 25934.0  * Y1 * Y * (1.0 + Y)
            + 3614.4   * Y1 * DL
            - 194868.0 * std::pow(DL, 2)
            - 4172.2   * std::pow(DL, 3)
            + 3924.3   * Y1 * DL1
            - 1324.9   * Y1 * std::pow(DL1, 2)
            - 12520.0  * Y1 * Y1 * std::pow(DL1, 2);
    }
    else {
        throw std::invalid_argument("Error in P3PSA: invalid nf (must be 3–6)");
    }

	double res{};
    if (IMOD == 1)
        res = P3psApp1;
    else if (IMOD == 2)
        res = P3psApp2;
    else
        res = 0.5 * (P3psApp1 + P3psApp2);
	return res/16.0;
}

// -------------------------------------------------------------
// P3NSMA
// -------------------------------------------------------------
double P3NSMA(double Y, int nf, int IMOD)
{
    double Y1  = 1.0 - Y;
    double DM  = 1.0 / Y1;
    double DL  = std::log(Y);
    double DL1 = std::log(1.0 - Y);

    // Leading large-nc
    double P3NSA0 =
        2.5e4 * (
            Y1 * (3.5254 + 8.6935 * Y - 1.5051 * std::pow(Y,2)
                 + 1.8300 * std::pow(Y,3))
            + 11.883 * Y * DL
            - 0.09066 * Y * std::pow(DL,2)
            + 11.410 * Y1 * DL1
            + 13.376 * DL * DL1
        )
        + 5.167133e4 * DL
        + 1.712095e4 * std::pow(DL,2)
        + 2.863226e3 * std::pow(DL,3)
        + 2.978255e2 * std::pow(DL,4)
        + 1.6e1 * std::pow(DL,5)
        + 5.0e-1 * std::pow(DL,6)
        - 2.973385e4
        + 1.906980e4 * DL1;

    double P3NSA1 =
        2.5e4 * (
            Y1 * (-0.74077 + 1.4860 * Y - 0.23631 * std::pow(Y,2)
                 + 0.31584 * std::pow(Y,3))
            + 2.5251 * Y1 * DL1
            + 2.5203 * DL * DL1
            + 2.2242 * Y * DL
            - 0.02460 * Y * std::pow(DL,2)
            + 0.00310 * Y * std::pow(DL,3)
        )
        - 9.239374e3 * DL
        - 2.917312e3 * std::pow(DL,2)
        - 4.305308e2 * std::pow(DL,3)
        - 3.6e1 * std::pow(DL,4)
        - (4.0/3.0) * std::pow(DL,5)
        + 8.115605e3
        - 3.079761e3 * DL1;

    // Nonleading large-nc approximations
    double P3NMA01 =
        (5992.88 * (1.0 + 2.0*Y) + 31321.44 * Y*Y) * Y1
        + 511.228
        - 1618.07 * DL
        + 2.25480 * std::pow(DL,3)
        + 31897.82 * DL1 * Y1
        + 4653.76 * std::pow(DL1,2) * Y1
        + 4.964335e-1 * (std::pow(DL,6) + 6.0 * std::pow(DL,5))
        - 2.601749e3
        - 2.118867e3 * DL1;

    double P3NMA02 =
        (4043.59 - 15386.6 * Y) * Y * Y1
        + 502.481
        + 1532.96 * std::pow(DL,2)
        + 31.6023 * std::pow(DL,3)
        - 3997.39 * DL1 * Y1
        + 511.567 * std::pow(DL1,3) * Y1
        + 4.964335e-1 * (std::pow(DL,6) + 18.0 * std::pow(DL,5))
        - 2.601749e3
        - 2.118867e3 * DL1;

    double P3NMA11 =
        (114.457 * (1.0 + 2.0*Y) + 2570.73 * Y*Y) * Y1
        - 7.08645
        - 127.012 * std::pow(DL,2)
        + 2.69618 * std::pow(DL,4)
        + 1856.63 * DL1 * Y1
        + 440.17 * std::pow(DL1,2) * Y1
        + 3.121643e2
        + 3.379310e2 * DL1;

    double P3NMA12 =
        (-335.995 * (2.0 + Y) - 1605.91 * Y*Y) * Y1
        - 7.82077
        - 9.76627 * std::pow(DL,2)
        + 0.14218 * std::pow(DL,5)
        - 1360.04 * DL1 * Y1
        + 38.7337 * std::pow(DL1,3) * Y1
        + 3.121643e2
        + 3.379310e2 * DL1;

    // nf^2 and nf^3
    double P3NSMA2 =
        2.5e2 * (
            Y1 * (3.2206 + 1.7507 * Y + 0.13281 * std::pow(Y,2)
                 + 0.45969 * std::pow(Y,3))
            + 1.5641 * Y * DL
            - 0.37902 * Y * std::pow(DL,2)
            - 0.03248 * Y * std::pow(DL,3)
            + 2.7511 * Y1 * DL1
            + 3.2709 * DL * DL1
        )
        + 4.378810e2 * DL
        + 1.282948e2 * std::pow(DL,2)
        + 1.959945e1 * std::pow(DL,3)
        + 9.876543e-1 * std::pow(DL,4)
        - 3.760092e2
        + 2.668861e1 * DL1;

    double P3NSA3 =
        -2.426296
        - 0.8460488 * Y
        + (0.5267490 * DM - 3.687243 + 3.160494 * Y) * DL
        - (1.316872 * (DM + 0.1) - 1.448560 * Y) * std::pow(DL,2)
        - (0.2633744 * DM - 0.131687 * (1.0 + Y)) * std::pow(DL,3);

    double P3NSMAI =
        P3NSA0
        + nf * P3NSA1
        + std::pow(nf,2) * P3NSMA2
        + std::pow(nf,3) * P3NSA3;

	double res{};
    if (IMOD == 1)
        res = P3NSMAI + P3NMA01 + nf * P3NMA11;
    else if (IMOD == 2)
        res = P3NSMAI + P3NMA02 + nf * P3NMA12;
    else
        res = P3NSMAI
            + 0.5 * ((P3NMA01 + P3NMA02)
            + nf * (P3NMA11 + P3NMA12));
	return res/16.0;
}

double P3NSMB(double Y, int nf, int IMOD)
{
    double A4qI =
        2.120902e4
        - 5.179372e3 * nf
        + 1.955772e2 * std::pow(nf,2)
        + 3.272344e0 * std::pow(nf,3);

    double A4ap1 = -511.228 + 7.08645 * nf;
    double A4ap2 = -502.481 + 7.82077 * nf;

	double res{};
    if (IMOD == 1)
        res = (A4qI + A4ap1);
    else if (IMOD == 2)
        res = (A4qI + A4ap2);
    else
        res = (A4qI + 0.5 * (A4ap1 + A4ap2));
	return res/16.0;
}

double P3NSMC(double Y, int nf, int IMOD)
{
    double B4qI =
        2.579609e4 + 0.08
        - (5.818637e3 + 0.97) * nf
        + (1.938554e2 + 0.0037) * std::pow(nf,2)
        + 3.014982e0 * std::pow(nf,3);

    double B4ap1 = -2426.05 + 266.674 * nf - 0.05 * nf;
    double B4ap2 = -2380.255 + 270.518 * nf - 0.05 * nf;

	double res{};
    if (IMOD == 1)
        res = B4qI + B4ap1;
    else if (IMOD == 2)
        res = B4qI + B4ap2;
    else
        res = + B4qI + 0.5*(B4ap1 + B4ap2);
	return res/16.0;
}

// -------------------------------------------------------------
// P3NSPA
// -------------------------------------------------------------
double P3NSPA(double Y, int nf, int IMOD)
{
    double Y1  = 1.0 - Y;
    double DM  = 1.0 / Y1;
    double DL  = std::log(Y);
    double DL1 = std::log1p(-Y);

    // Leading large-nc
    double P3NSA0 =
        2.5e4 * (
            Y1 * (3.5254 + 8.6935*Y - 1.5051*std::pow(Y,2)
                 + 1.8300*std::pow(Y,3))
            + 11.883*Y*DL
            - 0.09066*Y*std::pow(DL,2)
            + 11.410*Y1*DL1
            + 13.376*DL*DL1
        )
        + 5.167133e4*DL
        + 1.712095e4*std::pow(DL,2)
        + 2.863226e3*std::pow(DL,3)
        + 2.978255e2*std::pow(DL,4)
        + 1.6e1*std::pow(DL,5)
        + 5.0e-1*std::pow(DL,6)
        - 2.973385e4
        + 1.906980e4*DL1;

    double P3NSA1 =
        2.5e4 * (
            Y1 * (-0.74077 + 1.4860*Y - 0.23631*std::pow(Y,2)
                 + 0.31584*std::pow(Y,3))
            + 2.5251*Y1*DL1
            + 2.5203*DL*DL1
            + 2.2242*Y*DL
            - 0.02460*Y*std::pow(DL,2)
            + 0.00310*Y*std::pow(DL,3)
        )
        - 9.239374e3*DL
        - 2.917312e3*std::pow(DL,2)
        - 4.305308e2*std::pow(DL,3)
        - 3.6e1*std::pow(DL,4)
        - (4.0/3.0)*std::pow(DL,5)
        + 8.115605e3
        - 3.079761e3*DL1;

    // Nonleading approximations
    double P3NPA01 =
        3948.16*Y1
        - 2464.61*(2.0*Y - Y*Y)*Y1
        - 1839.44*std::pow(DL,2)
        - 402.156*std::pow(DL,3)
        - 1777.27*std::pow(DL1,2)*Y1
        - 204.183*std::pow(DL1,3)*Y1
        + 507.152
        - 5.587553e1*std::pow(DL,4)
        - 2.831276e0*std::pow(DL,5)
        - 1.488340e-1*std::pow(DL,6)
        - 2.601749e3
        - 2.118867e3*DL1;

    double P3NPA02 =
        (8698.39 - 10490.47*Y)*Y*Y1
        + 1389.73*DL
        + 189.576*std::pow(DL,2)
        - 173.936*std::pow(DL1,2)*Y1
        + 223.078*std::pow(DL1,3)*Y1
        + 505.209
        - 5.587553e1*std::pow(DL,4)
        - 2.831276e0*std::pow(DL,5)
        - 1.488340e-1*std::pow(DL,6)
        - 2.601749e3
        - 2.118867e3*DL1;

    double P3NPA11 =
        (-1116.34 + 1071.24*Y)*Y*Y1
        - 59.3041*std::pow(DL,2)
        - 8.4620*std::pow(DL,3)
        - 143.813*DL1*Y1
        - 18.8803*std::pow(DL1,3)*Y1
        - 7.33927
        + 4.658436e0*std::pow(DL,4)
        + 2.798354e-1*std::pow(DL,5)
        + 3.121643e2
        + 3.379310e2*DL1;

    double P3NPA12 =
        (-690.151 - 656.386*Y*Y)*Y1
        + 133.702*std::pow(DL,2)
        + 34.0569*std::pow(DL,3)
        - 745.573*DL1*Y1
        + 8.61438*std::pow(DL1,3)*Y1
        - 7.53662
        + 4.658437e0*std::pow(DL,4)
        + 2.798354e-1*std::pow(DL,5)
        + 3.121643e2
        + 3.379310e2*DL1;

    double P3NSPA2 =
        2.5e2 * (
            Y1 * (3.0008 + 0.8619*Y - 0.12411*std::pow(Y,2)
                 + 0.31595*std::pow(Y,3))
            - 0.37529*Y*DL
            - 0.21684*Y*std::pow(DL,2)
            - 0.02295*Y*std::pow(DL,3)
            + 0.03394*Y1*DL1
            + 0.40431*DL*DL1
        )
        + 3.930056e2*DL
        + 1.125705e2*std::pow(DL,2)
        + 1.652675e1*std::pow(DL,3)
        + 7.901235e-1*std::pow(DL,4)
        - 3.760092e2
        + 2.668861e1*DL1;

    double P3NSA3 =
        -2.426296
        - 0.8460488*Y
        + (0.5267490*DM - 3.687243 + 3.160494*Y)*DL
        - (1.316872*(DM + 0.1) - 1.448560*Y)*std::pow(DL,2)
        - (0.2633745*DM - 0.131687*(1.0+Y))*std::pow(DL,3);

    double P3NSPAI =
        P3NSA0
        + nf*P3NSA1
        + std::pow(nf,2)*P3NSPA2
        + std::pow(nf,3)*P3NSA3;

	double res{};
    if (IMOD == 1)
        res = P3NSPAI + P3NPA01 + nf*P3NPA11;
    else if (IMOD == 2)
        res = P3NSPAI + P3NPA02 + nf*P3NPA12;
    else
        res = P3NSPAI
            + 0.5*((P3NPA01 + P3NPA02)
            + nf*(P3NPA11 + P3NPA12));
	return res;
}

double P3NSPB(double Y, int nf, int IMOD)
{
    double A4qI =
        2.120902e4
        - 5.179372e3*nf
        + 1.955772e2*std::pow(nf,2)
        + 3.272344e0*std::pow(nf,3);

    double A4ap1 = -507.152 + 7.33927*nf;
    double A4ap2 = -505.209 + 7.53662*nf;

	double res{};
    if (IMOD == 1)
        res = (A4qI + A4ap1);
    else if (IMOD == 2)
        res = (A4qI + A4ap2);
    else
        res = (A4qI + 0.5*(A4ap1 + A4ap2));
	return res;
}

double P3NSPC(double Y, int nf, int IMOD)
{
    double B4qI =
        2.579609e4 + 0.08
        - (5.818637e3 + 0.97)*nf
        + (1.938554e2 + 0.0037)*std::pow(nf,2)
        + 3.014982e0*std::pow(nf,3);

    double B4ap1 = -2405.03 + 267.965*nf;
    double B4ap2 = -2394.47 + 269.028*nf;

	double res{};
    if (IMOD == 1)
        res = B4qI + B4ap1;
    else if (IMOD == 2)
        res = B4qI + B4ap2;
    else
        res = B4qI + 0.5*(B4ap1 + B4ap2);
	return res;
}

// -------------------------------------------------------------
// P3NSSA
// -------------------------------------------------------------
double P3NSSA(double Y, int nf, int IMOD)
{
    double Y1  = 1.0 - Y;
    double DM  = 1.0 / Y1;          // kept for structural consistency
    double DL  = std::log(Y);
    double DL1 = std::log(1.0 - Y);

    // ---------------------------------------------------------
    // nf^1 : two approximations
    // ---------------------------------------------------------

    double P3NSA11 =
          Y1 * Y * (4989.2 - 1607.73 * Y)
        + 3687.6  * DL
        + 3296.6  * std::pow(DL, 2)
        + 1271.11 * std::pow(DL, 3)
        + 533.44  * std::pow(DL, 4)
        + 97.27   * std::pow(DL, 5)
        + 4.0     * std::pow(DL, 6)
        + 60.40   * Y1 * std::pow(DL1, 2)
        + 4.685   * Y1 * std::pow(DL1, 3);

    double P3NSA12 =
          1030.79 * Y1 * Y
        + 1266.77 * Y1 * (2.0 - Y * Y)
        + 2987.83 * DL
        + 273.05  * std::pow(DL, 2)
        - 923.48  * std::pow(DL, 3)
        - 236.76  * std::pow(DL, 4)
        - 33.886  * std::pow(DL, 5)
        - 4.0     * std::pow(DL, 6)
        - 254.63  * Y1 * DL1
        - 0.28953 * Y1 * std::pow(DL1, 3);

    // ---------------------------------------------------------
    // nf^2 : parametrized
    // ---------------------------------------------------------

    double P3NSSA2 =
        2.5e2 * (
            Y1 * (-4.7656 + 1.6908 * Y + 0.1703 * std::pow(Y, 2))
            - 0.41652 * Y * DL
            + 0.90777 * Y * std::pow(DL, 2)
            + 0.12478 * Y * std::pow(DL, 3)
            + 0.17155 * Y1 * DL1
            + 0.17191 * DL * DL1
        )
        - 6.473971e2 * DL
        - 6.641219e1 * std::pow(DL, 2)
        - 5.353347e0 * std::pow(DL, 3)
        - 5.925926e0 * std::pow(DL, 4)
        - 3.950617e-1 * std::pow(DL, 5)
        + 1.970002e1 * Y1 * DL1
        - 3.435474e0 * Y1 * std::pow(DL1, 2);

    // ---------------------------------------------------------
    // Assembly
    // ---------------------------------------------------------

	double res{};
    if (IMOD == 1)
        res = nf * P3NSA11 + std::pow(nf, 2) * P3NSSA2;
    else if (IMOD == 2)
        res = nf * P3NSA12 + std::pow(nf, 2) * P3NSSA2;
    else
        res = 0.5 * nf * (P3NSA11 + P3NSA12)
             + std::pow(nf, 2) * P3NSSA2;
	return res/16.0;
}

double p3psa(double y, int nf, int imod) {
	double ym = 1.0 / y;
	double y1 = 1.0 - y;
	double dl = std::log(y);
	double dl1 = std::log(1.0 - y);

	double nf_d = static_cast<double>(nf);
	double nf2 = nf_d * nf_d;
	double nf3 = nf_d * nf2;

	// Known large-x coefficients
	double x1l4cff = -5.6460905e1 * nf_d + 3.6213992 * nf2;
	double x1l3cff = -2.4755054e2 * nf_d + 4.0559671e1 * nf2 - 1.5802469 * nf3;
	double y1l4cff = -1.3168724e1 * nf_d;
	double y1l3cff = -1.9911111e2 * nf_d + 1.3695473e1 * nf2;

	// Known small-x coefficients
	double bfkl1 = 1.7492273e3 * nf_d;
	double x0l6cff = -7.5061728 * nf_d + 7.9012346e-1 * nf2;
	double x0l5cff =  2.8549794e1 * nf_d + 3.7925926 * nf2;
	double x0l4cff = -8.5480010e2 * nf_d + 7.7366255e1 * nf2 - 1.9753086e-1 * nf3;

	// The resulting part of the function
	double p3ps01 = bfkl1 * std::pow(dl, 2) * ym
		+ x0l6cff * std::pow(dl, 6)
		+ x0l5cff * std::pow(dl, 5)
		+ x0l4cff * std::pow(dl, 4)
		+ x1l3cff * y1 * std::pow(dl1, 3)
		+ x1l4cff * y1 * std::pow(dl1, 4)
		+ y1l3cff * (y1 * y1) * std::pow(dl1, 3)
		+ y1l4cff * (y1 * y1) * std::pow(dl1, 4);

	double p3psapp1 = 0.0;
	double p3psapp2 = 0.0;

	if (nf == 3) {
		p3psapp1 = p3ps01
			+ 67731.0   * y1 * dl * ym
			+ 274100.0  * y1 * ym
			- 104493.0  * y1 * (1.0 + 2.0 * y)
			+ 34403.0   * y1 * (y * y)
			+ 353656.0  * y1 * dl
			+ 10620.0   * std::pow(dl, 2)
			+ 40006.0   * std::pow(dl, 3)
			- 7412.1    * y1 * dl1
			- 2365.1    * y1 * std::pow(dl1, 2)
			+ 1533.0    * (y1 * y1) * std::pow(dl1, 2);

		p3psapp2 = p3ps01
			+ 54593.0   * y1 * dl * ym
			+ 179748.0  * y1 * ym
			- 195263.0  * y1
			+ 12789.0   * y1 * y * (1.0 + y)
			+ 4700.0    * y1 * dl
			- 103604.0  * std::pow(dl, 2)
			- 2758.3    * std::pow(dl, 3)
			- 2801.2    * y1 * dl1
			- 1986.9    * y1 * std::pow(dl1, 2)
			- 6005.9    * (y1 * y1) * std::pow(dl1, 2);

	} else if (nf == 4) {
		p3psapp1 = p3ps01
			+ 90154.0   * y1 * dl * ym
			+ 359084.0  * y1 * ym
			- 136319.0  * y1 * (1.0 + 2.0 * y)
			+ 45379.0   * y1 * (y * y)
			+ 461167.0  * y1 * dl
			+ 13869.0   * std::pow(dl, 2)
			+ 52525.0   * std::pow(dl, 3)
			- 7498.2    * y1 * dl1
			- 2491.5    * y1 * std::pow(dl1, 2)
			+ 1727.2    * (y1 * y1) * std::pow(dl1, 2);

		p3psapp2 = p3ps01
			+ 72987.0   * y1 * dl * ym
			+ 235802.0  * y1 * ym
			- 254921.0  * y1
			+ 17138.0   * y1 * y * (1.0 + y)
			+ 5212.9    * y1 * dl
			- 135378.0  * std::pow(dl, 2)
			- 3350.9    * std::pow(dl, 3)
			- 1472.7    * y1 * dl1
			- 1997.2    * y1 * std::pow(dl1, 2)
			- 8123.3    * (y1 * y1) * std::pow(dl1, 2);

	} else if (nf == 5) {
		p3psapp1 = p3ps01
			+ 112481.0  * y1 * dl * ym
			+ 440555.0  * y1 * ym
			- 166581.0  * y1 * (1.0 + 2.0 * y)
			+ 56087.0   * y1 * (y * y)
			+ 562992.0  * y1 * dl
			+ 16882.0   * std::pow(dl, 2)
			+ 64577.0   * std::pow(dl, 3)
			- 6570.1    * y1 * dl1
			- 2365.7    * y1 * std::pow(dl1, 2)
			+ 1761.7    * (y1 * y1) * std::pow(dl1, 2);

		p3psapp2 = p3ps01
			+ 91468.0   * y1 * dl * ym
			+ 289658.0  * y1 * ym
			- 311749.0  * y1
			+ 21521.0   * y1 * y * (1.0 + y)
			+ 4908.9    * y1 * dl
			- 165795.0  * std::pow(dl, 2)
			- 3814.9    * std::pow(dl, 3)
			+ 804.5     * y1 * dl1
			- 1760.8    * y1 * std::pow(dl1, 2)
			- 10295.0   * (y1 * y1) * std::pow(dl1, 2);

	} else if (nf == 6) {
		p3psapp1 = p3ps01
			+ 134701.0  * y1 * dl * ym
			+ 518318.0  * y1 * ym
			- 195241.0  * y1 * (1.0 + 2.0 * y)
			+ 66517.0   * y1 * (y * y)
			+ 658832.0  * y1 * dl
			+ 19605.0   * std::pow(dl, 2)
			+ 76125.0   * std::pow(dl, 3)
			- 4734.5    * y1 * dl1
			- 2035.2    * y1 * std::pow(dl1, 2)
			+ 1633.1    * (y1 * y1) * std::pow(dl1, 2);

		p3psapp2 = p3ps01
			+ 110032.0  * y1 * dl * ym
			+ 341158.0  * y1 * ym
			- 365676.0  * y1
			+ 25934.0   * y1 * y * (1.0 + y)
			+ 3614.4    * y1 * dl
			- 194868.0  * std::pow(dl, 2)
			- 4172.2    * std::pow(dl, 3)
			+ 3924.3    * y1 * dl1
			- 1324.9    * y1 * std::pow(dl1, 2)
			- 12520.0   * (y1 * y1) * std::pow(dl1, 2);

	} else {
		std::cerr << " Error in function p3psa: choice of nf " << std::endl;
		return 0.0;
	}

	double result;
	if (imod == 1) {
		result = p3psapp1;
	} else if (imod == 2) {
		result = p3psapp2;
	} else {
		result = 0.5 * (p3psapp1 + p3psapp2);
	}

	return result / 16.0;
}
}


int main()
{
	auto compute_diff = [](double x, double y){ return std::abs((x-y)/((x+y)/2.0));};
	SplittingFunction::setN3LOApproxType(3);
	SplittingFunction::update(4, 0, 0.0);
	
	std::unique_ptr<SplittingFunction> p3gg = std::make_unique<P3gg>();
	std::unique_ptr<SplittingFunction> p3gq = std::make_unique<P3gq>();
	std::unique_ptr<SplittingFunction> p3qg = std::make_unique<P3qg>();
	std::unique_ptr<SplittingFunction> p3qq = std::make_unique<P3qq>();
	std::unique_ptr<SplittingFunction> p3nsm = std::make_unique<P3nsm>();
	std::unique_ptr<SplittingFunction> p3nsp = std::make_unique<P3nsp>();
	std::unique_ptr<SplittingFunction> p3nsv = std::make_unique<P3nsv>();

	std::unique_ptr<SplittingFunction> mvv_p3gg = std::make_unique<mvv_p3::P3gg>();
	std::unique_ptr<SplittingFunction> mvv_p3gq = std::make_unique<mvv_p3::P3gq>();
	std::unique_ptr<SplittingFunction> mvv_p3qg = std::make_unique<mvv_p3::P3qg>();
	std::unique_ptr<SplittingFunction> mvv_p3qq = std::make_unique<mvv_p3::P3qq>();
	std::unique_ptr<SplittingFunction> mvv_p3nsm = std::make_unique<mvv_p3::P3nsm>();
	std::unique_ptr<SplittingFunction> mvv_p3nsp = std::make_unique<mvv_p3::P3nsp>();
	std::unique_ptr<SplittingFunction> mvv_p3nsv = std::make_unique<mvv_p3::P3nsv>();

	int four = 4;
	int five = 5;
	int three = 3;
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
	for (double x = 1e-5; x<1.0; x+= (1.0-1e-5)/200.0) {
		double diff = compute_diff(p3gg->calcDelta(x), mvv_p3gg->calcDelta(x));
		std::cout << "x=" << x << "  ->  " << diff << '\n';
	}
}

// p3nsp->regular
//

