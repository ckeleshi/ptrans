#include <math.h>
#include <vector>
#include "ConversionLab.h"

std::vector<int> RGB_to_lab(std::vector<int> rgb)
{
    int unsigned i = 0;
    for (double value : rgb)
    {
        value /= 255;
        value = (value > 0.04045) ? pow((value + 0.055) / 1.055, 2.4) : value /= 12.92;
        rgb[i] = value * 100;
        i++;
    }

    std::vector<double> XYZ{0, 0, 0};

    XYZ[0] = rgb[0] * 0.4124 + rgb[1] * 0.3576 + rgb[2] * 0.1805;
    XYZ[1] = rgb[0] * 0.2126 + rgb[1] * 0.7152 + rgb[2] * 0.0722;
    XYZ[2] = rgb[0] * 0.0193 + rgb[1] * 0.1192 + rgb[2] * 0.9505;

    XYZ[0] /= 95.047;
    XYZ[1] /= 100.0;
    XYZ[2] /= 108.883;

    i = 0;
    for (double value : XYZ)
    {
        value = (value > 0.008856) ? pow(value, 0.33333333333) : (7.787 * value) + (16. / 116.);
        XYZ[i] = value;
        i++;
    }

    std::vector<int> LAB{0, 0, 0};
    LAB[0] = (116 * XYZ[1]) - 16.;
    LAB[1] = 500 * (XYZ[0] - XYZ[1]);
    LAB[2] = 200 * (XYZ[1] - XYZ[2]);

    return LAB;
}