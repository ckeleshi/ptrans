#pragma once

//##########################################################################
//#                                                                        #
//#            CLOUDCOMPARE PLUGIN: ColorimetricSegmenter                  #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 of the License.               #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#    COPYRIGHT:	Tri-Thien TRUONG, Ronan COLLIER, Mathieu LETRONE       #
//#                                                                        #
//##########################################################################

//qCC_db
#include <ccColorTypes.h>

//! HSV color
struct Hsv
{
	//! Default constrctor
	Hsv()
		: h(0)
		, s(0)
		, v(0)
	{
	}

	//! Constrctor from a RGB color
	Hsv(const ccColor::Rgb& rgb)
	{
		float r = rgb.r / 255.0f;
		float g = rgb.g / 255.0f;
		float b = rgb.b / 255.0f;
		float maxComp = std::max(std::max(r, g), b);
		float minComp = std::min(std::min(r, g), b);
		float deltaComp = maxComp - minComp;

		h = 0;
		if (deltaComp != 0)
		{
			if (r == maxComp)
			{
				h = (g - b) / deltaComp;
			}
			else
			{
				if (g == maxComp)
				{
					h = 2 + (b - r) / deltaComp;
				}
				else
				{
					h = 4 + (r - g) / deltaComp;
				}
			}
			h *= 60;
			if (h < 0)
				h += 360;
		}

		s = (maxComp == 0 ? 0 : (deltaComp / maxComp) * 100);
		v = maxComp * 100;
	}

	// HSV components
	float h, s, v;
};
