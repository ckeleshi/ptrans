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
#include "HSVDialog.h"

//local
#include "mainwindow.h"

//Qt
#include <QVariant>
#include <QApplication>
#include <QClipboard>
#include <QFileDialog>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>

//common
#include <ccPickingHub.h>

#include <ccGenericPointCloud.h>
#include <ccPointCloud.h>

//qCC_gl
#include <ccGLWidget.h>
#include <ccGLWindow.h>

/*
	Constructor
*/
HSVDialog::HSVDialog(ccPickingHub* pickingHub, QWidget* parent)
	: QDialog(parent)
	, Ui::HSVDialog()
	, m_pickingWin(0)
	, m_pickingHub(pickingHub)
{
	assert(pickingHub);

	setModal(false);
	setupUi(this);

	//restore semi-persistent parameters
	red->setValue(0);
	green->setValue(0);
	blue->setValue(0);

	//Link between Ui and actions
	connect(pointPickingButton_first, &QCheckBox::toggled, this, &HSVDialog::pickPoint);
	connect(red, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &HSVDialog::updateValues);
	connect(green, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &HSVDialog::updateValues);
	connect(blue, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &HSVDialog::updateValues);
		
	//auto disable picking mode on quit
	connect(this, &QDialog::finished, [&]()
	{
		if (pointPickingButton_first->isChecked()) pointPickingButton_first->setChecked(false); 
	}
	);
}

/*
	Method for the picking point functionnality
*/
void HSVDialog::pickPoint(bool state)
{
	if (!m_pickingHub)
	{
		return;
	}
	if (state)
	{
		if (!m_pickingHub->addListener(this, true))
		{
			ccLog::Error("Can't start the picking process (another tool is using it)");
			state = false;
		}
	}
	else
	{
		m_pickingHub->removeListener(this);
	}
		pointPickingButton_first->blockSignals(true);
		pointPickingButton_first->setChecked(state);
		pointPickingButton_first->blockSignals(false);
}

/*
	Method applied after a point is picked by picking point functionnality
*/
void HSVDialog::onItemPicked(const PickedItem& pi)
{
	assert(pi.entity);
	m_pickingWin = m_pickingHub->activeWindow();

	if (pi.entity->isKindOf(CC_TYPES::POINT_CLOUD))
	{
		//Get RGB values of the picked point
		ccGenericPointCloud* cloud = static_cast<ccGenericPointCloud*>(pi.entity);
		const ccColor::Rgb& rgb = cloud->getPointColor(pi.itemIndex);
		if (pointPickingButton_first->isChecked()) {
			ccLog::Print("Point picked");

			//blocking signals to avoid updating 2 times hsv values for nothing
			red->blockSignals(true);
			green->blockSignals(true);

			red->setValue(rgb.r);
			green->setValue(rgb.g);
			blue->setValue(rgb.b);

			red->blockSignals(false);
			green->blockSignals(false);

			pointPickingButton_first->setChecked(false);
		}
	}
	
}

/*
	Method applied after entering a value in RGB text fields
*/
void HSVDialog::updateValues()
{
	ccColor::Rgb rgb(red->value(), green->value(), blue->value());

	hsv hsv_values = rgb2hsv(rgb);
	hue_first->setValue(hsv_values.h);
	sat_first->setValue(hsv_values.s);
	val_first->setValue(hsv_values.v);
}

/*
	Method to convert from rgb values to hsv values
	return : hsv struct
*/
hsv HSVDialog::rgb2hsv(ccColor::Rgb rgb) {
	hsv res;
	float r = rgb.r / 255.0f;
	float g = rgb.g / 255.0f;
	float b = rgb.b / 255.0f;
	float max = std::max(std::max(r, g), b);
	float min = std::min(std::min(r, g), b);
	float delta = max - min;

	res.v = max;
	if (delta != 0) {
		float hue;
		if (r == max) {
			hue = (g - b) / delta;
		}
		else {
			if (g == max) {
				hue = 2 + (b - r) / delta;
			}
			else {
				hue = 4 + (r - g) / delta;
			}
		}
		hue *= 60;
		if (hue < 0) hue += 360;
		res.h = hue;
	}
	else {
		res.h = 0;
	}
	res.s = max == 0 ? 0 : ((max - min) / max) * 100;
	res.v = max * 100;
	
	return res;
}