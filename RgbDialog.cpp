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
#include "RgbDialog.h"

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

//semi-persistent parameters
static  ccColor::Rgb rgb(0,0,0);

RgbDialog::RgbDialog(ccPickingHub* pickingHub, QWidget* parent)
	: QDialog(parent)
	, Ui::RgbDialog()
	, m_pickingWin(0)
	, m_pickingHub(pickingHub)
{
	assert(pickingHub);

	setModal(false);
	setupUi(this);


	//restore semi-persistent parameters
	area_red->setValue(rgb.r);
	area_green->setValue(rgb.g);
	area_blue->setValue(rgb.b);

	//Link between Ui and actions
	connect(pointPickingButton, &QCheckBox::toggled, this, &RgbDialog::pickPoint);
		
	//auto disable picking mode on quit
	connect(this, &QDialog::finished, [&]()
	{
		if (pointPickingButton->isChecked()) pointPickingButton->setChecked(false); }
	);
}

void RgbDialog::pickPoint(bool state)
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
	pointPickingButton->blockSignals(true);
	pointPickingButton->setChecked(state);
	pointPickingButton->blockSignals(false);
}


void RgbDialog::onItemPicked(const PickedItem& pi)
{
	ccLog::Print("Point Picked");
	assert(pi.entity);
	m_pickingWin = m_pickingHub->activeWindow();

	if (pi.entity->isKindOf(CC_TYPES::POINT_CLOUD))
	{
		//Get RGB values of the picked point
		ccGenericPointCloud* cloud = static_cast<ccGenericPointCloud*>(pi.entity);
		const ccColor::Rgb& rgb = cloud->getPointColor(pi.itemIndex);

		area_red->setValue(rgb.r);
		area_green->setValue(rgb.g);
		area_blue->setValue(rgb.b);
	}
	pointPickingButton->setChecked(false);
}