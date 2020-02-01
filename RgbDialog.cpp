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

//qCC_gl
#include <ccGLWidget.h>
#include <ccGLWindow.h>

//semi-persistent parameters
static  ccColor::Rgb rgb(0,0,0);

RgbDialog::RgbDialog(ccPickingHub* pickingHub, QWidget* parent)
	: QDialog(parent)
	, Ui::RgbDialog()
	, m_pickingWin(0)
	, m_associatedPlane(0)
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
	//connect(pointPickingButton, &QCheckBox::toggled, this, &RgbDialog::pickPoint);
	connect(pointPickingButton, &QAbstractButton::toggled, this, &RgbDialog::pickPoint);
		
	//auto disable picking mode on quit
	/*connect(this, &QDialog::finished, [&]()
	{
		if (pointPickingButton->isChecked()) pointPickingButton->setChecked(false); }
	);*/
}

void RgbDialog::pickPoint(bool state)
{
	if (m_pickingHub)
	{
		ccLog::Print("pickingHub");
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
	}
	else if (m_pickingWin)
	{
		ccLog::Print("pickinWin");
		if (state)
		{
			m_pickingWin->setPickingMode(ccGLWindow::POINT_OR_TRIANGLE_PICKING);
			connect(m_pickingWin, &ccGLWindow::itemPicked, this, &RgbDialog::processPickedItem);
		}
		else
		{
			m_pickingWin->setPickingMode(ccGLWindow::DEFAULT_PICKING);
			disconnect(m_pickingWin, &ccGLWindow::itemPicked, this, &RgbDialog::processPickedItem);
		}
	}

	pointPickingButton->blockSignals(true);
	pointPickingButton->setChecked(state);
	pointPickingButton->blockSignals(false);
}

void RgbDialog::onItemPicked(const PickedItem& pi)
{
	if (!pi.entity)
	{
		return;
	}

	m_pickingWin = m_pickingHub->activeWindow();


	/*if (pi.entity->isKindOf(CC_TYPES::POINT_CLOUD))
	{
		rgb = pi.entity->getTempColor;
	}*/

	area_red->setValue(10);
	area_green->setValue(20);
	area_blue->setValue(30);
	ccLog::Print("Test");

	pointPickingButton->setChecked(false);
}

void RgbDialog::processPickedItem(ccHObject* entity, unsigned, int, int, const CCVector3& P, const CCVector3d& uvw)
{
	//without picking hub (ccViewer)
	if (!m_pickingWin)
	{
		assert(false);
		return;
	}

	if (!entity)
	{
		return;
	}

	pickPoint(false);
}

