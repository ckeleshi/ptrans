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

#include <iostream>

#include <QtGui>

#include "ColorimetricSegmenter.h"

/* 
	Default constructor:
	- pass the Qt resource path to the info.json file (from <yourPluginName>.qrc file) 
	- constructor should mainly be used to initialize actions and other members
*/
/*
	Class that will apply all the algorithm of our plugin
*/
ColorimetricSegmenter::ColorimetricSegmenter( QObject *parent )
	: QObject( parent )
	, ccStdPluginInterface( ":/CC/plugin/ColorimetricSegmenter/info.json" )
	, m_action_filterRgb( nullptr )
	, m_action_filterHSV( nullptr )
{
}

void ColorimetricSegmenter::handleNewEntity(ccHObject* entity)
{
	assert(entity && m_app);
	m_app->addToDB(entity);
}

void ColorimetricSegmenter::handleEntityChange(ccHObject* entity)
{
	assert(entity && m_app);
	entity->prepareDisplayForRefresh_recursive();
	m_app->refreshAll();
	m_app->updateUI();
}

void ColorimetricSegmenter::handleErrorMessage(QString message)
{
	if (m_app)
		m_app->dispToConsole(message, ccMainAppInterface::ERR_CONSOLE_MESSAGE);
}

// This method should enable or disable your plugin actions
// depending on the currently selected entities ('selectedEntities').
void ColorimetricSegmenter::onNewSelection( const ccHObject::Container &selectedEntities )
{
	if (m_action_filterRgb == nullptr)
	{
		return;
	}
	if (m_action_filterHSV == nullptr)
	{
		return;
	}
	
	// Only enable our action if something is selected.
	bool activateColorFilters = false;
	bool activateScalarFilter = false;
	for (ccHObject *entity : selectedEntities)
	{
		if (entity->isKindOf(CC_TYPES::POINT_CLOUD))
		{
			if (entity->hasColors()) {
				activateColorFilters = true;
			}
			else if (entity->hasDisplayedScalarField()) {
				activateScalarFilter = true;
			}
		}
	}

	m_action_filterRgb->setEnabled(false);
	m_action_filterHSV->setEnabled(false);

	//Activate only if only one of them is activated
	if ((activateColorFilters != activateScalarFilter) && !selectedEntities.empty()) {
		m_action_filterRgb->setEnabled(activateColorFilters);
		m_action_filterHSV->setEnabled(activateColorFilters);
	}

}

// This method returns all the 'actions' your plugin can perform.
// getActions() will be called only once, when plugin is loaded.
QList<QAction*> ColorimetricSegmenter::getActions()
{

	// RGB Filter
	if (!m_action_filterRgb)
	{
		m_action_filterRgb = new QAction("Filter RGB", this);
		m_action_filterRgb->setToolTip("Filter the points on the selected cloud by RGB color");
		m_action_filterRgb->setIcon(getIcon());

		// Connect appropriate signal
		connect(m_action_filterRgb, &QAction::triggered, this, &ColorimetricSegmenter::filterRgb);

		connect(m_action_filterRgb, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}

	// HSV Filter
	if (!m_action_filterHSV)
	{
		m_action_filterHSV = new QAction("Filter HSV", this);
		m_action_filterHSV->setToolTip("Filter the points on the selected cloud by HSV color");
		m_action_filterHSV->setIcon(getIcon());

		// Connect appropriate signal
		connect(m_action_filterHSV, &QAction::triggered, this, &ColorimetricSegmenter::filterHSV);

		connect(m_action_filterHSV, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterHSV, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterHSV, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}

	return { m_action_filterRgb, m_action_filterHSV };
}

// Get all point clouds that are selected in CC
// return a vector with ccPointCloud objects
std::vector<ccPointCloud*> ColorimetricSegmenter::getSelectedPointClouds()
{
	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);
		return std::vector<ccPointCloud*>{};
	}

	ccHObject::Container selectedEntities = m_app->getSelectedEntities();
	std::vector<ccPointCloud*> clouds;
	for (size_t i = 0; i < selectedEntities.size(); ++i)
	{
		if (selectedEntities[i]->isKindOf(CC_TYPES::POINT_CLOUD)) {
			clouds.push_back(static_cast<ccPointCloud*> (selectedEntities[i]));
		}
	}

	return clouds;
}


// Algorithm for the RGB filter
// It uses a color range with RGB values, and keeps the points with a color within that range.
void ColorimetricSegmenter::filterRgb()
{	
	if ( m_app == nullptr )
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT( false );
		
		return;
	}

	//check valid window
	if (!m_app->getActiveGLWindow())
	{
		m_app->dispToConsole("[ccCompass] Could not find valid 3D window.", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	// Retrieve parameters from dialog
	if (m_app->pickingHub()) {
		m_pickingHub = m_app->pickingHub();
	}
	
	rgbDlg = new RgbDialog(m_pickingHub,(QWidget*)m_app->getMainWindow());
	rgbDlg->show();

	if (!rgbDlg->exec())
		return;

	// Start timer
	auto start = std::chrono::high_resolution_clock::now();

	// Get margin value (percent)
	double marginError = static_cast<double>(rgbDlg->margin->value()) / 100.0;

	// Get all values to make the color range with RGB values
	int redInf = rgbDlg->red_first->value() - (marginError * rgbDlg->red_first->value());
	int redSup = rgbDlg->red_second->value() + marginError * rgbDlg->red_second->value();
	int greenInf = rgbDlg->green_first->value() - marginError * rgbDlg->green_first->value();
	int greenSup = rgbDlg->green_second->value() + marginError * rgbDlg->green_second->value();
	int blueInf = rgbDlg->blue_first->value() - marginError * rgbDlg->blue_first->value();
	int blueSup = rgbDlg->blue_second->value() + marginError * rgbDlg->blue_second->value();

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds) {
		if (cloud->hasColors())
		{
			// Use only references for speed reasons
			CCLib::ReferenceCloud* filteredCloud = new CCLib::ReferenceCloud(cloud);

			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgb& rgb = cloud->getPointColor(j);

				if (rgb.r > redInf && rgb.r < redSup &&
					rgb.g > greenInf && rgb.g < greenSup &&
					rgb.b > blueInf && rgb.b < blueSup) { // Points to select
					if (!filteredCloud->addPointIndex(j))
					{
						//not enough memory
						delete filteredCloud;
						filteredCloud = nullptr;
						m_app->dispToConsole("[ColorimetricSegmenter] Error, filter canceled.");
						break;
					}
				}

			}
			ccPointCloud* newCloud = cloud->partialClone(filteredCloud);
			newCloud->setName(QString::fromStdString("Rmin:" + std::to_string(redInf) + "/Gmin:" + std::to_string(greenInf) + "/Bmin:" + std::to_string(blueInf) +
													"/Rmax:" + std::to_string(redSup) + "/Gmax:" + std::to_string(greenSup) + "/Bmax:" + std::to_string(blueSup)));
			
			cloud->setEnabled(false);
			if (cloud->getParent()) {
				cloud->getParent()->addChild(newCloud);
			}

			m_app->addToDB(newCloud, false, true, false, false);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
			

		}
	}

	// Stop timer
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	QString s = QString::number(duration);

	//Print time of execution
	ccLog::Print("Time to execute : " + s + " milliseconds.");
}

// Algorithm for the HSV filter
// It uses the Hue-Saturation-Value (HSV) color space to filter the point cloud
void ColorimetricSegmenter::filterHSV()
{
	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);

		return;
	}

	//check valid window
	if (!m_app->getActiveGLWindow())
	{
		m_app->dispToConsole("[ccCompass] Could not find valid 3D window.", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	// Retrieve parameters from dialog
	if (m_app->pickingHub()) {
		m_pickingHub = m_app->pickingHub();
	}

	hsvDlg = new HSVDialog(m_pickingHub, (QWidget*)m_app->getMainWindow());
	hsvDlg->show();

	if (!hsvDlg->exec())
		return;

	// Start timer
	auto start = std::chrono::high_resolution_clock::now();

	// Get HSV values
	hsv hsv_first;
	hsv_first.h = hsvDlg->hue_first->value();
	hsv_first.s = hsvDlg->sat_first->value();
	hsv_first.v = hsvDlg->val_first->value();

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds) {
		if (cloud->hasColors()) {
			// Use only references for speed reasons
			CCLib::ReferenceCloud* filteredCloud = new CCLib::ReferenceCloud(cloud);

			// We manually add color ranges with HSV values
			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgb& rgb = cloud->getPointColor(j);
				hsv hsv_current = hsvDlg->rgb2hsv(rgb);
				
				// Hue is useless here because the saturation is not high enough
				if ( 0 <= hsv_first.s && hsv_first.s <= 25 && 0 <= hsv_current.s && hsv_current.s <= 25)
				{
					// We only check value
					if      (hsv_first.v >= 0 && hsv_first.v <= 25 && 0 <= hsv_current.v && hsv_current.v <= 25)   addPoint(filteredCloud, j); // black
					else if (hsv_first.v > 25 && hsv_first.v <= 60 && hsv_current.v > 25 && hsv_current.v <= 60)   addPoint(filteredCloud, j); // grey
					else if (hsv_first.v > 60 && hsv_first.v <= 100 && hsv_current.v > 60 && hsv_current.v <= 100) addPoint(filteredCloud, j); // white
				}
				else if (hsv_first.s > 25 && hsv_first.s <= 100 && hsv_current.s > 25 && hsv_current.s <= 100) {
					// We need to check value first
					if (0 <= hsv_first.v && hsv_first.v <= 25 && 0 <= hsv_current.v && hsv_current.v <= 25) addPoint(filteredCloud, j); // black
					// Then, we can check value
					else if (hsv_first.v > 25 && hsv_first.v <= 100 && hsv_current.v > 25 && hsv_current.v <= 100) 
					{
						if     (((hsv_first.h >= 0 && hsv_first.h <= 30) || (hsv_first.h >= 330 && hsv_first.h <= 360)) &&
							    ((hsv_current.h >= 0 && hsv_current.h <= 30) || (hsv_current.h >= 330 && hsv_current.h <= 360))) addPoint(filteredCloud, j); // red
						else if (hsv_first.h > 30 && hsv_first.h <= 90 && hsv_current.h > 30 && hsv_current.h <= 90)             addPoint(filteredCloud, j); // yellow
						else if (hsv_first.h > 90 && hsv_first.h <= 150 && hsv_current.h > 90 && hsv_current.h <= 150)           addPoint(filteredCloud, j); // green
						else if (hsv_first.h > 150 && hsv_first.h <= 210 && hsv_current.h > 150 && hsv_current.h <= 210)         addPoint(filteredCloud, j); // cyan
						else if (hsv_first.h > 210 && hsv_first.h <= 270 && hsv_current.h > 210 && hsv_current.h <= 270)         addPoint(filteredCloud, j); // blue
						else if (hsv_first.h > 270 && hsv_first.h <= 330 && hsv_current.h > 270 && hsv_current.h <= 330)         addPoint(filteredCloud, j); // magenta
					}
				}

			}

			// Copy new cloud
			ccPointCloud* newCloud = cloud->partialClone(filteredCloud);
			newCloud->setName(QString::fromStdString("h:" + std::to_string((int)hsv_first.h) + "/s:" + std::to_string((int)hsv_first.s) + "/v:" + std::to_string((int)hsv_first.v)));

			cloud->setEnabled(false);
			if (cloud->getParent()) {
				cloud->getParent()->addChild(newCloud);
			}

			m_app->addToDB(newCloud, false, true, false, false);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);


		}
	}

	// Stop timer
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	QString s = QString::number(duration);

	//Print time of execution
	ccLog::Print("Time to execute : " + s + " milliseconds");
}

void ColorimetricSegmenter::addPoint(CCLib::ReferenceCloud* filteredCloud, unsigned int j)
{
	if (!filteredCloud->addPointIndex(j))
	{
		//not enough memory
		delete filteredCloud;
		filteredCloud = nullptr;
		m_app->dispToConsole("[ColorimetricSegmenter] Error, filter canceled.");
	}
}