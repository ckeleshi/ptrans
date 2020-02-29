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

// First:
//	Replace all occurrences of 'ExamplePlugin' by your own plugin class name in this file.
//	This includes the resource path to info.json in the constructor.

// Second:
//	Open ExamplePlugin.qrc, change the "prefix" and the icon filename for your plugin.
//	Change the name of the file to <yourPluginName>.qrc

// Third:
//	Open the info.json file and fill in the information about the plugin.
//	 "type" should be one of: "Standard", "GL", or "I/O" (required)
//	 "name" is the name of the plugin (required)
//	 "icon" is the Qt resource path to the plugin's icon (from the .qrc file)
//	 "description" is used as a tootip if the plugin has actions and is displayed in the plugin dialog
//	 "authors", "maintainers", and "references" show up in the plugin dialog as well

#include <iostream>

#include <QtGui>

#include "ColorimetricSegmenter.h"



// Default constructor:
//	- pass the Qt resource path to the info.json file (from <yourPluginName>.qrc file) 
//  - constructor should mainly be used to initialize actions and other members
ColorimetricSegmenter::ColorimetricSegmenter( QObject *parent )
	: QObject( parent )
	, ccStdPluginInterface( ":/CC/plugin/ColorimetricSegmenter/info.json" )
	, m_action_filterScalar(nullptr )
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
	if (m_action_filterScalar == nullptr )
	{
		return;
	}
	if (m_action_filterRgb == nullptr)
	{
		return;
	}
	if (m_action_filterHSV == nullptr)
	{
		return;
	}
	
	// If you need to check for a specific type of object, you can use the methods
	// in ccHObjectCaster.h or loop and check the objects' classIDs like this:
	//
	//	for ( ccHObject *object : selectedEntities )
	//	{
	//		if ( object->getClassID() == CC_TYPES::VIEWPORT_2D_OBJECT )
	//		{
	//			// ... do something with the viewports
	//		}
	//	}
	
	// For example - only enable our action if something is selected.
	m_action_filterScalar->setEnabled( !selectedEntities.empty() );
	m_action_filterRgb->setEnabled(!selectedEntities.empty());
	m_action_filterHSV->setEnabled(!selectedEntities.empty());
}

// This method returns all the 'actions' your plugin can perform.
// getActions() will be called only once, when plugin is loaded.
QList<QAction *> ColorimetricSegmenter::getActions()
{

	// default action (if it has not been already created, this is the moment to do it)
	if ( !m_action_filterScalar )
	{
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action_filterScalar = new QAction( "Filter with scalar field", this );
		m_action_filterScalar->setToolTip( "Filter the points on the selected cloud by scalar value" );
		m_action_filterScalar->setIcon( getIcon() );
		
		// Connect appropriate signal
		connect(m_action_filterScalar, &QAction::triggered, this, &ColorimetricSegmenter::filterScalarField );

		connect(m_action_filterScalar, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterScalar, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterScalar, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));
	}

	if (!m_action_filterRgb)
	{
		m_action_filterRgb = new QAction( "Filter RGB", this);
		m_action_filterRgb->setToolTip( "Filter the points on the selected cloud by RGB color" );
		m_action_filterRgb->setIcon(getIcon());

		// Connect appropriate signal
		connect(m_action_filterRgb, &QAction::triggered, this, &ColorimetricSegmenter::filterRgb);

		connect(m_action_filterRgb, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}

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

	return { m_action_filterScalar, m_action_filterRgb, m_action_filterHSV };
}


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

void ColorimetricSegmenter::filterScalarField()
{
	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);
		return;
	}
	int minVal;
	int maxVal;

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();
	for (ccPointCloud* cloud : clouds) {
		ccPointCloud* filteredCloud = cloud->filterPointsByScalarValue(minVal, maxVal, false);
		cloud->getParent()->addChild(filteredCloud);
		cloud->setRGBColorWithCurrentScalarField();
		m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);

		emit newEntity(filteredCloud);
		emit entityHasChanged(cloud);
	}
}


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

	double marginError = static_cast<double>(rgbDlg->margin->value()) / 100.0;

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
	//std::cout << duration.count() << std::endl;
	ccLog::Print("Time to execute : " + s);
}

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

	hsv hsv_first;
	hsv_first.h = hsvDlg->hue_first->value();
	hsv_first.s = hsvDlg->sat_first->value();
	hsv_first.v = hsvDlg->val_first->value();

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds) {
		if (cloud->hasColors())
		{
			// Use only references for speed reasons
			CCLib::ReferenceCloud* filteredCloud = new CCLib::ReferenceCloud(cloud);

			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgb& rgb = cloud->getPointColor(j);
				hsv hsv_current = hsvDlg->rgb2hsv(rgb);
				
				if ( 0 < hsv_first.s && hsv_first.s <= 25 && 0 < hsv_current.s && hsv_current.s <= 25) // hue is useless here, so we will check value
				{
					if (0 < hsv_first.v && hsv_first.v <= 25 && 0 < hsv_current.v && hsv_current.v <= 25) { // black
						addPoint(filteredCloud, j);
					}
					else if (hsv_first.v > 25 && hsv_first.v <= 75 && hsv_current.v > 25 && hsv_current.v <= 75) { // grey
						addPoint(filteredCloud, j);
					}
					else if (hsv_first.v > 75 && hsv_first.v <= 100 && hsv_current.v > 75 && hsv_current.v <= 100) { // white
						addPoint(filteredCloud, j);
					}
				}
				else if (hsv_first.s > 25 && hsv_first.s <= 100 && hsv_current.s > 25 && hsv_current.s <= 100) { // we need to check value, then hue
					if (0 < hsv_first.v && hsv_first.v <= 25 && 0 < hsv_current.v && hsv_current.v <= 25) { // black
						addPoint(filteredCloud, j);
					}
					else if (hsv_first.v > 50 && hsv_first.v <= 100 && hsv_current.v > 50 && hsv_current.v <= 100) { // particular color
						if ((hsv_first.h >= 0 && hsv_first.h <= 30) || (hsv_first.h >= 330 && hsv_first.h <= 360) && (hsv_current.h >= 0 && hsv_current.h <= 30) || (hsv_current.h >= 330 && hsv_current.h <= 360)) // red
						{
							addPoint(filteredCloud, j);
						}
						else if (hsv_first.h > 30 && hsv_first.h <= 90 && hsv_current.h > 30 && hsv_current.h <= 90) // yellow
						{
							addPoint(filteredCloud, j);
						}
						else if (hsv_first.h > 90 && hsv_first.h <= 150 && hsv_current.h > 90 && hsv_current.h <= 150) // green
						{
							addPoint(filteredCloud, j);
						}
						else if (hsv_first.h > 150 && hsv_first.h <= 210 && hsv_current.h > 150 && hsv_current.h <= 210) // cyan
						{
							addPoint(filteredCloud, j);
						}
						else if (hsv_first.h > 210 && hsv_first.h <= 270 && hsv_current.h > 210 && hsv_current.h <= 270) // blue
						{
							addPoint(filteredCloud, j);
						}
						else if (hsv_first.h > 270 && hsv_first.h <= 330 && hsv_current.h > 270 && hsv_current.h <= 330) // magenta
						{
							addPoint(filteredCloud, j);
						}
					}
				}

			}
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
	//std::cout << duration.count() << std::endl;
	ccLog::Print("Time to execute : " + s);
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