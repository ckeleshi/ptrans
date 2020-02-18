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

#include <QtGui>

#include "ColorimetricSegmenter.h"
#include "RgbDialog.h"

#include "ccPointCloud.h"
#include "ccScalarField.h"


// Default constructor:
//	- pass the Qt resource path to the info.json file (from <yourPluginName>.qrc file) 
//  - constructor should mainly be used to initialize actions and other members
ColorimetricSegmenter::ColorimetricSegmenter( QObject *parent )
	: QObject( parent )
	, ccStdPluginInterface( ":/CC/plugin/ColorimetricSegmenter/info.json" )
	, m_action_filterScalar(nullptr )
	, m_action_filterRgb( nullptr )
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
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action_filterRgb = new QAction( "Filter RGB", this);
		m_action_filterRgb->setToolTip( "Filter the points on the selected cloud by RGB color" );
		m_action_filterRgb->setIcon(getIcon());

		// Connect appropriate signal
		connect(m_action_filterRgb, &QAction::triggered, this, &ColorimetricSegmenter::filterRgb);

		connect(m_action_filterRgb, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));
	}

	return { m_action_filterScalar, m_action_filterRgb };
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

	// Retrieve parameters from dialog
	RgbDialog rgbDlg((QWidget*)m_app->getMainWindow());

	if (!rgbDlg.exec())
		return;

	double marginError = static_cast<double>(rgbDlg.margin->value()) / 100.0;

	int redInf = rgbDlg.area_red->value() - (marginError * rgbDlg.area_red->value());
	int redSup = rgbDlg.area_red->value() + marginError * rgbDlg.area_red->value();
	int greenInf = rgbDlg.area_green->value() - marginError * rgbDlg.area_green->value();
	int greenSup = rgbDlg.area_green->value() + marginError * rgbDlg.area_green->value();
	int blueInf = rgbDlg.area_blue->value() - marginError * rgbDlg.area_blue->value();
	int blueSup = rgbDlg.area_blue->value() + marginError * rgbDlg.area_blue->value();

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
			
			cloud->setEnabled(false);
			if (cloud->getParent()) {
				cloud->getParent()->addChild(newCloud);
			}

			m_app->addToDB(newCloud, false, true, false, false);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
			

		}
	} 
}

void knn(std::vector<ccPointCloud*> clouds, CCVector3& point, int k, CCLib::ReferenceCloud points) {
    for (ccPointCloud* cloud : clouds) {
        ccOctree::Shared octree = cloud->computeOctree()->findPointNeighbourhood(point, points)
    }
}
