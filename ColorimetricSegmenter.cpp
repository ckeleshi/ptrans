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
#include <algorithm>

#include "ColorimetricSegmenter.h"
#include "RgbDialog.h"

#include "ccPointCloud.h"
#include "ccScalarField.h"
#include "DistanceComputationTools.h"


// Default constructor:
//	- pass the Qt resource path to the info.json file (from <yourPluginName>.qrc file) 
//  - constructor should mainly be used to initialize actions and other members
ColorimetricSegmenter::ColorimetricSegmenter( QObject *parent )
	: QObject( parent )
	, ccStdPluginInterface( ":/CC/plugin/ColorimetricSegmenter/info.json" )
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
	if (m_action_filterRgb == nullptr)
	{
		return;
	}

    if (m_action_filterRgbWithSegmentation == nullptr)
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
	m_action_filterRgb->setEnabled(!selectedEntities.empty());
    m_action_filterRgbWithSegmentation->setEnabled(!selectedEntities.empty());

}

// This method returns all the 'actions' your plugin can perform.
// getActions() will be called only once, when plugin is loaded.
QList<QAction *> ColorimetricSegmenter::getActions()
{
	// default action (if it has not been already created, this is the moment to do it)
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

    if (!m_action_filterRgbWithSegmentation)
    {
        // Here we use the default plugin name, description, and icon,
        // but each action should have its own.
        m_action_filterRgbWithSegmentation = new QAction( "Filter RGB", this);
        m_action_filterRgbWithSegmentation->setToolTip( "Filter the points on the selected cloud by RGB color" );
        m_action_filterRgbWithSegmentation->setIcon(getIcon());

        // Connect appropriate signal
        connect(m_action_filterRgbWithSegmentation, &QAction::triggered, this, &ColorimetricSegmenter::filterRgb);

        connect(m_action_filterRgbWithSegmentation, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
        connect(m_action_filterRgbWithSegmentation, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
        connect(m_action_filterRgbWithSegmentation, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));
    }

    return { m_action_filterRgb, m_action_filterRgbWithSegmentation };
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

void knn(ccPointCloud* cloud, const CCVector3* point, unsigned k, CCLib::ReferenceCloud* neighbours, unsigned thresholdDistance) {
    double maxSquareDist;
    neighbours = new CCLib::ReferenceCloud(cloud);
    unsigned count = cloud->computeOctree()->findPointNeighbourhood(point, neighbours, k, 1, maxSquareDist, thresholdDistance);
}

void knnRegions(ccPointCloud* basePointCloud, std::vector<CCLib::ReferenceCloud*>* regions, const CCLib::ReferenceCloud* region, unsigned k, std::vector<CCLib::ReferenceCloud*>* neighbours, unsigned thresholdDistance) {
    std::vector<int> distances;
    ccPointCloud* computedRegion = basePointCloud->partialClone(region);
    // compute distances
    CCLib::DistanceComputationTools::Cloud2CloudDistanceComputationParams params = CCLib::DistanceComputationTools::Cloud2CloudDistanceComputationParams();
    params.kNNForLocalModel = k;
    params.maxSearchDist = thresholdDistance;
    neighbours = new std::vector<CCLib::ReferenceCloud*>();
    for(CCLib::ReferenceCloud* r : *regions)
    {
        CCLib::DistanceComputationTools::computeCloud2CloudDistance(computedRegion, basePointCloud->partialClone(r), params);
        neighbours->push_back(r);
    }
    // sort the vectors
    std::sort(
        neighbours->begin(), neighbours->end(),
        [&](std::size_t a, std::size_t b) { return distances[a] < distances[b]; });
    while(neighbours->size() > k)
    {
        neighbours->pop_back();
    }
}

double colorimetricalDifference(ccColor::Rgb c1, ccColor::Rgb c2) {
    return sqrt(pow(c1.r-c2.r, 2) + pow(c1.g-c2.g, 2) + pow(c1.b-c2.b, 2));
}

ccColor::Rgb* meanRgb(ccPointCloud* basePointCloud, CCLib::ReferenceCloud* c)
{
    unsigned red = 0;
    unsigned green = 0;
    unsigned blue = 0;
    for (unsigned j = 0; j < c->size(); ++j)
    {
        const ccColor::Rgb& pRgb = basePointCloud->getPointColor(c->getPointGlobalIndex(j));
        red += pRgb.r;
        green += pRgb.g;
        blue += pRgb.b;
    }
    ccColor::Rgb* rgb = new ccColor::Rgb();
    rgb->r = red / c->size();
    rgb->g = green / c->size();
    rgb->b = blue / c->size();
    return rgb;
}

/**
 * @brief colorimetricalDifference compute mean coloremtrical difference between two reference clouds.
 * @param basePointCloud
 * @param c1
 * @param c2
 * @return
 */
double colorimetricalDifference(ccPointCloud* basePointCloud, CCLib::ReferenceCloud* c1, CCLib::ReferenceCloud* c2) {
    ccColor::Rgb* rgb1 = meanRgb(basePointCloud, c1);

    ccColor::Rgb* rgb2 = meanRgb(basePointCloud, c2);

    return colorimetricalDifference(*rgb1, *rgb2);
}


std::vector<CCLib::ReferenceCloud*>* regionGrowing(ccPointCloud* pointCloud, const unsigned TNN, const double TPP, const double TD)
{
    std::vector<unsigned> unlabeledPoints;
    for (unsigned j = 0; j < pointCloud->size(); ++j)
    {
        unlabeledPoints.push_back(j);
    }
    std::vector<CCLib::ReferenceCloud*>* regions = new std::vector<CCLib::ReferenceCloud*>();
    std::vector<unsigned>* points = new std::vector<unsigned>();
    // while there is any point in {P} that hasn’t been labeled
    while(unlabeledPoints.size() > 0)
    {
        // push an unlabeled point into stack Points
        points->push_back(unlabeledPoints.back());
        // initialize a new region Rc and add current point to R
        CCLib::ReferenceCloud* rc = new CCLib::ReferenceCloud(pointCloud);
        rc->addPointIndex(unlabeledPoints.back());
        // while stack Points is not empty
        while(points->size() > 0)
        {
            // pop Points’ top element Tpoint
            unsigned tPointIndex = points->back();
            points->pop_back();

            // for each point p in {KNNTNN(Tpoint)}
            CCLib::ReferenceCloud* knnResult;
            knn(pointCloud, pointCloud->getPoint(tPointIndex), TNN, knnResult, TD);
            for(int i=0; i < knnResult->size(); i++)
            {
                unsigned p = knnResult->getPointGlobalIndex(i);
                // if p is labelled
                if(std::find(unlabeledPoints.begin(), unlabeledPoints.end(), p) != unlabeledPoints.end())
                {
                    continue;
                }
                // if CD(Tpoint,p)<TPP
                if(colorimetricalDifference(pointCloud->getPointColor(p), pointCloud->getPointColor(tPointIndex)) < TPP)
                {
                    points->push_back(p);
                    rc->addPointIndex(p);
                }
            }

        }
        regions->push_back(rc);
    }
    return regions;
}

std::vector<CCLib::ReferenceCloud*>* findRegion(std::vector<std::vector<CCLib::ReferenceCloud*>*>* container, CCLib::ReferenceCloud* region)
{
    if(container->size() == 0)
    {
        return nullptr;
    }
    for(std::vector<CCLib::ReferenceCloud*>* l : *container)
    {
        if(std::find(l->begin(), l->end(), region) != l->end())
        {
           return l;
        }
    }
    return nullptr;
}

std::vector<CCLib::ReferenceCloud*>* regionMergingAndRefinement(ccPointCloud* basePointCloud, std::vector<CCLib::ReferenceCloud*>* regions, const unsigned TNN, const double TRR, const double TD, const unsigned Min)
{
    std::vector<std::vector<CCLib::ReferenceCloud*>*>* homogeneous = new std::vector<std::vector<CCLib::ReferenceCloud*>*>();

    // for each region Ri in {R}
    for(CCLib::ReferenceCloud* ri : *regions)
    {
        // if Ri is not in {H}
       if(findRegion(homogeneous, ri) == nullptr)
        {
           // create a new list to record Ri
           std::vector<CCLib::ReferenceCloud*>* newRegionGroup = new std::vector<CCLib::ReferenceCloud*>();
           newRegionGroup->push_back(ri);
           homogeneous->push_back(newRegionGroup);
        }

        // for each region Rj in {KNNTNN2,TD2(Ri)}
        std::vector<CCLib::ReferenceCloud*>* knnResult;
        knnRegions(basePointCloud, regions, ri, TNN, knnResult, TD);
        for(CCLib::ReferenceCloud* rj : *knnResult)
        {
            // if CD(Ri,Rj)<TRR
            if(colorimetricalDifference(basePointCloud, ri, rj) < TNN)
            {
                // if Rj is in {H}
                std::vector<CCLib::ReferenceCloud*>* regionContainer = findRegion(homogeneous, rj);
                if(regionContainer != nullptr)
                {
                    continue;
                }
                else
                {
                    // add Rj to the list which contains Ri
                    regionContainer->push_back(rj);
                }
            }
        }
    }

    // merge all the regions in the same list in {H} and get {R’}
    std::vector<CCLib::ReferenceCloud*>* mergedRegionsRef = new std::vector<CCLib::ReferenceCloud*>();
    for(std::vector<CCLib::ReferenceCloud*>* l : *homogeneous)
    {
        CCLib::ReferenceCloud* merged = l->at(0);
        for(int i = 1; i<l->size(); i++)
        {
           merged->add(*l->at(i));
        }
        mergedRegionsRef->push_back(merged);
    }
    std::vector<CCLib::ReferenceCloud*>* knnResult;
    // for each region Ri in {R’}
    /*for (CCLib::ReferenceCloud* r : *mergedRegionsRef)
    {
        // if sizeof(Ri)<Min
        if(r->size() < Min)
        {
            // merge Ri to its nearest neighbors
            knnRegions(basePointCloud, mergedRegionsRef, r, 1, knnResult, 0);
            mergedRegionsRef.
        }
    }*/
    //Return the merged and refined {R’}
    return mergedRegionsRef;
}


void ColorimetricSegmenter::filterRgbWithSegmentation()
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

            std::vector<CCLib::ReferenceCloud*>* regions = regionGrowing(cloud, TNN, TPP, TD);
            regions = regionMergingAndRefinement(cloud, regions, TNN, TRR, TD, Min);

            // retrieve the nearest region (in color range)
            for(CCLib::ReferenceCloud* r : *regions)
            {
                ccColor::Rgb* mean = meanRgb(cloud, r);
                if(redInf >= mean->r && redSup <= mean->r &&
                        greenInf >= mean->g && greenSup <= mean->g &&
                        blueInf >= mean->b && blueSup <= mean->b)
                {
                    ccPointCloud* newCloud = cloud->partialClone(r);

                    cloud->setEnabled(false);
                    if (cloud->getParent()) {
                        cloud->getParent()->addChild(newCloud);
                    }

                    m_app->addToDB(newCloud, false, true, false, false);

                    m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
                }

            }


        }
    }
}

