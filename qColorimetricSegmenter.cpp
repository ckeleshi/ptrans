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

//Qt
#include <QtGui>

//System
#include <iostream>
#include <algorithm>
#include <map>

//CloudCompare
#include <ccLog.h>
#include <ccPointCloud.h>
#include <ccScalarField.h>

//CCCoreLib
#include <DistanceComputationTools.h>

#include "qColorimetricSegmenter.h"

ColorimetricSegmenter::ColorimetricSegmenter(QObject* parent)
	: QObject(parent)
	, ccStdPluginInterface(":/CC/plugin/ColorimetricSegmenter/info.json")
	, m_action_filterRgb(nullptr)
	/*, m_action_filterRgbWithSegmentation(nullptr)*/
	, m_action_filterHSV(nullptr)
	, m_action_filterScalar(nullptr)
	, m_action_ToonMapping_Hist(nullptr)
	, m_action_ToonMapping_KMeans(nullptr)
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
void ColorimetricSegmenter::onNewSelection(const ccHObject::Container& selectedEntities)
{
	if (m_action_filterRgb == nullptr)
	{
		return;
	}

	if (m_action_filterHSV == nullptr)
	{
		return;
	}

	/*if (m_action_filterRgbWithSegmentation == nullptr)
	{
		return;
	}*/

	if (m_action_filterScalar == nullptr)
	{
		return;
	}
	if (m_action_ToonMapping_Hist == nullptr)
	{
		return;
	}
	if (m_action_ToonMapping_KMeans == nullptr)
	{
		return;
	}


	// For example - only enable our action if something is selected.
	// Only enable our action if something is selected.
	bool activateColorFilters = false;
	bool activateScalarFilter = false;
	for (ccHObject* entity : selectedEntities)
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
	//m_action_filterRgbWithSegmentation->setEnabled(false);
	m_action_filterScalar->setEnabled(false);
	m_action_ToonMapping_Hist->setEnabled(false);
	m_action_ToonMapping_KMeans->setEnabled(false);


	//Activate only if only one of them is activated
	if ((activateColorFilters != activateScalarFilter) && !selectedEntities.empty()) {
		m_action_filterRgb->setEnabled(activateColorFilters);
		m_action_filterHSV->setEnabled(activateColorFilters);
		//m_action_filterRgbWithSegmentation->setEnabled(activateColorFilters);
		m_action_filterScalar->setEnabled(activateScalarFilter);
		m_action_ToonMapping_Hist->setEnabled(activateColorFilters);
		m_action_ToonMapping_KMeans->setEnabled(activateColorFilters);

	}

}

QList<QAction*> ColorimetricSegmenter::getActions()
{
	// RGB Filter
	if (!m_action_filterRgb)
	{
		m_action_filterRgb = new QAction("Filter RGB", this);
		m_action_filterRgb->setToolTip("Filter the points on the selected cloud by RGB color");
        m_action_filterRgb->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_rgb.png"));

		// Connect appropriate signal
		connect(m_action_filterRgb, &QAction::triggered, this, &ColorimetricSegmenter::filterRgb);

		connect(m_action_filterRgb, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterRgb, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}

	
	/*if (!m_action_filterRgbWithSegmentation)
	{
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action_filterRgbWithSegmentation = new QAction("Filter RGB using segmentation", this);
		m_action_filterRgbWithSegmentation->setToolTip("Filter the points on the selected cloud by RGB color using segmentation");
		m_action_filterRgbWithSegmentation->setIcon(getIcon());

		// Connect appropriate signal
		connect(m_action_filterRgbWithSegmentation, &QAction::triggered, this, &ColorimetricSegmenter::filterRgbWithSegmentation);

		connect(m_action_filterRgbWithSegmentation, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterRgbWithSegmentation, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterRgbWithSegmentation, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));
	}*/

	// HSV Filter
	if (!m_action_filterHSV)
	{
		m_action_filterHSV = new QAction("Filter HSV", this);
		m_action_filterHSV->setToolTip("Filter the points on the selected cloud by HSV color");
        m_action_filterHSV->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_hsv.png"));

		// Connect appropriate signal
		connect(m_action_filterHSV, &QAction::triggered, this, &ColorimetricSegmenter::filterHSV);

		connect(m_action_filterHSV, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterHSV, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterHSV, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}

	// Scalar filter
	if (!m_action_filterScalar)
	{
		m_action_filterScalar = new QAction("Filter scalar", this);
		m_action_filterScalar->setToolTip("Filter the points on the selected cloud using scalar field");
        m_action_filterScalar->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_scalar.png"));

		// Connect appropriate signal
		connect(m_action_filterScalar, &QAction::triggered, this, &ColorimetricSegmenter::filterScalar);

		connect(m_action_filterScalar, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterScalar, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterScalar, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}
	if (!m_action_ToonMapping_Hist)
	{
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action_ToonMapping_Hist = new QAction("Histogram Clustering", this);
        m_action_ToonMapping_Hist->setToolTip("Quantify the number of colors using Histogram Clustering");
        m_action_ToonMapping_Hist->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_quantif_h.png"));

		// Connect appropriate signal
		connect(m_action_ToonMapping_Hist, &QAction::triggered, this, &ColorimetricSegmenter::HistogramClustering);

		connect(m_action_ToonMapping_Hist, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_ToonMapping_Hist, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_ToonMapping_Hist, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}
	if (!m_action_ToonMapping_KMeans)
	{
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action_ToonMapping_KMeans = new QAction("Kmeans Clustering", this);
        m_action_ToonMapping_KMeans->setToolTip("Quantify the number of colors using Kmeans Clustering");
        m_action_ToonMapping_KMeans->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_quantif_k.png"));

		// Connect appropriate signal
		connect(m_action_ToonMapping_KMeans, &QAction::triggered, this, &ColorimetricSegmenter::KmeansClustering);

		connect(m_action_ToonMapping_KMeans, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_ToonMapping_KMeans, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_ToonMapping_KMeans, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}


	return { m_action_filterRgb, m_action_filterHSV, /*m_action_filterRgbWithSegmentation,*/ m_action_filterScalar,m_action_ToonMapping_Hist, m_action_ToonMapping_KMeans };
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

	rgbDlg = new RgbDialog(m_pickingHub, (QWidget*)m_app->getMainWindow());
	rgbDlg->show();

	if (!rgbDlg->exec())
		return;

	// Start timer
	auto start = std::chrono::high_resolution_clock::now();

	// Get all values to make the color range with RGB values
	int redInf = std::min(rgbDlg->red_first->value(), rgbDlg->red_second->value());
	int redSup = std::max(rgbDlg->red_first->value(), rgbDlg->red_second->value());
	int greenInf = std::min(rgbDlg->green_first->value(), rgbDlg->green_second->value());
	int greenSup = std::max(rgbDlg->green_first->value(), rgbDlg->green_second->value());
	int blueInf = std::min(rgbDlg->blue_first->value(), rgbDlg->blue_second->value());
	int blueSup = std::max(rgbDlg->blue_first->value(), rgbDlg->blue_second->value());

	if (rgbDlg->margin->value() > 0) {
		// Get margin value (percent)
		double marginError = static_cast<double>(rgbDlg->margin->value()) / 100.0;

		redInf   -= marginError * redInf;
		redSup   += marginError * redSup;
		greenInf -= marginError * greenInf;
		greenSup += marginError * greenSup;
		blueInf  -= marginError * blueInf;
		blueSup  += marginError * blueSup;
	}

	// Set to min or max value (0-255)
	redInf = (redInf < MIN_VALUE ? MIN_VALUE : redInf);
	greenInf = (greenInf < MIN_VALUE ? MIN_VALUE : greenInf);
	blueInf = (blueInf < MIN_VALUE ? MIN_VALUE : blueInf);

	redSup = (redSup > MAX_VALUE ? MAX_VALUE : redSup);
	greenSup = (greenSup > MAX_VALUE ? MAX_VALUE : greenSup);
	blueSup = (blueSup > MAX_VALUE ? MAX_VALUE : blueSup);

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds) {
		if (cloud->hasColors())
		{
			// Use only references for speed reasons
			CCCoreLib::ReferenceCloud* filteredCloudInside = new CCCoreLib::ReferenceCloud(cloud);
			CCCoreLib::ReferenceCloud* filteredCloudOutside = new CCCoreLib::ReferenceCloud(cloud);

			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgb& rgb = cloud->getPointColor(j);
				(rgb.r >= redInf&& rgb.r <= redSup &&
					rgb.g >= greenInf&& rgb.g <= greenSup &&
					rgb.b >= blueInf&& rgb.b <= blueSup) ? addPoint(filteredCloudInside, j) : addPoint(filteredCloudOutside, j);
			}
			std::string name = "Rmin:" + std::to_string(redInf) + "/Gmin:" + std::to_string(greenInf) + "/Bmin:" + std::to_string(blueInf) +
				"/Rmax:" + std::to_string(redSup) + "/Gmax:" + std::to_string(greenSup) + "/Bmax:" + std::to_string(blueSup);

			createClouds<RgbDialog*>(rgbDlg, cloud, filteredCloudInside, filteredCloudOutside, name);

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


void ColorimetricSegmenter::filterScalar()
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

	scalarDlg = new ScalarDialog(m_pickingHub, (QWidget*)m_app->getMainWindow());
	scalarDlg->show();

	if (!scalarDlg->exec())
		return;

	auto start = std::chrono::high_resolution_clock::now();

	double marginError = static_cast<double>(scalarDlg->margin->value()) / 100.0;
	ScalarType min = std::min(scalarDlg->first->value(), scalarDlg->second->value());
	ScalarType max = std::max(scalarDlg->first->value(), scalarDlg->second->value());
	min -= (marginError * min); 
	max += (marginError * max); 

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();
	for (ccPointCloud* cloud : clouds) {

		// Use only references for speed reasons
		CCCoreLib::ReferenceCloud* filteredCloudInside = new CCCoreLib::ReferenceCloud(cloud);
		CCCoreLib::ReferenceCloud* filteredCloudOutside = new CCCoreLib::ReferenceCloud(cloud);

		for (unsigned j = 0; j < cloud->size(); ++j)
		{
			const ScalarType val = cloud->getPointScalarValue(j);
			(val > min&& val < max)
				? addPoint(filteredCloudInside, j) : addPoint(filteredCloudOutside, j);
		}
		std::string name = "min:" + std::to_string(min) + "/max:" + std::to_string(max);

		createClouds<ScalarDialog*>(scalarDlg, cloud, filteredCloudInside, filteredCloudOutside, name);

		m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
	}
	// Stop timer
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	QString s = QString::number(duration);

	//Print time of execution
	ccLog::Print("Time to execute : " + s + " milliseconds.");
}

/**
 * @brief knnRegions Determines the neighboring regions of a region.
 * @param basePointCloud The base cloud containing the points.
 * @param regions The list containing all the regions in the base cloud point.
 * @param region The region to compare with others.
 * @param k Max number of neighbours to find.
 * @param neighbours The resulting nearest regions.
 * @param thresholdDistance The maximum distance to search for neighbors.
 */
void knnRegions(ccPointCloud* basePointCloud, std::vector<CCCoreLib::ReferenceCloud*>* regions, const CCCoreLib::ReferenceCloud* region, unsigned k, std::vector<CCCoreLib::ReferenceCloud*>* neighbours, unsigned thresholdDistance) {
	ccPointCloud* computedRegion = basePointCloud->partialClone(region);
	// compute distances
	CCCoreLib::DistanceComputationTools::Cloud2CloudDistanceComputationParams params = CCCoreLib::DistanceComputationTools::Cloud2CloudDistanceComputationParams();
	params.kNNForLocalModel = k;
	params.maxSearchDist = thresholdDistance;
	// create to array, one containing regions, and another containing the distances to these regions.
	std::vector<CCCoreLib::ReferenceCloud*>* tempNeighbours = new std::vector<CCCoreLib::ReferenceCloud*>();
	std::vector<int>* distances = new std::vector<int>();
	for (CCCoreLib::ReferenceCloud* r : *regions)
	{
		distances->push_back(CCCoreLib::DistanceComputationTools::computeCloud2CloudDistance(computedRegion, basePointCloud->partialClone(r), params));
		tempNeighbours->push_back(r);
	}
	// sort the vectors
	std::vector<int>* index = new std::vector<int>(tempNeighbours->size());
	int n = 0;
	std::generate(index->begin(), index->end(),
		[n]
	()
		mutable
	{
		return n++;
	});

	std::sort(
		index->begin(), index->end(),
		[&](int a, int b) { return distances[a] < distances[b]; });

	// then extract the 'k' nearest neighbors.
	for (int i : *index)
	{
		if (neighbours->size() < k)
		{
			neighbours->push_back(tempNeighbours->at(i));
		}
	}

}

/**
 * @brief colorimetricalDifference Compute colorimetrical difference between two RGB color values.
 * @param c1 First color value.
 * @param c2 Second color value.
 * @return Colorimetrical difference.
 */
double colorimetricalDifference(ccColor::Rgb c1, ccColor::Rgb c2) {
	return sqrt(pow(c1.r - c2.r, 2) + pow(c1.g - c2.g, 2) + pow(c1.b - c2.b, 2));
}

ccColor::Rgb* meanRgb(ccPointCloud* basePointCloud, CCCoreLib::ReferenceCloud* c)
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
 * @brief colorimetricalDifference compute mean colorimetrical difference between two reference clouds.
 * The points in both clouds must be represented in RGB value.
 * @param basePointCloud The base cloud on which the reference clouds are based.
 * @param c1 The first reference cloud.
 * @param c2 The second reference cloud.
 * @return Colorimetrical difference.
 */
double colorimetricalDifference(ccPointCloud* basePointCloud, CCCoreLib::ReferenceCloud* c1, CCCoreLib::ReferenceCloud* c2) {
	ccColor::Rgb* rgb1 = meanRgb(basePointCloud, c1);

	ccColor::Rgb* rgb2 = meanRgb(basePointCloud, c2);

	return colorimetricalDifference(*rgb1, *rgb2);
}


std::vector<CCCoreLib::ReferenceCloud*>* ColorimetricSegmenter::regionGrowing(ccPointCloud* pointCloud, const unsigned TNN, const double TPP, const double TD)
{
	std::vector<unsigned> unlabeledPoints;
	for (unsigned j = 0; j < pointCloud->size(); ++j)
	{
		unlabeledPoints.push_back(j);
	}
	std::vector<CCCoreLib::ReferenceCloud*>* regions = new std::vector<CCCoreLib::ReferenceCloud*>();
	std::vector<unsigned>* points = new std::vector<unsigned>();
	CCCoreLib::DgmOctree* octree = new CCCoreLib::DgmOctree(pointCloud);// used to search nearest neighbors
	octree->build();
	// while there is any point in {P} that hasn’t been labeled
	while (unlabeledPoints.size() > 0)
	{
		// push an unlabeled point into stack Points
		points->push_back(unlabeledPoints.back());
		unlabeledPoints.pop_back();
		// initialize a new region Rc and add current point to R
		CCCoreLib::ReferenceCloud* rc = new CCCoreLib::ReferenceCloud(pointCloud);
		rc->addPointIndex(unlabeledPoints.back());
		// while stack Points is not empty
		while (points->size() > 0)
		{
			// pop Points’ top element Tpoint
			unsigned tPointIndex = points->back();
			points->pop_back();

			// for each point p in {KNNTNN(Tpoint)}
			CCCoreLib::DgmOctree::NearestNeighboursSearchStruct nNSS = CCCoreLib::DgmOctree::NearestNeighboursSearchStruct();
			nNSS.level = 1;
			nNSS.queryPoint = *(pointCloud->getPoint(tPointIndex));
			Tuple3i cellPos = Tuple3i();
			octree->getCellPos(octree->getCellCode(tPointIndex), 1, cellPos, false);
			nNSS.cellPos = cellPos;
			CCVector3 cellCenter;
			octree->computeCellCenter(octree->getCellCode(tPointIndex), 1, cellCenter);
			nNSS.cellCenter = cellCenter;
			nNSS.maxSearchSquareDistd = TD;
			nNSS.minNumberOfNeighbors = TNN;

			octree->findNearestNeighborsStartingFromCell(nNSS);
			CCCoreLib::DgmOctree::NeighboursSet knnResult = nNSS.pointsInNeighbourhood;

			for (int i = 0; i < knnResult.size(); i++)
			{
				unsigned p = knnResult.at(i).pointIndex;
				// if p is labelled
				if (std::find(unlabeledPoints.begin(), unlabeledPoints.end(), p) != unlabeledPoints.end())
				{
					continue;
				}
				// if CD(Tpoint,p)<TPP
				if (colorimetricalDifference(pointCloud->getPointColor(p), pointCloud->getPointColor(tPointIndex)) < TPP)
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

/**
 * @brief findRegion Find a given region in a vector of reference clouds.
 * @param container Container containing all the regions.
 * @param region Region to search for in the vector.
 * @return The pointer to the region if found, nullptr in the other case.
 */
std::vector<CCCoreLib::ReferenceCloud*>* findRegion(std::vector<std::vector<CCCoreLib::ReferenceCloud*>*>* container, CCCoreLib::ReferenceCloud* region)
{
	if (container->size() == 0)
	{
		return nullptr;
	}
	for (std::vector<CCCoreLib::ReferenceCloud*>* l : *container)
	{
		if (std::find(l->begin(), l->end(), region) != l->end())
		{
			return l;
		}
	}
	return nullptr;
}

std::vector<CCCoreLib::ReferenceCloud*>* ColorimetricSegmenter::regionMergingAndRefinement(ccPointCloud* basePointCloud, std::vector<CCCoreLib::ReferenceCloud*>* regions, const unsigned TNN, const double TRR, const double TD, const unsigned Min)
{
	std::vector<std::vector<CCCoreLib::ReferenceCloud*>*>* homogeneous = new std::vector<std::vector<CCCoreLib::ReferenceCloud*>*>();

	// for each region Ri in {R}
	for (CCCoreLib::ReferenceCloud* ri : *regions)
	{
		// if Ri is not in {H}
		if (findRegion(homogeneous, ri) == nullptr)
		{
			// create a new list to record Ri
			std::vector<CCCoreLib::ReferenceCloud*>* newRegionGroup = new std::vector<CCCoreLib::ReferenceCloud*>();
			newRegionGroup->push_back(ri);
			homogeneous->push_back(newRegionGroup);
		}

		// for each region Rj in {KNNTNN2,TD2(Ri)}
		std::vector<CCCoreLib::ReferenceCloud*>* knnResult = new std::vector<CCCoreLib::ReferenceCloud*>();
		knnRegions(basePointCloud, regions, ri, TNN, knnResult, TD);
		for (CCCoreLib::ReferenceCloud* rj : *knnResult)
		{
			// if CD(Ri,Rj)<TRR
			if (colorimetricalDifference(basePointCloud, ri, rj) < TNN)
			{
				// if Rj is in {H}
				std::vector<CCCoreLib::ReferenceCloud*>* regionContainer = findRegion(homogeneous, rj);
				if (regionContainer != nullptr)
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
	std::vector<CCCoreLib::ReferenceCloud*>* mergedRegionsRef = new std::vector<CCCoreLib::ReferenceCloud*>();
	for (std::vector<CCCoreLib::ReferenceCloud*>* l : *homogeneous)
	{
		CCCoreLib::ReferenceCloud* merged = l->at(0);
		for (int i = 1; i < l->size(); i++)
		{
			merged->add(*l->at(i));
		}
		mergedRegionsRef->push_back(merged);
	}
	
	//std::vector<CCCoreLib::ReferenceCloud*>* knnResult;
	// for each region Ri in {R’}
	/*for (CCCoreLib::ReferenceCloud* r : *mergedRegionsRef)
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
	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);

		return;
	}

    // Retrieve parameters from dialog
    if (m_app->pickingHub()) {
        m_pickingHub = m_app->pickingHub();
    }

	// Retrieve parameters from dialog
    rgbDlg = new RgbDialog(m_pickingHub, (QWidget*)m_app->getMainWindow());
    rgbDlg->show();

	auto start = std::chrono::high_resolution_clock::now();

	if (!rgbDlg->exec())
		return;
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
			std::vector<CCCoreLib::ReferenceCloud*>* regions = regionGrowing(cloud, TNN, TPP, TD);
			regions = regionMergingAndRefinement(cloud, regions, TNN, TRR, TD, Min);
			//m_app->dispToConsole(QString("[ColorimetricSegmenter] regions %1").arg(regions->size()), ccMainAppInterface::STD_CONSOLE_MESSAGE);

			// retrieve the nearest region (in color range)
			for (CCCoreLib::ReferenceCloud* r : *regions)
			{
				ccColor::Rgb* mean = meanRgb(cloud, r);
				if (redInf >= mean->r && redSup <= mean->r &&
					greenInf >= mean->g && greenSup <= mean->g &&
					blueInf >= mean->b && blueSup <= mean->b)
				{
					ccPointCloud* newCloud = cloud->partialClone(r);

					cloud->setEnabled(false);
					if (cloud->getParent()) {
						cloud->getParent()->addChild(newCloud);
					}

					m_app->addToDB(newCloud, false, true, false, false);

					m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered with segmentation ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
				}

			}


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
			CCCoreLib::ReferenceCloud* filteredCloudInside = new CCCoreLib::ReferenceCloud(cloud);
			CCCoreLib::ReferenceCloud* filteredCloudOutside = new CCCoreLib::ReferenceCloud(cloud);

			// We manually add color ranges with HSV values
			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgb& rgb = cloud->getPointColor(j);
				hsv hsv_current = hsvDlg->rgb2hsv(rgb);

				// Hue is useless here because the saturation is not high enough
				if (0 <= hsv_first.s && hsv_first.s <= 25 && 0 <= hsv_current.s && hsv_current.s <= 25)
				{
					// We only check value
					if (hsv_first.v >= 0 && hsv_first.v <= 25 && 0 <= hsv_current.v && hsv_current.v <= 25)  addPoint(filteredCloudInside, j); // black
					else if (hsv_first.v > 25 && hsv_first.v <= 60 && hsv_current.v > 25 && hsv_current.v <= 60)  addPoint(filteredCloudInside, j); // grey
					else if (hsv_first.v > 60 && hsv_first.v <= 100 && hsv_current.v > 60 && hsv_current.v <= 100) addPoint(filteredCloudInside, j); // white
					else addPoint(filteredCloudOutside, j);
				}
				else if (hsv_first.s > 25 && hsv_first.s <= 100 && hsv_current.s > 25 && hsv_current.s <= 100) {
					// We need to check value first
					if (0 <= hsv_first.v && hsv_first.v <= 25 && 0 <= hsv_current.v && hsv_current.v <= 25) addPoint(filteredCloudInside, j); // black
					// Then, we can check value
					else if (hsv_first.v > 25 && hsv_first.v <= 100 && hsv_current.v > 25 && hsv_current.v <= 100)
					{
						if (((hsv_first.h >= 0 && hsv_first.h <= 30) || (hsv_first.h >= 330 && hsv_first.h <= 360)) &&
							((hsv_current.h >= 0 && hsv_current.h <= 30) || (hsv_current.h >= 330 && hsv_current.h <= 360))) addPoint(filteredCloudInside, j); // red
						else if (hsv_first.h > 30 && hsv_first.h <= 90 && hsv_current.h > 30 && hsv_current.h <= 90)         addPoint(filteredCloudInside, j); // yellow
						else if (hsv_first.h > 90 && hsv_first.h <= 150 && hsv_current.h > 90 && hsv_current.h <= 150)         addPoint(filteredCloudInside, j); // green
						else if (hsv_first.h > 150 && hsv_first.h <= 210 && hsv_current.h > 150 && hsv_current.h <= 210)         addPoint(filteredCloudInside, j); // cyan
						else if (hsv_first.h > 210 && hsv_first.h <= 270 && hsv_current.h > 210 && hsv_current.h <= 270)         addPoint(filteredCloudInside, j); // blue
						else if (hsv_first.h > 270 && hsv_first.h <= 330 && hsv_current.h > 270 && hsv_current.h <= 330)         addPoint(filteredCloudInside, j); // magenta
						else addPoint(filteredCloudOutside, j);
					}
					else addPoint(filteredCloudOutside, j);
				}
				else addPoint(filteredCloudOutside, j);
			}

			std::string name = "h:" + std::to_string((int)hsv_first.h) + "/s:" + std::to_string((int)hsv_first.s) + "/v:" + std::to_string((int)hsv_first.v);
			createClouds<HSVDialog*>(hsvDlg, cloud, filteredCloudInside, filteredCloudOutside, name);

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

// Method to add point to a ReferenceCloud*
void ColorimetricSegmenter::addPoint(CCCoreLib::ReferenceCloud* filteredCloud, unsigned int j)
{
	if (!filteredCloud->addPointIndex(j))
	{
		//not enough memory
		delete filteredCloud;
		filteredCloud = nullptr;
		m_app->dispToConsole("[ColorimetricSegmenter] Error, filter canceled.");
	}
}

// Method to interact with the component "Which points to keep"
template <typename T>
void ColorimetricSegmenter::createClouds(T& dlg, ccPointCloud* cloud, CCCoreLib::ReferenceCloud* filteredCloudInside, CCCoreLib::ReferenceCloud* filteredCloudOutside, std::string name)
{

	if (dlg->retain->isChecked()) {
		createCloud(cloud, filteredCloudInside, name, true);
	}
	else if (dlg->exclude->isChecked()) {
		createCloud(cloud, filteredCloudOutside, name, false);
	}
	else if (dlg->both->isChecked()) {
		createCloud(cloud, filteredCloudInside, name, true);
		createCloud(cloud, filteredCloudOutside, name, false);
	}

}

// Method to create a new cloud
void ColorimetricSegmenter::createCloud(ccPointCloud* cloud, CCCoreLib::ReferenceCloud* referenceCloud, std::string name, bool inside) {
	ccPointCloud* newCloud = cloud->partialClone(referenceCloud);
	if (inside) {
		name += ".inside";
	}
	else {
		name += ".outside";
	}

	newCloud->setName(QString::fromStdString(name));
	cloud->setEnabled(false);
	if (cloud->getParent()) {
		cloud->getParent()->addChild(newCloud);
	}

	m_app->addToDB(newCloud, false, true, false, false);
}

/**
Generate nxnxn clusters of points according to their color value (RGB)
@param cloud : the cloud which we work with
@param clusterPerDim : coefficient uses to split each RGB component
Returns a map of nxnxn keys, for each key a vector of the points index in the partition
*/
std::map<int, std::vector<unsigned>> getKeyCluster(const ccPointCloud& cloud, int clusterPerDim) {

	float clusterSize = 256 / clusterPerDim;

	std::map<int, std::vector<unsigned>> keyMap;
	std::map<int, std::vector<unsigned>>::iterator it;


	for (unsigned i = 0; i < cloud.size(); i++) {

		const ccColor::Rgb& rgb = cloud.getPointColor(i);

		int redCluster = rgb.r / clusterSize;
		int greenCluster = rgb.g / clusterSize;
		int blueCluster = rgb.b / clusterSize;

		int index = redCluster + greenCluster * 10 + blueCluster * 100;
		it = keyMap.find(index);
		//check if the entry with this index already exists
		if (it == keyMap.end()) {
			//if no, we create it
			std::vector<unsigned> points = { i };
			keyMap.insert(std::pair<int, std::vector<unsigned>>(index, points));
		}
		else {
			//else we add the point in the container
			it->second.push_back(i);
		}


	}

	return keyMap;
}
/**
Compute the average color (RGB)
@param cloud : cloud which contains the points
@param bucket : vector of indexes of points
Returns average color (RGB)
*/
ccColor::Rgb computeAverageColor(const ccPointCloud& cloud, const std::vector<unsigned>& bucket) {
	unsigned  red = 0, green = 0, blue = 0;
	unsigned length = bucket.size();

	//other formula to compute the average can be used
	for (unsigned point : bucket) {
		const ccColor::Rgb rgb = cloud.getPointColor(point);
		red += rgb.r;
		green += rgb.g;
		blue += rgb.b;

	}
	red = floor(red / length);
	blue = floor(blue / length);
	green = floor(green / length);
	ccColor::Rgb res;

	res.r = static_cast<unsigned char>(red); res.b = static_cast<unsigned char>((blue)); res.g = static_cast<unsigned char>((green));
	return res;
}

/**
Compute the distance between two colors
/!\ the formula can be modified, here it is simple to be as quick as possible
*/
double ColorDistance(const ccColor::Rgb& c1, const ccColor::Rgb& c2) {
	return (c1.r - c2.r) + (c1.b - c2.b) + (c1.g - c2.g);
}
/**
Generate a pointcloud quantified using an histogram clustering
The purpose is to counter luminance variation due to the merge of different scans
*/
void ColorimetricSegmenter::HistogramClustering() {

	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);

		return;
	}
	// creation of the window
	quantiDlg = new QuantiDialog((QWidget*)m_app->getMainWindow());
	if (!quantiDlg->exec())
		return;

	// Start timer
	auto start = std::chrono::high_resolution_clock::now();

	int nbClusterByComponent = quantiDlg->area_quanti->value();



	std::vector<ccPointCloud*> clouds = ColorimetricSegmenter::getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds) {

		if (cloud->hasColors()) {

			ccPointCloud* histCloud = cloud->cloneThis();
			histCloud->setName(QString::fromStdString("HistogramClustering : Indice Q : " + std::to_string(nbClusterByComponent) + " //Couleurs : " + std::to_string(nbClusterByComponent * nbClusterByComponent * nbClusterByComponent)));

			std::map<int, std::vector<unsigned>> clusterMap;

			clusterMap = getKeyCluster(*histCloud, nbClusterByComponent);

			for (std::map<int, std::vector<unsigned>>::iterator it = clusterMap.begin(), end = clusterMap.end(); it != end; it++)
			{

				ccColor::Rgb averageColor = computeAverageColor(*histCloud, it->second);

				for (auto point : it->second) {
					(*histCloud).setPointColor(point, averageColor);
				}


			}

			cloud->setEnabled(false);
			if (cloud->getParent()) {
				cloud->getParent()->addChild(histCloud);
			}

			m_app->addToDB(histCloud, false, true, false, false);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully clustering ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);

		}
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	QString s = QString::number(duration);

	//Print time of execution
	ccLog::Print("Time to execute : " + s + " milliseconds.");

}
/**
K-means algorithm
@param k : k clusters
@param it : limit of iterations before returns a result
Returns a cloud quantified
*/
ccPointCloud* computeKmeansClustering(ccPointCloud* theCloud, unsigned char K, int it)
{
	//valid parameters?
	if (!theCloud || K == 0)
	{
		assert(false);
		return nullptr;
	}

	unsigned n = theCloud->size();
	if (n == 0)
		return nullptr;

	//on a besoin de memoire ici !
	std::vector<ccColor::Rgb> theKMeans;	//K clusters centers
	std::vector<unsigned char> belongings;	//index of the cluster the point belongs to
	std::vector<double> minDistsToMean; 	//distance to the nearest cluster center
	std::vector<unsigned> theKNums;			//number of points per clusters
	std::vector<unsigned> theOldKNums;		//number of points per clusters (prior to iteration)
	try
	{
		theKMeans.resize(n);
		belongings.resize(n);
		minDistsToMean.resize(n);
		theKNums.resize(K);
		theOldKNums.resize(K);
	}
	catch (const std::bad_alloc&)
	{
		//not enough memory
		return nullptr;
	}

	//init classes centers (regularly sampled
	unsigned step = n / K;
	for (unsigned char j = 0; j < K; ++j)
		theKMeans[j] = theCloud->getPointColor(step * j);


	//let's start
	bool meansHaveMoved = false;
	int iteration = 0;
	do
	{
		meansHaveMoved = false;
		++iteration;
		//
		std::map<unsigned char, std::vector<unsigned>> KGroups;
		{
			for (unsigned i = 0; i < n; ++i)
			{
				unsigned char minK = 0;

				ccColor::Rgb color = theCloud->getPointColor(i);
				minDistsToMean[i] = std::abs(ColorDistance(color, theKMeans[minK]));

				//we look for the nearest cluster center
				for (unsigned char j = 1; j < K; ++j)
				{
					double distToMean = std::abs(ColorDistance(color, theKMeans[j]));
					if (distToMean < minDistsToMean[i])
					{
						minDistsToMean[i] = distToMean;
						minK = j;
					}
				}


				belongings[i] = minK;
				//minDistsToMean[i] = V;
			}
		}

		//compute the clusters centers

		theOldKNums = theKNums;
		std::fill(theKNums.begin(), theKNums.end(), static_cast<unsigned>(0));
		for (unsigned i = 0; i < n; ++i)
		{
			auto it = KGroups.find(belongings[i]);
			if (it == KGroups.end()) {
				std::vector<unsigned> points = { i };
				KGroups.insert(std::pair<unsigned char, std::vector<unsigned>>(belongings[i], points));
			}
			else {
				it->second.push_back(i);
			}
			++theKNums[belongings[i]];
		}



		for (unsigned char j = 0; j < K; ++j)
		{
			ccColor::Rgb newMean = (KGroups[j].size() > 0 ? computeAverageColor(*theCloud, KGroups[j]) : theKMeans[j]);

			if (theOldKNums[j] != theKNums[j]) {
				meansHaveMoved = true;
			}

			theKMeans[j] = newMean;
		}



	} while (iteration < it);

	ccPointCloud* KCloud = theCloud->cloneThis();
	KCloud->setName(QString::fromStdString("Kmeans clustering : K : " + std::to_string(K)));

	//set color for each cluster
	for (unsigned i = 0; i < n; i++) {
		(*KCloud).setPointColor(i, theKMeans[belongings[i]]);
	}

	return KCloud;
}

/**
Algorithm based on k-means for clustering points cloud by its colors
*/
void ColorimetricSegmenter::KmeansClustering() {

	kmeansDlg = new KmeansDlg((QWidget*)m_app->getMainWindow());
	if (!kmeansDlg->exec())
		return;

	// Start timer
	auto start = std::chrono::high_resolution_clock::now();

	std::vector<ccPointCloud*> clouds = ColorimetricSegmenter::getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds)
	{

		ccPointCloud* kcloud = computeKmeansClustering(cloud, kmeansDlg->spinBox_k->value(), kmeansDlg->spinBox_it->value());

		cloud->setEnabled(false);
		if (cloud->getParent())
		{
			cloud->getParent()->addChild(kcloud);
		}

		m_app->addToDB(kcloud, false, true, false, false);
		m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully clustering ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);

	}
	// Stop timer
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	QString s = QString::number(duration);

	//Print time of execution
	ccLog::Print("Time to execute : " + s + " milliseconds");
}