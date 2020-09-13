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

//Local
#include "qColorimetricSegmenter.h"
#include "HSV.h"
#include "RgbDialog.h"
#include "HSVDialog.h"
#include "ScalarDialog.h"
#include "QuantiDialog.h"
#include "KmeansDlg.h"

//CloudCompare
#include <ccLog.h>
#include <ccPointCloud.h>

//CCCoreLib
#include <DistanceComputationTools.h>

//System
#include <algorithm>
#include <map>

//Qt
#include <QMainWindow>

static void ShowDurationNow(const std::chrono::high_resolution_clock::time_point& startTime)
{
	auto stopTime = std::chrono::high_resolution_clock::now();
	auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count();

	//Print duration of execution
	ccLog::Print("Time to execute: " + QString::number(duration_ms) + " milliseconds");
}

ColorimetricSegmenter::ColorimetricSegmenter(QObject* parent)
	: QObject(parent)
	, ccStdPluginInterface(":/CC/plugin/ColorimetricSegmenter/info.json")
	, m_action_filterRgb(nullptr)
	/*, m_action_filterRgbWithSegmentation(nullptr)*/
	, m_action_filterHSV(nullptr)
	, m_action_filterScalar(nullptr)
	, m_action_histogramClustering(nullptr)
	, m_action_kMeansClustering(nullptr)
	, m_addPointError(false)
{
}

void ColorimetricSegmenter::handleNewEntity(ccHObject* entity)
{
	Q_ASSERT(entity && m_app);
	m_app->addToDB(entity);
}

void ColorimetricSegmenter::handleEntityChange(ccHObject* entity)
{
	Q_ASSERT(entity && m_app);
	entity->prepareDisplayForRefresh_recursive();
	m_app->refreshAll();
	m_app->updateUI();
}

void ColorimetricSegmenter::handleErrorMessage(QString message)
{
	if (m_app)
		m_app->dispToConsole(message, ccMainAppInterface::ERR_CONSOLE_MESSAGE);
}

void ColorimetricSegmenter::onNewSelection(const ccHObject::Container& selectedEntities)
{
	// Only enable our action if something is selected.
	bool activateColorFilters = false;
	bool activateScalarFilter = false;
	for (ccHObject* entity : selectedEntities)
	{
		if (entity->isKindOf(CC_TYPES::POINT_CLOUD))
		{
			if (entity->hasColors())
			{
				activateColorFilters = true;
			}
			else if (entity->hasDisplayedScalarField())
			{
				activateScalarFilter = true;
			}
		}
	}

	//Activate only if only one of them is activated
	if (activateColorFilters && activateScalarFilter)
	{
		activateColorFilters = activateScalarFilter = false;
	}

	if (m_action_filterRgb)
		m_action_filterRgb->setEnabled(activateColorFilters);
	if (m_action_filterHSV)
		m_action_filterHSV->setEnabled(activateColorFilters);
	//if (m_action_filterRgbWithSegmentation)
	//	m_action_filterRgbWithSegmentation->setEnabled(activateColorFilters);
	if (m_action_filterScalar)
		m_action_filterScalar->setEnabled(activateScalarFilter);
	if (m_action_histogramClustering)
		m_action_histogramClustering->setEnabled(activateColorFilters);
	if (m_action_kMeansClustering)
		m_action_kMeansClustering->setEnabled(activateColorFilters);
}

QList<QAction*> ColorimetricSegmenter::getActions()
{
	// RGB Filter
	if (!m_action_filterRgb)
	{
		m_action_filterRgb = new QAction("Filter RGB", this);
		m_action_filterRgb->setToolTip("Filter the points of the selected cloud by RGB color");
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
		m_action_filterHSV->setToolTip("Filter the points of the selected cloud by HSV color");
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
		m_action_filterScalar->setToolTip("Filter the points of the selected cloud using scalar field");
        m_action_filterScalar->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_scalar.png"));

		// Connect appropriate signal
		connect(m_action_filterScalar, &QAction::triggered, this, &ColorimetricSegmenter::filterScalar);

		connect(m_action_filterScalar, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_filterScalar, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_filterScalar, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}
	
	if (!m_action_histogramClustering)
	{
		m_action_histogramClustering = new QAction("Histogram Clustering", this);
        m_action_histogramClustering->setToolTip("Quantify the number of colors using Histogram Clustering");
        m_action_histogramClustering->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_quantif_h.png"));

		// Connect appropriate signal
		connect(m_action_histogramClustering, &QAction::triggered, this, &ColorimetricSegmenter::HistogramClustering);

		connect(m_action_histogramClustering, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_histogramClustering, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_histogramClustering, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));

	}
	
	if (!m_action_kMeansClustering)
	{
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action_kMeansClustering = new QAction("Kmeans Clustering", this);
        m_action_kMeansClustering->setToolTip("Quantify the number of colors using Kmeans Clustering");
        m_action_kMeansClustering->setIcon(QIcon(":/CC/plugin/ColorimetricSegmenter/images/icon_quantif_k.png"));

		// Connect appropriate signal
		connect(m_action_kMeansClustering, &QAction::triggered, this, &ColorimetricSegmenter::KmeansClustering);

		connect(m_action_kMeansClustering, SIGNAL(newEntity(ccHObject*)), this, SLOT(handleNewEntity(ccHObject*)));
		connect(m_action_kMeansClustering, SIGNAL(entityHasChanged(ccHObject*)), this, SLOT(handleEntityChange(ccHObject*)));
		connect(m_action_kMeansClustering, SIGNAL(newErrorMessage(QString)), this, SLOT(handleErrorMessage(QString)));
	}

	return {	m_action_filterRgb,
				m_action_filterHSV,
				//m_action_filterRgbWithSegmentation,
				m_action_filterScalar,
				m_action_histogramClustering,
				m_action_kMeansClustering
	};
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
		if (selectedEntities[i]->isKindOf(CC_TYPES::POINT_CLOUD))
		{
			clouds.push_back(static_cast<ccPointCloud*>(selectedEntities[i]));
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
		m_app->dispToConsole("[ColorimetricSegmenter] No active 3D view", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();
	if (clouds.empty())
	{
		Q_ASSERT(false);
		return;
	}

	// Retrieve parameters from dialog
	RgbDialog rgbDlg(m_app->pickingHub(), m_app->getMainWindow());
	
	rgbDlg.show(); //necessary for setModal to be retained
	
	if (!rgbDlg.exec())
		return;

	// Start timer
	auto startTime = std::chrono::high_resolution_clock::now();

	// Get all values to make the color range with RGB values
	int redInf   = std::min( rgbDlg.red_first->value(),   rgbDlg.red_second->value()   );
	int redSup   = std::max( rgbDlg.red_first->value(),   rgbDlg.red_second->value()   );
	int greenInf = std::min( rgbDlg.green_first->value(), rgbDlg.green_second->value() );
	int greenSup = std::max( rgbDlg.green_first->value(), rgbDlg.green_second->value() );
	int blueInf  = std::min( rgbDlg.blue_first->value(),  rgbDlg.blue_second->value()  );
	int blueSup  = std::max( rgbDlg.blue_first->value(),  rgbDlg.blue_second->value()  );

	if (rgbDlg.margin->value() > 0)
	{
		// Get margin value (percent)
		double marginError = rgbDlg.margin->value() / 100.0;

		redInf   -= marginError * redInf;
		redSup   += marginError * redSup;
		greenInf -= marginError * greenInf;
		greenSup += marginError * greenSup;
		blueInf  -= marginError * blueInf;
		blueSup  += marginError * blueSup;
	}

	// Set to min or max value (0-255)
	{
		const int MIN_VALUE = 0;
		redInf   = std::max(redInf,   MIN_VALUE);
		greenInf = std::max(greenInf, MIN_VALUE);
		blueInf  = std::max(blueInf,  MIN_VALUE);

		const int MAX_VALUE = 255;
		redSup   = std::min(redSup,   MAX_VALUE);
		greenSup = std::min(greenSup, MAX_VALUE);
		blueSup  = std::min(blueSup,  MAX_VALUE);
	}

	for (ccPointCloud* cloud : clouds)
	{
		if (cloud && cloud->hasColors())
		{
			// Use only references for speed reasons
			CCCoreLib::ReferenceCloud filteredCloudInside(cloud);
			CCCoreLib::ReferenceCloud filteredCloudOutside(cloud);

			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgba& rgb = cloud->getPointColor(j);
				if (	rgb.r >= redInf   && rgb.r <= redSup
					&&	rgb.g >= greenInf && rgb.g <= greenSup
					&&	rgb.b >= blueInf  && rgb.b <= blueSup)
				{
					addPoint(filteredCloudInside, j);
				}
				else
				{
					addPoint(filteredCloudOutside, j);
				}

				if (m_addPointError)
				{
					return;
				}
			}
			QString name = "Rmin:" + QString::number(redInf) + "/Gmin:" + QString::number(greenInf) + "/Bmin:" + QString::number(blueInf) +
				"/Rmax:" + QString::number(redSup) + "/Gmax:" + QString::number(greenSup) + "/Bmax:" + QString::number(blueSup);

			createClouds<const RgbDialog&>(rgbDlg, cloud, filteredCloudInside, filteredCloudOutside, name);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
		}
	}

	ShowDurationNow(startTime);
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
		m_app->dispToConsole("[ColorimetricSegmenter] No active 3D view", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	// Retrieve parameters from dialog
	ScalarDialog scalarDlg(m_app->pickingHub(), m_app->getMainWindow());

	scalarDlg.show(); //necessary for setModal to be retained

	if (!scalarDlg.exec())
		return;

	// Start timer
	auto startTime = std::chrono::high_resolution_clock::now();

	double marginError = static_cast<double>(scalarDlg.margin->value()) / 100.0;
	ScalarType min = std::min(scalarDlg.first->value(), scalarDlg.second->value());
	ScalarType max = std::max(scalarDlg.first->value(), scalarDlg.second->value());
	min -= (marginError * min); 
	max += (marginError * max); 

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();
	for (ccPointCloud* cloud : clouds)
	{

		// Use only references for speed reasons
		CCCoreLib::ReferenceCloud filteredCloudInside(cloud);
		CCCoreLib::ReferenceCloud filteredCloudOutside(cloud);

		for (unsigned j = 0; j < cloud->size(); ++j)
		{
			const ScalarType val = cloud->getPointScalarValue(j);
			addPoint(val > min && val < max ? filteredCloudInside : filteredCloudOutside, j);

			if (m_addPointError)
			{
				return;
			}
		}
		QString name = "min:" + QString::number(min) + "/max:" + QString::number(max);

		createClouds<const ScalarDialog&>(scalarDlg, cloud, filteredCloudInside, filteredCloudOutside, name);

		m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
	}

	ShowDurationNow(startTime);
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
void knnRegions(ccPointCloud* basePointCloud,
				std::vector<CCCoreLib::ReferenceCloud*>* regions,
				const CCCoreLib::ReferenceCloud* region,
				unsigned k,
				std::vector<CCCoreLib::ReferenceCloud*>* neighbours,
				unsigned thresholdDistance)
{
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

// filterRgbWithSegmentation parameters
static const unsigned TNN = 1;
static const double TPP = 2.0;
static const double TD = 2.0;
static const double TRR = 2.0;
static const unsigned Min = 2;

void ColorimetricSegmenter::filterRgbWithSegmentation()
{
	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);

		return;
	}

	// Retrieve parameters from dialog
	RgbDialog rgbDlg(m_app->pickingHub(), m_app->getMainWindow());

	rgbDlg.show(); //necessary for setModal to be retained

	if (!rgbDlg.exec())
		return;

	// Start timer
	auto startTime = std::chrono::high_resolution_clock::now();

	// Get margin value (percent)
	double marginError = rgbDlg.margin->value() / 100.0;

	// Get all values to make the color range with RGB values
	int redInf = rgbDlg.red_first->value() - (marginError * rgbDlg.red_first->value());
	int redSup = rgbDlg.red_second->value() + marginError * rgbDlg.red_second->value();
	int greenInf = rgbDlg.green_first->value() - marginError * rgbDlg.green_first->value();
	int greenSup = rgbDlg.green_second->value() + marginError * rgbDlg.green_second->value();
	int blueInf = rgbDlg.blue_first->value() - marginError * rgbDlg.blue_first->value();
	int blueSup = rgbDlg.blue_second->value() + marginError * rgbDlg.blue_second->value();
	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();

	for (ccPointCloud* cloud : clouds)
	{
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
					if (cloud->getParent())
					{
						cloud->getParent()->addChild(newCloud);
					}

					m_app->addToDB(newCloud, false, true, false, false);

					m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered with segmentation ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
				}
			}
		}
	}

	ShowDurationNow(startTime);
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

	// Check valid window
	if (!m_app->getActiveGLWindow())
	{
		m_app->dispToConsole("[ColorimetricSegmenter] No active 3D view", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return;
	}

	std::vector<ccPointCloud*> clouds = getSelectedPointClouds();
	if (clouds.empty())
	{
		Q_ASSERT(false);
		return;
	}

	// Retrieve parameters from dialog
	HSVDialog hsvDlg(m_app->pickingHub(), m_app->getMainWindow());
	
	hsvDlg.show(); //necessary for setModal to be retained
	
	if (!hsvDlg.exec())
		return;

	// Start timer
	auto startTime = std::chrono::high_resolution_clock::now();

	// Get HSV values
	Hsv hsv_first;
	hsv_first.h = hsvDlg.hue_first->value();
	hsv_first.s = hsvDlg.sat_first->value();
	hsv_first.v = hsvDlg.val_first->value();

	for (ccPointCloud* cloud : clouds)
	{
		if (cloud->hasColors())
		{
			// Use only references for speed reasons
			CCCoreLib::ReferenceCloud filteredCloudInside(cloud);
			CCCoreLib::ReferenceCloud filteredCloudOutside(cloud);

			// We manually add color ranges with HSV values
			for (unsigned j = 0; j < cloud->size(); ++j)
			{
				const ccColor::Rgb& rgb = cloud->getPointColor(j);
				Hsv hsv_current(rgb);

				// If Saturation is too small, considering Hue is useless
				if (0 <= hsv_first.s && hsv_first.s <= 25 && 0 <= hsv_current.s && hsv_current.s <= 25)
				{
					// We only check Value
					if (	(hsv_first.v >=  0 && hsv_first.v <=  25 && hsv_current.v >=  0 && hsv_current.v <=  25) //black
						||	(hsv_first.v >  25 && hsv_first.v <=  60 && hsv_current.v >  25 && hsv_current.v <=  60) //grey
						||	(hsv_first.v >  60 && hsv_first.v <= 100 && hsv_current.v >  60 && hsv_current.v <= 100) //white
						)
					{
						addPoint(filteredCloudInside, j);
					}
					else
					{
						addPoint(filteredCloudOutside, j);
					}
				}
				else if (hsv_first.s > 25 && hsv_first.s <= 100 && hsv_current.s > 25 && hsv_current.s <= 100)
				{
					if (0 <= hsv_first.v && hsv_first.v <= 25 && 0 <= hsv_current.v && hsv_current.v <= 25)
					{
						addPoint(filteredCloudInside, j); // black
					}
					else if (hsv_first.v > 25 && hsv_first.v <= 100 && hsv_current.v > 25 && hsv_current.v <= 100)
					{
						if (((hsv_first.h   >= 0 && hsv_first.h   <= 30) || (hsv_first.h   >= 330 && hsv_first.h   <= 360)) &&
							((hsv_current.h >= 0 && hsv_current.h <= 30) || (hsv_current.h >= 330 && hsv_current.h <= 360))
							)
						{
							addPoint(filteredCloudInside, j); // red
						}
						else if (	(hsv_first.h >  30 && hsv_first.h <=  90 && hsv_current.h >  30 && hsv_current.h <=  90) // yellow
								||	(hsv_first.h >  90 && hsv_first.h <= 150 && hsv_current.h >  90 && hsv_current.h <= 150) // green
								||	(hsv_first.h > 150 && hsv_first.h <= 210 && hsv_current.h > 150 && hsv_current.h <= 210) // cyan
								||	(hsv_first.h > 210 && hsv_first.h <= 270 && hsv_current.h > 210 && hsv_current.h <= 270) // blue
								||	(hsv_first.h > 270 && hsv_first.h <= 330 && hsv_current.h > 270 && hsv_current.h <= 330) // magenta
							)
						{
							addPoint(filteredCloudInside, j);
						}
						else
						{
							addPoint(filteredCloudOutside, j);
						}
					}
					else
					{
						addPoint(filteredCloudOutside, j);
					}
				}
				else
				{
					addPoint(filteredCloudOutside, j);
				}

				if (m_addPointError)
				{
					return;
				}
			}

			QString name = "h:" + QString::number(hsv_first.h, 'f', 0) + "/s:" + QString::number(hsv_first.s, 'f', 0) + "/v:" + QString::number(hsv_first.v, 'f', 0);
			createClouds<const HSVDialog&>(hsvDlg, cloud, filteredCloudInside, filteredCloudOutside, name);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully filtered ! ", ccMainAppInterface::STD_CONSOLE_MESSAGE);
		}
	}

	ShowDurationNow(startTime);
}

// Method to add point to a ReferenceCloud*
bool ColorimetricSegmenter::addPoint(CCCoreLib::ReferenceCloud& filteredCloud, unsigned int j)
{
	m_addPointError = !filteredCloud.addPointIndex(j);
	
	if (m_addPointError)
	{
		//not enough memory
		m_app->dispToConsole("[ColorimetricSegmenter] Error, filter canceled.");
	}

	return m_addPointError;
}

// Method to interact with the component "Which points to keep"
template <typename T>
void ColorimetricSegmenter::createClouds(	T& dlg,
											ccPointCloud* cloud,
											const CCCoreLib::ReferenceCloud& filteredCloudInside,
											const CCCoreLib::ReferenceCloud& filteredCloudOutside,
											QString name )
{
	if (dlg.retain->isChecked())
	{
		createCloud(cloud, filteredCloudInside, name + ".inside");
	}
	else if (dlg.exclude->isChecked())
	{
		createCloud(cloud, filteredCloudOutside, name + ".outside");
	}
	else if (dlg.both->isChecked())
	{
		createCloud(cloud, filteredCloudInside, name + ".inside");
		createCloud(cloud, filteredCloudOutside, name + ".outside");
	}

}

// Method to create a new cloud
void ColorimetricSegmenter::createCloud(ccPointCloud* cloud,
										const CCCoreLib::ReferenceCloud& referenceCloud,
										QString name)
{
	if (!cloud)
	{
		Q_ASSERT(false);
		return;
	}
	
	ccPointCloud* newCloud = cloud->partialClone(&referenceCloud);
	if (!newCloud)
	{
		m_app->dispToConsole("Not enough memory");
		return;
	}
	
	newCloud->setName(name);
	cloud->setEnabled(false);
	if (cloud->getParent())
	{
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
typedef std::map< size_t, std::vector<unsigned> > ClusterMap;
static bool GetKeyCluster(const ccPointCloud& cloud, size_t clusterPerDim, ClusterMap& clusterMap)
{
	Q_ASSERT(ccColor::MAX == 255);

	try
	{
		for (unsigned i = 0; i < cloud.size(); i++)
		{
			const ccColor::Rgb& rgb = cloud.getPointColor(i);

			size_t redCluster   = (static_cast<size_t>(rgb.r) * clusterPerDim) >> 8; // shift 8 bits (= division by 256)
			size_t greenCluster = (static_cast<size_t>(rgb.g) * clusterPerDim) >> 8; // shift 8 bits (= division by 256)
			size_t blueCluster  = (static_cast<size_t>(rgb.b) * clusterPerDim) >> 8; // shift 8 bits (= division by 256)

			size_t index = redCluster + (greenCluster + blueCluster * clusterPerDim) * clusterPerDim;

			//we add the point to the right container
			clusterMap[index].push_back(i);
		}
	}
	catch (const std::bad_alloc&)
	{
		return false;
	}

	return true;
}
/**
Compute the average color (RGB)
@param cloud : cloud which contains the points
@param bucket : vector of indexes of points
Returns average color (RGB)
*/
static ccColor::Rgba ComputeAverageColor(const ccPointCloud& cloud, const std::vector<unsigned>& bucket)
{
	size_t count = bucket.size();
	if (count == 0)
	{
		return ccColor::white;
	}
	else if (count == 1)
	{
		return cloud.getPointColor(bucket.front());
	}

	//other formula to compute the average can be used
	size_t redSum = 0, greenSum = 0, blueSum = 0, alphaSum = 0;
	for (unsigned pointIndex : bucket)
	{
		const ccColor::Rgba& rgba = cloud.getPointColor(pointIndex);
		redSum   += rgba.r;
		greenSum += rgba.g;
		blueSum  += rgba.b;
		alphaSum += rgba.a;
	}

	ccColor::Rgba res(	static_cast<ColorCompType>(std::min(redSum   / count, static_cast<size_t>(ccColor::MAX))),
						static_cast<ColorCompType>(std::min(greenSum / count, static_cast<size_t>(ccColor::MAX))),
						static_cast<ColorCompType>(std::min(blueSum  / count, static_cast<size_t>(ccColor::MAX))),
						static_cast<ColorCompType>(std::min(alphaSum / count, static_cast<size_t>(ccColor::MAX))));

	return res;
}

/**
Compute the distance between two colors
/!\ the formula can be modified, here it is simple to be as quick as possible
*/
static int ColorDistance(const ccColor::Rgb& c1, const ccColor::Rgb& c2)
{
	return (static_cast<int>(c1.r) - c2.r) + (static_cast<int>(c1.b) - c2.b) + (static_cast<int>(c1.g) - c2.g);
}

/**
Generate a pointcloud quantified using an histogram clustering
The purpose is to counter luminance variation due to the merge of different scans
*/
void ColorimetricSegmenter::HistogramClustering()
{
	if (m_app == nullptr)
	{
		// m_app should have already been initialized by CC when plugin is loaded
		Q_ASSERT(false);

		return;
	}

	std::vector<ccPointCloud*> clouds = ColorimetricSegmenter::getSelectedPointClouds();
	if (clouds.empty())
	{
		Q_ASSERT(false);
		return;
	}

	// creation of the window
	QuantiDialog quantiDlg(m_app->getMainWindow());
	if (!quantiDlg.exec())
		return;

	// Start timer
	auto startTime = std::chrono::high_resolution_clock::now();

	int nbClusterByComponent = quantiDlg.area_quanti->value();
	if (nbClusterByComponent < 0)
	{
		Q_ASSERT(false);
		return;
	}

	for (ccPointCloud* cloud : clouds)
	{
		if (cloud->hasColors())
		{
			ClusterMap clusterMap;
			if (!GetKeyCluster(*cloud, static_cast<size_t>(nbClusterByComponent), clusterMap))
			{
				m_app->dispToConsole("Not enough memory", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
				break;
			}

			ccPointCloud* histCloud = cloud->cloneThis();
			if (!histCloud)
			{
				m_app->dispToConsole("Not enough memory", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
				break;
			}

			histCloud->setName(QString("HistogramClustering: Indice Q = %1 // colors = %2").arg(nbClusterByComponent).arg(nbClusterByComponent * nbClusterByComponent * nbClusterByComponent));

			for (auto it = clusterMap.begin(); it != clusterMap.end(); it++)
			{
				ccColor::Rgba averageColor = ComputeAverageColor(*histCloud, it->second);

				for (unsigned pointIndex : it->second)
				{
					(*histCloud).setPointColor(pointIndex, averageColor);
				}
			}

			cloud->setEnabled(false);
			if (cloud->getParent())
			{
				cloud->getParent()->addChild(histCloud);
			}

			m_app->addToDB(histCloud, false, true, false, false);

			m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully clustered!", ccMainAppInterface::STD_CONSOLE_MESSAGE);

		}
	}

	ShowDurationNow(startTime);
}

/**
K-means algorithm
@param k : k clusters
@param it : limit of iterations before returns a result
Returns a cloud quantified
*/
static ccPointCloud* ComputeKmeansClustering(ccPointCloud* theCloud, unsigned K, int maxIterationCount)
{
	//valid parameters?
	if (!theCloud || K == 0)
	{
		Q_ASSERT(false);
		return nullptr;
	}

	unsigned pointCount = theCloud->size();
	if (pointCount == 0)
		return nullptr;

	if (K >= pointCount)
	{
		ccLog::Warning("Cloud %1 has less point than the expected number of classes.");
		return nullptr;
	}

	ccPointCloud* KCloud = nullptr;

	try
	{
		std::vector<ccColor::Rgba> clusterCenters;	//K clusters centers
		std::vector<int> clusterIndex;				//index of the cluster the point belongs to

		clusterIndex.resize(pointCount);
		clusterCenters.resize(K);

		//init (regularly sampled) classes centers
		double step = static_cast<double>(pointCount) / K;
		for (unsigned j = 0; j < K; ++j)
		{
			//TODO: this initialization is pretty biased... To be improved?
			clusterCenters[j] = theCloud->getPointColor(static_cast<unsigned>(std::ceil(step * j)));
		}

		//let's start
		int iteration = 0;
		for (; iteration < maxIterationCount; ++iteration)
		{
			bool meansHaveMoved = false;

			// assign each point (color) to the nearest cluster
			for (unsigned i = 0; i < pointCount; ++i)
			{
				const ccColor::Rgba& color = theCloud->getPointColor(i);

				int minK = 0;
				int minDistsToMean = std::abs(ColorDistance(color, clusterCenters[minK]));

				//we look for the nearest cluster center
				for (unsigned j = 1; j < K; ++j)
				{
					int distToMean = std::abs(ColorDistance(color, clusterCenters[j]));
					if (distToMean < minDistsToMean)
					{
						minDistsToMean = distToMean;
						minK = j;
					}
				}

				clusterIndex[i] = minK;
			}

			//update the clusters centers
			std::vector< std::vector<unsigned> > clusters;
			clusters.resize(K);
			for (unsigned i = 0; i < pointCount; ++i)
			{
				unsigned index = clusterIndex[i];
				clusters[index].push_back(i);
			}

			ccLog::Print("Iteration " + QString::number(iteration));
			for (unsigned j = 0; j < K; ++j)
			{
				const std::vector<unsigned>& cluster = clusters[j];
				if (cluster.empty())
				{
					continue;
				}
				
				ccColor::Rgba newMean = ComputeAverageColor(*theCloud, cluster);

				if (!meansHaveMoved && ColorDistance(clusterCenters[j], newMean) != 0)
				{
					meansHaveMoved = true;
				}

				clusterCenters[j] = newMean;
			}

			if (!meansHaveMoved)
			{
				break;
			}
		}

		KCloud = theCloud->cloneThis();
		if (!KCloud)
		{
			//not enough memory
			return nullptr;
		}
		KCloud->setName("Kmeans clustering: K = " + QString::number(K) + " / it = " + QString::number(iteration));

		//set color for each cluster
		for (unsigned i = 0; i < pointCount; i++)
		{
			KCloud->setPointColor(i, clusterCenters[clusterIndex[i]]);
		}

	}
	catch (const std::bad_alloc&)
	{
		//not enough memory
		return nullptr;
	}

	return KCloud;
}

/**
Algorithm based on k-means for clustering points cloud by its colors
*/
void ColorimetricSegmenter::KmeansClustering()
{
	std::vector<ccPointCloud*> clouds = ColorimetricSegmenter::getSelectedPointClouds();
	if (clouds.empty())
	{
		Q_ASSERT(false);
		return;
	}

	KmeansDlg kmeansDlg(m_app->getMainWindow());
	if (!kmeansDlg.exec())
		return;

	assert(kmeansDlg.spinBox_k->value() >= 0);
	unsigned K = static_cast<unsigned>(kmeansDlg.spinBox_k->value());
	int iterationCount = kmeansDlg.spinBox_it->value();

	// Start timer
	auto startTime = std::chrono::high_resolution_clock::now();

	for (ccPointCloud* cloud : clouds)
	{
		ccPointCloud* kcloud = ComputeKmeansClustering(cloud, K, iterationCount);
		if (!kcloud)
		{
			m_app->dispToConsole(QString("[ColorimetricSegmenter] Failed to cluster cloud %1").arg(cloud->getName()), ccMainAppInterface::WRN_CONSOLE_MESSAGE);
			continue;
		}

		cloud->setEnabled(false);
		if (cloud->getParent())
		{
			cloud->getParent()->addChild(kcloud);
		}

		m_app->addToDB(kcloud, false, true, false, false);
		m_app->dispToConsole("[ColorimetricSegmenter] Cloud successfully clustered!", ccMainAppInterface::STD_CONSOLE_MESSAGE);

	}

	ShowDurationNow(startTime);
}
