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

#include "ccStdPluginInterface.h"
#include "ccPointCloud.h"
#include <QObject>
#include <QtGui>

#include <ccPickingListener.h>
#include <ccPickingHub.h>
#include <ccGLWindow.h>
#include "ccPointCloud.h"
#include "ccScalarField.h"

#include "RgbDialog.h"
#include "HSVDialog.h"
#include "ScalarDialog.h"
#include "QuantiDialog.h"
#include "KmeansDlg.h"

const int MIN_VALUE = 0;
const int MAX_VALUE = 255;

class ColorimetricSegmenter : public QObject, public ccStdPluginInterface
{
	Q_OBJECT
	Q_INTERFACES(ccPluginInterface ccStdPluginInterface)
	Q_PLUGIN_METADATA(IID "cccorp.cloudcompare.plugin.ColorimetricSegmenter" FILE "info.json")

public:
	explicit ColorimetricSegmenter(QObject* parent = nullptr);
	~ColorimetricSegmenter() override = default;

	// inherited from ccStdPluginInterface
	void onNewSelection(const ccHObject::Container& selectedEntities) override;
	QList<QAction*> getActions() override;

public slots:
	//! Handles new entity
	void handleNewEntity(ccHObject*);

	//! Handles entity (visual) modification
	void handleEntityChange(ccHObject*);

	//! Handles new error message
	void handleErrorMessage(QString);

signals:

	//! Signal emitted when an entity is (visually) modified
	void entityHasChanged(ccHObject*);

	//! Signal emitted when a new entity is created by the filter
	void newEntity(ccHObject*);

	//! Signal emitted when a new error message is produced
	void newErrorMessage(QString);

private:
	std::vector<ccPointCloud*> getSelectedPointClouds();

	//! Filter a cloud with RGB color
	void filterRgb();

	void filterHSV();

	void filterScalar();

	void HistogramClustering();

	void KmeansClustering();

	void addPoint(CCCoreLib::ReferenceCloud* filteredCloud, unsigned int j);

	template <typename T>
	void createClouds(T& dlg, ccPointCloud* cloud, CCCoreLib::ReferenceCloud* filteredCloudInside, CCCoreLib::ReferenceCloud* filteredCloudOutside, std::string name);

	void createCloud(ccPointCloud* cloud, CCCoreLib::ReferenceCloud* referenceCloud, std::string name, bool inside);

	//picked point callbacks
	//void pointPicked(ccHObject* entity, unsigned itemIdx, int x, int y, const CCVector3& P);
	//virtual void onItemPicked(const ccPickingListener::PickedItem& pi); //inherited from ccPickingListener

	//! Segment a cloud with RGB color
	void filterRgbWithSegmentation();

	/**
	 * @brief regionGrowing Segmentation method grouping the points into regions of similar colors.
	 * Method described in Qingming Zhan, Yubin Liang, Yinghui Xiao, 2009 "Color-based segmentation of point clouds".
	 * @param pointCloud The point cloud to segment.
	 * @param TNN Point-point colorimetrical similarity threshold.
	 * @param TPP Number of neighbours to search using KNN.
	 * @param TD Threshold distance between neighbouring points.
	 * @return Vector containing the resulting regions.
	 */
	std::vector<CCCoreLib::ReferenceCloud*>* regionGrowing(ccPointCloud* pointCloud, const unsigned TNN, const double TPP, const double TD);

	/**
	 * @brief regionMergingAndRefinement Merge previously created regions in 'regionGrowing' method.
	 * @param basePointCloud The base segmented point cloud used to create the regions.
	 * @param regions Vector containing the regions.
	 * @param TNN Point-point colorimetrical similarity threshold.
	 * @param TRR Region-region colorimetrical similarity threshold.
	 * @param TD Threshold distance between neighbouring regions. Used to merge close regions.
	 * @param Min Minimal size for a region.
	 * @return Vector containing the resulting merged and refined regions.
	 */
	std::vector<CCCoreLib::ReferenceCloud*>* regionMergingAndRefinement(ccPointCloud* basePointCloud, std::vector<CCCoreLib::ReferenceCloud*>* regions, const unsigned TNN, const double TRR, const double TD, const unsigned Min);


	//! Default action
	/** You can add as many actions as you want in a plugin.
		Each action will correspond to an icon in the dedicated
		toolbar and an entry in the plugin menu.
	**/
	QAction* m_action_filterRgb;
	//QAction* m_action_filterRgbWithSegmentation;
	QAction* m_action_filterHSV;
	QAction* m_action_filterScalar;
	QAction* m_action_ToonMapping_Hist;
	QAction* m_action_ToonMapping_KMeans;


	//! Picking hub
	ccPickingHub* m_pickingHub = nullptr;

	RgbDialog* rgbDlg;
	HSVDialog* hsvDlg;
	ScalarDialog* scalarDlg;
	QuantiDialog* quantiDlg;
	KmeansDlg* kmeansDlg;


	//link to application windows
	//ccGLWindow* m_window;
	//QMainWindow* m_main_window;

	const unsigned TNN = 1;
	const double TPP = 2.0;
	const double TD = 2.0;
	const double TRR = 2.0;
	const unsigned Min = 2;
};