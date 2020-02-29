#ifndef EXAMPLE_PLUGIN_HEADER
#define EXAMPLE_PLUGIN_HEADER

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
#include <QObject>
#include <QtGui>

#include <ccPickingListener.h>
#include <ccPickingHub.h>
#include <ccGLWindow.h>
#include "ccPointCloud.h"
#include "ccScalarField.h"

#include "RgbDialog.h"
#include "HSVDialog.h"

//! Example qCC plugin
/** Replace 'ExamplePlugin' by your own plugin class name throughout and then
	check 'ExamplePlugin.cpp' for more directions.
	
	Each plugin requires an info.json file to provide information about itself -
	the name, authors, maintainers, icon, etc..
	
	The one method you are required to implement is 'getActions'. This should
	return all actions (QAction objects) for the plugin. CloudCompare will
	automatically add these with their icons in the plugin toolbar and to the
	plugin menu. If	your plugin returns	several actions, CC will create a
	dedicated toolbar and a	sub-menu for your plugin. You are responsible for
	connecting these actions to	methods in your plugin.
	
	Use the ccStdPluginInterface::m_app variable for access to most of the CC
	components (database, 3D views, console, etc.) - see the ccMainAppInterface
	class in ccMainAppInterface.h.
**/
class ColorimetricSegmenter : public QObject, public ccStdPluginInterface
{
	Q_OBJECT
	Q_INTERFACES(ccStdPluginInterface)
	
	// Replace "Example" by your plugin name (IID should be unique - let's hope your plugin name is unique ;)
	// The info.json file provides information about the plugin to the loading system and
	// it is displayed in the plugin information dialog.
	Q_PLUGIN_METADATA(IID "cccorp.cloudcompare.plugin.ColorimetricSegmenter" FILE "info.json")
	
public:
	explicit ColorimetricSegmenter( QObject *parent = nullptr );
	~ColorimetricSegmenter() override = default;
	
	// inherited from ccStdPluginInterface
	void onNewSelection( const ccHObject::Container &selectedEntities ) override;
	QList<QAction *> getActions() override;

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

	//! Filter a cloud with a scalar field (add one if there is none)
	void filterScalarField();

	//! Filter a cloud with RGB color
	void filterRgb();

	void filterHSV();

	void addPoint(CCLib::ReferenceCloud* filteredCloud, unsigned int j);

	//picked point callbacks
	//void pointPicked(ccHObject* entity, unsigned itemIdx, int x, int y, const CCVector3& P);
	//virtual void onItemPicked(const ccPickingListener::PickedItem& pi); //inherited from ccPickingListener
	
	//! Default action
	/** You can add as many actions as you want in a plugin.
		Each action will correspond to an icon in the dedicated
		toolbar and an entry in the plugin menu.
	**/
	QAction* m_action_filterScalar;
	QAction* m_action_filterRgb;
	QAction* m_action_filterHSV;

	//! Picking hub
	ccPickingHub* m_pickingHub = nullptr;

	RgbDialog* rgbDlg;
	HSVDialog* hsvDlg;

	//link to application windows
	//ccGLWindow* m_window;
	//QMainWindow* m_main_window;

};

#endif
