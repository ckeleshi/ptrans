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

#ifndef RgbDialog_H
#define RgbDialog_H

#include <ui_RgbDialog.h>
#include "ccPickingListener.h"

//Qt
#include <ccHObject.h>
#include <QDialog>
#include <qcheckbox.h>

class ccGLWindow;
class ccPlane;
class ccHObject;
class ccPickingHub;

class RgbDialog : public QDialog, public ccPickingListener, public Ui::RgbDialog
{
	Q_OBJECT
public:
	explicit RgbDialog(ccPickingHub* pickingHub, QWidget* parent = 0);

	//! Inherited from ccPickingListener
	virtual void onItemPicked(const PickedItem& pi);

public slots:
	void pickPoint_first(bool);
	void pickPoint_second(bool);

protected: //members

	//! Picking window (if any)
	ccGLWindow* m_pickingWin;

	//! Picking hub
	ccPickingHub* m_pickingHub;

};

#endif // RgbDialog_H