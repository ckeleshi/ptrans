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

#ifndef HSVDialog_H
#define HSVDialog_H

#include <ui_HSVDialog.h>
#include "ccPickingListener.h"

//Qt
#include <ccHObject.h>
#include <QDialog>
#include <qcheckbox.h>

class ccGLWindow;
class ccPlane;
class ccHObject;
class ccPickingHub;

/*
	Struct for HSV
*/
typedef struct {
	double h;
	double s;
	double v;
} hsv;

/*
	Get the values of the HSV interface, and interactions
*/
class HSVDialog : public QDialog, public ccPickingListener, public Ui::HSVDialog
{
	Q_OBJECT
public:
	explicit HSVDialog(ccPickingHub* pickingHub, QWidget* parent = 0);


	//! Inherited from ccPickingListener
	virtual void onItemPicked(const PickedItem& pi);

	hsv rgb2hsv(ccColor::Rgb rgb);

public slots:
	void pickPoint(bool);

protected: //members

	//! Picking window (if any)
	ccGLWindow* m_pickingWin;

	//! Picking hub
	ccPickingHub* m_pickingHub;

};

#endif // HSVDialog_H
