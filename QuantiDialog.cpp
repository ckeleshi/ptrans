#include "QuantiDialog.h"


#include <QVariant>
QuantiDialog::QuantiDialog(QWidget* parent)
	: QDialog(parent)
	, Ui::QuantiDialog()
{
	
	setupUi(this);
	//connect(area_quanti, static_cast<void(QSpinBox::*)(double)>(&QSpinBox::valueChanged), this, SLOT(QuantiDialog::updateLabe()));
	
	//connect(area_quanti, SIGNAL(valueChanged(int)), this, &QuantiDialog::updateLabelValue);
	//connect(area_quanti, static_cast<void(QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &QuantiDialog::updateLabelValue);
	
		
	
}
