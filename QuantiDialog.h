#ifndef QUANTIDIALOG_H
#define QUANTIDIALOG_H


#include <ui_QuantiDialog.h>

//Qt
#include <QDialog>

class QuantiDialog : public QDialog, public Ui::QuantiDialog
{
Q_OBJECT
public:
	explicit QuantiDialog(QWidget* parent = 0);


public slots:
	void updateLabelValue(double value) {
		nb_color_label->setText(QString::fromStdString(std::to_string(value * value)));
	}
};

#endif 