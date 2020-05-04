#ifndef KMEANSDIALOG_H
#define KMEANSDIALOG_H


#include <ui_KmeansDlg.h>

//Qt
#include <QDialog>

class KmeansDlg : public QDialog, public Ui::KmeansDialog
{
Q_OBJECT
public:
	explicit KmeansDlg(QWidget* parent = 0);
	
};

#endif 