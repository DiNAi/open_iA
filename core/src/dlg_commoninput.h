/*************************************  open_iA  ************************************ *
* **********  A tool for scientific visualisation and 3D image processing  ********** *
* *********************************************************************************** *
* Copyright (C) 2016-2017  C. Heinzl, M. Reiter, A. Reh, W. Li, M. Arikan,            *
*                          J. Weissenböck, Artem & Alexander Amirkhanov, B. Fröhler   *
* *********************************************************************************** *
* This program is free software: you can redistribute it and/or modify it under the   *
* terms of the GNU General Public License as published by the Free Software           *
* Foundation, either version 3 of the License, or (at your option) any later version. *
*                                                                                     *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY     *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     *
* PARTICULAR PURPOSE.  See the GNU General Public License for more details.           *
*                                                                                     *
* You should have received a copy of the GNU General Public License along with this   *
* program.  If not, see http://www.gnu.org/licenses/                                  *
* *********************************************************************************** *
* Contact: FH OÖ Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*          Stelzhamerstraße 23, 4600 Wels / Austria, Email: c.heinzl@fh-wels.at       *
* ************************************************************************************/
#pragma once

#include "open_iA_Core_export.h"
#include "ui_CommonInput.h"

#include <QDialog>
#include <QString>
#include <QStringList>
#include <QList>
#include <QTextDocument>

class MdiChild;
class QWidget;
class QErrorMessage;
class QLabel;
class QScrollArea;

class open_iA_Core_API dlg_commoninput : public QDialog, public Ui_CommonInput
{
	Q_OBJECT
public:
	dlg_commoninput ( QWidget *parent, QString winTitel, QStringList inList, QList<QVariant> inPara, QTextDocument *fDescr = new QTextDocument( 0 ));
	int getCheckValue(int index) const;
	QString getComboBoxValue(int index) const;
	int getComboBoxIndex(int index) const;
	QString getText(int index) const;
	int getIntValue(int index) const;
	double getDblValue(int index) const;
	void showROI(MdiChild *child);
	int exec() override;
private:
	QWidget *container;
	int m_roi[6];
	MdiChild *m_roiMdiChild;
	bool m_roiMdiChildClosed;
	void updateValues(QList<QVariant>);
	void UpdateROIPart(QString const & partName, QString const & value);
private slots:
	void ROIUpdated(QString text);
	void ROIChildClosed();
protected:
	QStringList widgetList;
};
