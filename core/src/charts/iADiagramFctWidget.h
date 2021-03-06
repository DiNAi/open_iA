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

#include "iAChartWidget.h"
#include "iAPlotData.h"
#include "open_iA_Core_export.h"

#include <QSharedPointer>
#include <QPoint>
#include <QString>

#include <vector>

class QDomNode;
class QMenu;
class QPaintEvent;
class QPainter;

class vtkPiecewiseFunction;
class vtkColorTransferFunction;

class dlg_function;
class dlg_TFTable;
class MdiChild;

class open_iA_Core_API iADiagramFctWidget : public iAChartWidget
{
	Q_OBJECT

public:
	enum AdditionalMode { MOVE_NEW_POINT_MODE=Y_ZOOM_MODE+1, MOVE_POINT_MODE};
	
	static const int SELECTED_POINT_RADIUS = 10;
	static const int SELECTED_POINT_SIZE = 2*SELECTED_POINT_RADIUS;
	static const int POINT_RADIUS = 4;
	static const int POINT_SIZE = 2*POINT_RADIUS;
	
	static const int SELECTED_PIE_RADIUS = 16;
	static const int SELECTED_PIE_SIZE = 2 * SELECTED_PIE_RADIUS;
	static const int PIE_RADIUS = 16;
	static const int PIE_SIZE = 2 * PIE_RADIUS;

	iADiagramFctWidget(QWidget *parent, MdiChild *mdiChild,
		QString const & label = "Greyvalue", QString const & yLabel = "");
	virtual ~iADiagramFctWidget();

	int getSelectedFuncPoint() const;
	bool isFuncEndPoint(int index) const;
	bool isUpdateAutomatically() const;
	int ChartHeight() const;

	void SetTransferFunctions(vtkColorTransferFunction* ctf, vtkPiecewiseFunction* pwf);

	dlg_function *getSelectedFunction();
	std::vector<dlg_function*> &getFunctions();

	void SetAllowTrfReset(bool allow);
	void SetEnableAdditionalFunctions(bool enable);

	bool isTFTableCreated() const;
	void closeTFTable();
	QPoint getTFTablePos() const;
	void setTFTablePos(QPoint pos);

protected:
	void mousePressEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event); 
	void mouseMoveEvent(QMouseEvent *event);
	void mouseDoubleClickEvent(QMouseEvent *event);
	void enterEvent(QEvent *event);
	void keyPressEvent(QKeyEvent *event);
	void AddContextMenuEntries(QMenu* contextMenu) override;
	virtual void changeMode(int newMode, QMouseEvent *event) override;

	MdiChild*      activeChild;
	std::vector<dlg_function*> functions;
	bool updateAutomatically;
	unsigned int selectedFunction;
signals:
	void updateViews();
	void pointSelected();
	void noPointSelected();
	void endPointSelected();
	void active();
	void autoUpdateChanged(bool toggled);
	void applyTFForAll();
	void DblClicked();
	void updateTFTable();

public slots:
	int deletePoint();
	void changeColor(QMouseEvent *event = NULL);
	void autoUpdate(bool toggled);
	void resetTrf();
	void updateTrf();
	void loadTransferFunction();
	void loadTransferFunction(QDomNode &functionsNode);
	void saveTransferFunction();
	void applyTransferFunctionForAll();
	void addBezierFunction();
	void addGaussianFunction();
	bool loadFunctions();
	bool saveFunctions();
	void removeFunction();
	void showTFTable();
	void TFTableIsFinished();

protected:
	virtual void drawFunctions(QPainter &painter);
	bool m_showFunctions;
private:
	bool m_allowTrfReset;
	bool m_enableAdditionalFunctions;
	dlg_TFTable* TFTable;
	void NewTransferFunction();
	void DrawAfterPlots(QPainter& painter) override;
};
