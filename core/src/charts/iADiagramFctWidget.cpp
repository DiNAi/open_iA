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
#include "pch.h"
#include "iADiagramFctWidget.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include "dlg_bezier.h"
#include "dlg_function.h"
#include "dlg_gaussian.h"
#include "dlg_TFTable.h"
#include "dlg_transfer.h"
#include "iAPlotData.h"
#include "iAPlot.h"
#include "iAPlotTypes.h"
#include "iAMathUtility.h"
#include "iASettings.h"
#include "mainwindow.h"		// TODO: get rid of this inclusion!
#include "mdichild.h"

#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>

#include <QFileDialog>
#include <QMenu>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPainter>
#include <QToolTip>
#include <QtXml/QDomDocument>
#include <QXmlStreamWriter>

#include <cassert>

namespace
{
	const int CATEGORICAL_TEXT_ROTATION = 15;
	const int CATEGORICAL_FONT_SIZE = 7;

	const int X_AXIS_STEPS = 10;
	const int Y_AXIS_STEPS =  5;

class LinearConverter : public CoordinateConverter
{
public:
	LinearConverter(double yZoom, double yMax, int height)
	{
		LinearConverter::update(yZoom, yMax, 0, height);
	}
	double Diagram2ScreenY(double y) const override
	{
		return y * yScaleFactor;
	}
	double Screen2DiagramY(double y) const override
	{
		return y / yScaleFactor;
	}
	bool equals(QSharedPointer<CoordinateConverter> other) const override
	{
		LinearConverter* linearOther = dynamic_cast<LinearConverter*>(other.data());
		return (linearOther != 0 && yScaleFactor == linearOther->yScaleFactor);
	}
	QSharedPointer<CoordinateConverter> clone() override
	{
		return QSharedPointer<CoordinateConverter>(new LinearConverter(*this));
	}
	void update(double yZoom, double yMax, double yMinValueBiggerThanZero, int height) override
	{
		if (yMax)
			yScaleFactor = (double)(height-1) / yMax *yZoom;
		else
			yScaleFactor = 1;
	}
private:
	LinearConverter(LinearConverter const & other):
		yScaleFactor(other.yScaleFactor)
	{
	}
	double yScaleFactor;
};

class LogarithmicConverter : public CoordinateConverter
{
public:
	LogarithmicConverter(double yZoom, double yMax, double yMinValueBiggerThanZero, int height)
	{
		LogarithmicConverter::update(yZoom, yMax, yMinValueBiggerThanZero, height);
	}

	double Diagram2ScreenY(double y) const override
	{
		if (y <= 0)
			return 0;

		double yLog = LogFunc(y);

		yLog = clamp(yMinLog, yMaxLog, yLog);

		return mapValue(
			yMinLog, yMaxLog,
			0.0, static_cast<double>(height * yZoom),
			yLog
		);
	}
	double Screen2DiagramY(double y) const override
	{
		double yLog = mapValue(
			0.0, static_cast<double>(height * yZoom),
			yMinLog, yMaxLog,
			y
		);
		return std::pow(LogBase, yLog);
	}
	bool equals(QSharedPointer<CoordinateConverter> other) const  override
	{
		LogarithmicConverter* logOther = dynamic_cast<LogarithmicConverter*>(other.data());
		return (logOther && yZoom == logOther->yZoom &&
			yMaxLog == logOther->yMaxLog && yMinLog == logOther->yMinLog &&
			height == logOther->height);
	}
	QSharedPointer<CoordinateConverter> clone() override
	{
		return QSharedPointer<CoordinateConverter>(new LogarithmicConverter(*this));
	}
	void update(double yZoom, double yMax, double yMinValueBiggerThanZero, int height) override
	{
		this->yZoom = yZoom;
		yMaxLog = LogFunc(yMax);
		yMinLog = LogFunc(yMinValueBiggerThanZero)-1;
		this->height = height;
	}
private:
	LogarithmicConverter(LogarithmicConverter const & other):
		yZoom(other.yZoom),
		yMaxLog(other.yMaxLog),
		yMinLog(other.yMinLog),
		height(other.height)
	{}
	double yZoom;
	double yMaxLog, yMinLog;
	int height;
};
}

iADiagramFctWidget::iADiagramFctWidget(QWidget *parent, MdiChild *mdiChild,
	QString const & xLabel, QString const & yLabel) :
	iADiagramWidget(parent),
	TFTable(0),
	contextMenu(new QMenu(this)),
	xCaption(xLabel),
	yCaption(yLabel),
	m_yAxisSteps(Y_AXIS_STEPS),
	m_xAxisSteps(X_AXIS_STEPS),
	m_requiredPlacesAfterComma(0),
	m_yDrawMode(Linear),
	m_allowTrfReset(true),
	m_enableAdditionalFunctions(true),
	m_showXAxisLabel(true),
	m_captionPosition(Qt::AlignCenter | Qt::AlignBottom),
	m_maxYAxisValue(std::numeric_limits<iAPlotData::DataType>::lowest()),
	contextMenuVisible(false),
	m_showFunctions(false),
	m_customYAxisValue(false),
	selectedFunction(0),
	activeChild(mdiChild),
	updateAutomatically(true)
{
	dlg_transfer *transferFunction = new dlg_transfer(this, QColor(0, 0, 0, 255));
	functions.push_back(transferFunction);
	leftMargin   = (yLabel == "") ? 0 : 60;
}


iADiagramFctWidget::~iADiagramFctWidget()
{
	delete contextMenu;
	for (auto fct: functions)
		delete fct;
}


int iADiagramFctWidget::getSelectedFuncPoint() const
{
	auto it = functions.begin();
	if (it == functions.end())
		return -1;
	dlg_function *func = *(it + selectedFunction);
	return func->getSelectedPoint();
}

bool iADiagramFctWidget::isFuncEndPoint(int index) const
{
	auto it = functions.begin();
	if (it == functions.end())
		return false;
	dlg_function *func = *(it + selectedFunction);
	return func->isEndPoint(index);
}

bool iADiagramFctWidget::isUpdateAutomatically() const
{
	return updateAutomatically;
}

void iADiagramFctWidget::paintEvent(QPaintEvent * e)
{
	if (draw) drawEverything();

	QPainter painter(this);
	painter.drawImage(QRectF(0, 0, geometry().width(), geometry().height()), image);
}


void iADiagramFctWidget::CreateYConverter()
{
	if (m_maxYAxisValue == std::numeric_limits<iAPlotData::DataType>::lowest())
	{
		m_maxYAxisValue = GetMaxYValue();
	}
	if (m_yDrawMode == Linear)
	{
		m_yConverter = QSharedPointer<CoordinateConverter>(new LinearConverter(yZoom, m_maxYAxisValue, ActiveHeight() - 1));
	}
	else
	{																	// 1 - smallest value larger than 0. TODO: find that from data!
		m_yConverter = QSharedPointer<CoordinateConverter>(new LogarithmicConverter(yZoom, m_maxYAxisValue, 1, ActiveHeight() - 1));
	}
}


void iADiagramFctWidget::drawEverything()
{
	if (!m_yConverter)
	{
		CreateYConverter();
	}
	// TODO: update converter every time one of these values changes
	//		 alternative: give converter direct access to values (via some interface)
	m_yConverter->update(yZoom, m_maxYAxisValue, 1, ActiveHeight());

	QPainter painter(&image);
	painter.setRenderHint(QPainter::Antialiasing);

	drawBackground(painter);
	painter.translate(translationX+LeftMargin(), -BottomMargin());
	drawImageOverlays(painter);
	//change the origin of the window to left bottom
	painter.translate(0, height);
	painter.scale(1, -1);

	drawDatasets(painter);

	if (m_showFunctions)
	{
		drawFunctions(painter);
	}

	painter.scale(1, -1);
	painter.setRenderHint(QPainter::Antialiasing, false);

	drawAxes(painter);

	draw = false;
}

void iADiagramFctWidget::drawFunctions(QPainter &painter)
{
	int counter = 0;
	std::vector<dlg_function*>::iterator it = functions.begin();
	while (it != functions.end())
	{
		dlg_function *func = (*it);
		
		if (counter == selectedFunction)
		{
			func->draw(painter, QColor(255,128,0,255), 2);
		}
		else
		{
			func->draw(painter);
		}
		++it;
		counter++;
	}

	it = functions.begin();
	while (it != functions.end())
	{
		dlg_function *func = (*it);
		func->drawOnTop(painter);
		++it;
	}
}

void iADiagramFctWidget::redraw()
{
	draw = true;
	update();
}

void iADiagramFctWidget::mousePressEvent(QMouseEvent *event)  
{
	switch(event->button())
	{
		case Qt::LeftButton:
			if (((event->modifiers() & Qt::ShiftModifier) == Qt::ShiftModifier) &&
				((event->modifiers() & Qt::ControlModifier) == Qt::ControlModifier) &&
				((event->modifiers() & Qt::AltModifier) == Qt::AltModifier))
			{
				zoomY = event->y();
				changeMode(Y_ZOOM_MODE, event);
			}
			else if (((event->modifiers() & Qt::ShiftModifier) == Qt::ShiftModifier) &&
				((event->modifiers() & Qt::ControlModifier) == Qt::ControlModifier))
			{
				zoomX = event->x();
				zoomY = event->y();
				xZoomStart = xZoom;
				changeMode(X_ZOOM_MODE, event);
			}
			else if ((event->modifiers() & Qt::ShiftModifier) == Qt::ShiftModifier)
			{
				translationStartX = translationX;
				changeMode(MOVE_VIEW_MODE, event);
			}
			else if (!contextMenuVisible || contextMenuResult != NULL)
				changeMode(MOVE_POINT_MODE, event);
			break;
		case Qt::RightButton:
		{
			if (!m_showFunctions)
			{
				return;
			}
			std::vector<dlg_function*>::iterator it = functions.begin();
			dlg_function *func = *(it + selectedFunction);
			int selectedPoint = func->selectPoint(event);
			if (selectedPoint == -1)
				emit noPointSelected();
			redraw();
			break;
		}
		default:
			iADiagramWidget::mousePressEvent(event);
			break;
	}
}

void iADiagramFctWidget::mouseReleaseEvent(QMouseEvent *event)  
{
	if (!m_showFunctions)
	{
		return;
	}
	std::vector<dlg_function*>::iterator it = functions.begin();
	dlg_function *func = *(it + selectedFunction);
	switch (event->button())
	{
		case Qt::LeftButton:
			if (mode == MOVE_NEW_POINT_MODE)
				func->mouseReleaseEventAfterNewPoint(event);
			redraw();
			emit updateTFTable();

			if (isUpdateAutomatically())
				emit updateViews();
			break;
		default:
			break;
	}
	mode = NO_MODE;
	contextMenuVisible = false;
	func->mouseReleaseEvent(event);
}

void iADiagramFctWidget::mouseDoubleClickEvent(QMouseEvent *event)  
{
	if (event->button() == Qt::LeftButton)
	{
		changeColor(event);
	}
	emit DblClicked();
}

void iADiagramFctWidget::mouseMoveEvent(QMouseEvent *event)
{
	switch (mode)
	{
		case NO_MODE: /* do nothing */ break;
		case MOVE_POINT_MODE:
		case MOVE_NEW_POINT_MODE:
		{
			int viewX = event->x() - LeftMargin();
			int viewY = geometry().height() -event->y() -BottomMargin();
			if (m_showFunctions)
			{
				std::vector<dlg_function*>::iterator it = functions.begin();
				dlg_function *func = *(it + selectedFunction);
				func->moveSelectedPoint(viewX, viewY);
				redraw();
				emit updateTFTable();
			}
			showDataTooltip(event);
		}
		break;
		default:
			iADiagramWidget::mouseMoveEvent(event);
	}
}	
	
void iADiagramFctWidget::showDataTooltip(QMouseEvent *event)
{
	if (m_plots.empty())
		return;
	int xPos = clamp(0, width-1, event->x());
	size_t numBin = m_plots[0]->GetData()->GetNumBin();
	assert(numBin > 0);
	int nthBin = static_cast<int>((((xPos - translationX - LeftMargin()) * numBin) / (ActiveWidth())) / xZoom);
	nthBin = clamp(0, static_cast<int>(numBin), nthBin);
	if (xPos == width - 1)
		nthBin = static_cast<int>(numBin) - 1;
	QString toolTip(QString("%1: %2\n%3:").arg(xCaption).arg(m_plots[0]->GetData()->GetBinStart(nthBin)).arg(yCaption));
	for (auto plot : m_plots)
	{
		auto data = plot->GetData();
		if (!data || !data->GetRawData())
			continue;
		toolTip += QString::number(data->GetRawData()[nthBin]);
	}
	QToolTip::showText( event->globalPos(), toolTip, this);
}

void iADiagramFctWidget::enterEvent(QEvent*)
{
	emit active();
}

void iADiagramFctWidget::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Down)
	{
		selectedFunction++;
		if (selectedFunction >= functions.size())
			selectedFunction = 0;

		redraw();
	}
	else if (event->key() == Qt::Key_Up)
	{
		if(selectedFunction > 0)
			selectedFunction--;
		else
			selectedFunction = (int)(functions.size()-1);

		redraw();
	}
}

void iADiagramFctWidget::keyReleaseEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Alt || 
		event->key() == Qt::Key_AltGr ||
		event->key() == Qt::Key_Escape)
	contextMenuVisible = false;
}

void iADiagramFctWidget::contextMenuEvent(QContextMenuEvent *event)
{
	contextPos = event->pos();
	contextMenu->clear();
	
	if (m_showFunctions)
	{
		std::vector<dlg_function*>::iterator it = functions.begin();
		dlg_function *func = *(it + selectedFunction);

		if (func->getSelectedPoint() != -1)
		{
			if (func->isColored())
			{
				QAction *changeColorAction = new QAction(QIcon(":/images/changeColor.png"), tr("Change Color"), this);
				contextMenu->setDefaultAction(changeColorAction);
				connect(changeColorAction, SIGNAL(triggered()), this, SLOT(changeColor()));
				contextMenu->addAction(changeColorAction);
			}

			if (func->isDeletable(func->getSelectedPoint()))
				contextMenu->addAction(QIcon(":/images/deletePoint.png"), tr("Delete"), this, SLOT(deletePoint()));
			contextMenu->addSeparator();
		}
		contextMenu->addAction(QIcon(":/images/TFTableView.png"), tr("Transfer Function Table View"), this, SLOT(showTFTable()));
		contextMenu->addAction(QIcon(":/images/loadtrf.png"), tr("Load transfer function"), this, SLOT(loadTransferFunction()));
		contextMenu->addAction(QIcon(":/images/savetrf.png"), tr("Save transfer function"), this, SLOT(saveTransferFunction()));
		contextMenu->addAction(QIcon(":/images/savetrf.png"), tr("Apply transfer function for all"), this, SLOT(applyTransferFunctionForAll()));

		QAction *autoUpdateAction = new QAction(QIcon(":/images/autoUpdate.png"), tr("Update automatically"), this);
		autoUpdateAction->setCheckable(true);
		autoUpdateAction->setChecked(updateAutomatically);
		connect(autoUpdateAction, SIGNAL(toggled(bool)), this, SLOT(autoUpdate(bool)));
		contextMenu->addAction(autoUpdateAction);

		contextMenu->addAction(QIcon(":/images/update.png"), tr("Update views"), this, SIGNAL(updateViews()));
	}

	contextMenu->addAction(QIcon(":/images/save.png"), tr("Export data"), this, SLOT(ExportData()));

	if (m_showFunctions && m_allowTrfReset)
	{
		contextMenu->addAction(QIcon(":/images/resetTrf.png"), tr("Reset transfer function"), this, SLOT(resetTrf()));
	}
	contextMenu->addAction(QIcon(":/images/resetView.png"), tr("Reset view"), this, SLOT(resetView()));
	contextMenu->addSeparator();

	if (m_enableAdditionalFunctions)
	{
		contextMenu->addAction(QIcon(":/images/addBezier.png"), tr("Add bezier function"), this, SLOT(addBezierFunction()));
		contextMenu->addAction(QIcon(":/images/addGaussian.png"), tr("Add gaussian function"), this, SLOT(addGaussianFunction()));
		contextMenu->addAction(QIcon(":/images/openFkt.png"), tr("Load functions"), this, SLOT(loadFunctions()));
		contextMenu->addAction(QIcon(":/images/saveFkt.png"), tr("Save functions"), this, SLOT(saveFunctions()));

		if (selectedFunction != 0)
			contextMenu->addAction(QIcon(":/images/removeFkt.png"), tr("Remove selected function"), this, SLOT(removeFunction()));
	}

	contextMenuVisible = true;
	contextMenuResult = contextMenu->exec(event->globalPos());
}

QSharedPointer<CoordinateConverter> const iADiagramFctWidget::GetCoordinateConverter() const
{
	return m_yConverter;
}

void iADiagramFctWidget::drawDatasets(QPainter &painter)
{
	if (m_plots.empty())
	{
		painter.scale(1, -1);
		painter.drawText(QRect(0, 0, ActiveWidth(), -ActiveHeight()), Qt::AlignCenter, "Chart not (yet) available.");
		painter.scale(1, -1);
		return;
	}
	double binWidth = ActiveWidth() * xZoom / m_plots[0]->GetData()->GetNumBin();
	for (auto it = m_plots.constBegin(); it != m_plots.constEnd();	++it)
	{
		if ((*it)->Visible())
			(*it)->draw(painter, binWidth, m_yConverter);
	}
}


void iADiagramFctWidget::drawImageOverlays( QPainter &painter )
{
	QRect targetRect = geometry();
	int yTranslate = -(yZoom-1) * (targetRect.height());
	targetRect.setHeight(targetRect.height()-targetRect.top()-1);
	targetRect.setWidth( (targetRect.width() - LeftMargin()) * xZoom);
	targetRect.setTop(targetRect.top() + yTranslate);
	targetRect.setLeft(0);
	for (int i = 0; i < m_overlays.size(); ++i)
	{
		painter.drawImage(targetRect, *(m_overlays[i].data()), m_overlays[i]->rect());
	}
}

void iADiagramFctWidget::drawAxes(QPainter &painter)
{
	drawXAxis(painter);
	drawYAxis(painter);
}

int requiredDigits(double value)
{
	return (value >= -1.0 && value < 1.0) ?
		1 :	std::floor(std::log10(std::abs(value))) + 1;
}

double deg2rad(double const & number)
{
	return number * M_PI / 180;
}

int markerPos(int step, int stepNr, int width, bool binMiddle, size_t numBin)
{
	int x = static_cast<int>( static_cast<double>(step)/stepNr * width);
	if (binMiddle)
	{
		// for discrete x axis variables, the marker and caption should be "in the middle" of the bar
		x += static_cast<int>(0.5 * width/numBin);
	}
	if (step == stepNr) x--;
	return x;
}

int textPos(int markerX, int step, int stepNr, int textWidth)
{
	return (step == 0)
			? markerX					// right aligned to indicator line
			: (step < stepNr)
				? markerX-textWidth/2	// centered to the indicator line
				: markerX-textWidth;	// left aligned to the indicator line
}

void iADiagramFctWidget::drawXAxis(QPainter &painter)
{
	painter.setPen(QWidget::palette().color(QPalette::Text));

	const int MINIMUM_MARGIN = 8;
	const int TextAxisDistance = 2;
	QFontMetrics fm = painter.fontMetrics();
	// TODO: find better solution than just use first plot here:
	if (m_plots.empty())
		return;
	auto data = m_plots[0]->GetData();
	
	int stepNumber = static_cast<int>(data->GetNumBin());
	size_t numBin = data->GetNumBin();
	int stepSize = 1;
	
	if (data->GetRangeType() != Categorical)
	{
		bool overlap;
		do
		{
			overlap =  false;
			for (int i=0; i<stepNumber && !overlap; ++i)
			{
				int nthBin = static_cast<int>(static_cast<double>(i) / stepNumber * numBin);
				double value = data->GetBinStart(nthBin);
				int placesBeforeComma = requiredDigits(value);
				double stepToNextBin = data->GetBinStart(i < numBin -1? i+1: i)
					- data->GetBinStart(i < numBin -1? i: i-1);
				m_requiredPlacesAfterComma = (stepToNextBin < 10) ? requiredDigits(10 / stepToNextBin) : 0;
				QString text = GetXAxisCaption(value, placesBeforeComma, m_requiredPlacesAfterComma);

				int markerX = markerPos(i, stepNumber, ActiveWidth()*xZoom, IsDrawnDiscrete(), numBin);
				int textX = textPos(markerX, i, stepNumber, fm.width(text));
				int next_markerX = markerPos(i+1, stepNumber, ActiveWidth()*xZoom, IsDrawnDiscrete(), numBin);
				int next_textX = textPos(next_markerX, i+1, stepNumber, fm.width(text));
				int textWidth = fm.width(text) + fm.width("M");
				overlap = (textX + textWidth) >= next_textX;
			}
			if (overlap)
			{
				do
				{
					++stepSize;
				}
				while (stepSize < numBin && numBin % stepSize != 0);
				stepNumber = static_cast<int>(numBin / stepSize);
			}
		} while (overlap && stepSize < numBin);
	}
	else
	{
		QFont font = painter.font() ;
		font.setPointSize(CATEGORICAL_FONT_SIZE);
		painter.setFont(font);
		fm = painter.fontMetrics();
	}

	stepNumber = std::max(1, stepNumber); // at least two steps
	for (int i=0; i <= stepNumber; ++i)
	{
		int nthBin = (int)(static_cast<double>(i) / stepNumber * numBin);
		double value = data->GetBinStart(nthBin);
		int placesBeforeComma = requiredDigits(value);
		double stepToNextBin = data->GetBinStart(nthBin < numBin-1? nthBin+1: nthBin)
			- data->GetBinStart(nthBin < numBin -1? nthBin: nthBin-1);
		m_requiredPlacesAfterComma = (stepToNextBin < 10) ? requiredDigits(10 / stepToNextBin) : 0;

		QString text = GetXAxisCaption(value, placesBeforeComma, m_requiredPlacesAfterComma);

		// dirty hack to avoid last tick for discrete ranges:
		if (data->GetRangeType() == Discrete && i == stepNumber && text.length() < 3)
		{
			break;
		}
		int markerX = markerPos(i, stepNumber, ActiveWidth()*xZoom, IsDrawnDiscrete(), numBin);
		painter.drawLine(markerX, (int)(BottomMargin()*0.1), markerX, -1);

		int textX = textPos(markerX, i, stepNumber, fm.width(text));
		int textY = fm.height() + TextAxisDistance;
		painter.translate(textX, textY);
		painter.drawText(0, 0, text);
		painter.translate(-textX, -textY);
	}
	
	//draw the x axis
	painter.setPen(QWidget::palette().color(QPalette::Text));
	painter.drawLine(0, -1, (int)((ActiveWidth())*xZoom), -1);
	
	if (m_showXAxisLabel)
	{
		//write the x axis label
		QPointF textPos(
			m_captionPosition.testFlag(Qt::AlignCenter) ? (int)(ActiveWidth() * 0.45 - translationX): 0 /* left-aligned */ ,
			m_captionPosition.testFlag(Qt::AlignBottom) ? BottomMargin()-fm.descent()-1 : -Height() + BottomMargin() + fm.height()
		);
		painter.drawText(textPos, xCaption);
	}
}

void iADiagramFctWidget::drawYAxis(QPainter &painter)
{
	if ( LeftMargin() <= 0 )
	{
		return;
	}
	painter.save();
	painter.translate(-translationX, 0);
	QFontMetrics fm = painter.fontMetrics();
	int fontHeight = fm.height();

	int activeHeight = ActiveHeight()-1;
	painter.fillRect(QRect(0, BottomMargin(), -LeftMargin(), -(activeHeight+BottomMargin()+1)),
		QBrush(QWidget::palette().color(QWidget::backgroundRole())));
	painter.setPen(QWidget::palette().color(QPalette::Text));

	// at most, make Y_AXIS_STEPS, but reduce to number actually fitting in current height:
	int stepNumber = std::min(Y_AXIS_STEPS, static_cast<int>(activeHeight / (fontHeight*1.1)) );
	stepNumber = std::max(1, stepNumber);	// make sure there's at least 2 steps
	const double step = 1.0 / (stepNumber * yZoom);
	double logMax = LogFunc(static_cast<double>(m_maxYAxisValue));

	for (int i = 0; i <= stepNumber; ++i)
	{
		//calculate the nth bin located at a given pixel, actual formula is (i/100 * width) * (rayLength / width)
		double pos = step * i;
		//get the intensity value into a string

		double yValue =
			(m_yDrawMode == Linear) ?
			pos * m_maxYAxisValue :
			/* Logarithmic: */
			std::pow(LogBase, logMax/yZoom - (Y_AXIS_STEPS - i));

		QString text;
		if (yValue < 1.0)
			text = QString::number(yValue, 'g', 3);
		else
			if (yValue > 1000000000)
			{
				text = QString::number(yValue/1000000000.0, 'f', (yValue < 10000000000) ? 2 : 1)+"G";
			}
			else if (yValue > 1000000)
			{
				text = QString::number(yValue / 1000000, 'f', (yValue < 10000000) ? 2 : 1) + "M";
			}
			else if (yValue > 1000)
			{
				text = QString::number(yValue / 1000, 'f', (yValue < 10000) ? 2 : 1) + "K";
			}
			else
			{
				text = QString::number((int)yValue, 10);
			}

		//calculate the y coordinate
		int y = -(int)(pos * activeHeight * yZoom)-1;
		//draw a small indicator line
		painter.drawLine((int)(-LeftMargin()*0.1), y, 0, y);

		if(i == stepNumber)
			painter.drawText(TEXT_X-LeftMargin(), y+0.75*fontHeight, text); //write the text top aligned to the indicator line
		else
			painter.drawText(TEXT_X-LeftMargin(), y+0.25*fontHeight, text); //write the text centered to the indicator line

	}
	painter.drawLine(0, -1, 0, -(int)(activeHeight*yZoom));
	//write the y axis label
	painter.save();
	painter.rotate(-90);
	QPointF textPos(
		activeHeight*0.5 - 0.5*fm.width(yCaption),
		-LeftMargin() + fontHeight - 5);
	painter.drawText(textPos, yCaption);
	painter.restore();
	painter.restore();
}

int iADiagramFctWidget::BottomMargin() const
{
	if (!m_showXAxisLabel)
	{
		return BOTTOM_MARGIN/2.0 /* TODO: estimation for font height only. use real values! */;
	}
	return BOTTOM_MARGIN;
}

void iADiagramFctWidget::changeMode(int newMode, QMouseEvent *event)
{
	switch(newMode)
	{
		case MOVE_POINT_MODE:
		{
			if (!m_showFunctions)
			{
				return;
			}
			std::vector<dlg_function*>::iterator it = functions.begin();
			dlg_function *func = *(it + selectedFunction);
			int x = event->x() - LeftMargin();
			int y = geometry().height() - event->y() -BottomMargin();
			int selectedPoint = func->selectPoint(event, &x);

			// don't do anything if outside of diagram region:
			if (selectedPoint == -1 && x < 0)
			{
				return;
			}
			// disallow removal and reinsertion of first point; instead, insert a point after it:
			if (selectedPoint == -1 && x == 0)
			{
				x = 1;
			}
			if (selectedPoint == -1)
			{
				if (y < 0) y = 0;
				selectedPoint = func->addPoint(x, y);
				func->addColorPoint(x);
				
				mode = MOVE_NEW_POINT_MODE;
			}
			else
				mode = MOVE_POINT_MODE;

			redraw();

			bool endPoint = func->isEndPoint(selectedPoint);//selectedPoint == 0 || selectedPoint == opacityTF->GetSize()-1;
			if (endPoint)
				emit endPointSelected();
			else
				emit pointSelected();
		}
			break;
		default:
			iADiagramWidget::changeMode(newMode, event);
			break;
	}
}

int iADiagramFctWidget::deletePoint()
{
	std::vector<dlg_function*>::iterator it = functions.begin();
	dlg_function *func = *(it + selectedFunction);

	int selectedPoint = func->getSelectedPoint();
	func->removePoint(selectedPoint);
	redraw();
	emit updateTFTable();
	emit updateViews();

	return selectedPoint;
}

void iADiagramFctWidget::changeColor(QMouseEvent *event)
{
	if (!m_showFunctions)
	{
		return;
	}
	std::vector<dlg_function*>::iterator it = functions.begin();
	dlg_function *func = *(it + selectedFunction);

	func->changeColor(event);

	redraw();
	emit updateTFTable();
	emit updateViews();
}

void iADiagramFctWidget::autoUpdate(bool toggled)
{
	updateAutomatically = toggled;
	emit autoUpdateChanged(toggled);
}

void iADiagramFctWidget::resetView()
{
	iADiagramWidget::resetView();
}

void iADiagramFctWidget::resetTrf()
{
	std::vector<dlg_function*>::iterator it = functions.begin();
	dlg_function *func = *(it + selectedFunction);

	func->reset();
	redraw();
	emit updateTFTable();
	emit updateViews();
}

void iADiagramFctWidget::updateTrf()
{
	((dlg_transfer*)functions[0])->TranslateToNewRange(XBounds());
	redraw();
}

bool iADiagramFctWidget::loadTransferFunction()
{
	QString filePath = (activeChild) ? activeChild->getFilePath() : "";
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), filePath ,tr("XML (*.xml)"));
	if (!fileName.isEmpty())
	{
		Settings s(fileName);
		s.LoadTransferFunction((dlg_transfer*)functions[0], XBounds());

		emit noPointSelected();

		redraw();
		emit updateTFTable();
		emit updateViews();
	}
	return true;
}

void iADiagramFctWidget::loadTransferFunction(QDomNode &functionsNode)
{
	((dlg_transfer*)functions[0])->loadTransferFunction(functionsNode, XBounds());
}

bool iADiagramFctWidget::saveTransferFunction()
{
	dlg_transfer *transferFunction = (dlg_transfer*)functions[0];
	QString filePath = (activeChild) ? activeChild->getFilePath() : "";
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"), filePath ,tr("XML (*.xml)"));
	if (!fileName.isEmpty()) {
		Settings s;
		s.StoreTransferFunction(transferFunction);
		s.Save(fileName);
	}

	return true;
}

void iADiagramFctWidget::applyTransferFunctionForAll()
{
	emit applyTFForAll();
}

void iADiagramFctWidget::addBezierFunction()
{
	dlg_bezier *bezier = new dlg_bezier(this, PredefinedColors()[functions.size() % 7]);

	bezier->addPoint(contextPos.x(), ActiveHeight()-contextPos.y());

	selectedFunction = (unsigned int)functions.size();
	functions.push_back(bezier);

	redraw();
	emit updateTFTable();
	emit updateViews();
}

void iADiagramFctWidget::addGaussianFunction()
{
	dlg_gaussian *gaussian = new dlg_gaussian(this, PredefinedColors()[functions.size() % 7]);

	gaussian->setMean(contextPos.x());
	gaussian->setSigma(width/6);
	gaussian->setMultiplier((int)((ActiveHeight()-contextPos.y())*YZoom()));
	
	selectedFunction = (unsigned int)functions.size();
	functions.push_back(gaussian);

	redraw();
	emit updateTFTable();
	emit updateViews();
}

bool iADiagramFctWidget::loadFunctions()
{
	if (!activeChild)
	{
		return false;
	}
	QString filePath = activeChild->getFilePath();
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), filePath ,tr("XML (*.xml)"));
	if (!fileName.isEmpty())
	{
		MdiChild* child = activeChild;
		MainWindow *mw = (MainWindow*)child->window();
		QDomDocument doc = mw->loadSettingsFile(fileName);
		QDomElement root = doc.documentElement();

		for (unsigned int i = 1; i < functions.size(); i++)
		{
			delete functions.back();
			functions.pop_back();
		}

		QDomNode functionsNode = root.namedItem("functions");
		if (functionsNode.isElement())
				mw->loadProbabilityFunctions(functionsNode);
		
		emit noPointSelected();
		
		redraw();

		emit updateViews();
	}

	return true;
}

bool iADiagramFctWidget::saveFunctions()
{
	if (!activeChild)
	{
		return false;
	}
	QString filePath = activeChild->getFilePath();
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"), filePath ,tr("XML (*.xml)"));
	if (!fileName.isEmpty())
	{
		MdiChild* child = activeChild;
		MainWindow *mw = (MainWindow*)child->window();

		QDomDocument doc = mw->loadSettingsFile(fileName);
		mw->saveProbabilityFunctions(doc);
		mw->saveSettingsFile(doc, fileName);
	}

	return true;
}

void iADiagramFctWidget::removeFunction()
{
	std::vector<dlg_function*>::iterator it = functions.begin()+selectedFunction;
	dlg_function *function = *it;
	functions.erase(it);

	delete function;

	selectedFunction--;

	redraw();

	emit updateViews();
}

void iADiagramFctWidget::SetTransferFunctions(vtkColorTransferFunction* ctf, vtkPiecewiseFunction* pwf)
{
	m_showFunctions = ctf && pwf;
	if (!m_showFunctions)
		return;
	((dlg_transfer*)functions[0])->setColorFunction(ctf);
	((dlg_transfer*)functions[0])->setOpacityFunction(pwf);
	redraw();
	emit updateTFTable();
	emit updateViews();
}

void iADiagramFctWidget::ExportData()
{
	if (m_plots.empty())
		return;
	QString filePath = (activeChild) ? activeChild->getFilePath() : "";
	QString fileName = QFileDialog::getSaveFileName(
		this,
		tr("Save File"),
		filePath, tr("CSV (*.csv)")
	);
	if (fileName.isEmpty())
	{
		return;
	}
	std::ofstream out(fileName.toStdString());
	out << tr("Start of Bin").toStdString();
	for (int p = 0; p < m_plots.size(); ++p)
	{
		out << "," << QString("%1%2").arg(yCaption).arg(p).toStdString();
	}
	out << endl;
	for (int b = 0; b < m_plots[0]->GetData()->GetNumBin(); ++b)
	{
		out << m_plots[0]->GetData()->GetBinStart(b);
		for (int p = 0; p < m_plots.size(); ++p)
		{
			out << "," << m_plots[p]->GetData()->GetRawData()[b];
		}
		out << std::endl;
	}
	out.close();
}

double const * iADiagramFctWidget::XBounds() const
{
	// TODO: use all plots here?
	assert(!m_plots.empty());
	return m_plots[0]->GetData()->XBounds();
}

double iADiagramFctWidget::XRange() const
{
	// TODO: use all plots here?
	assert(!m_plots.empty());
	return m_plots[0]->GetData()->XBounds()[1] - m_plots[0]->GetData()->XBounds()[0];
}

dlg_function *iADiagramFctWidget::getSelectedFunction()
{
	return functions[selectedFunction];
}

int iADiagramFctWidget::ChartHeight() const
{
	return height - BottomMargin();
}

std::vector<dlg_function*> &iADiagramFctWidget::getFunctions()
{
	return functions;
}

iAPlotData::DataType iADiagramFctWidget::GetMaxYValue() const
{
	iAPlotData::DataType maxVal = std::numeric_limits<iAPlotData::DataType>::lowest();
	for (auto plot : m_plots)
		maxVal = std::max(plot->GetData()->YBounds()[1], maxVal);
	return maxVal;
}

iAPlotData::DataType iADiagramFctWidget::GetMaxYAxisValue() const
{
	return m_maxYAxisValue;
}

void iADiagramFctWidget::SetMaxYAxisValue(iAPlotData::DataType val)
{
	m_customYAxisValue = true;
	m_maxYAxisValue = val;
	for (auto it = m_plots.constBegin(); it != m_plots.constEnd(); ++it)
		(*it)->update();
}

void iADiagramFctWidget::ResetMaxYAxisValue()
{
	m_customYAxisValue = false;
	m_maxYAxisValue = GetMaxYValue();
}

void iADiagramFctWidget::AddPlot(QSharedPointer<iAPlot> plot)
{
	assert(plot);
	if (!plot)
		return;
	m_plots.push_back(plot);
	if (!m_customYAxisValue)
		m_maxYAxisValue = GetMaxYValue();
}

void iADiagramFctWidget::RemovePlot(QSharedPointer<iAPlot> plot)
{
	if (!plot)
		return;
	int idx = m_plots.indexOf(plot);
	if (idx != -1)
		m_plots.remove(idx);
}

QVector<QSharedPointer<iAPlot> > const & iADiagramFctWidget::Plots()
{
	return m_plots;
}

void iADiagramFctWidget::AddImageOverlay( QSharedPointer<QImage> imgOverlay )
{
	m_overlays.push_back(imgOverlay);
}

void iADiagramFctWidget::RemoveImageOverlay( QImage * imgOverlay )
{
	for (int i = 0; i < m_overlays.size(); ++i)
	{
		if( m_overlays.at(i).data() == imgOverlay)
		{
			m_overlays.removeAt(i);
			return;
		}
	}
}

int iADiagramFctWidget::diagram2PaintX(double x)
{
	if (m_plots.empty())
		return x;
	double screenX = (x - XBounds()[0]) * ActiveWidth() * xZoom / XRange();
	screenX = clamp(0.0, ActiveWidth()*xZoom, screenX);
	return static_cast<int>(round(screenX));
}

long iADiagramFctWidget::screenX2DataBin(int x)
{
	if (m_plots.empty())
		return x;
	double numBin = m_plots[0]->GetData()->GetNumBin();
	double diagX = static_cast<double>(x-translationX-LeftMargin()) * numBin / (ActiveWidth() * xZoom);
	diagX = clamp(0.0, numBin, diagX);
	return static_cast<long>(round(diagX));
}

int iADiagramFctWidget::dataBin2ScreenX(long x)
{
	if (m_plots.empty())
		return x;
	double numBin = m_plots[0]->GetData()->GetNumBin();
	double screenX = static_cast<double>(x) * ActiveWidth() * xZoom / (numBin);
	screenX = clamp(0.0, ActiveWidth()*xZoom, screenX);
	return static_cast<int>(round(screenX));
}

double iADiagramFctWidget::getMaxXZoom() const
{
	if (m_plots.empty())
		return 1;
	double numBin = m_plots[0]->GetData()->GetNumBin();
	return (std::max)((std::min)( iADiagramWidget::getMaxXZoom(), numBin), 1.0);
}

void iADiagramFctWidget::SetYDrawMode(DrawModeType drawMode)
{
	if (m_yDrawMode == drawMode)
		return;
	m_yDrawMode = drawMode;
	CreateYConverter();
}

void iADiagramFctWidget::SetXAxisSteps(int xSteps)
{
	m_xAxisSteps = xSteps;
}

void iADiagramFctWidget::SetYAxisSteps( int ySteps )
{
	m_yAxisSteps = ySteps;
}

void iADiagramFctWidget::SetRequiredPlacesAfterComma( int requiredPlaces )
{
	m_requiredPlacesAfterComma = requiredPlaces;
}

void iADiagramFctWidget::SetAllowTrfReset(bool allow)
{
	m_allowTrfReset = allow;
}

void iADiagramFctWidget::SetEnableAdditionalFunctions(bool enable)
{
	m_enableAdditionalFunctions = enable;
}

int iADiagramFctWidget::GetTFGradientHeight() const
{
	return BottomMargin();
}

QString iADiagramFctWidget::GetXAxisCaption(double value, int placesBeforeComma, int requiredPlacesAfterComma)
{
	if (!m_plots.empty() && m_plots[0]->GetData()->GetRangeType() == Continuous && requiredPlacesAfterComma > 0 )
	{
		QString result =  QString::number(value, 'g', ((value > 0) ? placesBeforeComma + requiredPlacesAfterComma : requiredPlacesAfterComma ));
		if (result.contains("e")) // only 2 digits for scientific notation:
		{
			result =  QString::number(value, 'g', 2);
		}
		return result;
	} else
		return QString::number((int)value, 10);
}

bool iADiagramFctWidget::IsDrawnDiscrete() const
{
	return (!m_plots.empty() && (
		(m_plots[0]->GetData()->GetRangeType() == Discrete &&
		(XRange() <= m_plots[0]->GetData()->GetNumBin()))
		|| m_plots[0]->GetData()->GetRangeType() == Categorical));
}

void iADiagramFctWidget::SetXCaption(QString const & caption)
{
	xCaption = caption;
}

void iADiagramFctWidget::showTFTable()
{
	if ( !TFTable )
	{
		std::vector<dlg_function*>::iterator it = functions.begin();
		dlg_function* func = *( it + selectedFunction );

		TFTable = new dlg_TFTable( this, func );
		TFTable->setWindowTitle( "Transfer Function Table View" );
		TFTable->setWindowFlags( Qt::Dialog |Qt::WindowMinimizeButtonHint |Qt::WindowCloseButtonHint );
		TFTable->setAttribute( Qt::WA_DeleteOnClose, true);
		TFTable->show();

		connect( TFTable, SIGNAL( destroyed( QObject* ) ), this, SLOT( TFTableIsFinished() ) );
		connect( this, SIGNAL( updateTFTable() ), TFTable, SLOT( updateTable() ) );
	}
}

QPoint iADiagramFctWidget::getTFTablePos()
{
	return TFTable->pos();
}

void iADiagramFctWidget::setTFTablePos( QPoint pos )
{
	TFTable->move( pos );
}

bool iADiagramFctWidget::isTFTableCreated()
{
	bool created;
	TFTable ? created = true : created = false;
	return created;
}

void iADiagramFctWidget::closeTFTable()
{
	TFTable->close();
}

void iADiagramFctWidget::TFTableIsFinished()
{
	TFTable = NULL;
}
