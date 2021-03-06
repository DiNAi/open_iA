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

#include "iAColorable.h"
#include "open_iA_Core_export.h"

#include <QSharedPointer>

class iAPlotData;
class iAMapper;

class QColor;
class QPainter;

/** Interface for a function which is drawable in a diagram
 *	encapsulates both the data of the function and the drawing method */
class open_iA_Core_API iAPlot: public iAColorable
{
public:
	iAPlot(QColor const & color);
	virtual ~iAPlot();
	/** in case your drawer implements caching, clear the cache when this is called */
	virtual void update();
	/**
	* \brief method which draws the function
	* it is allowed to cache the result; when the data has changed, update() needs to be called
	*/
	virtual void draw(QPainter& painter, double binWidth, QSharedPointer<iAMapper> converter) const =0;
	/** retrieves the data used for drawing */
	virtual QSharedPointer<iAPlotData> GetData();
	virtual bool Visible() const;
	virtual void SetVisible(bool visible);
private:
	bool m_visible;
};
