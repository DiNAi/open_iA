﻿/*********************************  open_iA 2016 06  ******************************** *
* **********  A tool for scientific visualisation and 3D image processing  ********** *
* *********************************************************************************** *
* Copyright (C) 2016  C. Heinzl, M. Reiter, A. Reh, W. Li, M. Arikan, J. Weissenböck, *
*                     Artem & Alexander Amirkhanov, B. Fröhler                        *
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

#include "iAFoamCharacterizationItemFilter.h"

#include <QApplication>
#include <QPainter>

#include "iAFoamCharacterizationDialogFilter.h"

iAFoamCharacterizationItemFilter::iAFoamCharacterizationItemFilter(vtkImageData* _pImageData)
	                                                            : iAFoamCharacterizationItem(iAFoamCharacterizationItem::itFilter)
																, m_pImageData (_pImageData)
{
	setItemIcon();

	setText("Filter");
}

iAFoamCharacterizationItemFilter::iAFoamCharacterizationItemFilter(iAFoamCharacterizationItemFilter* _pFilter)
	                                                            : iAFoamCharacterizationItem(iAFoamCharacterizationItem::itFilter)
{
	setItemIcon();

	setText(_pFilter->text());

	m_pImageData = _pFilter->imageData();
}

void iAFoamCharacterizationItemFilter::dialog()
{
	QScopedPointer<iAFoamCharacterizationDialogFilter> pDialog(new iAFoamCharacterizationDialogFilter(this, qApp->focusWidget()));
	pDialog->exec();
}

void iAFoamCharacterizationItemFilter::execute()
{

}

vtkImageData* iAFoamCharacterizationItemFilter::imageData() const
{
	return m_pImageData;
}

void iAFoamCharacterizationItemFilter::setItemIcon()
{
	const int iImageLength(qMax(font().pixelSize(), font().pointSize()));

	QScopedPointer<QImage> pImage(new QImage(iImageLength, iImageLength, QImage::Format_ARGB32));
	pImage->fill(0);

	QScopedPointer<QPainter> pPainter(new QPainter(pImage.data()));
	pPainter->setBrush((m_bItemEnabled) ? QBrush(Qt::red) : Qt::NoBrush);
	pPainter->setPen(Qt::red);
	pPainter->drawEllipse(pImage->rect().adjusted(0, 0, -1, -1));

	setIcon(QIcon(QPixmap::fromImage(*pImage.data())));
}
