/*********************************  open_iA 2016 06  ******************************** *
* **********  A tool for scientific visualisation and 3D image processing  ********** *
* *********************************************************************************** *
* Copyright (C) 2016  C. Heinzl, M. Reiter, A. Reh, W. Li, M. Arikan, J. Weissenb�ck, *
*                     Artem & Alexander Amirkhanov, B. Fr�hler                        *
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
* Contact: FH O� Forschungs & Entwicklungs GmbH, Campus Wels, CT-Gruppe,              *
*          Stelzhamerstra�e 23, 4600 Wels / Austria, Email: c.heinzl@fh-wels.at       *
* ************************************************************************************/
 

/*
	Autor: AMA2
*/

#pragma once
#ifndef dlg_blobVisualization_h
#define dlg_blobVisualization_h

#include <QDialog>
#include <QDialogButtonBox>

class QCheckBox;
class QDoubleSpinBox;
class QSpinBox;
class QLabel;

class dlg_blobVisualization : public QDialog
{
	Q_OBJECT
public:
	dlg_blobVisualization (QWidget* parent = 0);
	~dlg_blobVisualization (void);

	bool	GetOverlapping (void) const;
	double	GetRange (void) const ;
	double	GetOverlapThreshold (void) const;
	double	GetBlurVariance (void) const;
	int		GetResolution(void) const;

private:
	QDialogButtonBox*	m_buttons;
	QLabel*	m_labelOverlapping;
	QLabel*	m_labelRange;
	QLabel*	m_labelOverlapThreshold;
	QLabel* m_labelVariance;
	QLabel* m_labelResolution;
	QDoubleSpinBox*	m_spinBoxRange;
	QDoubleSpinBox*	m_gaussianBlurVariance;
	QDoubleSpinBox*	m_spinBoxOverlapThreshold;
	QSpinBox*	m_spinBoxResolution;
	QCheckBox*	m_chekBoxOverlapping;
};

#endif // dlg_blobVisualization_h
