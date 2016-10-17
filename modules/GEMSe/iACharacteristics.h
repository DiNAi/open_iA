/*********************************  open_iA 2016 06  ******************************** *
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
*          Stelzhamerstraße 23, 4600 Wels / Austria, Email:                           *
* ************************************************************************************/
 
#ifndef IACHARACTERISTICS_H
#define IACHARACTERISTICS_H

#include <QSharedPointer>
#include <QThread>
#include <time.h>
#include <itkCastImageFilter.h>
#include <itkLabelOverlapMeasuresImageFilter.h>
#include <QDir>

#include "iAGEMSeConstants.h"

class iAAttributes;
class iASingleResult;

class QString;

class iACharacteristics
{
public:
	iACharacteristics();
	int ObjectCount() const;
	void SetObjectCount(int objCount);
	static iACharacteristics Create(QString const & Descriptor);
	QString GetDescriptor() const;

	void SetDuration(double value);
	//! retrieve duration of ERW calculation
	double Duration() const;

	void SetMeasures(double value[MeasureCount]);
	double Measure(int idx) const;
private:
	int m_objectCount;
	double m_duration;
	double m_measures[MeasureCount];

	double m_median_accuracy, m_mean_accuracy, m_dice_overall, m_dice_0, m_dice_1, m_dice_2, m_dice_3, m_dice_4;
	double m_mean_uncertainty, m_median_uncertainty, m_uncertainty_0, m_uncertainty_1, m_uncertainty_2, m_uncertainty_3, m_uncertainty_4;
	double m_median_confusion_uncertainy_false, m_median_confusion_uncertainy_true, m_median_confusion_bad_to_good;
};

class CharacteristicsCalculator : public QThread
{
public:
	CharacteristicsCalculator(QSharedPointer<iASingleResult> result, QSharedPointer<iAAttributes> range, int objCountIdx);
	void GetUncertaintyValues(QSharedPointer<iASingleResult> result, QString groundTruthPath, std::vector< double > &uncert);
	double vecMed(std::vector<double> vec);
	double vecMAD(std::vector<double> vec);
private:
	QSharedPointer<iASingleResult> m_result;
	QSharedPointer<iAAttributes>   m_range;
	int m_objCountIdx;
	virtual void run();

};

#endif // IACHARACTERISTICS_H