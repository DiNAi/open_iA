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
#include <QFile>
#include <QTime>

#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkPatchBasedDenoisingImageFilter.h>

#include "iAConnector.h"

#include "iAFoamCharacterizationDialogFilter.h"

iAFoamCharacterizationItemFilter::iAFoamCharacterizationItemFilter(vtkImageData* _pImageData)
	                                               : iAFoamCharacterizationItem(_pImageData, iAFoamCharacterizationItem::itFilter)
{

}

iAFoamCharacterizationItemFilter::iAFoamCharacterizationItemFilter(iAFoamCharacterizationItemFilter* _pFilter)
	                                                                                        : iAFoamCharacterizationItem(_pFilter)
{
	m_eItemFilterType = _pFilter->itemFilterType();

	setName(_pFilter->name());

	m_dAnisotropicConductance = _pFilter->anisotropicConductance();
	m_uiAnisotropicIteration = _pFilter->anisotropicIteration();
	m_dAnisotropicTimeStep = _pFilter->anisotropicTimeStep();

	m_dGaussianVariance = _pFilter->gaussianVariance();

	m_uiMedianRadius = _pFilter->medianRadius();
	m_uiNonLocalMeansRadius = _pFilter->nonLocalMeansRadius();
}

double iAFoamCharacterizationItemFilter::anisotropicConductance() const
{
	return m_dAnisotropicConductance;
}

unsigned int iAFoamCharacterizationItemFilter::anisotropicIteration() const
{
	return m_uiAnisotropicIteration;
}

double iAFoamCharacterizationItemFilter::anisotropicTimeStep() const
{
	return m_dAnisotropicTimeStep;
}

void iAFoamCharacterizationItemFilter::dialog()
{
	QScopedPointer<iAFoamCharacterizationDialogFilter> pDialog(new iAFoamCharacterizationDialogFilter(this, qApp->focusWidget()));
	pDialog->exec();
	pDialog.reset();
}

void iAFoamCharacterizationItemFilter::execute()
{
	QTime t;
	t.start();

	switch (m_eItemFilterType)
	{
		case iftAnisotropic:
		executeAnisotropic();
		break;

		case iftGauss:
		executeGaussian();
		break;

		case iftMedian:
		executeMedian();
		break;

		default:
		executeNonLocalMeans();
		break;
	}

	setTime(t.elapsed());
}

void iAFoamCharacterizationItemFilter::executeAnisotropic()
{
	typedef itk::GradientAnisotropicDiffusionImageFilter <itk::Image<unsigned short, 3>, itk::Image<float, 3>> itkFilter;

	itkFilter::Pointer pFilter(itkFilter::New());

	iAConnector connector1;
	connector1.SetImage(m_pImageData);

	pFilter->SetInput(dynamic_cast<itk::Image<unsigned short, 3>*> (connector1.GetITKImage()));
	pFilter->SetConductanceParameter(m_dAnisotropicConductance);
	pFilter->SetNumberOfIterations(m_uiAnisotropicIteration);
	pFilter->SetTimeStep(m_dAnisotropicTimeStep);
	pFilter->Update();

	iAConnector connector2;
	connector2.SetImage(pFilter->GetOutput());

	m_pImageData->DeepCopy(connector2.GetVTKImage());
}

void iAFoamCharacterizationItemFilter::executeGaussian()
{
	typedef itk::DiscreteGaussianImageFilter <itk::Image<unsigned short, 3>, itk::Image<unsigned short, 3>> itkFilter;

	itkFilter::Pointer pFilter(itkFilter::New());

	iAConnector connector1;
	connector1.SetImage(m_pImageData);

	pFilter->SetInput(dynamic_cast<itk::Image<unsigned short, 3>*> (connector1.GetITKImage()));
	pFilter->SetVariance(m_dGaussianVariance);
	pFilter->Update();

	iAConnector connector2;
	connector2.SetImage(pFilter->GetOutput());

	m_pImageData->DeepCopy(connector2.GetVTKImage());
}

void iAFoamCharacterizationItemFilter::executeMedian()
{
	typedef itk::MedianImageFilter <itk::Image<unsigned short, 3>, itk::Image<unsigned short, 3>> itkFilter;

	itkFilter::Pointer pFilter(itkFilter::New());

	itkFilter::InputSizeType radius;
	radius.Fill(m_uiMedianRadius);
	pFilter->SetRadius(radius);

	iAConnector connector1;
	connector1.SetImage(m_pImageData);

	pFilter->SetInput(dynamic_cast<itk::Image<unsigned short, 3>*> (connector1.GetITKImage()));
	pFilter->Update();

	iAConnector connector2;
	connector2.SetImage(pFilter->GetOutput());

	m_pImageData->DeepCopy(connector2.GetVTKImage());
}

void iAFoamCharacterizationItemFilter::executeNonLocalMeans()
{
	typedef itk::PatchBasedDenoisingImageFilter <itk::Image<unsigned short, 3>, itk::Image<unsigned short, 3>> itkFilter;

	itkFilter::Pointer pFilter(itkFilter::New());

	iAConnector connector1;
	connector1.SetImage(m_pImageData);

	pFilter->SetInput(dynamic_cast<itk::Image<unsigned short, 3>*> (connector1.GetITKImage()));
	pFilter->SetNumberOfIterations(m_uiNonLocalMeansIteration);
	pFilter->SetPatchRadius(m_uiNonLocalMeansRadius);
	pFilter->Update();

	iAConnector connector2;
	connector2.SetImage(pFilter->GetOutput());

	m_pImageData->DeepCopy(connector2.GetVTKImage());
}

double iAFoamCharacterizationItemFilter::gaussianVariance() const
{
	return m_dGaussianVariance;
}

iAFoamCharacterizationItemFilter::EItemFilterType iAFoamCharacterizationItemFilter::itemFilterType() const
{
	return m_eItemFilterType;
}

QString iAFoamCharacterizationItemFilter::itemFilterTypeString() const
{
	switch (m_eItemFilterType)
	{
		case iftAnisotropic:
		return "A";
		break;

     	case iftGauss:
		return "G";
		break;

		case iftMedian:
		return "M";
		break;

		default:
		return "N";
		break;
	}
}

unsigned int iAFoamCharacterizationItemFilter::medianRadius() const
{
	return m_uiMedianRadius;
}

unsigned int iAFoamCharacterizationItemFilter::nonLocalMeansIteration() const
{
	return m_uiNonLocalMeansIteration;
}

unsigned int iAFoamCharacterizationItemFilter::nonLocalMeansRadius() const
{
	return m_uiNonLocalMeansRadius;
}

void iAFoamCharacterizationItemFilter::open(QFile* _pFileOpen)
{
	iAFoamCharacterizationItem::open(_pFileOpen);

	_pFileOpen->read((char*) &m_eItemFilterType, sizeof(m_eItemFilterType));
	_pFileOpen->read((char*) &m_dAnisotropicConductance, sizeof(m_dAnisotropicConductance));
	_pFileOpen->read((char*) &m_dAnisotropicTimeStep, sizeof(m_dAnisotropicTimeStep));
	_pFileOpen->read((char*) &m_uiAnisotropicIteration, sizeof(m_uiAnisotropicIteration));
	_pFileOpen->read((char*) &m_dGaussianVariance, sizeof(m_dGaussianVariance));
	_pFileOpen->read((char*)&m_uiMedianRadius, sizeof(m_uiMedianRadius));
	_pFileOpen->read((char*)&m_uiNonLocalMeansIteration, sizeof(m_uiNonLocalMeansIteration));
	_pFileOpen->read((char*)&m_uiNonLocalMeansRadius, sizeof(m_uiNonLocalMeansRadius));

	setItemText();
}

void iAFoamCharacterizationItemFilter::save(QFile* _pFileSave)
{
	iAFoamCharacterizationItem::save(_pFileSave);

	_pFileSave->write((char*) &m_eItemFilterType, sizeof(m_eItemFilterType));
	_pFileSave->write((char*) &m_dAnisotropicConductance, sizeof(m_dAnisotropicConductance));
	_pFileSave->write((char*) &m_dAnisotropicTimeStep, sizeof(m_dAnisotropicTimeStep));
	_pFileSave->write((char*) &m_uiAnisotropicIteration, sizeof(m_uiAnisotropicIteration));
	_pFileSave->write((char*) &m_dGaussianVariance, sizeof(m_dGaussianVariance));
	_pFileSave->write((char*) &m_uiMedianRadius, sizeof(m_uiMedianRadius));
	_pFileSave->write((char*)&m_uiNonLocalMeansIteration, sizeof(m_uiNonLocalMeansIteration));
	_pFileSave->write((char*)&m_uiNonLocalMeansRadius, sizeof(m_uiNonLocalMeansRadius));
}

void iAFoamCharacterizationItemFilter::setAnisotropicConductance(const double& _dAnisotropicConductance)
{
	m_dAnisotropicConductance = _dAnisotropicConductance;
}

void iAFoamCharacterizationItemFilter::setAnisotropicIteration(const unsigned int& _uiAnisotropicIteration)
{
	m_uiAnisotropicIteration = _uiAnisotropicIteration;
}

void iAFoamCharacterizationItemFilter::setAnisotropicTimeStep(const double& _dAnisotropicTimeStep)
{
	m_dAnisotropicTimeStep = _dAnisotropicTimeStep;
}

void iAFoamCharacterizationItemFilter::setGaussianVariance(const double& _dGaussianVariance)
{
	m_dGaussianVariance = _dGaussianVariance;
}

void iAFoamCharacterizationItemFilter::setItemFilterType(const iAFoamCharacterizationItemFilter::EItemFilterType& _eItemFilterType)
{
	m_eItemFilterType = _eItemFilterType;
}

void iAFoamCharacterizationItemFilter::setItemText()
{
	if (m_dExecute > 0.0)
	{
		setText(m_sName + QString(" [%1] (%2)").arg(itemFilterTypeString()).arg(executeTimeString()));
	}
	else
	{
		setText(m_sName + QString(" [%1]").arg(itemFilterTypeString()));
	}
}

void iAFoamCharacterizationItemFilter::setMedianRadius(const unsigned int& _uiMedianRadius)
{
	m_uiMedianRadius = _uiMedianRadius;
}

void iAFoamCharacterizationItemFilter::setNonLocalMeansIteration(const unsigned int& _uiNonLocalMeansIteration)
{
	m_uiNonLocalMeansIteration = _uiNonLocalMeansIteration;
}

void iAFoamCharacterizationItemFilter::setNonLocalMeansRadius(const unsigned int& _uiNonLocalMeansRadius)
{
	m_uiNonLocalMeansRadius = _uiNonLocalMeansRadius;
}
