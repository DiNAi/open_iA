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
*          Stelzhamerstra�e 23, 4600 Wels / Austria, Email:                           *
* ************************************************************************************/
#pragma once

#include "open_iA_Core_export.h"

#include <vtkSmartPointer.h>

#include <QSharedPointer>
#include <QVector>

#include "ui_modalities.h"
#include <iAQTtoUIConnector.h>
typedef iAQTtoUIConnector<QDockWidget, Ui_modalities> dlg_modalitiesUI;

class dlg_planeSlicer;
class iAFast3DMagicLensWidget;
class iAHistogramWidget;
class iAModality;
class iAModalityList;
class iAVolumeRenderer;
class iAVolumeSettings;
class MdiChild;
class iAModalityTransfer;

class vtkActor;
class vtkColorTransferFunction;
class vtkImageData;
class vtkPiecewiseFunction;
class vtkPlane;


// TODO: VOLUME: split off volume manager for the management of the actual volume rendering stuff
class open_iA_Core_API dlg_modalities : public dlg_modalitiesUI
{
	Q_OBJECT
public:
	dlg_modalities(iAFast3DMagicLensWidget* renderer, int numBin, QDockWidget* histogramContainer);
	void SetModalities(QSharedPointer<iAModalityList> modalities);
	QSharedPointer<iAModalityList const>  GetModalities() const;
	QSharedPointer<iAModalityList>  GetModalities();
	int						GetSelected() const;
	vtkSmartPointer<vtkColorTransferFunction> GetCTF(int modality);
	vtkSmartPointer<vtkPiecewiseFunction> GetOTF(int modality);
	void ChangeRenderSettings(iAVolumeSettings const & rs);
	void Store(QString const & filename);
	bool Load(QString const & filename);
	iAHistogramWidget* GetHistogram();
	void ShowVolumes(bool show);
	void ShowSlicePlanes(bool enabled);
	void SetSlicePlanes(vtkPlane* plane1, vtkPlane* plane2, vtkPlane* plane3);

	void AddModality(vtkSmartPointer<vtkImageData>, QString const & name);

public slots:
	void Load();
	void Store();
signals:
	void ModalityAvailable();
	void ShowImage(vtkSmartPointer<vtkImageData> img);

	//! @{ for histogram:
	void PointSelected();
	void NoPointSelected();
	void EndPointSelected();
	void Active();
	void AutoUpdateChanged(bool toogled);
	void UpdateViews();
	//! @}

private slots:
	void AddClicked();
	void RemoveClicked();
	void EditClicked();
	void ManualRegistration();
	void MagicLens();
	//void CuttingPlane();
	void RendererMouseMoved();

	void EnableButtons();
	void ListClicked(QListWidgetItem* item);
	void ModalityAdded(QSharedPointer<iAModality> mod);

private:
	MdiChild* m_mdiChild;
	QSharedPointer<iAModalityList> modalities;
	QString m_FileName;
	int m_selectedRow;
	iAFast3DMagicLensWidget* m_renderer;
	int m_numBin;
	QDockWidget* m_histogramContainer;
	iAHistogramWidget* m_currentHistogram;
	void SwitchHistogram(QSharedPointer<iAModalityTransfer> modTrans);
	void ShowVolume(QSharedPointer<iAVolumeRenderer> renderer, bool enabled);
	bool m_showVolumes;
	bool m_showSlicePlanes;
	vtkPlane *m_plane1, *m_plane2, *m_plane3;
};