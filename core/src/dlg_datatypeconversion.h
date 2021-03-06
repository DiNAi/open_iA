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

#include "charts/iAPlotData.h"
#include "ui_DataTypeConversion.h"

class iAConnector;

class QVTKWidget2;
class vtkImageData;
class vtkPlaneSource;

class QCheckBox;
class QComboBox;

#include <vtkSmartPointer.h>

#include <QDialog>
#include <QString>

class dlg_datatypeconversion : public QDialog, public Ui_DataTypeConversion
{
	Q_OBJECT

public:
	dlg_datatypeconversion ( QWidget *parent, vtkImageData* input, const char* filename, int intype, double* b, double* c, double* inPara );
	~dlg_datatypeconversion();

	void DataTypeConversion(QString const & m_filename, double* b);
	void DataTypeConversionROI(QString const & m_filename, double* b, double *roi);
	void histogramdrawing(iAPlotData::DataType* histbinlist, float min, float max, int m_bins, double discretization);

	void xyprojectslices();
	void xzprojectslices();
	QString coreconversionfunction(QString filename, QString & finalfilename, double* para, int indatatype, int outdatatype, double minrange, double maxrange, double minout, double maxout, int check);
	QString coreconversionfunctionforroi(QString filename, QString & finalfilename, double* para, int outdatatype, double minrange, double maxrange, double minout, double maxout, int check, double* roi);
	void updatevalues(double* inPara);
	void updateroi( );

	double getRangeLower();
	double getRangeUpper();
	double getOutputMin();
	double getOutputMax();
	double getXOrigin();
	double getXSize();
	double getYOrigin();
	double getYSize();
	double getZOrigin();
	double getZSize();

	QString getDataType();
	int getConvertROI();
	private slots:
		void update(QString a);

private:
	QString text11;

	iAConnector* m_roiconvertimage;

	double * m_bptr;
	int m_bins;
	vtkImageData* imageData;
	int m_intype;
	double m_sliceskip, m_insizex,	m_insizey, m_insizez, m_inspacex, m_inspacey, m_inspacez;
	QString m_filename;
	iAPlotData::DataType * m_histbinlist;
	float m_min, m_max, m_dis;
	vtkImageData* m_testxyimage, * m_testxzimage, * m_testyzimage, * m_roiimage;
	QVTKWidget2* vtkWidgetXY, *vtkWidgetXZ, *vtkWidgetYZ;

	iAConnector* xyconvertimage, * xzconvertimage, * yzconvertimage;

	int m_xstart, m_xend, m_ystart, m_yend, m_zstart, m_zend;
	QLineEdit* leRangeLower, *leRangeUpper, *leOutputMin,*leOutputMax, *leXOrigin, *leXSize, *leYOrigin, *leYSize, *leZOrigin, *leZSize;
	QComboBox* cbDataType;
	QCheckBox* chConvertROI
		// , *chUseMaxDatatypeRange
	;
	double m_roi[6];
	double m_spacing[3];

	vtkSmartPointer<vtkPlaneSource> xyroiSource, xzroiSource;
};
