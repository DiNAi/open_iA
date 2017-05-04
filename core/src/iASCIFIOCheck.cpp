/*************************************  open_iA  ************************************ *
* **********  A tool for scientific visualisation and 3D image processing  ********** *
* *********************************************************************************** *
* Copyright (C) 2016-2017  C. Heinzl, M. Reiter, A. Reh, W. Li, M. Arikan,            *
*                          J. Weissenb�ck, Artem & Alexander Amirkhanov, B. Fr�hler   *
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
#include "iASCIFIOCheck.h"

#include "iAConsole.h"

#include <QCoreApplication>
#include <QFileInfo>

void CheckSCIFIO()
{
#ifdef USE_SCIFIO
	// Workaround for ITK requiring SCIFIO_PATH to be set when compiled with SCIFIO
	const char* SCIFIO_PATH = "SCIFIO_PATH";
	QString envScifioPath(getenv(SCIFIO_PATH));

	if (envScifioPath.length() > 0)
	{
		return;
	}

	QFileInfo fi(QCoreApplication::applicationDirPath());
	QString scifioPath(fi.absoluteFilePath() + "/scifio_jars");
	if (!QFile::exists(scifioPath))
	{
		DEBUG_LOG(QString("ITK was built with SCIFIO, SCIFIO_PATH environment variable is not set, and scifio_jars directory (%1) was not found."
			"You might not be able to load files!").arg(scifioPath));
		return;
	}
	scifioPath.replace("/", "\\");
	QString scifioPathAssign(QString(SCIFIO_PATH) + "=" + scifioPath);
#ifdef _MSC_VER
	_putenv(scifioPathAssign.toStdString().c_str());
#else
	static char * scifioPathBuffer = new char[scifioPathAssign.length() + 1];
	strcpy(scifioPathBuffer, scifioPathAssign.toStdString().c_str());
	putenv(scifioPathBuffer);
#endif

#endif
}
