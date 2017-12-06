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
#include "iABatchFilter.h"

#include "iAAttributeDescriptor.h"
#include "iAConnector.h"
#include "iAConsole.h"
#include "iAFilterRegistry.h"
#include "io/iAITKIO.h"
#include "iAStringHelper.h"

#include <QDirIterator>
#include <QFile>
#include <QTextStream>

iABatchFilter::iABatchFilter():
	iAFilter("Batch...", "Image Ensembles",
		"Runs a filter on a selected set of images.<br/>"
		"Specify an <em>Image folder</em> which contains the images to be processed. "
		"<em>Recursive</em> toggles whether or not to also consider subdirectories. "
		"The <em>File mask</em> is applied to match which files in the given folder are processed "
		"(separate multiple masks via ';', e.g. '*.mhd;*.tif'. "
		"The specified <em>Filter</em> is applied to all files specified with above settings, "
		"every time executed with the same set of <em>Parameters</em>. "
		"When <em>Output file</em> is not empty, all output values produced by the filter "
		"will be written to the file name given here, one row per image and filter. "
		"If the output csv file exists, and <em>Append to output</em> is enabled, "
		"the output values are appended at the end of each line. "
		"If <em>Add filename</em> is enabled, then the name of the file processed for that "
		"line will be appended before the first output value from that file.", 0, 0)
{
	AddParameter("Image folder", String, "");
	AddParameter("Recursive", Boolean, false);
	AddParameter("File mask", String, "*.mhd");
	AddParameter("Filter", String, "Image Quality");
	AddParameter("Parameters", String, "");
	AddParameter("Output csv file", String, "");
	AddParameter("Append to output", Boolean, true);
	AddParameter("Add filename", Boolean, true);
}

void iABatchFilter::PerformWork(QMap<QString, QVariant> const & parameters)
{
	auto filter = iAFilterRegistry::Filter(parameters["Filter"].toString());
	if (!filter)
	{
		DEBUG_LOG(QString("Batch: Cannot run filter '%1', it does not exist!").arg(parameters["Filter"].toString()));
		return;
	}
	QMap<QString, QVariant> filterParams;
	QStringList filterParamStrs = SplitPossiblyQuotedString(parameters["Parameters"].toString());
	if (filter->Parameters().size() != filterParamStrs.size())
	{
		DEBUG_LOG(QString("Batch: Invalid number of parameters: %1 expected, %2 given!")
			.arg(filter->Parameters().size())
			.arg(filterParamStrs.size()));
		return;
	}

	iAConnector* con = new iAConnector();
	QVector<iAConnector*> inputImages;
	inputImages.push_back(con);

	for (int i=0; i<filterParamStrs.size(); ++i)
		filterParams.insert(filter->Parameters()[i]->Name(), filterParamStrs[i]);

	QString outputFile = parameters["Output csv file"].toString();
	QStringList outputBuffer;
	if (parameters["Append to output"].toBool() && QFile(outputFile).exists())
	{
		QFile file(outputFile);
		if (file.open(QIODevice::ReadOnly | QIODevice::Text))
		{
			QTextStream textStream(&file);
			while (!textStream.atEnd())
				outputBuffer << textStream.readLine();
			file.close();
		}
	}
	filter->SetUp(inputImages, m_log, m_progress);

	QStringList filters = parameters["File mask"].toString().split(";");
	QDirIterator it(parameters["Image folder"].toString(), filters, QDir::Files,
			parameters["Recursive"].toBool() ?
					QDirIterator::Subdirectories :
					QDirIterator::NoIteratorFlags);
	size_t curLine = 0;
	while (it.hasNext())
	{
		QString fileName = it.next();
		iAITKIO::ScalarPixelType pixelType;
		iAITKIO::ImagePointer img = iAITKIO::readFile(fileName, pixelType, false);
		inputImages[0]->SetImage(img);
		filter->Run(filterParams);
		if (curLine == 0)
		{
			QStringList captions;
			if (parameters["Add filename"].toBool())
				captions << "filename";
			for (auto outValue : filter->OutputValues())
				captions << outValue.first;
			if (outputBuffer.empty())
				outputBuffer.append("");
			outputBuffer[0] += (outputBuffer[0].isEmpty() ? "" : ",") + captions.join(",");
			++curLine;
		}
		if (curLine >= outputBuffer.size())
			outputBuffer.append("");
		QStringList values;
		if (parameters["Add filename"].toBool())
			values << fileName;
		for (auto outValue : filter->OutputValues())
			values.append(outValue.second.toString());
		QString textToAdd = (outputBuffer[curLine].isEmpty() ? "" : ",") + values.join(",");
		outputBuffer[curLine] += textToAdd;
		++curLine;
	}
	if (!outputFile.isEmpty())
	{
		QFile file(outputFile);
		if (file.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			QTextStream textStream(&file);
			for (QString line : outputBuffer)
			{
				textStream << line << endl;
			}
			file.close();
		}
	}
}

IAFILTER_CREATE(iABatchFilter);
