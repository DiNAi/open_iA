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
 
#include "pch.h"
#include "iACharacteristics.h"

#include "iAConsole.h"
#include "iAImageTypes.h"
#include "iASingleResult.h"
#include "iAAttributes.h"

#include <itkImageFileWriter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkScalarConnectedComponentImageFilter.h>

#include <QString>

/**
* funcion vecMed
*/
double CharacteristicsCalculator::vecMed(std::vector<double> vec) {
	if (vec.empty()) return 0;
	else {
		std::sort(vec.begin(), vec.end());
		if (vec.size() % 2 == 0)
			return (vec[vec.size() / 2 - 1] + vec[vec.size() / 2]) / 2;
		else
			return vec[vec.size() / 2];
	}
}

/**
* funcion vecMAD
*/
double CharacteristicsCalculator::vecMAD(std::vector<double> vec) {
	if (vec.empty()) return 0;
	else {
		double median = vecMed(vec);
		std::vector<double> MAD;
		std::sort(vec.begin(), vec.end());

		for (std::vector<double>::iterator it = vec.begin(); it != vec.end(); ++it)
			MAD.push_back(fabs(*it - median));

		if (MAD.size() % 2 == 0)
			return (MAD[MAD.size() / 2 - 1] + MAD[MAD.size() / 2]) / 2;
		else
			return MAD[MAD.size() / 2];
	}
}

void sublistFiles1(QDir directory, QString indent, std::vector<std::string> &allMHAs)
{
	indent += "\t";
	QDir dir(directory);
	QFileInfoList list = dir.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot);
	std::string partMHA = std::string(".mha");
	foreach(QFileInfo finfo, list) {
		std::string fullFilename = std::string(finfo.fileName().toStdString());
		if (fullFilename.find(partMHA) != std::string::npos) {
			allMHAs.push_back(finfo.canonicalFilePath().toStdString());
		}
	}
}

iACharacteristics::iACharacteristics():
m_objectCount(0),
m_duration(0.0),
m_median_accuracy(0.0), m_mean_accuracy(0.0), m_dice_overall(0.0), m_dice_0(0.0), m_dice_1(0.0), m_dice_2(0.0), m_dice_3(0.0), m_dice_4(0.0),
m_mean_uncertainty(0.0), m_median_uncertainty(0.0), m_uncertainty_0(0.0), m_uncertainty_1(0.0), m_uncertainty_2(0.0), m_uncertainty_3(0.0), m_uncertainty_4(0.0),
m_median_confusion_uncertainy_false(0.0), m_median_confusion_uncertainy_true(0.0), m_median_confusion_bad_to_good(0.0)
{
	std::fill(m_measures, m_measures+MeasureCount, 0.0);
}

int iACharacteristics::ObjectCount() const
{
	return m_objectCount;
}

void iACharacteristics::SetObjectCount(int objCount)
{
	m_objectCount = objCount;
}

double iACharacteristics::Duration() const
{
	return m_duration;
}

void iACharacteristics::SetDuration(double value)
{
	m_duration = value;
}

void iACharacteristics::SetMeasures(double value[MeasureCount])
{
	for (int i=0; i<MeasureCount; ++i)
	{
		m_measures[i] = value[i];
	}
}

double iACharacteristics::Measure(int idx) const
{
	return m_measures[idx];
}

iACharacteristics iACharacteristics::Create(QString const & descriptor)
{
	iACharacteristics result;
	QStringList tokens = descriptor.split(" ");
	if (tokens.size() == 0)
	{
		DEBUG_LOG(QString("Characteristics: Invalid line '%1'\n").arg(descriptor));
		return result;
	}
	bool ok;
	result.m_objectCount = tokens[0].toInt(&ok);
	if (!ok)
	{
		DEBUG_LOG(QString("Characteristics: Invalid objectCount in line '%1'\n").arg(descriptor));
	}
	ok = true;
	result.m_duration = (tokens.size() > 1) ? tokens[1].toDouble(&ok) : 0;
	result.m_median_accuracy = (tokens.size() > 2) ? tokens[2].toDouble(&ok) : 0;
	result.m_mean_accuracy = (tokens.size() > 3) ? tokens[3].toDouble(&ok) : 0;
	result.m_dice_overall = (tokens.size() > 4) ? tokens[4].toDouble(&ok) : 0;
	result.m_dice_0 = (tokens.size() > 5) ? tokens[5].toDouble(&ok) : 0;
	result.m_dice_1 = (tokens.size() > 6) ? tokens[6].toDouble(&ok) : 0;
	result.m_dice_2 = (tokens.size() > 7) ? tokens[7].toDouble(&ok) : 0;
	result.m_dice_3 = (tokens.size() > 8) ? tokens[8].toDouble(&ok) : 0;
	result.m_dice_4 = (tokens.size() > 9) ? tokens[9].toDouble(&ok) : 0;
	result.m_mean_uncertainty = (tokens.size() > 10) ? tokens[10].toDouble(&ok) : 0;
	result.m_median_uncertainty = (tokens.size() > 11) ? tokens[11].toDouble(&ok) : 0;
	result.m_uncertainty_0 = (tokens.size() > 12) ? tokens[12].toDouble(&ok) : 0;
	result.m_uncertainty_1 = (tokens.size() > 13) ? tokens[13].toDouble(&ok) : 0;
	result.m_uncertainty_2 = (tokens.size() > 14) ? tokens[14].toDouble(&ok) : 0;
	result.m_uncertainty_3 = (tokens.size() > 15) ? tokens[15].toDouble(&ok) : 0;
	result.m_uncertainty_4 = (tokens.size() > 16) ? tokens[16].toDouble(&ok) : 0;
	result.m_median_confusion_uncertainy_false = (tokens.size() > 17) ? tokens[17].toDouble(&ok) : 0;
	result.m_median_confusion_uncertainy_true = (tokens.size() > 18) ? tokens[18].toDouble(&ok) : 0;
	result.m_median_confusion_bad_to_good = (tokens.size() > 19) ? tokens[19].toDouble(&ok) : 0;

	if (!ok)
	{
		DEBUG_LOG(QString("Characteristics:Invalid duration in line '%1'\n").arg(descriptor));
	}
	if (tokens.size() <= 1)
	{
		DEBUG_LOG(QString("Characteristics: Missing duration in line '%1'\n").arg(descriptor));
		return result;
	}

	return result;
}

QString iACharacteristics::GetDescriptor() const
{
	return QString::number(m_objectCount) + " " + QString::number(m_duration) + " " + QString::number(m_median_accuracy) + " " + QString::number(m_mean_accuracy) + " " + QString::number(m_dice_overall) +
		" " + QString::number(m_dice_0) + " " + QString::number(m_dice_1) + " " + QString::number(m_dice_2) + " " + QString::number(m_dice_3) + " " + QString::number(m_dice_4) +
		" " + QString::number(m_mean_uncertainty) + " " + QString::number(m_median_uncertainty) + " " + QString::number(m_uncertainty_0) + " " + QString::number(m_uncertainty_1) + " " + QString::number(m_uncertainty_2) +
		" " + QString::number(m_uncertainty_3) + " " + QString::number(m_uncertainty_4) + " " + QString::number(m_median_confusion_uncertainy_false) + " " + QString::number(m_median_confusion_uncertainy_true) + " " + QString::number(m_median_confusion_bad_to_good);
}

CharacteristicsCalculator::CharacteristicsCalculator(QSharedPointer<iASingleResult> result, QSharedPointer<iAAttributes> range, int objCountIdx):
m_result(result), m_range(range),
m_objCountIdx(objCountIdx)
{}

void CharacteristicsCalculator::GetUncertaintyValues(QSharedPointer<iASingleResult> result, QString groundTruthPath, std::vector< double > &uncert)
{
	double limit = -std::log(1.0 / 5);
	double normalizeFactor = 1 / limit;

	typedef itk::Image< double, 3 >	DoubleImageType;
	typedef itk::Image< int, 3 >	IntImageType;

	typedef itk::ImageFileReader<DoubleImageType> DoubleReaderType;
	typedef itk::ImageFileReader<IntImageType	> IntReaderType;

	DoubleReaderType::Pointer readerGT = DoubleReaderType::New();
	QString gtPath = QString::fromStdString(groundTruthPath.toStdString());
	readerGT->SetFileName(gtPath.toStdString());
	readerGT->Update();

	std::vector< std::vector<int> > seeds;
	//createSeedVector(seeds);

	try
	{
		IntImageType::Pointer m_InputImage = dynamic_cast<IntImageType *>(result->GetLabelledImage().GetPointer());
		DoubleImageType::Pointer m_InputImage0 = dynamic_cast<DoubleImageType *>(result->GetProbabilityImg(0).GetPointer());
		DoubleImageType::Pointer m_InputImage1 = dynamic_cast<DoubleImageType *>(result->GetProbabilityImg(1).GetPointer());
		DoubleImageType::Pointer m_InputImage2 = dynamic_cast<DoubleImageType *>(result->GetProbabilityImg(2).GetPointer());
		DoubleImageType::Pointer m_InputImage3 = dynamic_cast<DoubleImageType *>(result->GetProbabilityImg(3).GetPointer());
		DoubleImageType::Pointer m_InputImage4 = dynamic_cast<DoubleImageType *>(result->GetProbabilityImg(4).GetPointer());
		DoubleImageType::Pointer m_InputGT = dynamic_cast<DoubleImageType *>(readerGT->GetOutput());

		DoubleImageType::Pointer m_OutputImage;

		m_OutputImage = DoubleImageType::New();
		DoubleImageType::SpacingType m_ImgSpace = m_InputImage->GetSpacing();

		// Initialize output image
		m_OutputImage->SetSpacing(m_InputImage->GetSpacing());
		m_OutputImage->SetOrigin(m_InputImage->GetOrigin());
		m_OutputImage->SetRegions(m_InputImage->GetRequestedRegion());
		m_OutputImage->Allocate();
		m_OutputImage->FillBuffer(0.0);

		//initiate image iterators
		itk::ImageRegionIterator<IntImageType> inIter(m_InputImage, m_InputImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter0(m_InputImage0, m_InputImage0->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter1(m_InputImage1, m_InputImage1->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter2(m_InputImage2, m_InputImage2->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter3(m_InputImage3, m_InputImage3->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter4(m_InputImage4, m_InputImage4->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIterGT(m_InputGT, m_InputGT->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> outIter(m_OutputImage, m_OutputImage->GetLargestPossibleRegion());

		inIter.GoToBegin();
		inIter0.GoToBegin();
		inIter1.GoToBegin();
		inIter2.GoToBegin();
		inIter3.GoToBegin();
		inIter4.GoToBegin();
		inIterGT.GoToBegin();
		outIter.GoToBegin();
		//outIterAll.GoToBegin();

		int first100 = 0;

		double min = 1.0;
		double max = 0.0;
		//std::cout << "Sample nr: " << i << endl;

		std::vector<double> entropies;
		std::vector<double> First_Minus_Second;

		int nrOferror = 0; int nrOferror_l[5] = { 0 }; int nrOferror_c[5] = { 0 };
		int nrOfSeedErrors = 0;
		int nrOfVoxels = 0;
		int nrOfSeeds = 0;

		while (!inIter0.IsAtEnd())
		{
			double val0 = inIter0.Get(); double val1 = inIter1.Get(); double val2 = inIter2.Get(); double val3 = inIter3.Get(); double val4 = inIter4.Get();
			double entropy = val0*std::log(val0) + val1*std::log(val1) + val2*std::log(val2) + val3*std::log(val3) + val4*std::log(val4);
			entropy = (-entropy) * normalizeFactor;
			outIter.Set(entropy);
			//outIterAll.Set(entropy + outIterAll.Get());

			std::vector<double> probs; probs.push_back(val0); probs.push_back(val1); probs.push_back(val2); probs.push_back(val3); probs.push_back(val4);
			std::sort(probs.begin(), probs.end());
			double FirstMinusSecond = probs[4] - probs[3];

			if (entropy < min) {
				min = entropy;
			}

			if (entropy > max) {
				max = entropy;
			}

			entropies.push_back(entropy);
			First_Minus_Second.push_back(FirstMinusSecond);

			if (inIter.Get() != inIterGT.Get()) {
				nrOferror++;
				nrOferror_l[int(inIterGT.Get())]++; nrOferror_c[int(inIterGT.Get())]++;
			}
			else {
				nrOferror_c[int(inIterGT.Get())]++;
			}

			++inIter;
			++inIter0;
			++inIter1;
			++inIter2;
			++inIter3;
			++inIter4;
			++inIterGT;
			++outIter;
			//++outIterAll;
			nrOfVoxels++;
		}

		for (std::vector< std::vector<int> >::iterator it = seeds.begin(); it != seeds.end(); ++it) {
			DoubleImageType::IndexType index;
			std::vector<int> pts = *it;
			index[0] = pts[0]; index[1] = pts[1]; index[2] = pts[2];

			inIter.SetIndex(index);
			inIterGT.SetIndex(index);

			if (inIter.Get() != pts[3]) {
				nrOfSeedErrors++;
			}

			nrOfSeeds++;
		}

		double mean = -1;
		double sum = std::accumulate(entropies.begin(), entropies.end(), 0.0);
		mean = sum / entropies.size();

		std::vector<double> diff(entropies.size());
		std::transform(entropies.begin(), entropies.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / entropies.size());
		double median = vecMed(entropies);
		double MAD = vecMAD(entropies);

		double successRate_l[5] = { 0.0 };
		double accuracyGT = 100.0 - (nrOferror / (nrOfVoxels / 100.0));
		double accuracySeeds = 100.0 - (nrOfSeedErrors / (nrOfSeeds / 100.0));

		for (int j = 0; j < 5; j++) {
			successRate_l[j] = 100.0 - (nrOferror_l[j] / (nrOferror_c[j] / 100.0));
		}

		std::cout << "Mean: " << mean << "Accuracy: " << accuracyGT / 100 << "Accuracy(Seeds): " << accuracySeeds / 100 << "Median: " << median << "Sigma: " << stdev << "MAD: " << MAD << "S0: " << successRate_l[0] << "S1: " << successRate_l[1] << "S2: " << successRate_l[2] << "S3: " << successRate_l[3] << "S4: " << successRate_l[4] << std::endl;

		//x.push_back(mean);
		//y.push_back(accuracyGT / 100);

		cout << mean << ";" << accuracyGT / 100 << ";" << accuracySeeds / 100 << ";" << median << ";" << stdev << ";" << MAD << ";" << successRate_l[0] << ";" << successRate_l[1] << ";" << successRate_l[2] << ";" << successRate_l[3] << ";" << successRate_l[4] << std::endl;
	}
	catch (itk::ExceptionObject &err)
	{
		throw err;
	}
}


void CharacteristicsCalculator::run()
{
	typedef itk::Image< unsigned int, 3 > OutputImageType;
	typedef itk::ScalarConnectedComponentImageFilter <LabelImageType, OutputImageType > ConnectedComponentImageFilterType;
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
	connected->SetDistanceThreshold(0.5);
	if (m_result->GetLabelledImage().IsNull())
	{
		DEBUG_LOG("Labelled Image is null\n");
		return;
	}
	LabelImageType* lblImg = dynamic_cast<LabelImageType*>(m_result->GetLabelledImage().GetPointer());
	connected->SetInput(lblImg);
	connected->Update();
	typedef itk::RelabelComponentImageFilter <OutputImageType, OutputImageType >
	RelabelFilterType;
	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
		relabel->SetInput(connected->GetOutput());
	//relabel->SetSortByObjectSize(false);
	relabel->Update();
	int objCount = relabel->GetNumberOfObjects();
	m_result->SetAttribute(m_objCountIdx, objCount);

	srand(static_cast <unsigned> (time(0)));

	std::vector<std::string> patientMHAs;
	QDir dir(m_result->GetFolder());
	dir.cdUp();
	sublistFiles1(dir, "", patientMHAs);

	int nrModals = 4;
	// check accuracy
	{
		typedef itk::Image<double, 3> DoubleImageType;

		typedef itk::Image< int, 3 >	IntImageType;
		typedef itk::Image< double, 3 >	DoubleImageType;

		typedef itk::ImageFileReader<DoubleImageType> DoubleReaderType;
		DoubleReaderType::Pointer DoubleReader = DoubleReaderType::New();
		DoubleReader->SetFileName(patientMHAs[nrModals].c_str());	// ground_truth
		DoubleReader->Update();
		DoubleImageType::Pointer m_gtImage = dynamic_cast<DoubleImageType *>(DoubleReader->GetOutput());

		std::string output;
		int nrOferrorOverall = 0, nrOfVoxelsOverall = 0;
		LabelImageType* m_LabelImage = dynamic_cast<LabelImageType*>(m_result->GetLabelledImage().GetPointer());

		itk::ImageRegionIterator<LabelImageType> sampleIter(m_LabelImage, m_LabelImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> gtIter(m_gtImage, m_gtImage->GetLargestPossibleRegion());
		sampleIter.GoToBegin();
		gtIter.GoToBegin();

		while (!gtIter.IsAtEnd())
		{
			if (gtIter.Get() != -1) {
				if (gtIter.Get() != sampleIter.Get()) {
					nrOferrorOverall++;
				}
				nrOfVoxelsOverall++;
			}
			++sampleIter;
			++gtIter;
		}

		double successRate = 100.0 - (nrOferrorOverall / (nrOfVoxelsOverall / 100.0));

		typedef itk::CastImageFilter< DoubleImageType, IntImageType > CastFilterType;
		CastFilterType::Pointer castFilterGroundTruth = CastFilterType::New();
		castFilterGroundTruth->SetInput(m_gtImage);

		typedef itk::LabelOverlapMeasuresImageFilter<LabelImageType> DiceType;
		DiceType::Pointer DiceCompare = DiceType::New();
		DiceCompare->SetSourceImage(m_LabelImage);
		DiceCompare->SetTargetImage(castFilterGroundTruth->GetOutput());
		DiceCompare->Update();

		double diceValue = DiceCompare->GetDiceCoefficient();

		std::vector< double > uncert(5);
		GetUncertaintyValues(m_result, QString::fromStdString(patientMHAs[nrModals].c_str()), uncert);

		m_result->SetAttribute(m_objCountIdx + 2, successRate);
		m_result->SetAttribute(m_objCountIdx + 3, successRate);
		m_result->SetAttribute(m_objCountIdx + 4, diceValue);
		m_result->SetAttribute(m_objCountIdx + 5, DiceCompare->GetDiceCoefficient(0));
		m_result->SetAttribute(m_objCountIdx + 6, DiceCompare->GetDiceCoefficient(1));
		m_result->SetAttribute(m_objCountIdx + 7, DiceCompare->GetDiceCoefficient(2));
		m_result->SetAttribute(m_objCountIdx + 8, DiceCompare->GetDiceCoefficient(3));
		m_result->SetAttribute(m_objCountIdx + 9, DiceCompare->GetDiceCoefficient(4));
		m_result->SetAttribute(m_objCountIdx + 10, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 11, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 12, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 13, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 14, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 15, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 16, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 17, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 18, ((double)rand() / (double)RAND_MAX));
		m_result->SetAttribute(m_objCountIdx + 19, ((double)rand() / (double)RAND_MAX));
	}

	/*
	itk::ImageFileWriter<OutputImageType>::Pointer writer = itk::ImageFileWriter<OutputImageType>::New();
	writer->SetFileName(debugCount.toStdString() );
	writer->SetUseCompression(true);
	writer->SetInput(relabel->GetOutput() );
	writer->Update();
	*/
}