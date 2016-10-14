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
#include "iAMMSegSampler.h"

#include "iAConsole.h"
#include "iACharacteristics.h"
#include "iAAttributes.h"
#include "iAAttributeDescriptor.h"
#include "iAImageCoordinate.h"
#include "iAMMSegParameter.h"
#include "iAMMSegParameterGenerator.h"
#include "iAMMSegParameterRange.h"
#include "iASingleResult.h"
#include "iASamplingResults.h"
#include "iAModality.h"
#include "iANormalizerImpl.h"
#include "iAPCA.h"
#include "iASpectraDistance.h"
#include "SVMImageFilter.h"

#include <QMap>
#include <QTextStream>
#include <QDir>

/**
* funcion sublistFiles2
*/
void sublistFiles2(QDir directory, QString indent, std::vector<std::string> &allMHAs)
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

/**
* funcion vecMed
*/
double iAMMSegSampler::vecMed(std::vector<double> vec) {
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
double iAMMSegSampler::vecMAD(std::vector<double> vec) {
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

void listFiles(QDir directory, QString indent, std::vector<std::string> &allMHAs)
{
	indent += "\t";
	QDir dir(directory);
	QFileInfoList list = dir.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot);
	std::string partFilename = std::string(".OT.");
	std::string partMHA = std::string(".mha");
	foreach(QFileInfo finfo, list) {
		std::string fullFilename = std::string(finfo.fileName().toStdString());
		if (fullFilename.find(partFilename) != std::string::npos && fullFilename.find(partMHA) != std::string::npos) {
			//cout << finfo.absoluteFilePath().toStdString() << endl;
			allMHAs.push_back(finfo.canonicalFilePath().toStdString());
		}
		if (finfo.isDir()) {
			listFiles(QDir(finfo.absoluteFilePath()), indent, allMHAs);
		}
	}
}

const int CONCURRENT_ERW_RUNS = 1;

iAMMSegSampler::iAMMSegSampler(
		QSharedPointer<iAModalityList const> modalities,
		QSharedPointer<iAMMSegParameterRange> range,
		QSharedPointer<iAMMSegParameterGenerator> sampleGenerator,
		SVMImageFilter::SeedsPointer seeds,
		QString const & outputBaseDir,
		QString const & parameterRangeFile,
		QString const & parameterSetFile,
		QString const & characteristicsFile,
		bool storeProbabilities) :
	m_modalities(modalities),
	m_range(range),
	m_sampleGenerator(sampleGenerator),
	m_seeds(seeds),
	m_curLoop(0),
	m_parameters(0),
	m_outputBaseDir(outputBaseDir),
	m_aborted(false),
	m_storeProbabilities(storeProbabilities),
	m_parameterRangeFile (parameterRangeFile),
	m_parameterSetFile   (parameterSetFile),
	m_characteristicsFile(characteristicsFile),
	m_runningOperations(0),
	m_calcERWDuration(0),
	m_calcChaDuration(0)
{
}

void iAMMSegSampler::StatusMsg(QString const & msg)
{
	emit Status(msg);
	DEBUG_LOG(msg+"\n");
}

void iAMMSegSampler::GetUncertaintyValues(QString resultDir, QString groundTruthPath, std::vector< double > &uncert)
{
	double limit = -std::log(1.0 / 4);
	double normalizeFactor = 1 / limit;

	typedef itk::Image< double, 3 >	DoubleImageType;
	typedef itk::ImageFileReader<DoubleImageType> DoubleReaderType;

	DoubleReaderType::Pointer reader = DoubleReaderType::New();
	DoubleReaderType::Pointer reader0 = DoubleReaderType::New();
	DoubleReaderType::Pointer reader1 = DoubleReaderType::New();
	DoubleReaderType::Pointer reader2 = DoubleReaderType::New();
	DoubleReaderType::Pointer reader3 = DoubleReaderType::New();
	DoubleReaderType::Pointer readerGT = DoubleReaderType::New();

	QString temp = resultDir + "\\label.mhd";
	QString temp0 = resultDir + "\\prob0.mhd";
	QString temp1 = resultDir + "\\prob1.mhd";
	QString temp2 = resultDir + "\\prob2.mhd";
	QString temp3 = resultDir + "\\prob3.mhd";
	QString entPath = resultDir + "\\firstsecond.mhd";
	QString gtPath = QString::fromStdString(groundTruthPath.toStdString());

	reader->SetFileName(temp.toStdString());
	reader0->SetFileName(temp0.toStdString());
	reader1->SetFileName(temp1.toStdString());
	reader2->SetFileName(temp2.toStdString());
	reader3->SetFileName(temp3.toStdString());
	readerGT->SetFileName(gtPath.toStdString());

	reader->Update();
	reader0->Update();
	reader1->Update();
	reader2->Update();
	reader3->Update();
	readerGT->Update();

	std::vector< std::vector<int> > seeds;
	//createSeedVector(seeds);

	try
	{
		DoubleImageType::Pointer m_InputImage = dynamic_cast<DoubleImageType *>(reader->GetOutput());
		DoubleImageType::Pointer m_InputImage0 = dynamic_cast<DoubleImageType *>(reader0->GetOutput());
		DoubleImageType::Pointer m_InputImage1 = dynamic_cast<DoubleImageType *>(reader1->GetOutput());
		DoubleImageType::Pointer m_InputImage2 = dynamic_cast<DoubleImageType *>(reader2->GetOutput());
		DoubleImageType::Pointer m_InputImage3 = dynamic_cast<DoubleImageType *>(reader3->GetOutput());
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
		itk::ImageRegionIterator<DoubleImageType> inIter(m_InputImage, m_InputImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter0(m_InputImage0, m_InputImage0->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter1(m_InputImage1, m_InputImage1->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter2(m_InputImage2, m_InputImage2->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIter3(m_InputImage3, m_InputImage3->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> inIterGT(m_InputGT, m_InputGT->GetLargestPossibleRegion());
		itk::ImageRegionIterator<DoubleImageType> outIter(m_OutputImage, m_OutputImage->GetLargestPossibleRegion());

		inIter.GoToBegin();
		inIter0.GoToBegin();
		inIter1.GoToBegin();
		inIter2.GoToBegin();
		inIter3.GoToBegin();
		inIterGT.GoToBegin();
		outIter.GoToBegin();
		//outIterAll.GoToBegin();

		int first100 = 0;

		double min = 1.0;
		double max = 0.0;
		//std::cout << "Sample nr: " << i << endl;

		std::vector<double> entropies;
		std::vector<double> First_Minus_Second;

		int nrOferror = 0; int nrOferror_l[4] = { 0 }; int nrOferror_c[4] = { 0 };
		int nrOfSeedErrors = 0;
		int nrOfVoxels = 0;
		int nrOfSeeds = 0;

		while (!inIter0.IsAtEnd())
		{
			double val0 = inIter0.Get(); double val1 = inIter1.Get(); double val2 = inIter2.Get(); double val3 = inIter3.Get();
			double entropy = val0*std::log(val0) + val1*std::log(val1) + val2*std::log(val2) + val3*std::log(val3);
			entropy = (-entropy) * normalizeFactor;
			outIter.Set(entropy);
			//outIterAll.Set(entropy + outIterAll.Get());

			std::vector<double> probs; probs.push_back(val0); probs.push_back(val1); probs.push_back(val2); probs.push_back(val3);
			std::sort(probs.begin(), probs.end());
			double FirstMinusSecond = probs[3] - probs[2];

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

		double successRate_l[4] = { 0.0 };
		double accuracyGT = 100.0 - (nrOferror / (nrOfVoxels / 100.0));
		double accuracySeeds = 100.0 - (nrOfSeedErrors / (nrOfSeeds / 100.0));

		for (int j = 0; j < 4; j++) {
			successRate_l[j] = 100.0 - (nrOferror_l[j] / (nrOferror_c[j] / 100.0));
		}

		std::cout << "Mean: " << mean << "Accuracy: " << accuracyGT / 100 << "Accuracy(Seeds): " << accuracySeeds / 100 << "Median: " << median << "Sigma: " << stdev << "MAD: " << MAD << "S0: " << successRate_l[0] << "S1: " << successRate_l[1] << "S2: " << successRate_l[2] << "S3: " << successRate_l[3] << std::endl;

		//x.push_back(mean);
		//y.push_back(accuracyGT / 100);

		cout << mean << ";" << accuracyGT / 100 << ";" << accuracySeeds / 100 << ";" << median << ";" << stdev << ";" << MAD << ";" << successRate_l[0] << ";" << successRate_l[1] << ";" << successRate_l[2] << ";" << successRate_l[3] << std::endl;

		//DoubleWriterType::Pointer writer = DoubleWriterType::New();
		//writer->ReleaseDataFlagOn();
		//writer->SetFileName(entPath.toStdString());
		//writer->SetInput(m_OutputImage);
		//writer->SetUseCompression(true);
		//writer->Update();

		//if (i == nrSegmen - 1){
		//	//set output image as the process output
		//	image->SetImage(m_OutputImage);
		//	image->Modified();
		//}
	}
	catch (itk::ExceptionObject &err)
	{
		throw err;
	}
}


void iAMMSegSampler::run()
{
	if (m_range->modalityParamRange.size() < 1)
	{
		DEBUG_LOG("At least one modality is required!\n");
		return;
	}
	StatusMsg("Generating sampling parameter sets...");
	m_parameters = m_sampleGenerator->GetParameterSets(m_range);
	if (!m_parameters)
	{
		DEBUG_LOG("No Parameters available!\n");
		return;
	}

	QMap<std::pair<int, int>, QSharedPointer<iASpectralVoxelData const> > m_pcaReducedSpectralData;

	// generate AttributeDescriptor:

	QString descriptor = QString("ERW Beta	Parameter	Continuous	%1	%2	%3\n")
		.arg(m_range->erw_beta_From)
		.arg(m_range->erw_beta_To)
		.arg(m_range->erw_beta_logScale ? "Logarithmic" : "Linear");
	descriptor += QString("ERW Gamma	Parameter	Continuous	%1	%2	%3\n")
		.arg(m_range->erw_gamma_From)
			.arg(m_range->erw_gamma_To)
			.arg(m_range->erw_gamma_logScale ? "Logarithmic" : "Linear");
	descriptor += QString("ERW Max.Iter.	Parameter	Discrete	%1	%2	%3\n")
		.arg(m_range->erw_maxIter_From)
			.arg(m_range->erw_maxIter_To)
			.arg(m_range->erw_maxIter_logScale ? "Logarithmic" : "Linear");
	descriptor += QString("SVM C	Parameter	Continuous	%1	%2	%3\n")
		.arg(m_range->svm_C_From)
			.arg(m_range->svm_C_To)
			.arg(m_range->svm_C_logScale ? "Logarithmic" : "Linear");
	descriptor += QString("SVM Gamma	Parameter	Continuous	%1	%2	%3\n")
			.arg(m_range->svm_gamma_From)
			.arg(m_range->svm_gamma_To)
			.arg(m_range->svm_gamma_logScale ? "Logarithmic" : "Linear");
	descriptor += QString("SVM Seed Set Probability	Parameter	Continuous	%1	%2	%3\n")
		.arg(m_range->svm_SeedProb_From)
		.arg(m_range->svm_SeedProb_To)
		.arg("Linear");
	descriptor += QString("GAD Iter	Parameter	Discrete	%1	%2	%3\n")
		.arg(m_range->gadIter_From)
		.arg(m_range->gadIter_To)
		.arg("Linear");
	descriptor += QString("GAD Step	Parameter	Continuous	%1	%2	%3\n")
		.arg(m_range->gadStep_From)
		.arg(m_range->gadStep_To)
		.arg("Linear");
	descriptor += QString("GAD Conductance	Parameter	Continuous	%1	%2	%3\n")
		.arg(m_range->gadCond_From)
		.arg(m_range->gadCond_To)
		.arg("Linear");
	descriptor += QString("SVM Channels	Parameter	Discrete	%1	%2	%3\n")
			.arg(m_range->svm_channels_From)
			.arg(m_range->svm_channels_To)
			.arg("Linear");

	for (int i = 0; i < m_range->modalityParamRange.size(); ++i)
	{
		QString distFuncStr;
		int distFncCnt = m_range->modalityParamRange[i].distanceFuncs.size();
		for (int j = 0; j < distFncCnt; ++j)
		{
			distFuncStr += m_range->modalityParamRange[i].distanceFuncs[j]->GetShortName();
			if (j < distFncCnt - 1)
			{
				distFuncStr += ",";
			}
		}
		descriptor += QString("Mod%1 PCA Components	Parameter	Discrete	%2	%3	Linear\n"
			"Mod%1 Weight	Parameter	Continuous	%4	%5	Linear\n"
			"Mod%1 Distance Functions	Parameter	Categorical	%6\n")
			.arg(i)
			.arg(m_range->modalityParamRange[i].pcaDimMin)
			.arg(m_range->modalityParamRange[i].pcaDimMax)
			.arg(m_range->modalityParamRange[i].weightFrom)
			.arg(m_range->modalityParamRange[i].weightTo)
			.arg(distFuncStr);
	}
	descriptor += QString("Object Count	Derived Output	Discrete	10000000000	1\n"
		"Duration	Derived Output	Continuous	10000000000	0\n"
		"MED A	Derived Output	Continuous	0	1\n"
		"M A	Derived Output	Continuous	0	1\n"
		"D	Derived Output	Continuous	0	1\n"
		"D 0	Derived Output	Continuous	0	1\n"
		"D 1	Derived Output	Continuous	0	1\n"
		"D 2	Derived Output	Continuous	0	1\n"
		"D 3	Derived Output	Continuous	0	1\n"
		"D 4	Derived Output	Continuous	0	1\n"
		"MED U	Derived Output	Continuous	0	1\n"
		"M U	Derived Output	Continuous	0	1\n"
		"U 0	Derived Output	Continuous	0	1\n"
		"U 1	Derived Output	Continuous	0	1\n"
		"U 2	Derived Output	Continuous	0	1\n"
		"U 3	Derived Output	Continuous	0	1\n"
		"U 4	Derived Output	Continuous	0	1\n"
		"MED C U F	Derived Output	Continuous	0	1\n"
		"MED C U T	Derived Output	Continuous	0	1\n"
		"MED C B G	Derived Output	Continuous	0	1\n");

	QTextStream textStream(&descriptor);

	m_attributes = iAAttributes::Create(textStream);

	m_derivedOutputStart = NonModalityParamCount + m_range->modalityParamRange.size() * ModalityParamCount;

	m_results = QSharedPointer<iASamplingResults>(new iASamplingResults(m_attributes, m_sampleGenerator->GetName()));
	for (m_curLoop=0; !m_aborted && m_curLoop<m_parameters->size(); ++m_curLoop)
	{
		QSharedPointer<iAMMSegParameter> param = (*m_parameters)[m_curLoop];
		param->setID(m_curLoop);

		typedef itk::Image< double, 3 > InputImageType;

		for (int iM = 0; iM < m_modalities->size(); iM++) {
			iAConnector* m_conn = new iAConnector;
			m_conn->SetImage(m_modalities->Get(iM)->GetImage());

			typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType, InputImageType > GGADIFType;
			typename GGADIFType::Pointer filter = GGADIFType::New();

			filter->SetNumberOfIterations(param->gaditer());
			filter->SetTimeStep(param->gadstep());
			filter->SetConductanceParameter(param->gadcond());
			filter->SetInput(dynamic_cast< InputImageType * >(m_conn->GetITKImage()));

			filter->Update();
			m_conn->SetImage(filter->GetOutput());
			m_modalities->Get(iM)->SetImage(m_conn->GetVTKImage());
			m_modalities->Get(iM)->GetImage()->Modified();
			filter->ReleaseDataFlagOn();
		}

		while (m_runningERW.size() >= CONCURRENT_ERW_RUNS && !m_aborted)
		{
			QThread::msleep(100);
		}
		if (m_aborted)
		{
			break;
		}
		StatusMsg(QString("SVM Run %1").arg(m_curLoop));
		iAPerformanceTimer svmPerfTimer;
		svmPerfTimer.start();

		SVMImageFilter svm(
			param->svm_c(),
			param->svm_gamma(),
			m_modalities,
			m_seeds,
			param->svm_channels());
		svm.Run();

		double svmTime = svmPerfTimer.elapsed();
		m_calcERWDuration += svmTime;
		
		
		StatusMsg(QString("SVM finished in %1 seconds; starting ERW Run %2")
			.arg(QString::number(svmTime))
			.arg(QString::number(m_curLoop)));

#if _DEBUG
		DEBUG_LOG(QString("%1: beta=%2, gamma=%3\n")
			.arg(m_curLoop)
			.arg(param->erw_beta())
			.arg(param->erw_gamma()));
		for (int i=0; i<param->modalityCount(); ++i)
		{
			DEBUG_LOG(QString(", modality %1 (weight=%2, metric=%3, pcaDim=%4)\n")
				.arg(i + 1)
				.arg(param->modalityWeight(i))
				.arg(param->distanceFunc(i)->GetShortName())
				.arg(param->pcaDim(i)));
		}
		DEBUG_LOG("\n");
#endif

		QSharedPointer<iAGaussianNormalizer> gauss(new iAGaussianNormalizer);
		gauss->SetBeta(param->erw_beta());
			
		QSharedPointer<QVector<iARWInputChannel> > inputChannels(new QVector<iARWInputChannel>());
		for (int i=0; i<param->modalityCount(); ++i)
		{
			iARWInputChannel input;
			QSharedPointer<iAModality const> modI = m_modalities->Get(i);
			if (modI->GetData()->channelCount() != param->pcaDim(i))
			{
				if (!m_pcaReducedSpectralData.contains(std::pair<int, int>(i, param->pcaDim(i))))
				{
					iAPCA pca(modI->GetData());
					m_pcaReducedSpectralData.insert(std::pair<int, int>(i, param->pcaDim(i)),
						pca.GetReduced(modI->GetConverter(), param->pcaDim(i)));
				}
				input.image = m_pcaReducedSpectralData[std::pair<int, int>(i, param->pcaDim(i))];
			}
			else
			{
				input.image = modI->GetData();
			}
			input.distanceFunc = param->distanceFunc(i);
			input.normalizeFunc = gauss;
			input.weight = param->modalityWeight(i);
			inputChannels->push_back(input);
		}
		QSharedPointer<iAModality const> mod0 = m_modalities->Get(0);
		iAExtendedRandomWalker* extendedRandomWalker = new iAExtendedRandomWalker(
			mod0->GetWidth(), mod0->GetHeight(), mod0->GetDepth(), mod0->GetSpacing(),
			inputChannels,
			svm.GetResult(),
			param->erw_gamma(),
			param->erw_maxIter()
		);
		m_erwPerfTimer.start();
		m_runningERW.insert(extendedRandomWalker, param);
		connect(extendedRandomWalker, SIGNAL(finished()), this, SLOT(erwFinished()) );
		
		m_mutex.lock();
		m_runningOperations++;
		m_mutex.unlock();
		extendedRandomWalker->start();

		cout << "m_outputBaseDir: " << m_outputBaseDir.toStdString() << endl;

		QString resultDir = m_outputBaseDir + "/sample" + QString::number(param->id());

		std::vector< double > uncert;

		std::vector<std::string> patientMHAs;
		QDir dir(m_outputBaseDir);
		//dir.cdUp();
		sublistFiles2(dir, "", patientMHAs);

		int nrModals = 4;

		GetUncertaintyValues(m_outputBaseDir, QString::fromStdString(patientMHAs[nrModals].c_str()), uncert);
		// use threadpool:
		//QThreadPool::globalInstance()->start(extendedRandomWalker);
	}
	if (m_aborted)
	{
		StatusMsg("Aborted by user!");
		return;
	}
	// wait for running operations to finish:
	while (m_runningOperations > 0)
	{
		msleep(100);
	}
}

void iAMMSegSampler::erwFinished()
{
	QObject* o = QObject::sender();
	iAExtendedRandomWalker* extendedRandomWalker = dynamic_cast<iAExtendedRandomWalker*>(o);
	if (!extendedRandomWalker)
	{
		DEBUG_LOG("FATAL: erwFinished - NULL or non-ERW-object sender!\n");
		m_aborted = true;
		return;
	}
	iAPerformanceTimer::DurationType time = m_erwPerfTimer.elapsed();
	m_calcERWDuration += time;
	if (extendedRandomWalker->GetError() != "")
	{
		DEBUG_LOG(QString("ERW Error: %1\n")
			.arg(extendedRandomWalker->GetError()));
		m_aborted = true;

		// we don't start characteristics calculation (at which's end we would do this otherwise):
		m_mutex.lock();
		m_runningOperations--;
		m_mutex.unlock();
		return;
	}
	QSharedPointer<iAMMSegParameter> param = m_runningERW[extendedRandomWalker];

	QSharedPointer<iARWResult> rwResult(extendedRandomWalker->GetResult());

	QString resultDir = m_outputBaseDir + "/sample" + QString::number(param->id());
	QSharedPointer<iASingleResult> result(new iASingleResult(param->id(), resultDir/*, time.count()*/));

	double value = -1;
	// set attributes from param:
	for (int i = 0; i < NonModalityParamCount; ++i)
	{
		switch (i)
		{
			case erwBeta: value = param->erw_beta(); break;
			case erwGamma: value = param->erw_gamma(); break;
			case erwMaxIter: value = param->erw_maxIter(); break;
			case svmC: value = param->svm_c(); break;
			case svmGamma: value = param->svm_gamma(); break;
			case svmSeedProb: value = param->svm_seedprob(); break;
			case gadIter: value = param->gaditer(); break;
			case gadStep: value = param->gadstep(); break;
			case gadCond: value = param->gadcond(); break;
			case svmChannelCount: value = param->svm_channels(); break;
		}
		result->SetAttribute(i, value);
	}
	int modalityCount = m_range->modalityParamRange.size();
	for (int m = 0; m < modalityCount; ++m)
	{
		for (int i = 0; i < ModalityParamCount; ++i)
		{
			switch (i)
			{
			case pca: value = param->pcaDim(m); break;
			case weight: value = param->modalityWeight(m); break;
			case distance: value = param->distanceFuncIdx(m); break;
			}
			int idx = NonModalityParamCount + (m*ModalityParamCount) + i;
			result->SetAttribute(idx, value);
		}
	}

	result->SetLabelImage(rwResult->labelledImage);
	result->AddProbabilityImages(rwResult->probabilityImages);
	result->SetAttribute(m_derivedOutputStart+1, time);
	m_results->GetAttributes()->Get(m_derivedOutputStart + 1)->AdjustMinMax(time);
	//m_range->AdaptDurationMinMax(time);

	CharacteristicsCalculator * newCharCalc = new CharacteristicsCalculator (result, m_results->GetAttributes(), m_derivedOutputStart);
	m_runningCharCalc.insert(newCharCalc, result);
	connect(newCharCalc, SIGNAL(finished()), this, SLOT(charactCalcFinished()) );
	m_ChaPerfTimer.start();
	newCharCalc->start();

	m_runningERW.remove(extendedRandomWalker);
	delete extendedRandomWalker;

}

#include "iAITKIO.h"
#include "iAToolsITK.h"

#include <QDir>

bool WriteResult(QSharedPointer<iASingleResult> singleResult, int labelCount, bool storeProbabilities)
{
	QString folder(singleResult->GetFolder());
	assert(singleResult);
	QDir dir(folder);
		
	if (dir.exists())
	{
		DEBUG_LOG(QString("ERROR: output directory '%1' already exists!\n").arg(folder));
		singleResult->DiscardDetails();
		return false;
	}
	dir.mkpath(".");
	QString labelFileName(folder + "\\label.mhd");
	itk::ImageIOBase::IOComponentType pixelType = GetITKScalarPixelType(singleResult->GetLabelledImage());
	iAITKIO::writeFile(labelFileName, singleResult->GetLabelledImage(), pixelType, true);
		
	for (int i = 0; storeProbabilities && (i<labelCount); ++i)
	{
		QString probFileName(folder + "\\prob" + QString::number(i) + ".mhd");
		itk::ImageIOBase::IOComponentType pixelType = GetITKScalarPixelType(singleResult->GetProbabilityImg(i));
		iAITKIO::writeFile(probFileName, singleResult->GetProbabilityImg(i), pixelType, true);
	}
	return true;
}


void iAMMSegSampler::charactCalcFinished()
{
	QObject* o = QObject::sender();
	CharacteristicsCalculator* charactCalc = dynamic_cast<CharacteristicsCalculator*>(o);
	if (!charactCalc)
	{
		DEBUG_LOG("ERROR: charactCalcFinished - NULL or on-CharacteristicsCalculator sender!\n");
		return;
	}
	QSharedPointer<iASingleResult> result = m_runningCharCalc[charactCalc];
	m_results->GetAttributes()->Get(m_derivedOutputStart)->AdjustMinMax(result->GetAttribute(m_derivedOutputStart));
	double time = m_ChaPerfTimer.elapsed();
	m_calcChaDuration += time;
	
	// TODO: pass in from somewhere! Or don't store here at all? but what in case of a power outage/error?
	QString sampleMetaFile      = m_outputBaseDir + "/" + m_parameterRangeFile;
	QString parameterSetFile    = m_outputBaseDir + "/" + m_parameterSetFile;
	QString characteristicsFile = m_outputBaseDir + "/" + m_characteristicsFile;
	if (WriteResult(result, m_seeds->size(), m_storeProbabilities))
	{
		result->DiscardDetails();
		m_results->AddResult(result);
		emit Progress((100*m_results->size()) / m_parameters->size());
		if (!m_results->Store(sampleMetaFile, parameterSetFile, characteristicsFile))
		{
			DEBUG_LOG("Error writing parameter file.\n");
		}
	}
	m_runningCharCalc.remove(charactCalc);
	delete charactCalc;
	m_mutex.lock();
	m_runningOperations--;
	m_mutex.unlock();
}

double iAMMSegSampler::elapsed() const
{
	return m_calcERWDuration+m_calcChaDuration;
}

double iAMMSegSampler::estimatedTimeRemaining() const
{
	return 
		((m_calcERWDuration+m_calcChaDuration)/(m_curLoop+1)) // average duration of one cycle
		* static_cast<double>(m_parameters->size()-m_curLoop-1) // remaining cycles
	;
}

QSharedPointer<iASamplingResults> iAMMSegSampler::GetResults()
{
	return m_results;
}

void iAMMSegSampler::Abort()
{
	m_aborted = true;
}

bool iAMMSegSampler::IsAborted()
{
	return m_aborted;
}

SVMImageFilter::SeedsPointer iAMMSegSampler::MethodForSelectingSeeds(QString pathToPool, float percentageSeeds) {
	SVMImageFilter::SeedsPointer seeds(new SVMImageFilter::SeedsType);
	std::vector<std::string> allMHAs;
	typedef itk::Image< double, 3 > InputImageType;

	{
		QDir dir(pathToPool);
		listFiles(dir, "", allMHAs);
	}

	cout << "nr of mhas: " << allMHAs.size() << endl;

	/* code snippet to select seeds */
	{
		for (std::vector<std::string>::iterator it = allMHAs.begin(); it != allMHAs.end(); ++it) {
			InputImageType::Pointer otImage = NULL;

			std::string curStr = *it;

			typedef itk::ImageFileReader< InputImageType > OTReaderType;
			OTReaderType::Pointer otreader = OTReaderType::New();
			otreader->SetFileName(curStr.c_str());
			otreader->Update();
			otreader->ReleaseDataFlagOn();

			otImage = otreader->GetOutput();

			/* Define seed points */

			typedef itk::ImageDuplicator< InputImageType > DuplicatorType;
			DuplicatorType::Pointer duplicator = DuplicatorType::New();
			duplicator->SetInputImage(otImage);
			duplicator->Update();
			InputImageType::Pointer clonedImage = duplicator->GetOutput();

			bool safe = false;
			float curProb = percentageSeeds / 100;
			std::string seedpathStr;

			/* Count number of voxel for every label */
			int counter[5] = { 0 };
			itk::ImageRegionConstIterator<InputImageType> imageIterator(clonedImage, clonedImage->GetLargestPossibleRegion());
			imageIterator.GoToBegin();
			while (!imageIterator.IsAtEnd()) {
				int val = imageIterator.Get();
				counter[val]++; ++imageIterator;
			}

			while (!safe) {
				typedef itk::SaltAndPepperNoiseImageFilter < InputImageType, InputImageType > SaltPepperFilterType;
				SaltPepperFilterType::Pointer saltpepper = SaltPepperFilterType::New();
				curProb += 0.00000001;
				saltpepper->SetProbability(curProb);
				saltpepper->SetInput(otImage);
				saltpepper->Update();

				std::size_t found = curStr.find_last_of("//");
				found = curStr.substr(0, found).find_last_of("//");
				std::string pathStr(curStr.substr(0, found).substr(0, found));
				seedpathStr = pathStr + "/seeds.seed";

				std::vector<std::vector<int>> labels[5];

				itk::ImageRegionConstIterator<InputImageType> imageIterator2(saltpepper->GetOutput(), saltpepper->GetOutput()->GetLargestPossibleRegion());

				imageIterator.GoToBegin();
				imageIterator2.GoToBegin();

				while (!imageIterator.IsAtEnd())
				{
					std::vector<int> pos;
					// Get the value of the current pixel
					int val = imageIterator.Get();
					int val2 = imageIterator2.Get();
					if (val != val2) {
						pos.push_back(imageIterator.GetIndex()[0]); pos.push_back(imageIterator.GetIndex()[1]); pos.push_back(imageIterator.GetIndex()[2]);
						if (val == 0) { labels[0].push_back(pos); }
						else if (val == 1) { labels[1].push_back(pos); }
						else if (val == 2) { labels[2].push_back(pos); }
						else if (val == 3) { labels[3].push_back(pos); }
						else if (val == 4) { labels[4].push_back(pos); }
					}
					++imageIterator;
					++imageIterator2;
				}

				cout << "Prob: " << curProb << "-";
				float percentages[5];
				for (int i = 0; i < 5; i++) { percentages[i] = float(labels[i].size()) / (float(counter[i]) / 100); cout << percentages[i] << ","; }
				cout << endl;

				//if (percentages[0] > 75 && percentages[1] > 75 && percentages[2] > 75 && percentages[3] > 75 && percentages[4] > 75){
				if (labels[0].size() > 0 && labels[1].size() > 0 && labels[2].size() > 0 && labels[3].size() > 0 && labels[4].size() > 0) {
					safe = true;

					ofstream myfile;
					myfile.open(seedpathStr);
					myfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
					myfile << "<Labels>" << endl;
					myfile << "\t<Label id=\"0\" name =\"0\">" << endl;

					QList<iAImageCoordinate> result0;
					for (std::vector<std::vector<int>>::iterator labelit = labels[0].begin(); labelit != labels[0].end(); ++labelit) {
						myfile << "\t\t<Seed x=\"" << (*labelit)[0] << "\" y=\"" << (*labelit)[1] << "\" z=\"" << (*labelit)[2] << "\"/>" << endl;
						iAImageCoordinate coord;
						coord.x = (*labelit)[0];
						coord.y = (*labelit)[1];
						coord.z = (*labelit)[2];
						result0.append(coord);
					}
					myfile << "\t</Label>" << endl;
					myfile << "\t<Label id=\"1\" name =\"1\">" << endl;
					QList<iAImageCoordinate> result1;
					for (std::vector<std::vector<int>>::iterator labelit = labels[1].begin(); labelit != labels[1].end(); ++labelit) {
						myfile << "\t\t<Seed x=\"" << (*labelit)[0] << "\" y=\"" << (*labelit)[1] << "\" z=\"" << (*labelit)[2] << "\"/>" << endl;
						iAImageCoordinate coord;
						coord.x = (*labelit)[0];
						coord.y = (*labelit)[1];
						coord.z = (*labelit)[2];
						result1.append(coord);
					}
					myfile << "\t</Label>" << endl;
					myfile << "\t<Label id=\"2\" name =\"2\">" << endl;

					QList<iAImageCoordinate> result2;
					for (std::vector<std::vector<int>>::iterator labelit = labels[2].begin(); labelit != labels[2].end(); ++labelit) {
						myfile << "\t\t<Seed x=\"" << (*labelit)[0] << "\" y=\"" << (*labelit)[1] << "\" z=\"" << (*labelit)[2] << "\"/>" << endl;
						iAImageCoordinate coord;
						coord.x = (*labelit)[0];
						coord.y = (*labelit)[1];
						coord.z = (*labelit)[2];
						result2.append(coord);
					}
					myfile << "\t</Label>" << endl;
					myfile << "\t<Label id=\"3\" name =\"3\">" << endl;

					QList<iAImageCoordinate> result3;
					for (std::vector<std::vector<int>>::iterator labelit = labels[3].begin(); labelit != labels[3].end(); ++labelit) {
						myfile << "\t\t<Seed x=\"" << (*labelit)[0] << "\" y=\"" << (*labelit)[1] << "\" z=\"" << (*labelit)[2] << "\"/>" << endl;
						iAImageCoordinate coord;
						coord.x = (*labelit)[0];
						coord.y = (*labelit)[1];
						coord.z = (*labelit)[2];
						result3.append(coord);
					}
					myfile << "\t</Label>" << endl;
					myfile << "\t<Label id=\"4\" name =\"4\">" << endl;

					QList<iAImageCoordinate> result4;
					for (std::vector<std::vector<int>>::iterator labelit = labels[4].begin(); labelit != labels[4].end(); ++labelit) {
						myfile << "\t\t<Seed x=\"" << (*labelit)[0] << "\" y=\"" << (*labelit)[1] << "\" z=\"" << (*labelit)[2] << "\"/>" << endl;

						iAImageCoordinate coord;
						coord.x = (*labelit)[0];
						coord.y = (*labelit)[1];
						coord.z = (*labelit)[2];
						result4.append(coord);
					}
					myfile << "\t</Label>" << endl;
					myfile << "</Labels>" << endl;
					cout << "Next ...." << endl;
					myfile.close();

					seeds->push_back(result0);
					seeds->push_back(result1);
					seeds->push_back(result2);
					seeds->push_back(result3);
					seeds->push_back(result4);
				}
			}

			cout << "Seed file: " << seedpathStr << " Prob: " << curProb << endl;
		}
	}

	return seeds;
}