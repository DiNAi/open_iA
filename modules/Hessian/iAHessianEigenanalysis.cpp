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
#include "iAHessianEigenanalysis.h"

#include "iAConnector.h"
#include "iAProgress.h"
#include "iATypedCallHelper.h"

#include <itkDerivativeImageFilter.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include <itkLaplacianSegmentationLevelSetFunction.h>
#include <itkLaplacianImageFilter.h>
#include <itkZeroCrossingImageFilter.h>

#include <vtkImageData.h>

#include <QLocale>

void find_and_replace(std::string& source, std::string const& find, std::string const& replace)
{
	for (std::string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
	{
		source.replace(i, find.length(), replace);
		i += replace.length();
	}
}

/**
* template computeHessian
*
* This template is used for calculating the hessian matrix
* \param	sigma			Sigma
* \param	hessianComputed Is the hessian matrix already computed.
* \param	p				Filter progress information.
* \param	image			Input image.
* \param	T				Input type
* \return	int Status-Code.
*/
template<class T> int computeHessian_template(int sigma, bool hessianComputed, int nr, iAProgress* p, iAConnector* image, std::string fileDir, std::string fileName)
{
	typedef itk::Image< T, 3 >	InputImageType;
	typedef itk::Image< int, 3 >	IntImageType;
	typedef itk::Image< float, 3 > FloatImageType;
	typedef itk::Image< double, 3 > DoubleImageType;

	std::vector<float> params;
	getParam(params);

	if (sigma == 101) {
		typedef  itk::CurvatureAnisotropicDiffusionImageFilter< InputImageType, InputImageType > SmoothingFilterType;
		SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
		smoothing->SetTimeStep(params[0]);
		smoothing->SetNumberOfIterations(params[1]);
		smoothing->SetConductanceParameter(params[2]);
		smoothing->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));

		image->SetImage(smoothing->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 102) {
		typedef  itk::GradientMagnitudeRecursiveGaussianImageFilter< InputImageType, InputImageType > GradientFilterType;
		GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
		gradientMagnitude->SetSigma(params[0]);
		gradientMagnitude->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));

		image->SetImage(gradientMagnitude->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 103)
	{
		typedef itk::SobelEdgeDetectionImageFilter <InputImageType, FloatImageType> SobelEdgeDetectionImageFilterType;
		SobelEdgeDetectionImageFilterType::Pointer sobelFilter = SobelEdgeDetectionImageFilterType::New();
		sobelFilter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));

		image->SetImage(sobelFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 104)
	{
		typedef itk::SigmoidImageFilter <InputImageType, FloatImageType> SigmoidImageFilterType;
		SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
		sigmoidFilter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		sigmoidFilter->SetOutputMinimum(params[2]);
		sigmoidFilter->SetOutputMaximum(params[3]);
		sigmoidFilter->SetAlpha(params[0]);
		sigmoidFilter->SetBeta(params[1]);

		image->SetImage(sigmoidFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 105)
	{
		int pos = fileDir.find_last_of("\\/");
		std::string modalDir = fileDir.substr(0, pos);
		std::vector<std::string> patientMHAs;
		QDir dir(modalDir.c_str());
		sublistFiles(dir, "", patientMHAs);
		int ot = patientMHAs.size();
		cout << "0: " << patientMHAs[ot-1] << endl;

		typedef itk::ConfidenceConnectedImageFilter<InputImageType, InputImageType> ConfidenceConnectedFilterType;
		ConfidenceConnectedFilterType::Pointer confidenceConnectedFilter = ConfidenceConnectedFilterType::New();
		confidenceConnectedFilter->SetInitialNeighborhoodRadius(params[0]);
		confidenceConnectedFilter->SetMultiplier(params[1]);
		confidenceConnectedFilter->SetNumberOfIterations(params[2]);
		confidenceConnectedFilter->SetReplaceValue(1);

		int label = params[3];

		// Set seed
		InputImageType::IndexType seed;
		seed[0] = params[4];
		seed[1] = params[5];
		seed[2] = params[6];
		confidenceConnectedFilter->SetSeed(seed);
		confidenceConnectedFilter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		confidenceConnectedFilter->Update();

		typedef itk::ImageFileReader<InputImageType> ReaderType;

		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(patientMHAs[ot-1]);
		reader->Update();

		itk::ImageRegionIterator<InputImageType> imageIteratorOutput(confidenceConnectedFilter->GetOutput(), confidenceConnectedFilter->GetOutput()->GetLargestPossibleRegion());
		itk::ImageRegionIteratorWithIndex<InputImageType> imageOT(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

		int error = 0;
		int all = 0;
		imageIteratorOutput.GoToBegin();
		imageOT.GoToBegin();

		while (!imageIteratorOutput.IsAtEnd()) {
			int value = imageIteratorOutput.Value();
			int valueOT = imageOT.Value();

			if (value == 1) {
				all++;
				if (valueOT != label) {
					error++;
				}
				else {
					cout << "\t\t<Seed x=\"" << imageOT.GetIndex()[0] << "\" y=\"" << imageOT.GetIndex()[1] << "\" z=\"" << imageOT.GetIndex()[2] << "\"/>" << endl;
				}
			}
			++imageIteratorOutput;
			++imageOT;
		}

		cout << "Anzahl fehlerhafter Pixel: " << error << " aus " << all << endl;

		image->SetImage(confidenceConnectedFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 106) {
		cout << "Begin 1" << endl;
		typedef  itk::CurvatureAnisotropicDiffusionImageFilter< InputImageType, FloatImageType > SmoothingFilterType;
		SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
		smoothing->SetTimeStep(0.06250);
		smoothing->SetNumberOfIterations(5.0);
		smoothing->SetConductanceParameter(15.0);
		smoothing->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		smoothing->Update();

		cout << "Begin 2" << endl;
		typedef  itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType > GradientFilterType;
		GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
		gradientMagnitude->SetSigma(params[0]);
		gradientMagnitude->SetInput(smoothing->GetOutput());
		gradientMagnitude->Update();

		cout << "Begin 3" << endl;
		typedef itk::SigmoidImageFilter <FloatImageType, FloatImageType> SigmoidImageFilterType;
		SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
		sigmoidFilter->SetInput(gradientMagnitude->GetOutput());
		sigmoidFilter->SetOutputMinimum(0);
		sigmoidFilter->SetOutputMaximum(1);
		sigmoidFilter->SetAlpha(params[1]);
		sigmoidFilter->SetBeta(params[2]);
		sigmoidFilter->Update();

		cout << "Begin 4" << endl;
		typedef  itk::FastMarchingImageFilter< FloatImageType, FloatImageType >    FastMarchingFilterType;
		typedef FastMarchingFilterType::NodeContainer  NodeContainer;
		typedef FastMarchingFilterType::NodeType       NodeType;
		NodeContainer::Pointer seeds = NodeContainer::New();
		InputImageType::IndexType  seedPosition;
		FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
		fastMarching->SetInput(sigmoidFilter->GetOutput());
		seedPosition[0] = params[3];
		seedPosition[1] = params[4];
		seedPosition[2] = params[5];
		const double initialDistance = params[6];
		NodeType node;
		const double seedValue = -initialDistance;
		node.SetValue(seedValue);
		node.SetIndex(seedPosition);
		seeds->Initialize();
		seeds->InsertElement(0, node);
		fastMarching->SetTrialPoints(seeds);
		fastMarching->SetSpeedConstant(1.0);
		fastMarching->SetOutputSize(sigmoidFilter->GetOutput()->GetBufferedRegion().GetSize());
		fastMarching->Update();

		cout << "Begin 5" << endl;
		typedef  itk::GeodesicActiveContourLevelSetImageFilter< FloatImageType, FloatImageType >    GeodesicActiveContourFilterType;
		GeodesicActiveContourFilterType::Pointer geodesicActiveContour = GeodesicActiveContourFilterType::New();
		geodesicActiveContour->SetPropagationScaling(params[7]);
		geodesicActiveContour->SetCurvatureScaling(params[8]);
		geodesicActiveContour->SetAdvectionScaling(params[9]);
		geodesicActiveContour->SetMaximumRMSError(0.02);
		geodesicActiveContour->SetNumberOfIterations(800);
		geodesicActiveContour->SetInput(fastMarching->GetOutput());
		geodesicActiveContour->SetFeatureImage(sigmoidFilter->GetOutput());

		image->SetImage(fastMarching->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 107) {
		typedef itk::SymmetricSecondRankTensor< double, 3 > HessianPixelType;
		typedef itk::Image< HessianPixelType, 3 >           HessianImageType;
		typedef itk::HessianToObjectnessMeasureImageFilter< HessianImageType, InputImageType > ObjectnessFilterType;
		ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
		objectnessFilter->SetBrightObject(false);
		objectnessFilter->SetScaleObjectnessMeasure(false);
		objectnessFilter->SetObjectDimension(params[6]);
		objectnessFilter->SetAlpha(params[0]);
		objectnessFilter->SetBeta(params[1]);
		objectnessFilter->SetGamma(params[2]);

		typedef itk::MultiScaleHessianBasedMeasureImageFilter< InputImageType, HessianImageType, InputImageType > MultiScaleEnhancementFilterType;
		MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
		multiScaleEnhancementFilter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		multiScaleEnhancementFilter->SetHessianToMeasureFilter(objectnessFilter);
		multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
		multiScaleEnhancementFilter->SetSigmaMinimum(params[3]);
		multiScaleEnhancementFilter->SetSigmaMaximum(params[4]);
		multiScaleEnhancementFilter->SetNumberOfSigmaSteps(params[5]);

		image->SetImage(multiScaleEnhancementFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 108) {
		typedef itk::HessianRecursiveGaussianImageFilter< InputImageType >  HessianFilterType;
		typedef HessianFilterType::OutputImageType                          HessianImageType;
		typedef HessianImageType::PixelType HessianPixelType;

		typedef  itk::FixedArray< double, HessianPixelType::Dimension >     EigenValueArrayType;
		typedef itk::Image< EigenValueArrayType, 3 > EigenValueImageType;

		typedef  itk::SymmetricEigenAnalysisImageFilter<HessianImageType, EigenValueImageType > EigenAnalysisFilterType;

		HessianFilterType::Pointer hessian = HessianFilterType::New();
		EigenAnalysisFilterType::Pointer eigen = EigenAnalysisFilterType::New();

		typedef itk::LocalStructureImageFilter< EigenValueImageType, FloatImageType >   LocalStructureFilterType;

		LocalStructureFilterType::Pointer localStructure = LocalStructureFilterType::New();

		hessian->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		eigen->SetInput(hessian->GetOutput());
		localStructure->SetInput(eigen->GetOutput());

		hessian->SetSigma(params[0]);
		localStructure->SetAlpha(params[1]);
		localStructure->SetGamma(params[2]);

		eigen->SetDimension(3);

		image->SetImage(localStructure->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 109) {
		typedef itk::DanielssonDistanceMapImageFilter< IntImageType, InputImageType, InputImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();

		typedef itk::ConnectedComponentImageFilter<InputImageType, IntImageType > LabelerType;
		LabelerType::Pointer labeler = LabelerType::New();

		labeler->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		filter->SetInput(labeler->GetOutput());

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 110) {
		typedef itk::FastChamferDistanceImageFilter< InputImageType, InputImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();

		filter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 111) {
		typedef itk::LaplacianImageFilter< DoubleImageType, DoubleImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();

		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 112) {
		typedef itk::RescaleIntensityImageFilter< DoubleImageType, DoubleImageType > RescaleFilterType;
		RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));
		rescaleFilter->SetOutputMinimum(params[0]);
		rescaleFilter->SetOutputMaximum(params[1]);

		image->SetImage(rescaleFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 113) {
		typedef itk::ExpNegativeImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));
		filter->SetFactor(params[0]);

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 114) {
		typedef itk::AtanImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 115) {
		typedef itk::SinImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 116) {
		typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 117) {
		typedef itk::SobelEdgeDetectionImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 118) {
		typedef itk::ZeroCrossingImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 119) {
		typedef itk::SimpleContourExtractorImageFilter <DoubleImageType, DoubleImageType> SimpleContourExtractorImageFilterType;
		SimpleContourExtractorImageFilterType::Pointer contourFilter = SimpleContourExtractorImageFilterType::New();
		contourFilter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));
		contourFilter->Update();

		image->SetImage(contourFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 120) {
		typedef itk::LaplacianSharpeningImageFilter< DoubleImageType, DoubleImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();

		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 121) {
		typedef itk::LaplacianRecursiveGaussianImageFilter< DoubleImageType, DoubleImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();

		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 122) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood
		IteratorType nit(radius, dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIterator<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		nit.GoToBegin();

		cout << "Neighbourhood size: " << nit.Size() << endl;

		while (!nit.IsAtEnd()) {
			std::vector<double> values;
			for (unsigned int i = 0; i < 27; i++) {
				bool IsInBounds;
				double value = nit.GetPixel(i, IsInBounds);
				if (IsInBounds) {
					values.push_back(value);
				}
			}

			double out = calcMedian(values);
			imageIteratorOutput.SetIndex(nit.GetIndex());
			imageIteratorOutput.Set(out);
			++nit;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 123) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood
		IteratorType nit(radius, dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIterator<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		nit.GoToBegin();

		cout << "Neighbourhood size: " << nit.Size() << endl;

		while (!nit.IsAtEnd()) {
			std::vector<double> values;
			for (unsigned int i = 0; i < 27; i++) {
				bool IsInBounds;
				double value = nit.GetPixel(i, IsInBounds);
				if (IsInBounds) {
					values.push_back(value);
				}
			}

			double out = calcMean(values);
			imageIteratorOutput.SetIndex(nit.GetIndex());
			imageIteratorOutput.Set(out);
			++nit;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 124) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood
		IteratorType nit(radius, dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIterator<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		nit.GoToBegin();

		cout << "Neighbourhood size: " << nit.Size() << endl;

		while (!nit.IsAtEnd()) {
			std::vector<double> values;
			for (unsigned int i = 0; i < 27; i++) {
				bool IsInBounds;
				double value = nit.GetPixel(i, IsInBounds);
				if (IsInBounds) {
					values.push_back(value);
				}
			}

			double out = calcVar(values);
			imageIteratorOutput.SetIndex(nit.GetIndex());
			imageIteratorOutput.Set(out);
			++nit;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 125) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood
		IteratorType nit(radius, dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIterator<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		nit.GoToBegin();

		cout << "Neighbourhood size: " << nit.Size() << endl;

		while (!nit.IsAtEnd()) {
			std::vector<double> values;
			for (unsigned int i = 0; i < 27; i++) {
				bool IsInBounds;
				double value = nit.GetPixel(i, IsInBounds);
				if (IsInBounds) {
					values.push_back(value);
				}
			}

			double out = calcSkewness(values);
			imageIteratorOutput.SetIndex(nit.GetIndex());
			imageIteratorOutput.Set(out);
			++nit;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 126) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood
		IteratorType nit(radius, dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIterator<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		nit.GoToBegin();

		cout << "Neighbourhood size: " << nit.Size() << endl;

		while (!nit.IsAtEnd()) {
			std::vector<double> values;
			for (unsigned int i = 0; i < 27; i++) {
				bool IsInBounds;
				double value = nit.GetPixel(i, IsInBounds);
				if (IsInBounds) {
					values.push_back(value);
				}
			}

			double out = calcKurtosis(values);
			imageIteratorOutput.SetIndex(nit.GetIndex());
			imageIteratorOutput.Set(out);
			++nit;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 127) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIterator(dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIterator.GoToBegin();
		imageIteratorOutput.GoToBegin();

		while (!imageIterator.IsAtEnd())
		{
			imageIteratorOutput.Set(imageIterator.GetIndex()[0]);

			++imageIterator;
			++imageIteratorOutput;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 128) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIterator(dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIterator.GoToBegin();
		imageIteratorOutput.GoToBegin();

		while (!imageIterator.IsAtEnd())
		{
			imageIteratorOutput.Set(imageIterator.GetIndex()[1]);

			++imageIterator;
			++imageIteratorOutput;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 129) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		typedef itk::NeighborhoodIterator<InputImageType> IteratorType;
		IteratorType::RadiusType radius;
		radius.Fill(1); // this is equivalent to 3x3x3 neigborhood

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIterator(dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIterator.GoToBegin();
		imageIteratorOutput.GoToBegin();

		while (!imageIterator.IsAtEnd())
		{
			imageIteratorOutput.Set(imageIterator.GetIndex()[2]);

			++imageIterator;
			++imageIteratorOutput;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 130) {
		typedef itk::ConnectedComponentImageFilter <DoubleImageType, IntImageType > ConnectedComponentImageFilterType;
		ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
		connected->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));
		connected->SetCoordinateTolerance(params[0] == 1);
		connected->SetDirectionTolerance(params[1] == 1);
		connected->SetFullyConnected(params[2] == 1);
		connected->Update();

		image->SetImage(connected->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 131) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIterator(dynamic_cast< InputImageType * >(image->GetITKImage()), dynamic_cast< InputImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIterator.GoToBegin();
		imageIteratorOutput.GoToBegin();

		double min = std::numeric_limits<float>::max(), max = std::numeric_limits<float>::min();

		while (!imageIterator.IsAtEnd()){
		
			if (imageIterator.Value() >= -std::numeric_limits<double>::max() && imageIterator.Value() <= std::numeric_limits<double>::max()){
				if (imageIterator.Value() < min) {
					min = imageIterator.Value();
				}
				if (imageIterator.Value() > max) {
					max = imageIterator.Value();
				}
			}

			++imageIterator;
		}

		cout << "Min: " << min << endl;
		cout << "Max: " << max << endl;

		imageIterator.GoToBegin();

		while (!imageIterator.IsAtEnd())
		{
			if (imageIterator.Value() >= -std::numeric_limits<double>::max() && imageIterator.Value() <= std::numeric_limits<double>::max()) {
				imageIteratorOutput.Set(imageIterator.Value()); 
			}
			else {
				imageIteratorOutput.Set(max);
			}

			++imageIterator;
			++imageIteratorOutput;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 132) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIteratorOutput.GoToBegin();

		typedef itk::SymmetricSecondRankTensor< double, 3 >  T_PixelType;
		typedef itk::HessianRecursiveGaussianImageFilter<InputImageType> HRGIFType;
		typedef typename HRGIFType::OutputImageType HessianImageType;
		typename HRGIFType::Pointer hessfilter = HRGIFType::New();
		hessfilter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		hessfilter->SetSigma(1);
		hessfilter->Update();

		itk::ImageRegionIteratorWithIndex<HessianImageType> hessIter(hessfilter->GetOutput(), hessfilter->GetOutput()->GetLargestPossibleRegion());

		hessIter.GoToBegin();
		while (!hessIter.IsAtEnd())
		{
			T_PixelType   hess_matrix;
			hess_matrix = hessIter.Get();

			imageIteratorOutput.Set(hess_matrix[(int)params[0]]);

			++hessIter;
			++imageIteratorOutput;
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 133) {
		
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIteratorOutput.GoToBegin();

		typedef itk::SymmetricSecondRankTensor< double, 3 >  T_PixelType;
		typedef itk::HessianRecursiveGaussianImageFilter<InputImageType> HRGIFType;
		typedef typename HRGIFType::OutputImageType HessianImageType;
		typename HRGIFType::Pointer hessfilter = HRGIFType::New();
		hessfilter->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
		hessfilter->SetSigma(1);
		hessfilter->Update();

		itk::ImageRegionIteratorWithIndex<HessianImageType> hessIter(hessfilter->GetOutput(), hessfilter->GetOutput()->GetLargestPossibleRegion());

		for (int i = 0; i < 6; i++) {
			hessIter.GoToBegin();
			imageIteratorOutput.GoToBegin();
			while (!hessIter.IsAtEnd())
			{
				T_PixelType   hess_matrix;
				hess_matrix = hessIter.Get();

				imageIteratorOutput.Set(hess_matrix[i]);

				++hessIter;
				++imageIteratorOutput;
			}

			QDir dir(QString::fromStdString(fileDir) + QString::number(i));
			if (!dir.exists()) {
				dir.mkpath(".");
			}
			QString writePath = QString::fromStdString(fileDir) + QString::number(i) + QString::fromStdString("/") + QString::fromStdString(fileName) + QString::number(i) + QString::fromStdString(".mha");

			typedef itk::RescaleIntensityImageFilter< InputImageType, InputImageType > RescaleFilterType;
			RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
			rescaleFilter->SetInput(newImage2);
			rescaleFilter->SetOutputMinimum(params[0]);
			rescaleFilter->SetOutputMaximum(params[1]);
			rescaleFilter->Update();

			cout << "Writing " << i << " to " << writePath.toStdString() << endl;
			StoreImage<InputImageType>(rescaleFilter->GetOutput(), writePath, false);
		}

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 134) {
		typedef itk::LogImageFilter<DoubleImageType, DoubleImageType> LogImageFilterType;
		LogImageFilterType::Pointer logfilter = LogImageFilterType::New();
		logfilter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));
		//logfilter->SetCoordinateTolerance(params[0] == 1);
		//logfilter->SetDirectionTolerance(params[1] == 1);
		//logfilter->SetFullyConnected(params[2] == 1);
		logfilter->Update();

		image->SetImage(logfilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 135) {
		typedef itk::LaplacianRecursiveGaussianImageFilter< DoubleImageType, DoubleImageType > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(dynamic_cast< DoubleImageType * >(image->GetITKImage()));
		filter->Update();

		int i = 6;

		QDir dir(QString::fromStdString(fileDir) + QString::number(i));
		if (!dir.exists()) {
			dir.mkpath(".");
		}
		QString writePath = QString::fromStdString(fileDir) + QString::number(i) + QString::fromStdString("/") + QString::fromStdString(fileName) + QString::number(i) + QString::fromStdString(".mha");

		typedef itk::RescaleIntensityImageFilter< DoubleImageType, InputImageType > RescaleFilterType;
		RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput(filter->GetOutput());
		rescaleFilter->SetOutputMinimum(params[0]);
		rescaleFilter->SetOutputMaximum(params[1]);
		rescaleFilter->Update();

		cout << "Writing " << i << " to " << writePath.toStdString() << endl;
		StoreImage<InputImageType>(rescaleFilter->GetOutput(), writePath, false);

		image->SetImage(filter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 136)
	{
		int pos = fileDir.find_last_of("\\/");
		std::string modalDir = fileDir.substr(0, pos);

		std::vector<std::string> patientMHAs;
		QDir dir(modalDir.c_str());
		sublistFiles(dir, "", patientMHAs);
		cout << "Size: " << patientMHAs.size() << endl;
		cout << "0: " << patientMHAs[0] << endl;

		typedef itk::ConfidenceConnectedImageFilter<InputImageType, FloatImageType> ConfidenceConnectedFilterType;
		ConfidenceConnectedFilterType::Pointer confidenceConnectedFilter = ConfidenceConnectedFilterType::New();
		for (int i = 0; i < 5; i++) {
			confidenceConnectedFilter = ConfidenceConnectedFilterType::New();
			confidenceConnectedFilter->SetInitialNeighborhoodRadius(params[3]);
			confidenceConnectedFilter->SetMultiplier(params[4]);
			confidenceConnectedFilter->SetNumberOfIterations(params[5]);
			confidenceConnectedFilter->SetReplaceValue(1);

			// Set seed
			InputImageType::IndexType seed;
			seed[0] = params[0];
			seed[1] = params[1];
			seed[2] = params[2];
			confidenceConnectedFilter->SetSeed(seed);
			confidenceConnectedFilter->SetInput(dynamic_cast<InputImageType *>(image->GetITKImage()));
		}

		image->SetImage(confidenceConnectedFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 137){
		int sumOt = 0;
		int pos = fileDir.find_last_of("\\/");
		std::string modalDir = fileDir.substr(0, pos);
		std::vector<std::string> patientMHAs;
		QDir dir(modalDir.c_str());
		sublistFiles(dir, "", patientMHAs);
		int ot = patientMHAs.size();
		//cout << "0: " << patientMHAs[ot - 1] << endl;

		typedef itk::ConfidenceConnectedImageFilter<InputImageType, InputImageType> ConfidenceConnectedFilterType;
		ConfidenceConnectedFilterType::Pointer confidenceConnectedFilter;

		typedef itk::ImageFileReader<InputImageType> ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(patientMHAs[ot - 1]);
		reader->Update();

		itk::ImageRegionIterator<InputImageType> imageOT(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
		imageOT.GoToBegin();

		while (!imageOT.IsAtEnd()) {
			int valueOT = imageOT.Value();
			if (valueOT == params[3]) {
				sumOt++;
			}
			++imageOT;
		}

		confidenceConnectedFilter = ConfidenceConnectedFilterType::New();
		confidenceConnectedFilter->SetInput(dynamic_cast<InputImageType *>(image->GetITKImage()));
		InputImageType::IndexType seed;
		seed[0] = params[4];
		seed[1] = params[5];
		seed[2] = params[6];
		confidenceConnectedFilter->SetSeed(seed);

		for (int nh = 0; nh < 11; nh++) {
			for (double mult = 0.1; mult < 3.1; mult = mult + 0.1) {
				for (int iter = 0; iter < 31; iter++) {
					confidenceConnectedFilter->SetInitialNeighborhoodRadius(nh);
					confidenceConnectedFilter->SetMultiplier(mult);
					confidenceConnectedFilter->SetNumberOfIterations(iter);
					confidenceConnectedFilter->SetReplaceValue(1);

					int label = params[3];

					// Set seed
					confidenceConnectedFilter->Update();

					itk::ImageRegionIterator<InputImageType> imageIteratorOutput(confidenceConnectedFilter->GetOutput(), confidenceConnectedFilter->GetOutput()->GetLargestPossibleRegion());
					itk::ImageRegionIterator<InputImageType> imageOT(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

					int error = 0;
					int all = 0;
					imageIteratorOutput.GoToBegin();
					imageOT.GoToBegin();

					while (!imageIteratorOutput.IsAtEnd()) {
						int value = imageIteratorOutput.Value();
						int valueOT = imageOT.Value();
						if (value == 1) {
							all++;
							if (valueOT != label) {
								error++;
							}
						}
						++imageIteratorOutput;
						++imageOT;
					}

					if (error < 1000) {
						cout << nh << "," << mult << "," << iter << "," << error << "," << all << "," << sumOt << endl;
					}
				}
			}
		}

		image->SetImage(confidenceConnectedFilter->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 138)
	{
		typedef itk::VectorImage<double, 3>  VectorImageType;
		typedef itk::VectorConfidenceConnectedImageFilter< VectorImageType, IntImageType > ConnectedFilterType;
		ConnectedFilterType::Pointer confidenceConnected = ConnectedFilterType::New();

		int pos = fileDir.find_last_of("\\/");
		std::string modalDir = fileDir.substr(0, pos);
		std::vector<std::string> patientMHAs;
		QDir dir(modalDir.c_str());
		sublistFiles(dir, "", patientMHAs);
		int ot = patientMHAs.size();

		typedef itk::ImageFileReader<InputImageType> ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(patientMHAs[ot - 1]);
		reader->Update();

		typedef itk::ComposeImageFilter<InputImageType, VectorImageType> ImageToVectorImageFilterType;
		ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

		typedef itk::ImageFileReader<InputImageType> ReaderType;
		ReaderType::Pointer reader2;

		for (int i = 0; i < ot - 1; i++) {
			reader2 = ReaderType::New();
			reader2->SetFileName(patientMHAs[i]);
			reader2->Update();

			imageToVectorImageFilter->SetInput(i, reader2->GetOutput());
		}
		imageToVectorImageFilter->Update();

		for (int iter = 1; iter < 3; iter++) {
			for (int nh = 1; nh < 3; nh++) {
				for (double mult = 0.01; mult < 3; mult = mult + 0.01) {
					confidenceConnected = ConnectedFilterType::New();
					confidenceConnected->SetInitialNeighborhoodRadius(nh);
					confidenceConnected->SetMultiplier(mult);
					confidenceConnected->SetNumberOfIterations(iter);
					confidenceConnected->SetReplaceValue(1);

					int label = params[3];

					// Set seed
					VectorImageType::IndexType seed;
					seed[0] = params[4];
					seed[1] = params[5];
					seed[2] = params[6];
					confidenceConnected->SetSeed(seed);
					confidenceConnected->SetInput(imageToVectorImageFilter->GetOutput());
					confidenceConnected->Update();

					itk::ImageRegionIterator<IntImageType> imageIteratorOutput(confidenceConnected->GetOutput(), confidenceConnected->GetOutput()->GetLargestPossibleRegion());
					itk::ImageRegionIterator<InputImageType> imageOT(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

					int error = 0;
					int all = 0;
					imageIteratorOutput.GoToBegin();
					imageOT.GoToBegin();

					while (!imageIteratorOutput.IsAtEnd()) {
						int value = imageIteratorOutput.Value();
						int valueOT = imageOT.Value();
						if (value == 1) {
							all++;
							if (valueOT != label) {
								error++;
							}
						}
						++imageIteratorOutput;
						++imageOT;
					}
					//if (error < 10000) {
						cout << iter << "," << nh << "," << mult << ":" << error << "," << all << endl;
					//}
				}
			}
		}

		image->SetImage(confidenceConnected->GetOutput());
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 139)
	{
		int pos = fileDir.find_last_of("\\/");

		std::string modalDir = fileDir.substr(0, pos);
		std::vector<std::string> patientMHAs;
		QDir dir(modalDir.c_str());
		sublistFiles(dir, "", patientMHAs);
		int ot = patientMHAs.size();

		typedef itk::ConfidenceConnectedImageFilter<InputImageType, InputImageType> ConfidenceConnectedFilterType;
		std::vector<ConfidenceConnectedFilterType::Pointer> confidenceConnectedFilter(4);

		typedef itk::ImageFileReader<InputImageType> ReaderType;
		ReaderType::Pointer reader2;

		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(patientMHAs[ot - 1]);
		reader->Update();

		InputImageType::Pointer newImage;
		newImage = InputImageType::New();
		newImage->SetRegions(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetLargestPossibleRegion());
		newImage->SetSpacing(dynamic_cast< DoubleImageType * >(image->GetITKImage())->GetSpacing());
		newImage->Allocate();

		int order[4] = {5, 0, 6, 5};

		for (int i = 0, j = 0; i < 28; i = i + 7, j++) {

			confidenceConnectedFilter[j] = ConfidenceConnectedFilterType::New();
			confidenceConnectedFilter[j]->SetInitialNeighborhoodRadius(	params[i + 0]);
			confidenceConnectedFilter[j]->SetMultiplier(				params[i + 1]);
			confidenceConnectedFilter[j]->SetNumberOfIterations(		params[i + 2]);
			confidenceConnectedFilter[j]->SetReplaceValue(				params[i + 3]);

			int label = params[i + 3];

			reader2 = ReaderType::New();
			reader2->SetFileName(patientMHAs[order[j]]);
			reader2->Update();

			// Set seed
			InputImageType::IndexType seed;
			seed[0] = params[i + 4];
			seed[1] = params[i + 5];
			seed[2] = params[i + 6];
			confidenceConnectedFilter[j]->SetSeed(seed);
			confidenceConnectedFilter[j]->SetInput(reader2->GetOutput());
			confidenceConnectedFilter[j]->Update();

			itk::ImageRegionIterator<InputImageType> imageIteratorOutput(confidenceConnectedFilter[j]->GetOutput(), confidenceConnectedFilter[j]->GetOutput()->GetLargestPossibleRegion());
			itk::ImageRegionIteratorWithIndex<InputImageType> imageOT(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
			itk::ImageRegionIterator<InputImageType> imageIteratorNewImage(newImage, newImage->GetLargestPossibleRegion());

			int error = 0;
			int all = 0;
			imageIteratorOutput.GoToBegin();
			imageOT.GoToBegin();
			imageIteratorNewImage.GoToBegin();

			while (!imageIteratorOutput.IsAtEnd()) {
				int value = imageIteratorOutput.Value();
				int valueOT = imageOT.Value();

				if (value == label) {
					all++;
					if (valueOT != label) {
						error++;
					}
					//else {
					//	cout << "\t\t<Seed x=\"" << imageOT.GetIndex()[0] << "\" y=\"" << imageOT.GetIndex()[1] << "\" z=\"" << imageOT.GetIndex()[2] << "\"/>" << endl;
					//}
					imageIteratorNewImage.Set(label);
				}
				++imageIteratorOutput;
				++imageOT;
				++imageIteratorNewImage;
			}

			cout << "\n Anzahl fehlerhafter Pixel: " << error << " aus " << all << " fuer " << label << endl;
		}

		image->SetImage(newImage);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 140) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIteratorOutput.GoToBegin();

		imageIteratorOutput.GoToBegin();
		while (!imageIteratorOutput.IsAtEnd())
		{
			imageIteratorOutput.Set(imageIteratorOutput.GetIndex()[0]);
			++imageIteratorOutput;
		}

		cout << fileName << endl;
		std::string fileName_X(fileName);
		find_and_replace(fileName_X, "Flair", "Flair_X");
		cout << endl << "fileName_X: " << fileName_X << endl;
		std::string fileName_X_Path = fileName_X.substr(0, fileName_X.find_last_of("\\/"));
		cout << "fileName_X_Path: " << fileName_X_Path << endl;
		QDir dir(QString::fromStdString(fileName_X_Path));
		if (!dir.exists()) {
			dir.mkpath(".");
		}
		QString writePath = QString::fromStdString(fileName_X);
		StoreImage<InputImageType>(newImage2, writePath, false);

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 141) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIteratorOutput.GoToBegin();

		imageIteratorOutput.GoToBegin();
		while (!imageIteratorOutput.IsAtEnd())
		{
			imageIteratorOutput.Set(imageIteratorOutput.GetIndex()[1]);
			++imageIteratorOutput;
		}

		QString writePath = "D:\\Workspace\\Workspace_BRATS\\brats\\Test_cases\\testcase2\\different_brains\\1_train_0\\brats_tcia_pat226_0090\\VSD.Brain.XX.O.MR_Flair_Y.40923\\VSD.Brain.XX.O.MR_Flair_Y.40923.mha";
		StoreImage<InputImageType>(newImage2, writePath, false);

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 142) {
		InputImageType::Pointer newImage2;
		newImage2 = InputImageType::New();
		newImage2->SetRegions(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetLargestPossibleRegion());
		newImage2->SetSpacing(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetSpacing());
		newImage2->Allocate();

		itk::ImageRegionIteratorWithIndex<InputImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
		imageIteratorOutput.GoToBegin();

		imageIteratorOutput.GoToBegin();
		while (!imageIteratorOutput.IsAtEnd())
		{
			imageIteratorOutput.Set(imageIteratorOutput.GetIndex()[2]);
			++imageIteratorOutput;
		}

		QString writePath = "D:\\Workspace\\Workspace_BRATS\\brats\\Test_cases\\testcase2\\different_brains\\1_train_0\\brats_tcia_pat226_0090\\VSD.Brain.XX.O.MR_Flair_Z.40923\\VSD.Brain.XX.O.MR_Flair_Z.40923.mha";
		StoreImage<InputImageType>(newImage2, writePath, false);

		image->SetImage(newImage2);
		image->Modified();

		return EXIT_SUCCESS;
	}
	if (sigma == 143) {
		std::vector<std::string> XYZ;
		XYZ.push_back("Flair_X"); XYZ.push_back("Flair_Y"); XYZ.push_back("Flair_Z");

		for (int i = 0; i < 3; i++) {
			DoubleImageType::Pointer newImage2;
			newImage2 = DoubleImageType::New();
			newImage2->SetRegions(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetLargestPossibleRegion());
			newImage2->SetSpacing(dynamic_cast<DoubleImageType *>(image->GetITKImage())->GetSpacing());
			newImage2->Allocate();

			itk::ImageRegionIteratorWithIndex<DoubleImageType> imageIteratorOutput(newImage2, newImage2->GetLargestPossibleRegion());
			imageIteratorOutput.GoToBegin();

			while (!imageIteratorOutput.IsAtEnd())
			{
				imageIteratorOutput.Set(imageIteratorOutput.GetIndex()[i]);
				++imageIteratorOutput;
			}

			cout << fileName << endl;
			std::string fileName_New(fileName);
			find_and_replace(fileName_New, "Flair", XYZ[i]);

			std::string fileName_New_Path = fileName_New.substr(0, fileName_New.find_last_of("\\/"));
			cout << "fileName_X_Path: " << fileName_New_Path << endl;
			QDir dir(QString::fromStdString(fileName_New_Path));
			if (!dir.exists()) {
				dir.mkpath(".");
			}
			QString writePath = QString::fromStdString(fileName_New);

			typedef itk::RescaleIntensityImageFilter< DoubleImageType, DoubleImageType > RescaleFilterType;
			RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
			rescaleFilter->SetInput(newImage2);
			rescaleFilter->SetOutputMinimum(params[0]);
			rescaleFilter->SetOutputMaximum(params[1]);

			StoreImage<DoubleImageType>(rescaleFilter->GetOutput(), writePath, false);
			image->SetImage(newImage2);
		}

		//image->SetImage(newImage2);
		//image->Modified();

		return EXIT_SUCCESS;
	}

	/************************************** Data and type definitions *********************/
	typedef itk::Vector<double, 3>	VectorPixelType;
	VectorPixelType eigenTempVector;
	typedef itk::Image< T, 3 >	InputImageType;
	typedef	itk::HessianRecursiveGaussianImageFilter<InputImageType >	HessianFilterType;
	typedef	typename HessianFilterType::OutputImageType	HessianImageType;

	typedef	typename HessianImageType::PixelType HessianPixelType;
	typedef	itk::FixedArray< double, HessianPixelType::Dimension > EigenValueArrayType;
	typedef	itk::Image< EigenValueArrayType, HessianImageType::ImageDimension > EigenValueImageType;
	typedef	itk::SymmetricEigenAnalysisImageFilter< HessianImageType, EigenValueImageType >     EigenAnalysisFilterType;

	itkStaticConstMacro(Dimension, unsigned int, 3);
	typedef	float	MeshPixelType;
	typedef itk::ImageAdaptor<  EigenValueImageType, EigenValueAccessor< EigenValueArrayType > > ImageAdaptorType;
	typedef itk::Image< MeshPixelType, Dimension > EachEigenValueImageType;
	typedef itk::CastImageFilter< ImageAdaptorType, EachEigenValueImageType >  CastImageFilterType;
	/************************************** Data and type definitions end *****************/

	/************************************** Hessian part **********************************/
	typename HessianFilterType::Pointer	m_Hessian;	// In m_Hessian werden die Matrizen f�r jedes Voxel gepsichert. MA
	m_Hessian = HessianFilterType::New();	// Objekt f�r die Matrizen deklariert. MA

	m_Hessian->SetInput(dynamic_cast< InputImageType * >(image->GetITKImage()));
	m_Hessian->SetSigma(sigma);
	if (!hessianComputed)
		m_Hessian->Update();
	/************************************** Hessian part end ******************************/

	/************************************** Eigen values part *****************************/
	// Compute eigen values.. order them in ascending order
	typename EigenAnalysisFilterType::Pointer m_EigenFilter;

	m_EigenFilter = EigenAnalysisFilterType::New();
	m_EigenFilter->SetDimension(HessianPixelType::Dimension);
	m_EigenFilter->SetInput(m_Hessian->GetOutput());
	typename HessianImageType::Pointer hessianImage = m_Hessian->GetOutput();
	m_EigenFilter->OrderEigenValuesBy(EigenAnalysisFilterType::FunctorType::OrderByValue);
	/************************************** Eigen values part end **********************************/

	/************************************** Eigen analysis **********************************/
	// Conversion or extraction of data from class itk::Image<class itk::FixedArray<double,3>,3> * __ptr64 to class itk::Image<float,3> * __ptr64
	typename ImageAdaptorType::Pointer               m_EigenAdaptor1;
	typename ImageAdaptorType::Pointer               m_EigenAdaptor2;
	typename ImageAdaptorType::Pointer               m_EigenAdaptor3;

	typename CastImageFilterType::Pointer            m_EigenCastfilter1;
	typename CastImageFilterType::Pointer            m_EigenCastfilter2;
	typename CastImageFilterType::Pointer            m_EigenCastfilter3;

	typename EigenValueImageType::Pointer eigenImage1;

	// Create an adaptor and plug the output to the parametric space
	m_EigenAdaptor1 = ImageAdaptorType::New();
	EigenValueAccessor< EigenValueArrayType > accessor1;
	accessor1.SetEigenIdx(0);
	m_EigenAdaptor1->SetImage(m_EigenFilter->GetOutput());
	m_EigenAdaptor1->SetPixelAccessor(accessor1);

	m_EigenAdaptor2 = ImageAdaptorType::New();
	EigenValueAccessor< EigenValueArrayType > accessor2;
	accessor2.SetEigenIdx(1);
	m_EigenAdaptor2->SetImage(m_EigenFilter->GetOutput());
	m_EigenAdaptor2->SetPixelAccessor(accessor2);

	m_EigenAdaptor3 = ImageAdaptorType::New();
	EigenValueAccessor< EigenValueArrayType > accessor3;
	accessor3.SetEigenIdx(2);
	m_EigenAdaptor3->SetImage(m_EigenFilter->GetOutput());
	m_EigenAdaptor3->SetPixelAccessor(accessor3);

	eigenImage1 = m_EigenFilter->GetOutput();
	typename EigenValueImageType::Pointer eigenRaRbS = m_EigenFilter->GetOutput();

	// m_EigenCastfilter1 will give the eigen values with the maximum eigen
	// value. m_EigenCastfilter3 will give the eigen values with the 
	// minimum eigen value.
	m_EigenCastfilter1 = CastImageFilterType::New();
	m_EigenCastfilter1->SetInput(m_EigenAdaptor3);
	m_EigenCastfilter2 = CastImageFilterType::New();
	m_EigenCastfilter2->SetInput(m_EigenAdaptor2);
	m_EigenCastfilter3 = CastImageFilterType::New();
	m_EigenCastfilter3->SetInput(m_EigenAdaptor1);
	/************************************** Eigen analysis part end **********************************/


	/************************************** Define output **********************************/
	if (nr == 1 || nr == 0) {
		p->Observe(m_EigenCastfilter1);
		m_EigenCastfilter1->Update();
		image->SetImage(m_EigenCastfilter1->GetOutput());
		m_EigenCastfilter1->ReleaseDataFlagOn();
	}
	else if (nr == 2) {
		p->Observe(m_EigenCastfilter2);
		m_EigenCastfilter2->Update();
		image->SetImage(m_EigenCastfilter2->GetOutput());
		m_EigenCastfilter2->ReleaseDataFlagOn();
	}
	else if (nr == 3) {
		p->Observe(m_EigenCastfilter3);
		m_EigenCastfilter3->Update();
		image->SetImage(m_EigenCastfilter3->GetOutput());
		m_EigenCastfilter3->ReleaseDataFlagOn();
	}
	/************************************** Define output end **********************************/
	image->Modified();

	/************************************** Iteration through eigenvalues to compute Ra, Rb and S values **********************************/
	typedef itk::ImageRegionIterator<EigenValueImageType> EigenIteratorType;
	EigenIteratorType eigenImageIterator(eigenRaRbS, eigenRaRbS->GetLargestPossibleRegion());
	EigenIteratorType eigenImageIt(eigenImage1, eigenImage1->GetLargestPossibleRegion());

	// iterate through image and get each eigen value
	int j = 0;
	eigenImageIterator.GoToBegin();
	for (eigenImageIt.GoToBegin(); !eigenImageIt.IsAtEnd() && !eigenImageIterator.IsAtEnd(); ++eigenImageIt)
	{
		std::cout << j << " -> ";
		EigenValueArrayType eigenArray = eigenImageIt.Get();

		eigenTempVector[0] = fabs(eigenArray[1]) / fabs(eigenArray[2]);
		eigenTempVector[1] = fabs(eigenArray[0]) / (sqrt(fabs(eigenArray[1] * eigenArray[2])));
		eigenTempVector[2] = sqrt(pow(eigenArray[0], 2) + pow(eigenArray[1], 2) + pow(eigenArray[2], 2));

		eigenImageIterator.Set(eigenTempVector);
		++eigenImageIterator;
		j++;
	}
	/************************************** Iteration through eigenvalues to compute Ra, Rb and S values end **********************************/

	return EXIT_SUCCESS;
}



/**
* template computeLaplacian
* This template is used for calculating the Laplacian of Gaussian (LoG)
* \param	sigma			Sigma
* \param	p				Filter progress information.
* \param	image			Input image.
* \param	T				Input type
* \return	int Status-Code.
*/

template<class T> int computeLaplacian_template(unsigned int sigma, iAProgress* p, iAConnector* image)
{
	typedef itk::Image< T, DIM > ImageType;
	typedef itk::Image<float, DIM> OutputImageType;
	typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, OutputImageType> LoGFilterType;

	typename LoGFilterType::RealType sigmaType(sigma);

	typename LoGFilterType::Pointer filter = LoGFilterType::New();
	filter->SetInput(dynamic_cast< ImageType * >(image->GetITKImage()));
	filter->SetSigma(sigmaType);

	image->SetImage(filter->GetOutput());
	image->Modified();

	return EXIT_SUCCESS;
}


iAHessianEigenanalysis::iAHessianEigenanalysis(QString fn, FilterID fid, vtkImageData* i, vtkPolyData* p, iALogger* logger, QObject* parent)
	: iAFilter(fn, fid, i, p, logger, parent)
{

}


iAHessianEigenanalysis::~iAHessianEigenanalysis()
{
}


void iAHessianEigenanalysis::run()
{
	switch (getFilterID())
	{
	case COMPUTEHESSIANEIGENANALYSIS:
		computeHessian(); break;
	case COMPUTE_LAPLACIAN:
		computeLaplacian(); break;
	case UNKNOWN_FILTER:
	default:
		addMsg(tr("  unknown filter type"));
	}
}


void iAHessianEigenanalysis::computeHessian()
{
	addMsg(tr("%1  %2 started.").arg(QLocale().toString(Start(), QLocale::ShortFormat))
		.arg(getFilterName()));
	getConnector()->SetImage(getVtkImageData()); getConnector()->Modified();

	try
	{
		VTK_TYPED_CALL(computeHessian_template, getVtkImageData()->GetScalarType(),
			sigma, hessianComputed, nr, getItkProgress(), getConnector(), fileDir, fileName);
	}
	catch (itk::ExceptionObject &excep)
	{
		addMsg(tr("%1  %2 terminated unexpectedly. Elapsed time: %3 ms. For learning only.").arg(QLocale().toString(QDateTime::currentDateTime(), QLocale::ShortFormat))
			.arg(getFilterName())
			.arg(Stop()));
		addMsg(tr("  %1 in File %2, Line %3. For learning only.").arg(excep.GetDescription())
			.arg(excep.GetFile())
			.arg(excep.GetLine()));
		return;
	}
	addMsg(tr("%1  %2 finished. Elapsed time: %3 ms").arg(QLocale().toString(QDateTime::currentDateTime(), QLocale::ShortFormat))
		.arg(getFilterName())
		.arg(Stop()));

	emit startUpdate();
}


void iAHessianEigenanalysis::computeLaplacian()
{
	addMsg(tr("%1  %2 started. Sigma %3").arg(QLocale().toString(Start(), QLocale::ShortFormat))
		.arg(getFilterName())
		.arg(this->sigma)
	);
	getConnector()->SetImage(getVtkImageData()); getConnector()->Modified();

	try
	{
		VTK_TYPED_CALL(computeLaplacian_template, getVtkImageData()->GetScalarType(),
			this->sigma, getItkProgress(), getConnector());
	}
	catch (itk::ExceptionObject &excep)
	{
		addMsg(tr("%1  %2 terminated unexpectedly. Elapsed time: %3 ms. ").arg(QLocale().toString(QDateTime::currentDateTime(), QLocale::ShortFormat))
			.arg(getFilterName())
			.arg(Stop()));
		addMsg(tr("  %1 in File %2, Line %3. ").arg(excep.GetDescription())
			.arg(excep.GetFile())
			.arg(excep.GetLine()));
		return;
	}
	addMsg(tr("%1  %2 finished. Elapsed time: %3 ms").arg(QLocale().toString(QDateTime::currentDateTime(), QLocale::ShortFormat))
		.arg(getFilterName())
		.arg(Stop()));

	emit startUpdate();
}

void getParam(std::vector<float> &params) {
	ifstream infile("D:\\params.csv");
	std::string sLine;
	if (infile.good()) {
		getline(infile, sLine);
		//cout << sLine << endl;
	}
	infile.close();

	std::vector<std::string> v;
	std::istringstream buf(sLine);
	for (std::string token; getline(buf, token, ','); )
		v.push_back(token);
	copy(v.begin(), v.end(), std::ostream_iterator<std::string>(std::cout, "."));
	for (int i = 0; i < v.size(); i++) {
		params.push_back(stof(v[i]));
	}
	//std::cout << '\n';
	//for (int i = 0; i < params.size(); i++) {
	//	std::cout << params[i] << '\n';
	//}
	//std::cout << '\n';
	
}

double calcMedian(std::vector<double> scores)
{
	double median = 0;
	size_t size = scores.size();

	sort(scores.begin(), scores.end());

	if (size % 2 == 0)
	{
		median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
	}
	else
	{
		median = scores[size / 2];
	}

	return median;
}

double calcMean(std::vector<double> scores)
{
	double mean = 0;

	sort(scores.begin(), scores.end());
	int sum = accumulate(scores.begin(), scores.end(), 0);
	mean = (double)sum / scores.size();

	return mean;
}

double calcVar(std::vector<double> scores)
{
	double mean = 0, var = 0;

	sort(scores.begin(), scores.end());
	int sum = accumulate(scores.begin(), scores.end(), 0);
	mean = (double)sum / scores.size();

	for (int i = 0; i < scores.size(); i++)
	{
		var += (scores[i] - mean)*(scores[i] - mean);
	}
	var = (double)(var) / (scores.size() - 1);

	return var;
}

double calcSkewness(std::vector<double> scores)
{
	double mean = 0, var = 0, S = 0, skewness = 0;

	sort(scores.begin(), scores.end());
	int sum = accumulate(scores.begin(), scores.end(), 0);
	mean = (double)sum / scores.size();

	for (int i = 0; i < scores.size(); i++)
	{
		var += (scores[i] - mean)*(scores[i] - mean);
	}
	var = (double)(var) / (scores.size() - 1);

	S = (double)sqrt(var);
	for (int i = 0; i < scores.size(); i++)
		skewness += (scores[i] - mean)*(scores[i] - mean)*(scores[i] - mean);
	skewness = skewness / (scores.size() * S * S * S);

	return skewness;
}

double calcKurtosis(std::vector<double> scores)
{
	double mean, var = 0, S = 0, skewness = 0, k = 0;

	sort(scores.begin(), scores.end());
	int sum = accumulate(scores.begin(), scores.end(), 0);
	mean = (double)sum / scores.size();

	for (int i = 0; i < scores.size(); i++)
	{
		var += (scores[i] - mean)*(scores[i] - mean);
	}
	var = (double)(var) / (scores.size() - 1);

	S = (double)sqrt(var);
	for (int i = 0; i < scores.size(); i++)
		skewness += (scores[i] - mean)*(scores[i] - mean)*(scores[i] - mean);
	skewness = skewness / (scores.size() * S * S * S);

	for (int i = 0; i < scores.size(); i++)
		k += (scores[i] - mean)*(scores[i] - mean)*(scores[i] - mean)*(scores[i] - mean);
	k = k / (scores.size()*S*S*S*S);
	k -= 3;

	return k;
}

void sublistFiles(QDir directory, QString indent, std::vector<std::string> &allMHAs)
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
		if (finfo.isDir()) {
			sublistFiles(QDir(finfo.absoluteFilePath()), indent, allMHAs);
		}
	}
}