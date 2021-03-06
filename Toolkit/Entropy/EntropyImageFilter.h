/*************************************  open_iA  ************************************ *
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
#pragma once

#include <itkImageToImageFilter.h>

 /** \class EntropyImageFilter
  *
  * \brief Calculates the pixelwise entropy out of a collection of input images
  *
  * \par INPUTS
  * All inputs to this filter must be floating point images representing a
  * distribution (i.e. their sum should be 1)
  * Input volumes must all contain the same size RequestedRegions.
  *
  * \par OUTPUTS
  * The filter produces a single output volume. Each output pixel contains the
  * entropy of the given pixel over all input volumes.
  * Higher values indicate higher "uncertainty", i.e. an equal distribution of
  * probabilities; lower values indicate that less evenly distributed
  * probabilities; a value of 0 indicates that the probability of one input is
  * 1, while the others are all 0.
  *
  * \par PARAMETERS
  * Normalize whether the output entropy should be normalized to the interval
  * [0..1]
  *
  * \author Bernhard Froehler, FH Wels
  *
  * \ ingroup ???
  */

namespace fhw
{

template< typename TInputImage, typename TOutputImage = TInputImage >
class EntropyImageFilter :
	public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef EntropyImageFilter                          Self;
	typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
	typedef itk::SmartPointer< Self >                            Pointer;
	typedef itk::SmartPointer< const Self >                      ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self);

	/** Run-time type information (and related methods) */
	itkTypeMacro(EntropyImageFilter, ImageToImageFilter);

	/** Extract some information from the image types.  Dimensionality
	 * of the two images is assumed to be the same. */
	typedef typename TOutputImage::PixelType OutputPixelType;
	typedef typename TInputImage::PixelType  InputPixelType;

	/** Extract some information from the image types.  Dimensionality
	 * of the two images is assumed to be the same. */
	itkStaticConstMacro(InputImageDimension, int, TInputImage::ImageDimension);
	itkStaticConstMacro(ImageDimension, int, TOutputImage::ImageDimension);

	/** Image typedef support */
	typedef TInputImage                           InputImageType;
	typedef TOutputImage                          OutputImageType;
	typedef typename InputImageType::ConstPointer InputImagePointer;
	typedef typename OutputImageType::Pointer     OutputImagePointer;

	typedef unsigned long                         LabelCountType;

	/** Superclass typedefs. */
	typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

	void SetNormalize(bool on)
	{
		m_normalize = on;
	}
	bool GetNormalize() const
	{
		return m_normalize;
	}

#ifdef ITK_USE_CONCEPT_CHECKING
	// Begin concept checking
	itkConceptMacro(SameDimensionCheck,
		(itk::Concept::SameDimension< InputImageDimension, ImageDimension >));
	// End concept checking
#endif

protected:
	EntropyImageFilter();
	virtual ~EntropyImageFilter() {}

	/** Determine maximum label value in all input images and initialize
	 * global data. */
	void BeforeThreadedGenerateData() ITK_OVERRIDE;

	void ThreadedGenerateData
	(const OutputImageRegionType & outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

	void PrintSelf(std::ostream &, itk::Indent) const ITK_OVERRIDE;

private:
	EntropyImageFilter(const Self &) ITK_DELETE_FUNCTION;
	void operator=(const Self &) ITK_DELETE_FUNCTION;

	bool m_normalize;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "EntropyImageFilter.hxx"
#endif
