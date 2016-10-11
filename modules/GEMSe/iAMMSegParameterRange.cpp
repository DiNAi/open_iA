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
#include "iAMMSegParameterRange.h"

#include "iAConsole.h"
#include "iASpectraDistanceImpl.h"

#include <QTextStream>

// TODO: use serialization library (e.g. boost:serialize?)

template <typename T>
struct Converter
{
	static T convert(QString str, bool * ok)
	{
		DEBUG_LOG("Unspecialized Converter called! This should not happen!\n");
		assert(false);
		return std::numeric_limits<T>::signaling_NaN();
	}
};

template <>
struct Converter<int>
{
	static int convert(QString str, bool * ok)
	{
		return str.toInt(ok);
	}
};

template <>
struct Converter<double>
{
	static double convert(QString str, bool * ok)
	{
		return str.toDouble(ok);
	}
};

template <typename T>
bool readLine(QTextStream & in, T & from, T & to, bool * log =0)
{
	if (in.atEnd()) return false;
	QString line = in.readLine();
	QStringList firstSplit = line.split(Output::NameSeparator);
	if (firstSplit.size() != 2) return false;
	QStringList important = firstSplit[1].split(Output::OptionalParamSeparator);
	if (important.size() == 0) return false;
	QStringList values = important[0].split(Output::ValueSeparator);
	if (values.size() != 2) return false;
	bool result = true;
	from = Converter<T>::convert(values[0], &result); if (!result) return false;
	to   = Converter<T>::convert(values[1], &result); if (!result) return false;
	if (log)
	{
		*log = (important.size() > 1 && important[1] == "log");
	}
	return true;
}

bool readDistanceFuncs(QTextStream & in, QVector<QSharedPointer<iASpectraDistance> > & result)
{
	if (in.atEnd()) return false;
	QString line = in.readLine();
	QStringList nameValues = line.split(Output::NameSeparator);
	if (nameValues.size() != 2) return false;
	QStringList distFuncNames = nameValues[1].split(Output::ValueSeparator);
	for (int i=0; i<distFuncNames.size(); ++i)
	{
		QSharedPointer<iASpectraDistance> distFunc = GetDistanceMeasureFromShortName(distFuncNames[i]);
		if (!distFunc)
		{
			return false;
		}
		result.push_back(distFunc);
	}
	return true;
}

iAMMSegParameterRange::iAMMSegParameterRange()
{
	for (int i = 0; i < MeasureCount; ++i)
	{
		measureMinMax.push_back(std::make_pair(0, 0));
	}
}

QSharedPointer<iAMMSegParameterRange> iAMMSegParameterRange::Create(QTextStream& in)
{
	QSharedPointer<iAMMSegParameterRange> result(new iAMMSegParameterRange);
	
	if (!readLine<double>(in, result->erw_beta_From, result->erw_beta_To, &result->erw_beta_logScale) ||
		!readLine<double>(in, result->erw_gamma_From, result->erw_gamma_To, &result->erw_gamma_logScale) ||
		!readLine<int>(in, result->erw_maxIter_From, result->erw_maxIter_To) ||
		!readLine<double>(in, result->svm_C_From, result->svm_C_To, &result->svm_C_logScale) ||
		!readLine<double>(in, result->svm_gamma_From, result->svm_gamma_To, &result->svm_gamma_logScale) ||
		!readLine<double>(in, result->svm_SeedProb_From, result->svm_SeedProb_To) ||
		!readLine<int>(in, result->svm_channels_From, result->svm_channels_To) ||
		!readLine<int>(in, result->objCountMin, result->objCountMax) ||
		!readLine<double>(in, result->durationMin, result->durationMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax) ||
		!readLine<double>(in, result->diceMin, result->diceMax))
	{
		DEBUG_LOG("Invalid Parameter Range text format!\n");
		return QSharedPointer<iAMMSegParameterRange>();
	}
	while (!in.atEnd())
	{
		iAMMSegModalityParamRange modParamRange;
		if (!readLine<double>(in, modParamRange.weightFrom, modParamRange.weightTo, &result->weightLogScale) ||
			!readLine<int>(in, modParamRange.pcaDimMin, modParamRange.pcaDimMax) ||
			!readDistanceFuncs(in, modParamRange.distanceFuncs))
		{
			DEBUG_LOG("Invalid Parameter Range text format in modality range!\n");
			return QSharedPointer<iAMMSegParameterRange>();
		}
		result->modalityParamRange.push_back(modParamRange);
	}
	return result;
}

QString iAMMSegParameterRange::GetDistanceFuncsStr(QVector<QSharedPointer<iASpectraDistance> > const & funcList)
{
	QString result;
	QTextStream out(&result);
	for (int i=0; i<funcList.size(); ++i)
	{
		out << funcList[i]->GetShortName();
		if (i < funcList.size()-1)
		{
			out << Output::ValueSeparator;
		}
	}
	return result;
}

bool iAMMSegParameterRange::Store(QTextStream& out)
{
	out << "ERW_Beta" << Output::NameSeparator << erw_beta_From << Output::ValueSeparator << erw_beta_To << Output::OptionalParamSeparator << (erw_beta_logScale?"log": "lin") << endl;
	out << "ERW_Gamma" << Output::NameSeparator << erw_gamma_From << Output::ValueSeparator << erw_gamma_To << Output::OptionalParamSeparator << (erw_gamma_logScale ? "log" : "lin") << endl;
	out << "ERW_MaxIter" << Output::NameSeparator << erw_maxIter_From << Output::ValueSeparator << erw_maxIter_To << Output::OptionalParamSeparator << (erw_maxIter_logScale ? "log" : "lin") << endl;
	out << "SVM_C" << Output::NameSeparator << svm_C_From << Output::ValueSeparator << svm_C_To << Output::OptionalParamSeparator << (svm_C_logScale ? "log" : "lin") << endl;
	out << "SVM_Gamma" << Output::NameSeparator << svm_gamma_From << Output::ValueSeparator << svm_gamma_To << Output::OptionalParamSeparator << (svm_gamma_logScale ? "log" : "lin") << endl;
	out << "SVM_Channels" << Output::NameSeparator << svm_channels_From << Output::ValueSeparator << svm_channels_To << endl;
	out << "SVM_SeedSetProb" << Output::NameSeparator << svm_SeedProb_From << Output::ValueSeparator << svm_SeedProb_To << endl;
	out << "GADIter" << Output::NameSeparator << gadIter_From << Output::ValueSeparator << gadIter_To << endl;
	out << "GADStep" << Output::NameSeparator << gadStep_From << Output::ValueSeparator << gadStep_To << endl;
	out << "GADCond" << Output::NameSeparator << gadCond_From << Output::ValueSeparator << gadCond_To << endl;
	out << "ObjectCount" << Output::NameSeparator << objCountMin << Output::ValueSeparator << objCountMax << endl;
	out << "Duration" << Output::NameSeparator << durationMin << Output::ValueSeparator << durationMax << endl;

	out << "MedianAccuracy" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "MeanAccuracy" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "DiceOverall" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Dice0" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Dice1" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Dice2" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Dice3" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Dice4" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "MedianUncertainty" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "MeanUncertainty" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Uncertainty0" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Uncertainty1" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Uncertainty2" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Uncertainty3" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "Uncertainty4" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "MedianConfusionUncertaintyFalse" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "MedianConfusionUncertaintyTrue" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;
	out << "MedianConfusionBadGood" << Output::NameSeparator << diceMin << Output::ValueSeparator << diceMax << endl;



	for (int i=0; i<modalityParamRange.size(); ++i)
	{
		out << "Weight" << Output::NameSeparator << modalityParamRange[i].weightFrom << Output::ValueSeparator << modalityParamRange[i].weightTo << Output::OptionalParamSeparator << (weightLogScale?"log": "lin") << endl;
		out << "PCADim" << Output::NameSeparator << modalityParamRange[i].pcaDimMin << Output::ValueSeparator << modalityParamRange[i].pcaDimMax << endl;
		out << "Dist" << Output::NameSeparator << GetDistanceFuncsStr(modalityParamRange[i].distanceFuncs) << endl;
	}
	return true;
}

bool iAMMSegParameterRange::CoversWholeRange(AttributeID id, double min, double max) const
{
	if (id < NonModalityParamCount)
	{
		switch (id)
		{
			case erwBeta:         return min <= erw_beta_From && erw_beta_To <= max;
			case erwGamma:	      return min <= erw_gamma_From && erw_gamma_To <= max;
			case erwMaxIter:      return min <= erw_maxIter_From && erw_maxIter_From <= max;
			case svmC:            return min <= svm_C_From && svm_C_To <= max;
			case svmGamma:        return min <= svm_gamma_From && svm_gamma_To <= max;
			case svmSeedProb:     return min <= svm_SeedProb_From && svm_SeedProb_To <= max;
			case gadIter:     return min <= gadIter_From && gadIter_To <= max;
			case gadStep:     return min <= gadStep_From && gadStep_To <= max;
			case gadCond:     return min <= gadCond_From && gadCond_To <= max;
			case svmChannelCount: return min <= svm_channels_From && svm_channels_To <= max;
		}
	}
	else if (id < GetInputParameterCount())
	{
		int modIdx  = (id-NonModalityParamCount) / ModalityParamCount;
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		switch (paramId)
		{
			case weight:          return min <=      modalityParamRange[modIdx].weightFrom && modalityParamRange[modIdx].weightTo <= max;
			case distance:        return min <= 0 && modalityParamRange[modIdx].distanceFuncs.size() <= max;
			case pca:             return min <=      modalityParamRange[modIdx].pcaDimMin && modalityParamRange[modIdx].pcaDimMax <= max;
		}
	}
	else
	{
		int paramId = id - GetInputParameterCount();
		switch (paramId)
		{
			case objectCount:     return min <= objCountMin && objCountMax <= max;
			case duration:        return min <= durationMin && durationMax <= max;
			case medianAcc:       return min <= diceMin && diceMax <= max;
			case meanAcc:         return min <= diceMin && diceMax <= max;
			case diceOverall:     return min <= diceMin && diceMax <= max;
			case dice0:           return min <= diceMin && diceMax <= max;
			case dice1:           return min <= diceMin && diceMax <= max;
			case dice2:           return min <= diceMin && diceMax <= max;
			case dice3:           return min <= diceMin && diceMax <= max;
			case dice4:           return min <= diceMin && diceMax <= max;
			case medianUncertainty:     return min <= diceMin && diceMax <= max;
			case meanUncertainty:       return min <= diceMin && diceMax <= max;
			case uncertainty0:          return min <= diceMin && diceMax <= max;
			case uncertainty1:          return min <= diceMin && diceMax <= max;
			case uncertainty2:          return min <= diceMin && diceMax <= max;
			case uncertainty3:          return min <= diceMin && diceMax <= max;
			case uncertainty4:          return min <= diceMin && diceMax <= max;
			case median_confusion_uncertainty_false:        return min <= diceMin && diceMax <= max;
			case median_confusion_uncertainty_true:         return min <= diceMin && diceMax <= max;
			case m_median_confusion_bad_to_good:            return min <= diceMin && diceMax <= max;
			case diceMetric:
			case kappa:
			case overallAcc:
			case precision:
			case recall:
				return min <= measureMinMax[paramId-diceMetric].first && measureMinMax[paramId-diceMetric].second <= max;
		}
	}
	return false;
}

bool iAMMSegParameterRange::IsLogScale(AttributeID id) const
{
	if (id < NonModalityParamCount)
	{
		switch (id)
		{
			case erwBeta:         return erw_beta_logScale;
			case erwGamma:	      return erw_gamma_logScale;
			case erwMaxIter:      return erw_maxIter_logScale;
			case svmC:            return svm_C_logScale;
			case svmGamma:        return svm_gamma_logScale;
			case svmSeedProb:	  return false;
			case gadIter:	  return false;
			case gadStep:	  return false;
			case gadCond:	  return false;
			case svmChannelCount: return false;
		}
	}
	else if (id < GetInputParameterCount())
	{
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		switch (paramId)
		{
			case weight:          return weightLogScale;
			default:              return false;
		}
	}
	return false;
}

double iAMMSegParameterRange::GetMin(AttributeID id) const
{
	if (id < NonModalityParamCount)
	{
		switch (id)
		{
			case erwBeta:         return erw_beta_From;
			case erwGamma:        return erw_gamma_From;
			case erwMaxIter:      return erw_maxIter_From;
			case svmC:            return svm_C_From;
			case svmGamma:        return svm_gamma_From;
			case svmSeedProb:	  return svm_SeedProb_From;
			case gadIter:	  return gadIter_From;
			case gadStep:	  return gadStep_From;
			case gadCond:	  return gadCond_From;
			case svmChannelCount: return svm_channels_From;
		}
	}
	else if (id < GetInputParameterCount())
	{
		int modIdx  = (id-NonModalityParamCount) / ModalityParamCount;
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		switch (paramId)
		{
			case weight:          return modalityParamRange[modIdx].weightFrom;
			case pca:             return modalityParamRange[modIdx].pcaDimMin;
			case distance:        return 0;
		}
	}
	else
	{
		int paramId = id - GetInputParameterCount();
		switch(paramId)
		{
			case objectCount:     return objCountMin;
			case duration:        return durationMin;
			case medianAcc:             return diceMin;
			case meanAcc:               return diceMin;
			case diceOverall:           return diceMin;
			case dice0:                 return diceMin;
			case dice1:                 return diceMin;
			case dice2:                 return diceMin;
			case dice3:                 return diceMin;
			case dice4:                 return diceMin;
			case medianUncertainty:     return diceMin;
			case meanUncertainty:       return diceMin;
			case uncertainty0:          return diceMin;
			case uncertainty1:          return diceMin;
			case uncertainty2:          return diceMin;
			case uncertainty3:          return diceMin;
			case uncertainty4:          return diceMin;
			case median_confusion_uncertainty_false:        return diceMin;
			case median_confusion_uncertainty_true:         return diceMin;
			case m_median_confusion_bad_to_good:            return diceMin;
			case diceMetric:
			case kappa:
			case overallAcc:
			case precision:
			case recall:	  return measureMinMax[paramId-diceMetric].first;
		}
	}
	return 0;
}

double iAMMSegParameterRange::GetMax(AttributeID id) const
{
	if (id < NonModalityParamCount)
	{
		switch (id)
		{
			case erwBeta:         return erw_beta_To;
			case erwGamma:        return erw_gamma_To;
			case erwMaxIter:      return erw_maxIter_To;
			case svmC:            return svm_C_To;
			case svmGamma:        return svm_gamma_To;
			case svmSeedProb:     return svm_SeedProb_To;
			case gadIter:	  return gadIter_To;
			case gadStep:	  return gadStep_To;
			case gadCond:	  return gadCond_To;
			case svmChannelCount: return svm_channels_To;
		}
	}
	else if (id < GetInputParameterCount())
	{
		int modIdx  = (id-NonModalityParamCount) / ModalityParamCount;
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		switch (paramId)
		{
			case weight:          return modalityParamRange[modIdx].weightTo;
			case pca:             return modalityParamRange[modIdx].pcaDimMax;
			case distance:        return modalityParamRange[modIdx].distanceFuncs.size()-1;
		}
	}
	else
	{
		int paramId = id - GetInputParameterCount();
		switch (paramId)
		{
			case objectCount:     return objCountMax;
			case duration:        return durationMax;
			case medianAcc:             return diceMax;
			case meanAcc:               return diceMax;
			case diceOverall:           return diceMax;
			case dice0:                 return diceMax;
			case dice1:                 return diceMax;
			case dice2:                 return diceMax;
			case dice3:                 return diceMax;
			case dice4:                 return diceMax;
			case medianUncertainty:     return diceMax;
			case meanUncertainty:       return diceMax;
			case uncertainty0:          return diceMax;
			case uncertainty1:          return diceMax;
			case uncertainty2:          return diceMax;
			case uncertainty3:          return diceMax;
			case uncertainty4:          return diceMax;
			case median_confusion_uncertainty_false:        return diceMax;
			case median_confusion_uncertainty_true:         return diceMax;
			case m_median_confusion_bad_to_good:            return diceMax;
			case diceMetric:
			case kappa:
			case overallAcc:
			case precision:
			case recall:	  return measureMinMax[paramId-diceMetric].second;
		}
	}
	return 1;
}


iAValueType iAMMSegParameterRange::GetRangeType(AttributeID id) const
{
	if (id < NonModalityParamCount)
	{
		switch (id)
		{
			case erwBeta:
			case erwGamma:        return Continuous;
			case erwMaxIter:      return Discrete;
			case svmSeedProb:     return Continuous;
			case gadIter:	        return Discrete;
			case gadStep:	  
			case gadCond:	  
			case svmC:
			case svmGamma:        return Continuous;
			case svmChannelCount: return Discrete;
		}
	}
	else if (id < GetInputParameterCount())
	{
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		switch (paramId)
		{
			case weight:          return Continuous;
			case pca:             return Discrete;
			case distance:        return Categorical;
		}
	}
	else
	{
		int paramId = id - GetInputParameterCount();
		switch (paramId)
		{
			case objectCount:     return Discrete;
			case duration:        return Continuous;
			case medianAcc:             return Continuous;
			case meanAcc:               return Continuous;
			case diceOverall:           return Continuous;
			case dice0:                 return Continuous;
			case dice1:                 return Continuous;
			case dice2:                 return Continuous;
			case dice3:                 return Continuous;
			case dice4:                 return Continuous;
			case medianUncertainty:     return Continuous;
			case meanUncertainty:       return Continuous;
			case uncertainty0:          return Continuous;
			case uncertainty1:          return Continuous;
			case uncertainty2:          return Continuous;
			case uncertainty3:          return Continuous;
			case uncertainty4:          return Continuous;
			case median_confusion_uncertainty_false:        return Continuous;
			case median_confusion_uncertainty_true:         return Continuous;
			case m_median_confusion_bad_to_good:            return Continuous;
			case diceMetric:
			case kappa:
			case overallAcc:
			case precision:
			case recall:	  return Continuous;
		}
	}
	return Continuous;
}

QString iAMMSegParameterRange::GetName(AttributeID id) const
{
	if (id < NonModalityParamCount)
	{
		switch (id)
		{
			case erwBeta:         return "ERW Beta";
			case erwGamma:        return "ERW Gamma";
			case erwMaxIter:      return "ERW Max. It.";
			case svmC:            return "SVM C";
			case svmGamma:        return "SVM Gamma";
			case svmSeedProb:	  return "SVM Seed Probability";
			case gadIter:	  return "G Iter";
			case gadStep:	  return "G Step";
			case gadCond:	  return "G Cond";
			case svmChannelCount: return "SVM Channel #";
		}
	}
	else if (id < GetInputParameterCount())
	{
		int modIdx  = (id-NonModalityParamCount) / ModalityParamCount;
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		switch (paramId)
		{	/* TODO: link to modalities to get name! modalityParamRange[modIdx].GetName */
			case weight:          return QString("Modality %1 Weight").arg(modIdx+1);
			case pca:             return QString("Modality %1 PCA Dim.").arg(modIdx+1);
			case distance:        return QString("Modality %1 Distance Metric").arg(modIdx+1);
		}
	}
	else
	{
		int paramId = id - GetInputParameterCount();
		switch (paramId)
		{
			case objectCount:           return "Object Count";
		 	case duration:              return "Performance";
		    case medianAcc:             return "MED A";
		    case meanAcc:               return "M A";
		    case diceOverall:           return "D";
		    case dice0:                 return "D 0";
		    case dice1:                 return "D 1";
		    case dice2:                 return "D 2";
		    case dice3:                 return "D 3";
		    case dice4:                 return "D 4";
		    case medianUncertainty:     return "MED U";
		    case meanUncertainty:       return "m U";
		    case uncertainty0:          return "U 0";
		    case uncertainty1:          return "U 1";
		    case uncertainty2:          return "U 2";
		    case uncertainty3:          return "U 3";
		    case uncertainty4:          return "U 4";
		    case median_confusion_uncertainty_false:        return "MED C U F";
		    case median_confusion_uncertainty_true:         return "MED C U T";
		    case m_median_confusion_bad_to_good:            return "MED C B G";
			case diceMetric:            return "Mean Overlap (Dice)";
			case kappa:			            return "Kappa Coefficient";
			case overallAcc:	          return "Overall Accuracy";
			case precision:		          return "Precision";
			case recall:		            return "Recall";
		}
	}
	return "Unknown";
}

#include "iANameMapper.h"

class iADistFuncNameMapper: public iANameMapper
{
public:
	iADistFuncNameMapper(QVector<QSharedPointer<iASpectraDistance> > const & distFuncs):
		m_distFuncs(distFuncs)
	{}
	virtual QString GetName(int idx) const
	{
		return m_distFuncs[idx]->GetName();
	}
	virtual int GetIdx(QString const & name, bool & ok) const
	{
		ok = false;
		return -1;
	}
	
	virtual int size() const
	{
		return m_distFuncs.size();
	}
private:
	QVector<QSharedPointer<iASpectraDistance> > const & m_distFuncs;
};


QSharedPointer<iANameMapper> iAMMSegParameterRange::GetNameMapper(AttributeID id) const
{
	if (id >= NonModalityParamCount && id < GetInputParameterCount())
	{
		int modIdx  = (id-NonModalityParamCount) / ModalityParamCount;
		int paramId = (id-NonModalityParamCount) % ModalityParamCount;
		if (paramId == distance)
		{
			return QSharedPointer<iANameMapper>(new iADistFuncNameMapper(modalityParamRange[modIdx].distanceFuncs));
		}
	}
	return QSharedPointer<iANameMapper>();
}


void iAMMSegParameterRange::AdaptDurationMinMax(double curDur)
{
	if (curDur < durationMin)
	{
		assert(curDur > 0);
		durationMin = curDur;
	}
	if (curDur > durationMax)
	{
		durationMax = curDur;
	}
}


void iAMMSegParameterRange::AdaptObjCountMinMax(int curCnt)
{
	if (curCnt < objCountMin)
	{
		assert(curCnt > 0);
		objCountMin = curCnt;
	}
	if (curCnt > objCountMax)
	{
		objCountMax = curCnt;
	}
}

int iAMMSegParameterRange::GetInputParameterCount() const
{
	return NonModalityParamCount+modalityParamRange.size()*ModalityParamCount;
}

int iAMMSegParameterRange::GetDerivedOutputCount() const
{
	return DerivedOutputCount;
}

int iAMMSegParameterRange::GetAttributeCount() const
{
	return NonModalityParamCount // beta + gamma
		+ modalityParamRange.size()*ModalityParamCount // per modality: weight, PCA, distance
		+ DerivedOutputCount // for object count
	;
}

int iAMMSegParameterRange::GetFirstDerivedOutputIndex() const
{
	return GetInputParameterCount();
}