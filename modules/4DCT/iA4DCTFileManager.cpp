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

#include "iA4DCTFileManager.h"
// vtk
#include <vtkMetaImageReader.h>
#include <vtkImageData.h>

iA4DCTFileManager& iA4DCTFileManager::getInstance( )
{
	static iA4DCTFileManager instance;
	return instance;
}

vtkImageData* iA4DCTFileManager::getImage( iA4DCTFileData file )
{
	return findOrCreateImage( file )->GetOutput( );
}

vtkAlgorithmOutput* iA4DCTFileManager::getOutputPort( iA4DCTFileData file )
{
	return findOrCreateImage( file )->GetOutputPort( );
}

iA4DCTFileManager::ReaderType iA4DCTFileManager::findOrCreateImage( iA4DCTFileData file )
{
	string key = file.Name.toStdString( );
	if( m_map.find( key ) == m_map.end( ) )
	{	// the key does not exist
		ReaderType reader = ReaderType::New( );
		reader->SetFileName( file.Path.toStdString( ).c_str( ) );
		reader->Update( );
		m_map[ key ] = reader;
	}
	return m_map[ key ];
}
