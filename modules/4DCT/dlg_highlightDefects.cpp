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
*          Stelzhamerstra�e 23, 4600 Wels / Austria, Email: c.heinzl@fh-wels.at       *
* ************************************************************************************/
 
#include "pch.h"
#include "dlg_highlightDefects.h"
// iA
#include "iA4DCTVisWin.h"

const QColor CrackColor(91, 155, 213);		// blue
const QColor DebondingColor(0, 176, 80);	// green
const QColor PulloutColor(255, 192, 0);		// yellow
const QColor BreakageColor(255, 0, 0);		// red

dlg_highlightDefects::dlg_highlightDefects( QWidget * parent )
	: QDialog( parent )
{
	setupUi( this );
	connect( pbPullouts,	SIGNAL( clicked( ) ), this, SLOT( pulloutsButtonClicked( ) ) );
	connect( pbBreakages,	SIGNAL( clicked( ) ), this, SLOT( breakagesButtonClicked( ) ) );
	connect( pbCracks,		SIGNAL( clicked( ) ), this, SLOT( cracksButtonClicked( ) ) );
	connect( pbDebondings,	SIGNAL( clicked( ) ), this, SLOT( debondingsButtonClicked( ) ) );
	connect( pbLabeledImg,	SIGNAL( clicked( ) ), this, SLOT( labeledImgButtonClicked( ) ) );

	cbPullouts->setColor( PulloutColor );
	cbBreakages->setColor( BreakageColor );
	cbCracks->setColor( CrackColor );
	cbDebondings->setColor( DebondingColor );
}

dlg_highlightDefects::~dlg_highlightDefects( )
{ /* not implemented */ }

void dlg_highlightDefects::setVisWin( iA4DCTVisWin * visWin )
{
	m_visWin = visWin;
}

void dlg_highlightDefects::pulloutsButtonClicked( )
{
	QString file;
	m_visWin->showDialog( file );
	lePullouts->setText( file );
}

void dlg_highlightDefects::breakagesButtonClicked( )
{
	QString file;
	m_visWin->showDialog( file );
	leBreakages->setText( file );
}

void dlg_highlightDefects::debondingsButtonClicked( )
{
	QString file;
	m_visWin->showDialog( file );
	leDebondings->setText( file );
}

void dlg_highlightDefects::cracksButtonClicked( )
{
	QString file;
	m_visWin->showDialog( file );
	leCracks->setText( file );
}

void dlg_highlightDefects::labeledImgButtonClicked( )
{
	QString file;
	m_visWin->showDialog( file );
	leLabeledImg->setText( file );
}
