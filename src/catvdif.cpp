/*****************************************************************************
*    <DataSim: VLBI data simulator>                                          * 
*    Copyright (C) <2015> <Zheng Meyer-Zhao>                                 *
*                                                                            *
*    This file is part of DataSim.                                           *
                                                                             *
*    DataSim is free software: you can redistribute it and/or modify         *
*    it under the terms of the GNU General Public License as published by    *
*    the Free Software Foundation, either version 3 of the License, or       *
*    (at your option) any later version.                                     *
*                                                                            *
*    DataSim is distributed in the hope that it will be useful,              *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
*    GNU General Public License for more details.                            *
*                                                                            *
*    You should have received a copy of the GNU General Public License       *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdint.h>
#include "architecture.h"
#include "configuration.h"
#include "vdifio.h"
#include "catvdif.h"
#include <mpi.h>

using namespace std;


// Read in vdif header produced after vdifzipper
// count the duration in seconds of each vdif
// set the seconds in the vdif header to the process number

// ????? Different antenna's????

void catvdif(size_t div, int myid, vector<string> vdiffiles)
{
	// One output vdif
	// Every process writes to the output vdif with offset
	ofstream outputvdif;

	// Only needs 'div' number of processes to work on it
	if(myid < div)
	{
		// Open all time-divided vdif files
		// Compute time duration of the vdif file
		// modify start time of vdif to myid*timeduration
		ifstream inputvdif;
		inputvdif.()

	}
}