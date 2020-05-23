/*

 $Id: success.cc,v 1.3 2007/04/27 06:01:51 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>
#include "success.h"
#include "timesyshms.h"


extern char *programname;
extern FILE *logFile;

void success( char hostnm[MAX_CHARS],
		Clock jobStart,
		struct tms tms_jobStart )

{
    char message[LINE_LEN];
    Clock jobEnd;
    struct tms tms_jobEnd;

    pr_2x( logFile, stderr, "\n" );
    pr_2x( logFile, stderr, UnderLine );
    prStr( message, "%s: Successful Completion on \"%s\"\n\n", programname, hostnm );
    pr_2x( logFile, stderr, message );

    jobEnd = times( &tms_jobEnd );

    timesyshms( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd );

    pr_2x( logFile, stderr, UnderLine );
}
/* EOF */
