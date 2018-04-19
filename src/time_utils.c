/***************************************************************************
 *   Copyright (C) 2009 by Nicola Di Mauro                                 *
 *   ndm@di.uniba.it                                                       *
 *                                                                         *
 *   Department of Computer Science                                        *
 *   University of Bari, Italy                                             *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "time_utils.h"

long long
timeval_diff (struct rusage *difference, struct rusage *end_time,
	      struct rusage *start_time)
{

  struct rusage temp_diff;

  if (difference == NULL)
    {
      difference = &temp_diff;
    }

  difference->ru_utime.tv_sec =
    end_time->ru_utime.tv_sec - start_time->ru_utime.tv_sec;
  difference->ru_utime.tv_usec =
    end_time->ru_utime.tv_usec - start_time->ru_utime.tv_usec;

  /* Using while instead of if below makes the code slightly more robust. */

  while (difference->ru_utime.tv_usec < 0)
    {
      difference->ru_utime.tv_usec += 1000000;
      difference->ru_utime.tv_sec -= 1;
    }

  return 1000000LL * difference->ru_utime.tv_sec +
    difference->ru_utime.tv_usec;

}				/* timeval_diff() */
