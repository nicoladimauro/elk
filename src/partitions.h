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

/*! \file partitions.h
    \brief 
*/

#ifndef _PARTITIONS_H
#define _PARTITIONS_H

void restricted_growth_first (const int *, const int, int *);
void reverse_int (int *, const int, const int);
int next_permutation (int *, const int);
int compare_ints (const void *, const void *);
int restricted_growth_next_permute (int *, const int, const int);
int restricted_growth_next (int *, const int, int *, const int);

typedef struct partition_node
{
  int m;
  int numparts;
  int *parts;
} partition;

int count_partitions (int);
void set_partitions (partition *, int);
void print_partitions (partition *, int);
void free_partitions (partition *, int);
void GenPartitions (int m);


#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)<(y)?(y):(x))


void firstKCoalitionStructure (int *, int *, int, int);
int nextKCoalitionStructure (int *, int *, int, int);
void pks_first (int *, int *, int, int, int);
int pks_next (int *, int *, int, int, int);

void output (partition);

void ip_initialize_first(int *, int *, int , int *, int *, int );
int ip_next_combination(int *, int *, int , int );
int ip_next_partition(int *, int *, int , int *, int *, int , int *, int *);
void ip_prune_partition(int *, int *, int, int *, int, int);

#endif
