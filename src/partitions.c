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

/**
   Exhaustive information about partitioning can be found in 
   Knuth's "pre-fascicle" of Volume 4x of TAOCP
*/

#include "partitions.h"
#include <stdio.h>
#include <stdlib.h>

#define swap_int( x, y ) {int __i = x; x = y; y = __i;}

/**
   Computes the first RG corresponding to the first partition of n elements
   into num_blocks blocks.

   Note that the block description must be ordered. In particular,
   integer partitions are ordered. Specify, [1,2,2] and not [2,2,1].
*/
void
restricted_growth_first (const int *block, const int num_blocks, int *rg)
{
  register int i, j, k;
  k = 0;
  for (i = 0; i < num_blocks; i++)
    for (j = 0; j < block[i]; j++)
      rg[k++] = i;
}

/**
   Reverse elements of an array of integers pointed by v from position start 
   to end 
*/
void
reverse_int (int *v, const int start, const int end)
{
  register int i, j;
  for (i = start, j = end; i < j; i++, j--)
    swap_int (v[i], v[j]);
}

/**
   Compute the next permutation of the elements in v. n is the length of v.
   It returns 1 if v contains the next permutation, 0 if there are not more 
   permutations of v.

   This is the C translation, with some minor modifications, of the Perl code
   contribution to www.perlmonks.org by blokhead reported in the article
   "Partitioning a set into parts of given sizes" availbale at 
   http://www.perlmonks.org/index.pl?node_id=533530
*/
int
next_permutation (int *v, const int n)
{
  register int i, j;

  if (n < 2)
    return (0);
  /* Find the rightmost position 'i' where the sequence increases */
  i = n - 2;
  while (i >= 0 && v[i] >= v[i + 1])
    i--;
  if (i < 0)
    return (0);
  /* Reverse everything to the right (now it's in increasing order) */
  reverse_int (v, i + 1, n - 1);
  /* Move to the right partition in order to find the first number larger 
     than 'v[i]' */
  j = i + 1;
  while (v[i] >= v[j])
    j++;
  swap_int (v[i], v[j]);
  return (1);
}

/**
   Function to compare integers useful for qsort
 */
int
compare_ints (const void *a, const void *b)
{
  return (*(int *) a - *(int *) b);
}


/**
   In order to obtain lexicographically next RG string, look for the rightmost
   position where we have an appropriate candidate available. The candidate 
   is the smallest number to the right of our current position such that:
     - candidate is larger than our current position
     - candidate is not >=2 larger than everything to the left
    (restricted growth property)

   This is the C translation, with some minor modifications, of the Perl code
   contribution to www.perlmonks.org by blokhead reported in the article
   "Partitioning a set into parts of given sizes" availbale at 
   http://www.perlmonks.org/index.pl?node_id=533530
*/
int
restricted_growth_next_permute (int *rg, const int n, const int num_blocks)
{
  int j, candidate, p, i, vc, found;
  int *avail = (int *) calloc (num_blocks, sizeof (int));
  for (j = 0; j < num_blocks; j++)
    avail[j] = -1;

  i = n;
  while (--i)  {
		candidate = -1;
		p = rg[i] + 1;
		while (p < num_blocks && candidate == -1)	{
			candidate = avail[p];
			p++;
		}
		if (candidate != -1)	{
			vc = rg[candidate] - 1;
			found = 0;
			for (j = 0; j < i && !found; j++)
				if (rg[j] >= vc)
					found = 1;
			if (found)
				break;
		}
		avail[rg[i]] = i;
	}
  if (i == 0)    {
		free (avail);
		return (0);
	}
  swap_int (rg[i], rg[candidate]);
  qsort (&rg[i + 1], n - (i + 1), sizeof (int), compare_ints);
  free (avail);
  return (1);

}

/**
   The function compute the next restricted growth of rg.
   
   "The basic way to iterate over partitions in general is to maintain
   a restricted-growth (RG) string.
   So it's settled that we must maintain some sort of RG string. 
   Indeed, it's not hard to see that we need to all the permutations
   of a certain RG string that are also RG strings themselves (because 
   permuting an RG string preserves the number of items sent to each 
   block of the corresponding partition). But there's a catch -- RG 
   strings are such that the resulting partitions are sorted by their 
   smallest element (in particular, the smallest element of the whole 
   set is always in partition 1). So we'll have to cycle through all 
   permutations of the block sizes as well, since the first block should 
   be any of the possible sizes.

   So here's roughly how we can get all partitions into blocks of the 
   given sizes:

    * for each permutation of block_size
          o initialize rg to the (lexicographically) first RG string, 
	   where the multiplicities of each element match block_size.  
	   For instance, (4,2,5) becomes the RG string 00001122222.
          o for each permutation of rg that is also an RG string, 
	   return the corresponding partition"
	                "Partitioning a set into parts of given sizes",
			http://www.perlmonks.org/index.pl?node_id=533530
*/
int
restricted_growth_next (int *block, const int num_blocks, int *rg, const int n){
  if (restricted_growth_next_permute (rg, n, num_blocks))
    return (1);
  else if (next_permutation (block, num_blocks)) {
		restricted_growth_first (block, num_blocks, rg);
		return (1);
	}
  return (0);
}


void
output (partition a)
/*
**  print out the partition a
*/
{
  int i, j;
  printf ("%d = ", a.m);
  for (i = 0; i < a.numparts; i = i + 1)
    {
      j = a.parts[i];
      printf ("%d", j);
      if (i != a.numparts - 1)
	printf ("+");
    }
}


static partition *partitions;

void
set_rec_partitions (partition * a, int m, int B, int N, partition * p)
{
  partitions = p;
  int i, j;
  if (m == 0)  {
		a->numparts = N;
		partitions->m = a->m;
		partitions->parts = (int *) calloc (a->numparts, sizeof (int));
		for (j = 0; j < a->numparts; j++)	{
			partitions->parts[j] = a->parts[a->numparts - j - 1];
		}
		partitions->numparts = a->numparts;
		partitions++;
	}
  else    {
		for (i = 1; i <= min (B, m); i = i + 1)			{
			a->parts[N] = i;
			set_rec_partitions (a, m - i, i, N + 1, partitions);
		}
	}
}


int
count_rec_partitions (partition * a, int m, int B, int N)
{
  int i, c = 0;
  if (m == 0)
    {
      return (1);
    }
  else
    {
      for (i = 1; i <= min (B, m); i = i + 1)
	{
	  a->parts[N] = i;
	  c += count_rec_partitions (a, m - i, i, N + 1);
	}
    }
  return (c);
}


int
count_partitions (int m)
{
  partition a;
  a.m = m;
  a.parts = (int *) calloc (m, sizeof (int));
  return (count_rec_partitions (&a, m, m, 0));
  free (a.parts);
}

void
set_partitions (partition * partitions, int m)
{
  if (m > 100)
    exit (0);
  partition a;
  a.m = m;
  a.parts = (int *) calloc (m, sizeof (int));
  set_rec_partitions (&a, m, m, 0, partitions);
  free (a.parts);
}


void
free_partitions (partition * p, int m)
{
  int i;
  partition *pp = p;
  for (i = 0; i < m; i++)
    {
      free (pp->parts);
      pp++;
    }
  free (p);
}

void
print_partitions (partition * p, int m)
{
  int i;
  partition a;
  for (i = 0; i < m; i++)
    {
      a = *p;
      output (a);
      printf ("\n");

      /***/
      int j, s = 1;
      for (j = 0; j < a.numparts; j = j + 1)
	s *= a.parts[j];

      int *ip1 = (int *) calloc (13, sizeof (int));
      int *ipm = (int *) calloc (13, sizeof (int));

      //pks_first(ip1, ipm, 8, a.numparts, s );
      firstKCoalitionStructure (ip1, ipm, 13, a.numparts);

      /*      int f;
         printf("\n");  
         for (f=0;f<5;f++)
         printf("%d ",ip1[f]);
         printf(" - ");
         for (f=0;f<5;f++)
         printf("%d ",ipm[f]); */

      while (nextKCoalitionStructure (ip1, ipm, 13, a.numparts))
	//      while (pks_next(ip1, ipm, 8, a.numparts, s ))
	{
	  /*      printf("\n");  
	     for (f=0;f<5;f++)
	     printf("%d ",ip1[f]);
	     printf(" - ");
	     for (f=0;f<5;f++)
	     printf("%d ",ipm[f]); */
	}

      /***/

      p++;
    }
}

void
firstKCoalitionStructure (int *CS, int *MCS, int nAgents, int K)
{
  int i;
  for (i = 0; i <= nAgents - K; i++)
    {
      CS[i] = 0;
      MCS[i] = 0;
    }
  for (i = nAgents - K + 1; i < nAgents; i++)
    {
      CS[i] = i - (nAgents - K);
      MCS[i] = i - (nAgents - K);
    }
  CS[nAgents] = K;
}

int
nextKCoalitionStructure (int *CS, int *MCS, int nAgents, int K)
{
  int i, j;
  for (i = nAgents - 1; i > 0; i--)
    {
      if (CS[i] < K - 1 && CS[i] <= MCS[i - 1])
	{
	  CS[i]++;
	  if (MCS[i] < CS[i])
	    MCS[i] = CS[i];
	  for (j = i + 1; j <= nAgents - (K - MCS[i]); j++)
	    {
	      CS[j] = 0;
	      MCS[j] = MCS[i];
	    }
	  for (j = nAgents - (K - MCS[i]) + 1; j < nAgents; j++)
	    {
	      CS[j] = K - (nAgents - j);
	      MCS[j] = K - (nAgents - j);
	    }
	  return 1;
	}
    }
  return 0;
}

int
p_get_id (int *CS, int nAgents, int K)
{
  int i, j, s1, s0;
  s1 = 1;
  for (j = 0; j < K; j++)
    {
      s0 = 0;
      for (i = 0; i < nAgents; i++)
	{
	  if (CS[i] == j)
	    s0++;
	}
      s1 *= s0;
    }
  return (s1);
}

/*
void pks_first(int * CS, int * MCS, int nAgents, int K, int s )
{
  firstKCoalitionStructure(CS, MCS, nAgents, K );
  while (p_get_id(CS,nAgents,K) != s)
    {
      nextKCoalitionStructure(CS, MCS, nAgents, K);
    }
}

int pks_next(int * CS, int * MCS, int nAgents, int K, int s )
{
  while (nextKCoalitionStructure(CS, MCS, nAgents, K))
    {
      if (p_get_id(CS,nAgents,K) == s)
	return(1);
    }
  return(0);
}

*/


/* ------------------------------------------------------------------------------- 
	 Purposely designed for IP
   ------------------------------------------------------------------------------- */

/*
	rg is the restricted growth representing the partition
	Min_k is the smallest element in M_k
	n is the number of elements
	parts_k is the size of the block k
	sum_k is the multiplicity of part_k
	np is the number of blocks
*/
void ip_initialize_first(int *rg, int *Min, int n, int *parts, int *sum, int np){
	int i,l;
	int j = 0;
	for (i=0; i<np; i++){
		Min[i] = 0;
		for (l = 0; l<parts[i]; l++){
			rg[j] = i;
			j++;
		}
	}

	for (i=0;i< np; i++){
		sum[i]=0;
		for (j = 0; j< np; j++)
			if (parts[j] == parts[i])
				sum[i] += parts[j];
	}

	int prec = 0;
	for (i=1;i<np;i++){
		if (parts[i] != parts[i-1]){
			prec = sum[i-1];
			sum[i] = sum[i] + sum[i-1];
		} else {
			sum[i] = sum[i] + prec;
		}
	}
}


int ip_next_combination(int *comb, int *elements, int h, int t){
	int i,j,k;

	i = h-1;
	while (i>=0){
		if (comb[i] < elements[t-h+i]){
			for (k=0;k<t && elements[k]!= comb[i];k++);
			k++;
			comb[i] = elements[k];
			k++;
			for (j=i+1;j<h;j++){
				comb[j] = elements[k];
				k++;
			}
			return(1);
		}
		i--;
	}
	return(0);
}

int ip_next_partition(int *rg, int *Min, int n, int *parts, int *sum, int np, int *comb, int *elem){
	int i,j,l,h,t,k;
	int found,tmp,ref,candidate;
	int alpha, smallest, check;

	i = np - 2;
	while (i>=0){
		if (i==0 || parts[i]!=parts[i-1])
			alpha = 0;
		else
			alpha = Min[i-1];

		/* controlla se può cambiare la combinazione della coalizione i */
		h =0; t=0;
		for (l=0;l<n;l++){
			if (rg[l]==i){
				comb[h]=l;
				h++;
			}
			if (rg[l]>=i){
				elem[t] = l;
				t++;
				//				rg[l]=-1;
			}
		}
		if (ip_next_combination(comb,elem,h,t)){
			for (smallest = 0; smallest < t && comb[0]!=elem[smallest]; smallest++);
			Min[i] = smallest;

			check = 0;
			if (i==0 && parts[i]==parts[i+1]) check = 1;
			if (i>0 && (parts[i]==parts[i-1] || parts[i]==parts[i+1])) check = 1;
			//									printf("\n--> %d - ",i);
			//									printa(comb,h);printf("; ");printa(elem,t);printf(" alpha %d, smallest %d, check %d",alpha,smallest,check);
			if ((check && comb[0] >= alpha && smallest < n - sum[i] + 1) || !check){
				for (l=0;l<n;l++)
					if (rg[l]>=i)
						rg[l]=-1;
				for (l=0;l<h;l++)
					rg[comb[l]]=i;
				j = i + 1;
				
				/* per ogni nuova parte trova la corretta posizione iniziale.
					 In particolare, se la dimensione è uguale alla precedente, 
					 la posizione deve essere maggiore della precedente */
				
				while (j<np){
					if (parts[j] != parts[j-1])
						l=0;
					else {
						/* TODO: memorizzare le posizioni dei minimi così evito di cercarla */
						found = 0;
						l=0;
						while (l<n && !found){
							if (rg[l] == j-1)
								found = 1;
							else
								l++;
						}
					}
					k=0;
					while (l<n && k<parts[j]){
						if (rg[l] == -1){
							rg[l]=j;
							if (!k)
								Min[j] = l;
							k++;
						}					
						l++;
					}
					j++;
				}
				
				/*					if (k>=parts[j]){
						j++;
						k=0;
					}
					l++;
				}

				
				/* reimposto i minimi */
				
				/*for (j=0;j<np;j++){
					found = 0;
					for (l=0;l<n && !found; l++)
						if (rg[l]==j) found = 1;
					Min[j] = l-1;
					}*/
				return(1);
			}
		}
		i--;
	}
	return(0);
}

void ip_prune_partition(int *rg, int *Min, int n, int *parts, int np, int pos){
	
	int i;
	int k = 0;
	int p = np - 1;
	for (i=0;i<n;i++){
		if (rg[i]>pos){
			rg[i] = p;
			k++;
		}
		if (k == parts[p]){
			k = 0;
			p--;
		}
	}
}
