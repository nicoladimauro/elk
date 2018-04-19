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

#include "ip.h"

/* The IP algorihtm 
   
   1. Scanning the input in order to compute the bounds (i.e. MAXG and AVGG ) for every 
      sub-space PG — while doing so, we can (at a very small cost):
    (a) find the best coalition structures within particular sub-spaces.
    (b) prune other sub-spaces based on their upper-bounds.
    (c) establish a worst-case bound on the quality of the best solution found so far.
   2. Searching within the remaining sub-spaces — the techniques we use allow us to:
    (a) avoid making unnecessary comparisons between coalitions to generate valid coalition
        structures (i.e. those that contain disjoint coalitions).
    (b) avoid computing the same coalition structure more than once.
    (c) apply branch-and-bound to further reduce the amount of search to be done.

*/

double
ip (int nAgents, double beta_star)
{
  /* scan_and_search */

  /* le coalizioni sono 2^n-1. In C da 0 a 2^n-2 */
  unsigned long long int numCoalitions = pow (2, nAgents) - 2;

  int i;
  cs best_cs = cs_new (nAgents), current_cs;

  double *max;
  double *avg;
  unsigned long int *counts;
  double cs_value, best_cs_value;


  printf ("\n== IP ==");
  fflush (stdout);


  max = (double *) calloc (nAgents, sizeof (double));
  avg = (double *) calloc (nAgents, sizeof (double));
  counts = (unsigned long int *) calloc (nAgents, sizeof (unsigned long int));


  /**
     scanning the input  
  */

	/* subspace P1 */
  double P1 = characteristic_function[numCoalitions];
  /* subspace PN */
	double PN = 0.0;
  for (i = 0; i < nAgents; i++)
    PN += characteristic_function[powers[i] - 1];

  if (__verbose)
    printf ("\n P1 %f PN %f", P1, PN);

  /** selecting the best CS between P1 and PN */
  if (PN > P1){
		best_cs_value = PN;
		for (i = 0; i < nAgents; i++)	{
			best_cs[i] = i + 1;
			best_cs[nAgents + i] = 1;
		}
		best_cs[nAgents * 2] = nAgents;
		best_cs[nAgents * 2 + 1] = nAgents;
	} else {
		best_cs_value = P1;
		for (i = 0; i < nAgents; i++)
			best_cs[i] = 1;
		best_cs[nAgents] = nAgents;
		best_cs[nAgents * 2] = nAgents;
		best_cs[nAgents * 2 + 1] = 1;
	}

  /**
     search through level P2 
  */

  for (i = 0; i < nAgents - 1; i++) {
		max[i] = LONG_MIN;
		avg[i] = 0.0;
		counts[i] = 0;
	}
  max[nAgents - 1] = PN;
  avg[nAgents - 1] = PN;
  counts[nAgents - 1] = 1;

  int *ip = (int *) calloc (nAgents + 1, sizeof (int));
  int *ipm = (int *) calloc (nAgents + 1, sizeof (int));


  int *MinPos = (int *) calloc (nAgents + 1, sizeof (int));
  int *SumParts = (int *) calloc (nAgents + 1, sizeof (int));

  int *SuppVet1 = (int *) calloc (nAgents + 1, sizeof (int));
  int *SuppVet2 = (int *) calloc (nAgents + 1, sizeof (int));

  firstKCoalitionStructure (ip, ipm, nAgents, 2);
  current_cs = cs_new (nAgents);
  int k, index, length;

  /*  printf("\n ..");
     for (i=0;i<nAgents;i++)
     printf("%d ",ip[i]); */

  /* compute the index of the CS */
  index = 0;
  length = 0;
  k = 0;
  for (i = 0; i < nAgents; i++)
    if (ip[i]) {
			index += powers[i];
			length++;
		}
  index--;

  /*  printf("Index %d, nc %d",index,numCoalitions-index-1); */
  double max_value = characteristic_function[index] + characteristic_function[numCoalitions - index - 1];
  int max_index = index;
  /*  printf("  %f",max_value); */
  if (characteristic_function[index] > max[length - 1])
    max[length - 1] = characteristic_function[index];
  if (characteristic_function[numCoalitions - index - 1] > max[nAgents - length - 1])
    max[nAgents - length - 1] = characteristic_function[numCoalitions - index - 1];
  avg[length - 1] += characteristic_function[index];
  avg[nAgents - length - 1] += characteristic_function[numCoalitions - index - 1];
  counts[length - 1]++;
  counts[nAgents - length - 1]++;

  while (nextKCoalitionStructure (ip, ipm, nAgents, 2))
    {

      /*      printf("\n ..");
							for (i=0;i<nAgents;i++)
							printf("%d ",ip[i]); */
			
      /* compute the index of the next CS */
      index = 0;
      length = 0;
      k = 0;
      for (i = 0; i < nAgents; i++)
				if (ip[i]) {
					index += powers[i];
					length++;
				}
      index--;

      cs_value = characteristic_function[index] + characteristic_function[numCoalitions - index - 1];

      /*      printf("  %f",cs_value); */

      if (cs_value > max_value)	{
				max_value = cs_value;
				max_index = index;
				if (__verbose)
					printf ("\n %f %d", max_value, max_index);
			}

      if (characteristic_function[index] > max[length - 1])
				max[length - 1] = characteristic_function[index];
      if (characteristic_function[numCoalitions - index - 1] > max[nAgents - length - 1])
				max[nAgents - length - 1] = characteristic_function[numCoalitions - index - 1];
      avg[length - 1] += characteristic_function[index];
      avg[nAgents - length - 1] +=characteristic_function[numCoalitions - index - 1];
      counts[length - 1]++;
      counts[nAgents - length - 1]++;
    }


  if (max_value > best_cs_value) {
		best_cs_value = max_value;
		int index_cs = max_index + 1, k = 0;
		while (index_cs > 0){
			best_cs[k] = index_cs % 2 + 1;
			best_cs[nAgents + best_cs[k] - 1]++;
			index_cs = index_cs / 2;
			k++;
		}
		for (i = k; i < nAgents; i++)
			best_cs[i] = 1;
		for (i = nAgents + 2; i < nAgents * 2; i++)
			best_cs[i] = 0;
		best_cs[nAgents * 2] = nAgents;
		best_cs[nAgents * 2 + 1] = 2;
	}

  //  printf("\nBest coalition: ");
  //  printf("{");
  //  for (i=0;i<nAgents;i++)
  //    printf("%d ",best_cs[i]);
  //  printf("}{");
  //  for (i=nAgents;i<nAgents*2;i++)
  //    printf("%d ",best_cs[i]);
  //  printf("} %d %d",best_cs[nAgents*2],best_cs[nAgents*2+1]);  
  //  printf(" Val: %f",best_cs_value);

  //  printf("\n");
  
	for (i = 0; i < nAgents; i++) {
		avg[i] /= counts[i];
		//      printf("\ni=%d: Max %f, AVG %f, C %d",i+1,max[i],avg[i],counts[i]); 
	}

	int h = count_partitions (nAgents);
  //  printf("\n%d Integer partitions for %d agents\n",h,nAgents);fflush(stdout);
  partition *LP = (partition *) calloc (h, sizeof (partition));
  set_partitions (LP, nAgents);
	/*  print_partitions(LP, h); */

  /** 
      compute upper and lower bounds for each sub-space (32-35) 
  */

  double *sub_space_ub = (double *) calloc (h, sizeof (double));
  double *sub_space_lb = (double *) calloc (h, sizeof (double));

  int sub_space;
  partition *sub_space_ptr = LP;
  double max_max_sub_space, max_avg_sub_space;
  int first = 1;
  for (sub_space = 0; sub_space < h; sub_space++) {
		if (__verbose){
			printf ("\n Sub-space ");
			output (*sub_space_ptr);
		}
		sub_space_ub[sub_space] = 0.0;
		sub_space_lb[sub_space] = 0.0;
		for (i = 0; i < sub_space_ptr->numparts; i++)	{
			sub_space_ub[sub_space] += max[sub_space_ptr->parts[i] - 1];
			sub_space_lb[sub_space] += avg[sub_space_ptr->parts[i] - 1];
		}

		if (sub_space_ptr->numparts > 2 && sub_space_ptr->numparts != nAgents){
			if (!first) {
	      if (sub_space_ub[sub_space] > max_max_sub_space)
					max_max_sub_space = sub_space_ub[sub_space];
	      if (sub_space_lb[sub_space] > max_avg_sub_space)
					max_avg_sub_space = sub_space_lb[sub_space];
	    } else {
	      max_max_sub_space = sub_space_ub[sub_space];
	      max_avg_sub_space = sub_space_lb[sub_space];
	      first = 0;
	    }
		}

		if (__verbose)
			printf (" -- UB: %f LB: %f", sub_space_ub[sub_space],sub_space_lb[sub_space]);
		sub_space_ptr++;
	}

  if (__verbose)
    printf ("\nMAX_max %f, MAX_avg %f", max_max_sub_space, max_avg_sub_space);

  /**
   compute best upper and lower bounds (36-37)
  */

  double ub_star, lb_star;
  if (best_cs_value > max_max_sub_space)
    ub_star = best_cs_value;
  else
    ub_star = max_max_sub_space;

  if (best_cs_value > max_avg_sub_space)
    lb_star = best_cs_value;
  else
    lb_star = max_avg_sub_space;

  if (__verbose)
    printf ("\nUB*: %f, LB*: %f", ub_star, lb_star);


  /* pruning sub-spaces */
  int *sub_space_pruned = (int *) calloc (h, sizeof (int));
  for (i = 0; i < h; i++)
    sub_space_pruned[i] = 0;
  sub_space_ptr = LP;
  int part;
  int sub_spaces_remaining = 0;
  for (part = 0; part < h; part++) {
		if ((sub_space_ptr->numparts <= 2) || (sub_space_ptr->numparts == nAgents) || (sub_space_ub[part] < lb_star))
			sub_space_pruned[part] = 1;
		else
			sub_spaces_remaining++;
		sub_space_ptr++;
	}
  /* END pruning sub-spaces */

  double beta = (double) nAgents / 2;
  if ((ub_star / best_cs_value) < beta)
    beta = ub_star / best_cs_value;

  if (__verbose)
    printf ("\nbeta=%f", beta);


  double max_selected = sub_space_ub[0];
  int part_selected = 0;
  partition *part_ptr;

  int *coal_ip = (int *) calloc (nAgents + 1, sizeof (int));

  double old_best;

	/***************************************
	 * selecting and searching a sub space *
	 ***************************************/

  partition *part_ptr_selected;
  while (sub_spaces_remaining && beta > beta_star) {
		old_best = best_cs_value;

		/* select a sub space */
		/*--------------------*/
		//      if (beta_star == 1.0)
		if (1){
			part_ptr = LP;
			part_ptr_selected = LP;
			int primo = 1;
			for (part = 0; part < h; part++){
	      if (!sub_space_pruned[part]){
					if (primo) {
						max_selected = sub_space_ub[part];
						part_selected = part;
						part_ptr_selected = part_ptr;
						primo = 0;
					}
					else if (sub_space_ub[part] > max_selected) {
						max_selected = sub_space_ub[part];
						part_selected = part;
						part_ptr_selected = part_ptr;
					}
				}
	      part_ptr++;
	    }
			sub_space_pruned[part_selected] = 1;
			sub_spaces_remaining--;

			//printf("\n Part %d %d",part_selected,sub_spaces_remaining);
			/* --------- */
			if (__verbose) {
	      printf ("\n");
	      output (*part_ptr_selected);
	    }

			ip_initialize_first(coal_ip, MinPos, nAgents, part_ptr_selected->parts, SumParts, part_ptr_selected->numparts);

			int coalitionIndex, c, k;
			double value = 0.0;
			for (c = 0; c < part_ptr_selected->numparts; c++) {
	      coalitionIndex = 0;
	      for (k = 0; k < nAgents; k++)
					if (coal_ip[k] == c)
						coalitionIndex += powers[k];
	      if (coalitionIndex != 0)
					value += characteristic_function[coalitionIndex - 1];
	    }
			if (value > best_cs_value) {
	      best_cs_value = value;
	      for (i = 0; i < nAgents; i++)
					best_cs[i] = coal_ip[i] + 1;
	      best_cs[nAgents * 2 + 1] = part_ptr_selected->numparts;
	    }
			
			int c1 = 0;
			int skip = 0;
			int pos_skip;
			//															printf("\n"); output (*part_ptr_selected); printf(" ; ");
			//															for (k = 0; k < nAgents; k++)
			//																printf("%d ",coal_ip[k]);
			while (ip_next_partition(coal_ip, MinPos, nAgents, part_ptr_selected->parts, SumParts, part_ptr_selected->numparts, SuppVet1, SuppVet2)) {
				//			while (restricted_growth_next(part_ptr_selected->parts, part_ptr_selected->numparts,coal_ip, nAgents)) {
				//																printf("\n"); output (*part_ptr_selected); printf(" ; ");
				//																				for (k = 0; k < nAgents; k++)
				//																					printf("%d ",coal_ip[k]);

				skip = 0;
				value = 0.0;
				for (c = 0; c < part_ptr_selected->numparts && !skip; c++)	{
					coalitionIndex = 0;
					for (k = 0; k < nAgents; k++)
						if (coal_ip[k] == c)
							coalitionIndex += powers[k];
					if (coalitionIndex != 0)
						value += characteristic_function[coalitionIndex - 1];
			
					/* branch and bound */
					double remaining_values = 0.0;
					for(c1 = c + 1; c1 <	part_ptr_selected->numparts; c1++)	{
						remaining_values += max[part_ptr_selected->parts[c1]-1];
					}
					if ((remaining_values + value) < best_cs_value){
						skip = 1;
						pos_skip = c;
					}
					//					if (skip) printf(" s %d",c);
				}
				if (skip && pos_skip < part_ptr_selected->numparts-2)
					ip_prune_partition(coal_ip, MinPos, nAgents, part_ptr_selected->parts, part_ptr_selected->numparts, pos_skip);

				if (value > best_cs_value){
					best_cs_value = value;
					for (i = 0; i < nAgents; i++)
						best_cs[i] = coal_ip[i] + 1;
					best_cs[nAgents * 2 + 1] = part_ptr_selected->numparts;
				}
			}
		}
      /* --------------- */
		else{
		}

		/* update ub and beta */

		if (best_cs_value > old_best)	{
			/* re pruning */
			for (part = 0; part < h; part++){
	      if (sub_space_ub[part] < best_cs_value && !sub_space_pruned[part]){
					sub_space_pruned[part] = 1;
					sub_spaces_remaining--;
				}
	    }
		}

		/* updating UB* */
		max_max_sub_space = sub_space_ub[0];
		for (part = 0; part < h; part++){
			if (sub_space_ub[part] > max_max_sub_space)
				max_max_sub_space = sub_space_ub[part];
		}
		if (best_cs_value > max_max_sub_space)
			ub_star = best_cs_value;
		else
			ub_star = max_max_sub_space;
		
		if (__verbose)
			printf ("\nUB*: %f", ub_star);

		if ((ub_star / best_cs_value) < beta)
			beta = ub_star / best_cs_value;

		if (__verbose)
			printf ("\nbeta=%f", beta);
		
	}

	
  /*part_ptr = LP;
     for (part=0;part<h;part++)
     {
     printf("\n");output(*part_ptr);
     if (!sub_space_pruned[part])
     {
     restricted_growth_first(part_ptr->parts,part_ptr->numparts,coal_ip);
     printf("\n");
     for (i=0;i<nAgents;i++){
     printf("%d ",coal_ip[i]);
     }
     // compute the value of the coalition 
     int coalitionIndex,c,k;
     double value = 0.0;
     for (c = 0; c < part_ptr->numparts; c++)
     {
     coalitionIndex = 0;
     for (k = 0; k < nAgents; k++)
     if (coal_ip[k] == c)
     coalitionIndex += powers[k];
     if (coalitionIndex != 0)
     value += characteristic_function[coalitionIndex - 1];
     }
     //   printf(" %f",value);
     if (value>best_cs_value)
     {
     best_cs_value = value;
     for (i=0;i<nAgents;i++)
     best_cs[i]=coal_ip[i]+1;
     best_cs[nAgents*2+1]=part_ptr->numparts;
     }

     while (restricted_growth_next(part_ptr->parts,part_ptr->numparts,coal_ip,nAgents))
     {

     printf("\n");
     for (i=0;i<nAgents;i++)
     printf("%d ",coal_ip[i]);

     // compute the value of the coalition 
     int coalitionIndex,c,k;
     double value = 0.0;
     for (c = 0; c < part_ptr->numparts; c++)
     {
     coalitionIndex = 0;
     for (k = 0; k < nAgents; k++)
     if (coal_ip[k] == c)
     coalitionIndex += powers[k];
     if (coalitionIndex != 0)
     value += characteristic_function[coalitionIndex - 1];
     }
     //       printf(" %f",value);
     if (value>best_cs_value)
     {
     best_cs_value = value;
     for (i=0;i<nAgents;i++)
     best_cs[i]=coal_ip[i]+1;
     best_cs[nAgents*2+1]=part_ptr->numparts;
     }
     }
     }
     part_ptr++;
     } */
  printf ("\nSolution: ");

  /*  for (i=0;i<nAgents;i++)
     printf("%d ",best_cs[i]);
     printf("   - %d\n",best_cs[nAgents*2+1]);
   */
  int p;
  for (p = 1; p <= best_cs[nAgents * 2 + 1]; p++)  {
		printf ("{");
		for (i = 0; i < nAgents; i++)
			if (best_cs[i] == p)
				printf ("%d ", i + 1);
		printf ("}");
	}
  printf ("  %f", best_cs_value);

  return (best_cs_value);


  free_partitions (LP, h);
}
