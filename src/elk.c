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

#include "dp.h"
#include "idp.h"
#include "grasp.h"
#include "ip.h"
#include "sandholm.h"

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <time.h>
//#include <sys/timex.h>

#define ELK_VERSION_No "1.0"


int __verbose = 0;
extern unsigned long int *powers;

/* n = 4

   l=1   [1234]

   l=2   [1123] [1213] [1231] [1223] [1232] [1233]

   l=3   [1112] [1121] [1211] [2111] [1122] [1212] [1221]

   l=4   [1111]


   [1123]
   [1213]
   [1231]
   [1223]
   [1232]
   [1233]
   

n = 5

   [11123]
   [11213]
   [11231]
   [12113]
   [12131]
   [12311]
   [12223]
   [12232]
   [12333]
   [12233]
   [12323]
   [12332]
   


 */

/* setpartitions 

Djoki B, Miyakawa M, Sekiguchi S, Semba I, Stojmenovi ́ I. A fast iterative algorithm
for generating set partitions. Comput J 1989;32(3):281–282.

 Procedura per individuare tutte le partizioni di un set di n elementi
*/

void
sp (int *codewords, int n)
{
  int r = 1;
  codewords[0] = 1;
  int j = 0;
  int *b = (int *) calloc (n, sizeof (int));
  b[0] = 1;
  int n1 = n - 1;
  int i, jj;

  do
    {
      while (r < n1)
	{
	  r++;
	  codewords[r - 1] = 1;
	  j++;
	  b[j] = r;
	}
      for (i = 1; i <= n - j; i++)
	{
	  codewords[n - 1] = i;
	  printf ("\n");
	  for (jj = 0; jj < n; jj++)
	    printf ("%d ", codewords[jj]);
	}
      r = b[j];
      codewords[r - 1]++;
      if (codewords[r - 1] > (r - j))
	j--;
    }
  while (r > 1);
}

/* setpartitions 

Djoki B, Miyakawa M, Sekiguchi S, Semba I, Stojmenovi ́ I. A fast iterative algorithm
for generating set partitions. Comput J 1989;32(3):281–282.

see also
 M.C. Er, A fast algorithm for generationg set partitions, 1988

 Procedura per individuare tutte le partizioni di un set di n elementi di k classi
*/

void
spk (int *codewords, int n, int k)
{
  int r = 1;
  codewords[0] = 1;
  int j = 0;
  int *b = (int *) calloc (n, sizeof (int));
  b[0] = 1;
  int n1 = n - 1;
  int i, jj;

  do
    {
      while (r < n1)
	{
	  r++;
	  codewords[r - 1] = 1;
	  j++;
	  b[j] = r;
	}
      for (i = 1; i <= n - j && i <= k; i++)
	{
	  codewords[n - 1] = i;
	  int kk;
	  int trovato;

	  trovato = 1;
	  for (jj = 1; jj <= k && trovato; jj++)
	    {
	      trovato = 0;
	      for (kk = 0; kk < n && !trovato; kk++)
		{
		  if (codewords[kk] == jj)
		    trovato++;
		}
	    }
	  /*      if (trovato)
	     {
	     printf("\n");
	     for (jj=0;jj<n;jj++)
	     printf("%d ",codewords[jj]);
	     } */
	}
      while (codewords[b[j] - 1] == k)
	{
	  r--;
	  j--;
	}
      r = b[j];
      codewords[r - 1]++;
      if (codewords[r - 1] > r - j)
	j--;
    }
  while (r > 1);
}



static void
display_usage (void)
{
  printf
    ("\nNAME\n\telk - A system for coalition formation algorithms (Ver. %s)",
     ELK_VERSION_No);
  printf ("\n\nSYNOPSIS\n\t elk [OPTION] ...\n\nDESCRIPTION");

  printf
    ("\n\tIt  executes  some  coalition  formation  algorithm on  games  with  a");
  printf
    ("\n\tcharacteristic  function  (CFG).  In  particular,  the  algorithms  it");
  printf
    ("\n\timplements are:  DP, IDP, IP  and GRASP. They  can be executed  on the");
  printf
    ("\n\tsame CF and  in sequence by listing in the  parameter oprions the name");
  printf ("\n\tof the algorithms.\n");

  printf ("\n--help");
  printf ("\n\tDisplay this text and exit");

  printf ("\n--agents nAgents");
  printf ("\n\tSolve a CFG with nAgents agents (default 10)");

  printf ("\n--ucf");
  printf
    ("\n\tAutomatically generates a uniform characteristic function (the default).");

  printf ("\n--uscf");
  printf
    ("\n\tAutomatically generates a uniform scaled characteristic function.");

  printf ("\n--ncf");
  printf
    ("\n\tAutomatically generates a normal characteristic function (the default).");

  printf ("\n--nscf");
  printf
    ("\n\tAutomatically generates a normal scaled characteristic function.");

  printf ("\n--ndcf");
  printf
    ("\n\tAutomatically generates a normally distributed characteristic function.");

  printf ("\n--dp");
  printf
    ("\n\tExecute the DP algorithm and compute the value of the best solution.");

  printf ("\n--sa");
  printf ("\n\tExecute the Sandholm algorithm.");

  printf ("\n--idp");
  printf
    ("\n\tExecute the IDP algorithm and compute the value of the best solution.");

  printf ("\n--ip");
  printf
    ("\n\tExecute the IP algorithm and compute the value of the best solution.");

  printf ("\n  --ipbeta");
  printf ("\n\tBeta star value for IP in [0,1] (default 1.0)");

  printf ("\n--grasp");
  printf ("\n\tExecute the grasp algorithm");

  printf ("\n  --gimin");
  printf ("\n\tMinumun number of iterations when set grasp (default 5)");

  printf ("\n  --gimax");
  printf ("\n\tMaximum number of iterations when set grasp (default 5)");

  //printf ("\n  --gamin");
  //printf ("\n\tMinumun alpha set grasp (default 0.9)");

  //printf ("\n  --gamax");
  //printf ("\n\tMaximum alpha when set grasp (default 0.9)");

  printf ("\n  --save");
  printf ("\n\tSave to the file the cgf (default out)");

  printf ("\n  --load");
  printf ("\n\tLoad the file containing the cgf (default out)");
  

  printf ("\n--verbose=VAL");
  printf ("\n\tSet verbosity to VAL");

  printf ("\n--version");
  printf ("\n\tPrint version number and exit\n");


  printf ("\n\nCOPYRIGHT\n\tCopyright © 2009 Nicola Di Mauro.");
  printf
    ("\n\tLicense GPL: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>");
  printf
    ("\n\tThis is free software: you are free to change and redistribute it.\n\tThere is NO WARRANTY, to the extent permitted by law.");

  printf ("\n\nAUTHOR");
  printf ("\n\tWritten by Nicola Di Mauro.");

  printf ("\n\nREPORTING BUGS");
  printf ("\n\tReport bugs to <ndm@di.uniba.it>.\n");

  exit (10);

  exit (EXIT_SUCCESS);
}


static void
show_version ()
{
  printf ("elk v. %s", ELK_VERSION_No);
  printf ("\nCopyright © 2009 Nicola Di Mauro.  License GPL.\n");
}


int
main (int argc, char *argv[])
{

  struct rusage earlier, later, interval;
  extern double grasp_nodes_succ;
  extern double grasp_nodes_neigh;
  /*int *ipp = (int *) calloc(6, sizeof(int));
     int f;
     ipp[0]=1;
     ipp[1]=1;
     while (idp_next_split_structure1 (ipp, 6, 2))
     {
     printf("\n");
     for (f=0;f<6;f++)
     printf("%d ",ipp[f]);
     }
     exit(1);

   */


  /*  int *ip1 = (int *) calloc(5, sizeof(int));
     int *ipm = (int *) calloc(5, sizeof(int));

     firstKCoalitionStructure(ip1,ipm,5,3);
     int f;
     printf("\n");  
     for (f=0;f<5;f++)
     printf("%d ",ip1[f]);
     printf(" - ");
     for (f=0;f<5;f++)
     printf("%d ",ipm[f]);
     while (nextKCoalitionStructure(ip1, ipm, 5, 3))  
     {
     printf("\n");  
     for (f=0;f<5;f++)
     printf("%d ",ip1[f]);
     printf(" - ");
     for (f=0;f<5;f++)
     printf("%d ",ipm[f]);
     }

     exit(1);
   */
  /*  int *ip1 = (int *) calloc(5, sizeof(int));
     int *ipm = (int *) calloc(5, sizeof(int));

     int q;
     for (q=1;q<=5;q++)
     {
     printf("\n%d",q);fflush(stdout);
     firstKCoalitionStructure(ip1,ipm,5,q);
     int f;
     printf("\n\t");  
     for (f=0;f<5;f++)
     printf("%d ",ip1[f]);
     printf(" - ");
     for (f=0;f<5;f++)
     printf("%d ",ipm[f]);
     while (nextKCoalitionStructure(ip1, ipm, 5, q))
     {
     printf("\n\t");  
     for (f=0;f<5;f++)
     printf("%d ",ip1[f]);
     printf(" - ");
     for (f=0;f<5;f++)
     printf("%d ",ipm[f]);
     }
     }
     exit(1);

   */

  unsigned int nAgents;

  unsigned int k;
  //double alpha, alphaMin, alphaMax;
  int Iterations, minIterations, maxIterations, runs;

  struct globalArgs_t
  {
    int nAgents;
    int ucf;
    int uscf;
    int ncf;
    int nscf;
    int ndcf;
    int dp;
    int idp;
    int ip;
    int sa;
    double ipbeta;
    int grasp;
    int gimin;
    int gimax;
    //    double gamin;
    //    double gamax;
    char * fileOut;
    char * fileIn;
  } global_args;

  global_args.nAgents = 10;
  global_args.ucf = 1;
  global_args.uscf = 0;
  global_args.ncf = 0;
  global_args.nscf = 0;
  global_args.ndcf = 0;
  global_args.dp = 0;
  global_args.idp = 0;
  global_args.ip = 0;
  global_args.sa = 0;
  global_args.ipbeta = 1.0;
  global_args.grasp = 0;
  global_args.gimin = 5;
  global_args.gimax = 5;
  global_args.fileOut = "out";
  global_args.fileIn = "out";

  //  global_args.gamin = 0.9;
  //  global_args.gamax = 0.9;


  static struct option long_options[] = {
    {"verbose", required_argument, 0, 1},
    {"help", no_argument, 0, 2},
    {"agents", required_argument, 0, 3},
    {"ucf", no_argument, 0, 4},
    {"uscf", no_argument, 0, 5},
    {"ncf", no_argument, 0, 6},
    {"nscf", no_argument, 0, 7},
    {"ndcf", no_argument, 0, 8},
    {"dp", no_argument, 0, 9},
    {"idp", no_argument, 0, 10},
    {"ip", no_argument, 0, 11},
    {"sa", no_argument, 0, 19},
    {"ipbeta", required_argument, 0, 18},
    {"grasp", no_argument, 0, 12},
    {"gimin", required_argument, 0, 13},
    {"gimax", required_argument, 0, 14},
    //    {"gamin", required_argument, 0, 15},
    //    {"gamax", required_argument, 0, 16},
    {"version", no_argument, 0, 17},
    {"load", required_argument, 0, 20},
    {"save", required_argument, 0, 21}
  };

  while (1)
    {


      /* getopt_long stores the option index here. */
      int option_index = 0;

      int c = getopt_long (argc, argv, "", long_options, &option_index);

      /* printf("\nc %d, optind %d, optarg %s, opterr %d, optopt %d\n",c,optind,optarg,opterr,optopt); */

      /* Detect the end of the options. */
      if (c == -1)
	break;

      switch (c)
	{
	case 0:
	  /* If this option set a flag, do nothing else now. */
	  if (long_options[option_index].flag != 0)
	    break;
	  printf ("option %s", long_options[option_index].name);
	  if (optarg)
	    printf (" with arg %s", optarg);
	  printf ("\n");
	  break;

	case 1:
	  __verbose = atoi (optarg);
	  break;

	case 2:
	  display_usage ();
	  break;

	case 3:
	  global_args.nAgents = atoi (optarg);
	  break;

	case 4:
	  global_args.ucf = 1;
	  global_args.uscf = 0;
	  global_args.ncf = 0;
	  global_args.nscf = 0;
	  global_args.ndcf = 0;
	  break;

	case 5:
	  global_args.ucf = 0;
	  global_args.uscf = 1;
	  global_args.ncf = 0;
	  global_args.nscf = 0;
	  global_args.ndcf = 0;
	  break;

	case 6:
	  global_args.ucf = 0;
	  global_args.uscf = 0;
	  global_args.ncf = 1;
	  global_args.nscf = 0;
	  global_args.ndcf = 0;
	  break;

	case 7:
	  global_args.ucf = 0;
	  global_args.uscf = 0;
	  global_args.ncf = 0;
	  global_args.nscf = 1;
	  global_args.ndcf = 0;
	  break;

	case 8:
	  global_args.ucf = 0;
	  global_args.uscf = 0;
	  global_args.ncf = 0;
	  global_args.nscf = 0;
	  global_args.ndcf = 1;
	  break;

	case 9:
	  global_args.dp = 1;
	  break;

	case 10:
	  global_args.idp = 1;
	  break;

	case 11:
	  global_args.ip = 1;
	  break;

	case 12:
	  global_args.grasp = 1;
	  break;

	case 13:
	  global_args.gimin = atoi (optarg);
	  break;

	case 14:
	  global_args.gimax = atoi (optarg);
	  break;

	  //	case 15:
	  //	  global_args.gamin = atof (optarg);
	  //	  break;

	  //	case 16:
	  //	  global_args.gamax = atof (optarg);
	  //	  break;

	case 18:
	  global_args.ipbeta = atof (optarg);
	  break;

	case 19:
	  global_args.sa = 1;
	  break;


	case 17:
	  show_version ();
	  exit (EXIT_SUCCESS);
	  break;

	case 20:
	  global_args.fileIn = optarg;
	  break;

	case 21:
	  global_args.fileOut = optarg;
	  break;

	case '?':
	  printf ("\nUse \"elk --help\" for more information.\n");
	  exit (EXIT_SUCCESS);
	  break;
	default:
	  display_usage ();
	  break;
	}
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
	printf ("%s ", argv[optind++]);
      putchar ('\n');
      exit (1);
    }

  show_version ();

  if (global_args.dp == 0 && global_args.idp == 0 && global_args.ip == 0
      && global_args.grasp == 0)
    {
      printf ("\nSpecify at least one algorithm");
      printf ("\nUse \"elk --help\" for more information.\n");
      exit (EXIT_SUCCESS);
    }

  /* print parameters */

  printf
    ("\n[agents:%d ucf:%d uscf:%d ncf:%d nscf:%d ndcf:%d dp:%d idp:%d sa:%d ip:%d ipbeta:%.2f grasp:%d gimin:%d gimax:%d fout:%s fin:%s]\n",
     global_args.nAgents, global_args.ucf, global_args.uscf, global_args.ncf,
     global_args.nscf, global_args.ndcf, global_args.dp, global_args.idp,
     global_args.sa, global_args.ip, global_args.ipbeta, global_args.grasp,
     global_args.gimin, global_args.gimax, global_args.fileOut, global_args.fileIn);

  nAgents = global_args.nAgents;

  powers = (unsigned long int *) calloc (nAgents, sizeof (unsigned long int));
  for (k = 0; k < nAgents; k++)
    powers[k] = pow (2, k);

  characteristic_function = chf_new (nAgents);
  if (global_args.ucf)
    chf_init_uniform (nAgents, &characteristic_function);
  else if (global_args.uscf)
    chf_init_uniform_scaled (nAgents, &characteristic_function);
  else if (global_args.ncf)
    chf_init_normal (nAgents, &characteristic_function);
  else if (global_args.nscf)
    chf_init_normal_scaled (nAgents, &characteristic_function);
  else if (global_args.ndcf)
    chf_init_normal_distributed (nAgents, &characteristic_function);

  chf_print(nAgents,&characteristic_function,global_args.fileOut);


  if (global_args.sa)
    {
      getrusage (RUSAGE_SELF, &earlier);
      sandholm (nAgents);
      getrusage (RUSAGE_SELF, &later);

      printf ("\n>Time: %.6f s, %lld micros",
	      (double) timeval_diff (&interval, &later, &earlier) / 1000000,
	      timeval_diff (&interval, &later, &earlier));
      printf (" (%ld seconds, %ld microseconds)\n",
	      interval.ru_utime.tv_sec, interval.ru_utime.tv_usec);
    }

  double bestVal;
  if (global_args.dp)
    {
      getrusage (RUSAGE_SELF, &earlier);
      bestVal = dp (nAgents);
      getrusage (RUSAGE_SELF, &later);

      printf ("\n>Time: %.6f s, %lld micros",
	      (double) timeval_diff (&interval, &later, &earlier) / 1000000,
	      timeval_diff (&interval, &later, &earlier));
      printf (" (%ld seconds, %ld microseconds)\n",
	      interval.ru_utime.tv_sec, interval.ru_utime.tv_usec);
    }

  if (global_args.idp)
    {
      getrusage (RUSAGE_SELF, &earlier);
      bestVal = idp (nAgents);
      getrusage (RUSAGE_SELF, &later);

      printf ("\n>Time: %.6f s, %lld micros",
	      (double) timeval_diff (&interval, &later, &earlier) / 1000000,
	      timeval_diff (&interval, &later, &earlier));
      printf (" (%ld seconds, %ld microseconds)\n",
	      interval.ru_utime.tv_sec, interval.ru_utime.tv_usec);
    }

  if (global_args.ip)
    {
      getrusage (RUSAGE_SELF, &earlier);
      double ipval = ip (nAgents, (1.0 / global_args.ipbeta));
      getrusage (RUSAGE_SELF, &later);

      printf ("\n>Time: %.6f s, %lld micros",
	      (double) timeval_diff (&interval, &later, &earlier) / 1000000,
	      timeval_diff (&interval, &later, &earlier));
      printf (" (%ld seconds, %ld microseconds)\n",
	      interval.ru_utime.tv_sec, interval.ru_utime.tv_usec);
      if (global_args.dp || global_args.idp)
	printf ("R: %.6f\n", ipval / bestVal);

    }

  if (global_args.grasp)
    {
      minIterations = global_args.gimin;
      maxIterations = global_args.gimax;
      //alphaMin = global_args.gamin;
      //alphaMax = global_args.gamax;

      double bestGrasp;
      printf ("\n== GRASP ==");
      long long int savedtime = time (NULL);
      for (Iterations = minIterations; Iterations <= maxIterations;
	   Iterations++)
	{
	  //for (alpha = alphaMin; alpha <= alphaMax; alpha += 0.05)
	  //  {
	  savedtime++;
	  runs = 0;
	  getrusage (RUSAGE_SELF, &earlier);
	  do
	    {
	      runs++;
	      srand (savedtime);
	      bestGrasp = grasp (nAgents,  Iterations);
	      getrusage (RUSAGE_SELF, &later);
	    }
	  while (((double) timeval_diff (&interval, &later, &earlier) /
		  1000000) < 0.125);
	  
	  printf ("\nI:%d, r:%d,", Iterations, runs);
	  //min = bestGrasp; max = bestGrasp; mean = bestGrasp;
	  
	  printf (" T: %.8f,",
		  (double) timeval_diff (&interval, &later,
					 &earlier) / 1000000 / runs);
	  //              printf(" Vmin: %.6f, Vmax: %.6f, Vavg: %.6f",min, max, mean / runs);
	  printf (" V: %.6f", bestGrasp);
	  if (global_args.dp || global_args.idp)
	    printf (", R: %.6f", bestGrasp / bestVal);
	  printf (", nodes %lld (s:%lld n:%lld)",
		      ((long long int) grasp_nodes_succ +
		       (long long int) grasp_nodes_neigh),
		  (long long int) grasp_nodes_succ,
		  (long long int) grasp_nodes_neigh);
	  printf (" srand(%lld)", savedtime);
	  //	    }
	}
    }
  printf ("\n");
  return (EXIT_SUCCESS);
}
