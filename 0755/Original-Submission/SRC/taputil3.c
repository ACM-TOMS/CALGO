/*
   --------------------------------------------------------------
   File taputil3.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   Provides the tape management routines tapestats, start_trace,
   stop_trace, etc.
*/

#include <string.h>
/* #include <malloc.h>  not always necessary */

#ifdef __cplusplus
extern "C" {
#endif

#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h"
#include "oplate.h"
#include "tayutil.h"
#include "taputil1.h"
#include "taputil2.h"


/* Max number of tapes currently in use */ 
static int maxtapes = 0;

/* File Names */
static char op_file[20];
static char int_file[20];
static char val_file[20];

/* Arrays of pointers to the various tape buffers */
static unsigned char **op_tape;
static locint **int_tape;
static double **val_tape;

/* Array of pointers to the stats arrays (of size 11) */
static int **stats;


void fail (int error)
{
  switch (error)
  {
  case -1:
    printf("ADOL-C error: Malloc of memory failed!");
    exit (0);
  }
}

/* 
   int2asc converts the integer num to a string, places it
   in the array string, and returns the pointer to the 
   string.  (I now that this is rather redundant, but I
   didn't write the original code for this.-- DWJ ;-)                         
*/
char* int2asc(int num, char string[])
{
  sprintf(string,"%d",num);
  return(string);
}


/********************************************************************/
/* The subroutine get_fstr appends to the filename FNAME            */
/* (found in usrparms.h) the number fnum, and puts the resulting    */
/* string in fstr.                                                  */
/********************************************************************/ 
void  get_fstr (char *fstr,short fnum)

/**** 
  The caller of this function is responsible for allocating the appropriate 
  amount of storage for fstr [strlen(FNAME)+1 <= strlen(fstr) 
                                              <= strlen(FNAME)+5] 
****/
{
  char tstr[10];

  if (fnum)
  {
    strcpy (fstr,FNAME);
    int2asc (fnum,tstr);
    strcat (fstr,tstr);
  }
  else
  {
    strcpy (fstr,FNAME);
    fstr[strlen(fstr)-1] = '\0';
  }

}
/********************************************************************/
/* The subroutine get_fstr appends to the filename FNAME            */
/* (found in usrparms.h) the number fnum, and puts the resulting    */
/* string in fstr.                                                  */
/********************************************************************/ 
void  get_fstr1 (char *fstr,short fnum)

/**** 
  The caller of this function is responsible for allocating the appropriate 
  amount of storage for fstr [strlen(FNAME)+1 <= strlen(fstr) 
                                              <= strlen(FNAME)+5] 
****/
{
  char tstr[10];

  if (fnum)
  {
    strcpy (fstr,FNAME1);
    int2asc (fnum,tstr);
    strcat (fstr,tstr);
  }
  else
  {
    strcpy (fstr,FNAME1);
    fstr[strlen(fstr)-1] = '\0';
  }

}
/********************************************************************/
/* The subroutine get_fstr appends to the filename FNAME            */
/* (found in usrparms.h) the number fnum, and puts the resulting    */
/* string in fstr.                                                  */
/********************************************************************/ 
void  get_fstr2 (char *fstr,short fnum)

/**** 
  The caller of this function is responsible for allocating the appropriate 
  amount of storage for fstr [strlen(FNAME)+1 <= strlen(fstr) 
                                              <= strlen(FNAME)+5] 
****/
{
  char tstr[10];

  if (fnum)
  {
    strcpy (fstr,FNAME2);
    int2asc (fnum,tstr);
    strcat (fstr,tstr);
  }
  else
  {
    strcpy (fstr,FNAME2);
    fstr[strlen(fstr)-1] = '\0';
  }

}
void get_op_stats(int tag, char **ret_op_file, int *ret_op_len, 
int *ret_op_access,unsigned char **ret_op_tape)
{
  get_fstr(op_file,tag);
  *ret_op_file = op_file;
  *ret_op_len = stats[tag][5];
  *ret_op_access = stats[tag][6];
  *ret_op_tape = op_tape[tag];

}
void get_int_stats(int tag, char **ret_int_file, int *ret_int_len, 
int *ret_int_access,locint **ret_int_tape)
{
  get_fstr1(int_file,tag);
  *ret_int_file = int_file;
  *ret_int_len = stats[tag][7];
  *ret_int_access = stats[tag][8];
  *ret_int_tape = int_tape[tag];

}
void get_val_stats(int tag, char **ret_val_file, int *ret_val_len, 
int *ret_val_access,double **ret_val_tape)
{
  get_fstr2(val_file,tag);
  *ret_val_file = val_file;
  *ret_val_len = stats[tag][9];
  *ret_val_access = stats[tag][10];
  *ret_val_tape = val_tape[tag];
}

static int tag;

static void init_stat_space(short tnum)
{
  unsigned char **t1;          /* t1,t2,t3 and t4 are temporaries */
  double **t2;
  locint **t3;
  int **t4;
  int jj,

  tag = tnum;
  /* Set up space for */ 
  if (maxtapes==0) /*this is only done at first call to start_trace or
                     init_stat_space */
    {
      maxtapes = 10;
      if ((op_tape = (unsigned char **)malloc(maxtapes*sizeof(unsigned char*)))==0) 
	fail(-1);
      if ((int_tape = (locint **)malloc(maxtapes*sizeof(locint *)))==0)
	fail(-1);
      if ((val_tape = (double **)malloc(maxtapes*sizeof(double *)))==0)
	fail(-1);
      
      if ((stats = (int**)malloc(maxtapes*sizeof(int*)))==0)
	fail(-1);
      for(jj=0;jj<maxtapes;jj++)
	{
	  op_tape[jj]=0;
	  int_tape[jj] = 0;
	  val_tape[jj]= 0;
	  stats[jj] = 0;
	}
    }
  
  if (tag>=maxtapes)
    {
      t1=op_tape;
      t3=int_tape;
      t2=val_tape;
      t4 = stats;
      if ((op_tape =(unsigned char**)malloc(2*tag*sizeof(unsigned char*)))==0) 
	fail(-1);
      if ((int_tape = (locint **)malloc(2*tag*sizeof(locint *)))==0)
	fail(-1);
      if ((val_tape = (double **)malloc(2*tag*sizeof(double *)))==0)
	fail(-1);
      if ((stats = (int**)malloc(tag*2*sizeof(int*)))==0)
	fail(-1);
      
      for (jj=0;jj<maxtapes;jj++)
	{
	  op_tape[jj]=t1[jj];
	  int_tape[jj]=t3[jj];
	  val_tape[jj]=t2[jj];
	  stats[jj]=t4[jj];
	}
      free((char *)t1);free((char *)t2);free((char *)t3);free((char *)t4);

      for(jj=maxtapes;jj<2*tag;jj++)
	{
	  op_tape[jj]=0;
	  int_tape[jj]=0;
	  val_tape[jj] = 0;
	  stats[jj] = 0;
	}
      maxtapes = 2*tag;
    }
}

static void set_up_buffers(short tag, int buffer_size)
{
  /* Return old memory ... if used */
  if(op_tape[tag]) free((char*)op_tape[tag]);
  if(int_tape[tag]) free((char*)int_tape[tag]);
  if(val_tape[tag]) free((char*)val_tape[tag]);
  if(stats[tag]) free((char*)stats[tag]);
  
  op_tape[tag] = (unsigned char *)malloc(buffer_size*sizeof(unsigned char));
  int_tape[tag]= (locint *)malloc(buffer_size*sizeof(locint));
  val_tape[tag]= (double *)malloc(buffer_size*sizeof(double));
  stats[tag] = (int*)malloc(11*sizeof(int));
}

static void read_tape_stats(short tag, int *stats)
{
  char int_file[20];

  FILE *int_file_in;

  get_fstr1(int_file,tag);
  
  if ((int_file_in = fopen(int_file,"r"))==0) 
    {
      printf("ADOL-C error: Error reading integer tape number %d \n",tag);
      printf("Fopen returned error number %d \n",tag);
      exit(-1);
    }
  if (fread((char *)stats,11*sizeof(int),1,int_file_in)!= 1)
    {
      printf("ADOL-C error: Error reading integer tape number %d \n",tag);
      printf("Fread returned error number %d \n",tag);
      exit(-1);
    }
  fclose(int_file_in);
}

/*-------------------------------------------------------------------------*/
/* Tapestats:                                                              */
/* Returns statistics on the tape tag.  The array tape_stat is assumed to  */
/* contain at least 11 elements.  The elements of the array are the        */
/* following.                                                              */
/* tape_stat[0] = # of independent variables.                              */
/* tape_stat[1] = # of dependent variables.                                */
/* tape_stat[2] = max # of live variables.                                 */
/* tape_stat[3] = # of death notices.                                      */
/* tape_stat[4] = buffer size (# of chars, # of doubles, # of locints)     */
/* tape_stat[5] = # of operations.                                         */ 
/* tape_stat[6] = operation file access flag (1 = file in use, 0 otherwise)*/
/* tape_stat[7] = # of saved locations.                                    */ 
/* tape_stat[8] = location file access flag (1 = file in use, 0 otherwise) */
/* tape_stat[9] = # of saved constant values.                              */ 
/* tape_stat[10]= value file access flag (1 = file in use, 0 otherwise)    */
/*                                                                         */
/*-------------------------------------------------------------------------*/

void tapestats(short tag,int *tape_stat)
{ 
  int i;

  /* Make sure that there is tape access */
  init_stat_space(tag);

  if(stats[tag] == 0) 
     { 
       /* Tape number does not exist , so read in tape data */

       read_tape_stats(tag,tape_stat);

       set_up_buffers(tag,tape_stat[4]);

       /* Copy data to stats for future use */

       for (i = 0;i <11;i++)
	 {
	   stats[tag][i]= tape_stat[i];
	 }
     }
  else
    {
      for (i=0; i<11;i++)
	{
	  tape_stat[i] = stats[tag][i];
	}
    }
}


  
/***************************************************************************/
/* start_trace: (part of trace_on)                                         */
/* Initialization for the taping process.  Sets up the arrays op_tape,     */
/* int_tape, val_tape, and stats.  Op_tape, int_tape, val_tape are arrays  */
/* of pointers to individual buffers for operations, integers (locints),   */
/* and values (doubles).  Also initializes buffers for this tape, sets     */
/* files names, and calls appropriate setup routines.                      */
/***************************************************************************/

void start_trace(short tnum,int revals)
{
  int revalso;
  int kk;

  double** dum = 0;
  int degree = 0;
  
  revalso = revals;
  tag = tnum;
  
  /* Set buffer size to be the default in usrparms.h */
  set_buf_size(bufsize);

  get_fstr(op_file,tag);
  get_fstr1(int_file,tag);
  get_fstr2(val_file,tag);
  
  init_stat_space(tag);

  /* Return old memory ... if used */
  if(op_tape[tag]) free((char*)op_tape[tag]);
  if(int_tape[tag]) free((char*)int_tape[tag]);
  if(val_tape[tag]) free((char*)val_tape[tag]);
  if(stats[tag]) free((char*)stats[tag]);
  
  op_tape[tag] = (unsigned char *)malloc(bufsize*sizeof(unsigned char));
  int_tape[tag]= (locint *)malloc(bufsize*sizeof(locint));
  val_tape[tag]= (double *)malloc(bufsize*sizeof(double));
  stats[tag] = (int*)malloc(11*sizeof(int));
  
  
  clear_stats(revalso);
  
  /* Initialize Tapes */
  set_buffers(op_file,op_tape[tag],
	      int_file,int_tape[tag],
	      val_file,val_tape[tag]);
  
  
  /* Put operation denoting the start_of_the tape */ 
  put_op(start_of_tape);
  /* Leave space for the stats */
  for (kk=0;kk<22;kk++)
    put_locint(0);
  
  
  if (revalso) taylor_begin(bufsize,dum,degree);
}


/*************************************************************************/
/* Stop Tracing.  Clean up, and turn off trace_flag.                    **/
/*************************************************************************/
void stop_trace(int locations, int flag)
{
  int tape_stats[11];
  int i,sizer;

  int loc_ptr,dep_ptr,ind_ptr,death_ptr,revalso;
  /* int op_cnt,loc_cnt,val_cnt,access_ptr; */

  loc_ptr = locations;     
  put_op(end_of_tape);        /* Mark end of tape. */

  get_write_stats(&dep_ptr,&ind_ptr,&death_ptr,&revalso);
  tape_stats[0]=ind_ptr;
  tape_stats[1]=dep_ptr;
  tape_stats[2]=loc_ptr;
  tape_stats[3]=death_ptr;
  tape_stats[4]=bufsize;
  close_tape(tape_stats,flag); /** closes the tape, files up stats, and
                                   writes the tape stats to the integer
                                   tape. **/
  
  sizer = sizeof(revreal)*(death_ptr);
  if (revalso) taylor_close(sizer/(1+sizer/bufsize),dep_ptr,ind_ptr);
  for (i=0;i<11;i++)
    {
      stats[tag][i]=tape_stats[i];
    }
}

#ifdef __cplusplus
}
#endif

