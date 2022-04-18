/*
   --------------------------------------------------------------
   File tayutil.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   Taylor series utilities - primarily called from the module
   hos_forward.c (--- a forward pass generates taylor coefficients which
   need to be saved when a variable dies, or is overwritten) and 
   the reverse modules (-- to retrieve these taylor coefficients to calculate
   the adjoints.
*/

#ifdef __cplusplus
extern "C" {
#endif

#include "dvlparms.h" /* Developers Parameters */

#include "usrparms.h"
#include "tayutil.h"


/* Global Variables to this routine */

static int numdep;
static int numind;
static int T_file_access = 0;
static FILE* temp2_file=0;
static int taylor_cnt;
static revreal * save_taylor = 0;
static int T_write_cnt;
static int T_blocks, T_tail, T_buf_size, T_length; 
static double** T;
static revreal ** Tr;
static revreal * Trs;
static int degsave;


/* Access routines */

/* Has a taylor file been written? */

int taylor_access()
{
  return T_file_access;
}

/* Close any open taylor file. */

void close_taylor()
{
  fclose(temp2_file);
}

/******************************************************************/
/** T_Put_Block puts a block of tape to the disk.  I assume this **/
/** is called only during a successive forward pass,             **/
/** computation.                                                 **/
/******************************************************************/
static void T_put_block(int nitems)
{
  int n;
  if(T_file_access == 0) temp2_file = fopen("adoltemp.xxx","w+");
  if(T_write_cnt == 0) fseek(temp2_file,0,0); 
  T_file_access = 1;
  taylor_cnt=0;
  if((n=fwrite((char *)save_taylor,sizeof(revreal)*nitems,
          1,temp2_file))!=1) {
  fprintf(stderr,"ADOL-C error: fatal error-doing a write %d--- error %d\n",n,errno);
      switch (errno) {
      case 28: /* ENOSPC */
	fprintf(stderr,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	fprintf(stderr,"File to big-- Taylor-tape space exhausted.\n");
	break;
      default:
	fprintf(stderr,"Unexpected error %d .\n",errno);
	break;
      }
      exit(-1);
    }
  T_write_cnt++;
}

/*-----------------------------------------------------------------------*/
/* Static function T_prev_block                                          */
/* called by taylor_back, taylor_back2, get_taylor, get_taylors          */
/* Gets the next (previous block) of size nitems                         */
/*-----------------------------------------------------------------------*/
static int T_prev_block(int nitems)
{
  int n;
#ifdef DEBUG
printf("ADOL-C debug: prev %d =nitems %d T_write_cnt \n", nitems, T_write_cnt);
#endif
  if (T_file_access)
    {
      if (T_write_cnt == 0)
	{
	  return 0;
	}
      T_write_cnt--;
      fseek(temp2_file,T_buf_size*T_write_cnt*sizeof(revreal),0);
      n=fread((char *)save_taylor,sizeof(revreal),nitems,temp2_file);
      if(n != nitems)
	{ 
	 printf("ADOL-C error: Read error on taylor file n= %d\n",n);
	 return 0;
	 }
      taylor_cnt = nitems;
      return 1;
    }
  return 0;
}

/*---------------------------------------------------------------------*/
/* Write_taylor writes the block of size depth of taylor coefficients  */
/* from point loc to the taylor buffer.  If the buffer is filled, then */
/* it is written to the taylor tape (T_put_block).                     */
/*---------------------------------------------------------------------*/
void write_taylor(locint loc,int depth)
{
  int i;
  double* Tloc = T[loc];
  for (i=0;i<depth;i++)
    {
      save_taylor[taylor_cnt++]=Tloc[i];  /* In this assignment the */

      if (taylor_cnt == T_buf_size) /* precision will be sacrificed if the */
	{                           /* type revreal is defined as float.   */
          T_put_block(T_buf_size);
	}
    }
}

/*-------------------------------------------------------------------------*/
/* Write_scaylor writes a single element (x) to the taylor buffer.  If full*/
/* the buffer is written out.                                              */
/*-------------------------------------------------------------------------*/
void write_scaylor(revreal x)
{
  save_taylor[taylor_cnt++]= x;
  if(taylor_cnt == T_buf_size) T_put_block(T_buf_size);
}


/*-------------------------------------------------------------------------*/
/* Write_scaylors writes # size elements from x to the taylor buffer.      */
/* If full, the buffer is written out.                                     */
/*-------------------------------------------------------------------------*/
void write_scaylors(revreal *x, int size)
{ 
  int i;
  for(i=0; i<size; i++)
    {
      save_taylor[taylor_cnt++]= x[i];
      if(taylor_cnt == T_buf_size) T_put_block(T_buf_size);
    }
}

/* Close the taylor file, reset data. */

void taylor_close(int buffer,int dep,int ind)
{
  int n;

  numdep = dep;
  numind = ind;
  T_tail = taylor_cnt;
  T_length = T_buf_size*T_write_cnt+taylor_cnt;
  if (T_write_cnt)
     {
     if(T_tail>0 ) T_put_block(T_tail);
     free((char *)save_taylor);
     save_taylor =0;
     }
  T_blocks = T_write_cnt;
  if ((T_blocks) && (T_length*sizeof(revreal) <= buffer))
    { 
      save_taylor = (revreal *) malloc(T_length*sizeof(revreal));
      fseek(temp2_file,0,0);
      n=fread((char *)save_taylor,sizeof(revreal),T_length,temp2_file);
      if ( n != T_length)
	{
	  printf("ADOL-C error: read error in taylor_close n= %d\n",n);
	  exit(-2);
	}
      T_tail = T_length;
      T_blocks = 0;
    }
#ifdef DEBUG
  if(T_blocks)
    printf("\n ADOL-C debug: taylor file of length %d bytes completed\n", T_length*sizeof(revreal));
  else
    printf("\n ADOL-C debug: taylor array of length %d bytes completed\n", T_length*sizeof(revreal));
#endif
}

/* Set up statics for writing taylor data */
 
void taylor_begin(int buffer,double** Tg,int degree)
{ 
  T = Tg;
  if(save_taylor) free((char *)save_taylor);
  T_buf_size = 1+buffer/sizeof(revreal);
  save_taylor = (revreal *)malloc(sizeof(revreal)*T_buf_size);
  T_write_cnt = 0;
  T_length = 0;
  taylor_cnt = 0; 
  degsave = degree;
}


void taylor_back2(revreal** Trg,int* dep,int* ind,int* degree)
{
  *dep = numdep;
  *ind = numind;
  *degree = degsave;
  Tr = Trg;
  T_write_cnt = T_blocks;
  taylor_cnt = T_tail;
  if (T_blocks == 0 && save_taylor == 0 )
  { 
    printf("ADOL-C error: no temp file or array for reverse sweep \n");
    exit(-2);
   }
  if(T_blocks)
     {
       if(save_taylor) free((char*) save_taylor);
       save_taylor = (revreal*) malloc(T_buf_size*sizeof(revreal));
       if(T_prev_block(T_tail) == 0) printf("ADOL-C error: problems in back \n");
     }
}

void taylor_back(revreal* Trg,int* dep,int* ind,int* degree)
{
  *dep = numdep;
  *ind = numind;
  *degree = degsave;
  Trs = Trg;
  T_write_cnt = T_blocks;
  taylor_cnt = T_tail;
  if (T_blocks == 0 && save_taylor == 0 )
  { 
    printf("ADOL-C error: no temp file or array for reverse sweep \n");
    exit(-2);
   }
  if(T_blocks)
     {
       if(save_taylor) free((char*) save_taylor);
       save_taylor = (revreal*) malloc(T_buf_size*sizeof(revreal));
       if(T_prev_block(T_tail) == 0) printf("ADOL-C error: problems in back \n");
     }
}

void get_taylors(locint loc,int depth)
{
  int i;
  revreal* Trloc = Tr[loc];
  for (i = depth-1;i >= 0;i--)
    {
      if (taylor_cnt ==0)
	{
	  if (!T_prev_block(T_buf_size))
	    {
	      fprintf(stderr,"ADOL-C error: Fatal Error in get_taylor ");
	      exit(-1);
	    }
	}
      Trloc[i] = save_taylor[--taylor_cnt];
    }
}

void get_taylor(locint loc)
{
  if (taylor_cnt ==0)
    {
      if (!T_prev_block(T_buf_size))
	{
	  fprintf(stderr,"ADOL-C error: Fatal Error in get_taylor ");
	  exit(-1);
	}
    }
  Trs[loc] = save_taylor[--taylor_cnt];
}


#ifdef __cplusplus
}
#endif

