/*
   --------------------------------------------------------------
   File taputil2.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   Provides the definition of the tape - interface functions,etc.
   The majority of these functions are called by the forward (hos_forward)
   and reverse (fos, fov, hos, hov -reverse) functions and the functions
   in taputil1.c.
*/

#ifdef __cplusplus
extern "C" {
#endif

#include "dvlparms.h" /* Developers Parameters */

/* Necessary Includes */

#include "usrparms.h"
#include "oplate.h"
#include "taputil3.h"

/* Buffers for the operation tape, location tape, real tape. */

unsigned char *op_codes;
locint *loc_tape;
double *real_tape;

/* File pointers to the operation tape, location tape, real tape */ 
static FILE *op_file_out;
static FILE *int_file_out;
static FILE *val_file_out;

/* Stats on operation tape, location tape, real tape */
static int op_access_ptr,int_access_ptr,val_access_ptr;
static int op_len_ptr,int_len_ptr,val_len_ptr;

/* Pointers into the operation tape, location tape, real tape */
int op_ptr;
int loc_ptr;
int real_ptr;

unsigned char *g_op_ptr;
locint *g_loc_ptr;
double *g_real_ptr;

/* Strings for the tape names (actual file names) */
static char *op_file,*int_file,*val_file;

/* Current buffer size */
static int buff_size;



void set_buf_size(int size)
{
  buff_size = size;
}

void set_buffers(char *file1, unsigned char *op_addr,
		 char *file2, locint *int_addr,
		 char *file3, double *real_addr)
{
  op_codes = op_addr;
  loc_tape = int_addr;
  real_tape = real_addr;
  op_file = file1;
  int_file = file2;
  val_file = file3;
  op_ptr=loc_ptr=real_ptr = 0;
  op_access_ptr = int_access_ptr = val_access_ptr = 0;
  op_len_ptr=int_len_ptr=val_len_ptr = 0;
}
/****************************************************************/
/** Put_Block puts a block of tape to the disk.  I assume this **/
/** is called only during a first forward pass or during the   **/
/** the taping itself. Its purpose is to record all of the     **/
/** computations.                                              **/
/****************************************************************/
void put_op_block(int buffer_size)
{
  int n;
  if (op_access_ptr == 0)
    {

    op_file_out = fopen(op_file,"r");
    if (op_file_out != 0)
      {
#ifdef DEBUG
	printf("ADOL-C debug: old tapefile %s exists and deleted\n",op_file);
#endif
	fclose(op_file_out);
        if(remove(op_file))  /*  Complies with ANSI C standard */
   /*   if(unlink(op_file))      works on some UNIX systems */
	  printf("ADOL-C error: unable to remove old tapefile\n");
	op_file_out = fopen(op_file,"w");
      }
    else
      {
	op_file_out = fopen(op_file,"w");
	errno =0; /* Clear Out the Error */
      }
    op_access_ptr = 1;
    }

  op_len_ptr += buffer_size;

  if ((n=fwrite((char *)op_codes,buffer_size,1,op_file_out))!=1)
    {
      printf("ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
      fprintf(stderr,"ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
      switch (errno) {
      case 28: /* ENOSPC */
	printf("No space left on device-contact sys. manager\n");
	fprintf(stderr,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	printf("File too big-- tape space exhausted.\n");
	fprintf(stderr,"File too big-- tape space exhausted.\n");
	break;
      default:
	printf("Unexpected unix file error-- %d.\n",errno);
	fprintf(stderr,"Unexpected unix file error-- %d.\n",errno);
	break;
      }
      exit(-3);
    }
  op_ptr=0;
  errno=0;
}
void put_locint_block(int buffer_size)
{
  int n;
  if (int_access_ptr == 0)
    {
    int_file_out = fopen(int_file,"r");
    if (int_file_out != 0)
      {
#ifdef DEBUG
	printf("ADOL-C debug: old tapefile %s exists and deleted\n",int_file);
#endif
	fclose(int_file_out);
        if(remove(int_file))  /*    Complies with ANSI C standard */
   /*   if(unlink(int_file))        works on some UNIX systems    */
	  printf("ADOL-C error: unable to remove old tapefile\n");
	int_file_out = fopen(int_file,"w");
      }
    else
      {
	int_file_out = fopen(int_file,"w");
	errno =0; /* Clear Out the Error */
      }
    int_access_ptr = 1;
    }

  int_len_ptr += buffer_size;

  if ((n=fwrite((locint *)loc_tape,buffer_size*sizeof(locint),1,int_file_out))!=1)
    {
      printf("ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
      fprintf(stderr,"ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
      switch (errno) {
      case 28: /* ENOSPC */
	printf("No space left on device-contact sys. manager\n");
	fprintf(stderr,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	printf("File too big-- tape space exhausted.\n");
	fprintf(stderr,"File too big-- tape space exhausted.\n");
	break;
      default:
	printf("Unexpected unix file error-- %d.\n",errno);
	fprintf(stderr,"Unexpected unix file error-- %d.\n",errno);
	break;
      }
      exit(-3);
    }
  loc_ptr=0;
  errno=0;
}
void put_val_block(int buffer_size)
{
  int n;
  if (val_access_ptr == 0)
    {

    val_file_out = fopen(val_file,"r");
    if (val_file_out != 0)
      {
#ifdef DEBUG
	printf("ADOL-C debug: old tapefile %s exists and deleted\n",val_file);
#endif
	fclose(val_file_out);
        if(remove(val_file))   /* Complies with ANSI C standard */
     /*	if(unlink(val_file))      works on some UNIX systems    */
	  printf("ADOL-C error: unable to remove old tapefile\n");
	val_file_out = fopen(val_file,"w");
      }
    else
      {
	val_file_out = fopen(val_file,"w");
	errno =0; /* Clear Out the Error */
      }
    val_access_ptr = 1;
    }

  val_len_ptr += buffer_size;

  if ((n=fwrite((double *)real_tape,buffer_size*sizeof(double),1,val_file_out))!=1)
    {
      printf("ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
      fprintf(stderr,"ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
      switch (errno) {
      case 28: /* ENOSPC */
	printf("No space left on device-contact sys. manager\n");
	fprintf(stderr,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	printf("File too big-- tape space exhausted.\n");
	fprintf(stderr,"File too big-- tape space exhausted.\n");
	break;
      default:
	printf("Unexpected unix file error-- %d.\n",errno);
	fprintf(stderr,"Unexpected unix file error-- %d.\n",errno);
	break;
      }
      exit(-3);
    }
  real_ptr=0;
  errno=0;
}

void close_tape(int *stats, int flag)
{
  int i;
  int access = (flag || op_access_ptr || int_access_ptr ||val_access_ptr);
  if (access)
    {
      if (op_ptr!=0)
	put_op_block(op_ptr);
      fclose(op_file_out);
    }
  else op_len_ptr=op_ptr;

  stats[5]=op_len_ptr;
  stats[6]=op_access_ptr;

  if (access)
    {
     if (real_ptr!=0)
       put_val_block(real_ptr);
      fclose(val_file_out);
    }
  else
    val_len_ptr = real_ptr;
  stats[9]=val_len_ptr;
  stats[10]=val_access_ptr;
  if (access)
    {
     if (loc_ptr!=0)
       put_locint_block(loc_ptr);
     stats[7]=int_len_ptr;
     stats[8]=int_access_ptr;
     fseek(int_file_out,0,0);
     fwrite(stats,11*sizeof(int),1,int_file_out);
     fclose(int_file_out);
   }
  else{
    int_len_ptr=loc_ptr;
    stats[7]=int_len_ptr;
    stats[8]=int_access_ptr;
    for(i=0;i<11;i++)
      {
	loc_tape[i]=stats[i];
      }
  }
}

static long op_file_cnt,int_file_cnt,val_file_cnt;

void end_sweep()
{
  if (op_access_ptr)
    { 
      fclose(op_file_out);
    } 
  if (int_access_ptr)
    { 
      fclose(int_file_out);
    }
  if (val_access_ptr)
    { 
      fclose(val_file_out);
    }

}


void init_rev_sweep(int tag)
{
  get_op_stats(tag,&op_file,&op_len_ptr,&op_access_ptr,&op_codes);
  if (op_access_ptr)
    { 
      op_file_out=fopen(op_file,"r");
      op_ptr=op_len_ptr % buff_size;
      fseek(op_file_out,0,2);
      op_file_cnt =  ftell(op_file_out);
      op_file_cnt -=op_ptr*sizeof(unsigned char);
      fseek(op_file_out,op_file_cnt,0);
      fread((char *)op_codes,op_ptr,1,op_file_out);
      op_file_cnt -= buff_size*sizeof(unsigned char);
      g_op_ptr=op_codes+op_ptr;
    } 
  else g_op_ptr = op_codes + op_len_ptr;
  get_int_stats(tag,&int_file,&int_len_ptr,&int_access_ptr,&loc_tape);
  if (int_access_ptr)
    { 
      int_file_out=fopen(int_file,"r");
      loc_ptr=int_len_ptr % buff_size;
      fseek(int_file_out,0,2);
      int_file_cnt =  ftell(int_file_out);
      int_file_cnt -=loc_ptr*sizeof(locint);
      fseek(int_file_out,int_file_cnt,0);
      fread((char *)loc_tape,loc_ptr*sizeof(locint),1,int_file_out);
      int_file_cnt -= buff_size*sizeof(locint);
      g_loc_ptr=loc_tape+loc_ptr;
    } 
  else g_loc_ptr= loc_tape+int_len_ptr;
  get_val_stats(tag,&val_file,&val_len_ptr,&val_access_ptr,&real_tape);
  if (val_access_ptr)
    { 
      val_file_out=fopen(val_file,"r");
      real_ptr= val_len_ptr % buff_size;
      fseek(val_file_out,0,2);
      val_file_cnt =  ftell(val_file_out);
      val_file_cnt -=real_ptr*sizeof(double);
      fseek(val_file_out,val_file_cnt,0);
      fread((char *)real_tape,real_ptr*sizeof(double),1,val_file_out);
      val_file_cnt -= buff_size*sizeof(double);
      g_real_ptr= real_tape+real_ptr;
    } 
  else g_real_ptr=real_tape+val_len_ptr;
}


#define min(a,b) ( (a)>(b)? (b):(a) )

void init_for_sweep(int tag)
{
  get_op_stats(tag,&op_file,&op_len_ptr,&op_access_ptr,&op_codes);
  if (op_access_ptr)
    { 
      op_file_out=fopen(op_file,"r");
      op_ptr=min(buff_size,op_len_ptr); 
      fread((char *)op_codes,op_ptr,1,op_file_out);
      op_len_ptr-=op_ptr;
    } 
  g_op_ptr=op_codes;
  get_int_stats(tag,&int_file,&int_len_ptr,&int_access_ptr,&loc_tape);
  if (int_access_ptr)
    { 
      int_file_out=fopen(int_file,"r");
      loc_ptr=min(buff_size,int_len_ptr);
      fread((locint *)loc_tape,sizeof(locint),loc_ptr,int_file_out);
      int_len_ptr-=loc_ptr;
    } 
  g_loc_ptr=loc_tape+22; 
  /* loc_tape = (loc_tape+loc_ptr); */
  get_val_stats(tag,&val_file,&val_len_ptr,&val_access_ptr,&real_tape);
  if (val_access_ptr)
    { 
      val_file_out=fopen(val_file,"r");
      real_ptr= min(val_len_ptr,buff_size);
      fread((char *)real_tape,real_ptr*sizeof(double),1,val_file_out);
      val_len_ptr-=real_ptr;
    } 
  g_real_ptr=real_tape;
}
void get_op_block_f()
{
  op_ptr=min(buff_size,op_len_ptr); 
  fread((char *)op_codes,op_ptr,1,op_file_out);
  op_len_ptr-=op_ptr;
  g_op_ptr=op_codes;
}

void get_loc_block_f()
{
  loc_ptr=min(buff_size,int_len_ptr);
  fread((char *)loc_tape,loc_ptr*sizeof(locint),1,int_file_out);
  int_len_ptr-=loc_ptr;
  g_loc_ptr=loc_tape;
}

void get_val_block_f()
{
  real_ptr= min(val_len_ptr,buff_size);
  fread((char *)real_tape,real_ptr*sizeof(double),1,val_file_out);
  val_len_ptr-=real_ptr;
  g_real_ptr=real_tape;
  g_loc_ptr++; /* get_locint_f(); value used in reverse only */
}
    
void get_op_block_r()
{
  fseek(op_file_out,op_file_cnt,0);
  fread((char *)op_codes,buff_size,1,op_file_out);
  op_file_cnt -= buff_size*sizeof(unsigned char);
  g_op_ptr=op_codes+buff_size;
}

void get_loc_block_r()
{
  fseek(int_file_out,int_file_cnt,0);
  fread((char *)loc_tape,buff_size*sizeof(locint),1,int_file_out);
  int_file_cnt -= buff_size*sizeof(locint);
  g_loc_ptr=loc_tape+buff_size-loc_tape[buff_size-1];
}
void get_val_block_r()
{
  locint temp;
  fseek(val_file_out,val_file_cnt,0);
  fread((char *)real_tape,buff_size*sizeof(double),1,val_file_out);
  val_file_cnt -= buff_size*sizeof(double);
  temp=*(--g_loc_ptr);   /*get_locint_r();*/ 
  g_real_ptr=real_tape+buff_size-temp;
}

void put_to_op(unsigned char op)
{
  if (op_ptr == buff_size-1){
    op_codes[op_ptr]=end_of_op;
    put_op_block(buff_size);
    op_codes[op_ptr++]=end_of_op;
  }
  op_codes[op_ptr++]=op;
}
void put_locint(locint);

void put_op(unsigned char op)
{
  if (loc_ptr > buff_size-5) 
    {
      loc_tape[buff_size-1]=buff_size-loc_ptr;
      put_locint_block(buff_size);
      put_to_op(end_of_int);
    }
  if (real_ptr > buff_size-5) 
    {
      put_locint(buff_size-real_ptr);
      put_val_block(buff_size);
      put_to_op(end_of_val);
    }
  put_to_op(op);
}




void put_val(double r_val)
{
  /* if (real_ptr == buff_size) put_val_block(buff_size); */
  real_tape[real_ptr++]=r_val;
}

void put_vals_p(double *r_val,int size)
{
  int j;
  for (j=0;j<size;j++)
    {
      real_tape[real_ptr++]=r_val[j];
    }
  put_locint(buff_size-real_ptr);
  put_val_block(buff_size);
  put_to_op(end_of_val);
}

void put_vals_r(double *r_val,int size)
{
  int j;
  for (j=0;j<size;j++)
    {
      real_tape[real_ptr++]=r_val[j];
    }
}

int get_val_space()
{
  if ((buff_size-real_ptr-5)<0) {
    put_locint(buff_size-real_ptr);
    put_val_block(buff_size);
    put_to_op(end_of_val);
  } /* endif */
  return (buff_size-real_ptr-5);
}



double *get_val_v_r(locint size)
{
  g_real_ptr -= size;
  return(g_real_ptr);
}
double *get_val_v_f(locint size)
{
  double *temp = g_real_ptr;
  g_real_ptr += size;
  return (temp);
}
void reset_val_r()
{
  if (g_real_ptr == real_tape)
    {
      get_val_block_r();
    }
}

void put_locint(locint loc1)
{
  /*if (loc_ptr == buff_size) put_locint_block(buff_size); */
  loc_tape[loc_ptr++]=loc1;
}

#ifdef DEBUG
unsigned char get_op_f()
{
  unsigned char temp;
  temp= *g_op_ptr++;
  printf("f_op: %i\n",temp-'\0');
  return temp;
}
#endif

#ifdef __cplusplus
}
#endif
