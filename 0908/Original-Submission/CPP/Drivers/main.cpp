// Author: Yong-Kang Zhu (yongkang.zhu@gmail.com)
// Code can be used only for academic purpose

// This file provides a driver which tests whether the summation
// algorithms are working correctly.

// How to run?
// Suppose after a successful compilation, the executable is "main".
// We can run it with three parameters: "main $1 $2 $3", where $1 is
// the number of summands, $2 is the exponent difference, and $3 is
// the data type (1, well-conditioned; 2, random; 3, Anderson's; and
// 4, exact sum equals zero). More details about $2 and $3 can be
// found in our paper. If it is running correctly, all the results
// should be same. Furthermore, the results are all zero when $3 is
// set to 4. Please check the compiler options in ExactSum.h.
//
// If this program does not work correctly on your machine, please
// email me the running arguments, the output, the error message if
// any, and your system information (CPU, Operating System, etc).

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ExactSum.h"

int pflag;  //program flag (Data type 1 to 4)
// random flop generator with exponent differece expo
int rflag=0;
double ret_saved;
double Rand(int expo)
{
	double ret;
	unsigned int low=0,high=0;
	int i;
	for(i=0;i<32;i++)
	{
		if(rand() & 0x1)
			low=(low>>1) | 0x80000000;
		else
			low=low>>1;
	}
	for(i=0;i<32;i++)
	{
		if(rand() & 0x1)
			high=(high>>1) | 0x80000000;
		else
			high=high>>1;
	}
	*(int*)(&ret)=low;
	*(((int*)&ret)+1)=high;
	if(expo!=0)
		((str_double*)&ret)->exponent=((rand()%(expo))-expo/2)+0x3ff;
	else 
		((str_double*)&ret)->exponent=0x3ff;
	if(pflag==2 || pflag==3)
		return ret;  //pflag=2, random data
	if(pflag==1)
		rflag=0;     //pflag=1, well-conditioned data
	if(rflag==0)  // exact zero
	{
		ret_saved=ret;
		((str_double*)(&ret))->sign=1;
		rflag=1;
	}
	else
	{
		ret=ret_saved;
		((str_double*)(&ret))->sign=0;
		rflag=0;
	}
	return ret;
}

int main(int argn, char* argc[])
{
	int deltaExp;
	int MAXNUM;

	if(argn!=4)  // n, exp, pflag 
	{
		printf("wrong arguments\n");
		return 0;
	}

	ExactSum mysum;

	MAXNUM=atoi(argc[1]);
	deltaExp=atoi(argc[2]);
	pflag=atoi(argc[3]);
	// pflag
	// 1 well-conditioned
	// 2 random
	// 3 anderson's
	// 4 sum = zero

	if(deltaExp<0)
	{
		printf("wrong exp\n");
		return 0;
	}	

	if(pflag>4 || pflag<1)
	{
		printf("wrong pflag\n");
		return 0;
	}

	if(pflag == 4 && MAXNUM % 2 == 1)
	{
		printf("N should be even for Data 4\n");
		return 0;
	}

	srand((unsigned)time(NULL));
	int i;

	double *num_list;
	double *original_list;

	double st=0;
	
	double result_ORS=-0.1, result_iFastSum=-0.1, result_OnlineExactSum=-0.1,
		result_OnlineExactSum_AddNumber=-0.1,
		result_OnlineExactSum_AddArray=-0.1;

	original_list=new double[MAXNUM+1];
	num_list=new double[MAXNUM+1];

	for(i=1;i<=MAXNUM;i++)
	{
		original_list[i]=(Rand(deltaExp));
		st+=original_list[i];
	}

	// Anderson's ill-conditioned data
	if(pflag==3) 
		for(i=1;i<=MAXNUM;i++)
			original_list[i]-=st/MAXNUM;

	//randomly change the order
	for(i=1;i<=MAXNUM*2;i++)
	{
		int x=rand()%MAXNUM+1;
		int y=rand()%MAXNUM+1;
		while(x==y)
			y=rand()%MAXNUM+1;
		double temp=original_list[x];
		original_list[x]=original_list[y];
		original_list[y]=temp;
	}

	// iFastSum
	for(i=1;i<=MAXNUM;i++)
			num_list[i]=original_list[i];
	result_iFastSum=mysum.iFastSum(num_list,MAXNUM);
	printf("iFastSum is done\n");

	// OnlineExactSum
	for(i=1;i<=MAXNUM;i++)
			num_list[i]=original_list[i];
	result_OnlineExactSum=mysum.OnlineExactSum(num_list,MAXNUM);
	printf("OnlineExactSum is done\n");

	for(i=1;i<=MAXNUM;i++)
			num_list[i]=original_list[i];
	mysum.Reset();
	for(i=1;i<=MAXNUM;i++)
		mysum.AddNumber(num_list[i]);
	result_OnlineExactSum_AddNumber=mysum.GetSum();
	printf("OnlineExactSum:AddNumber is done\n");

	for(i=1;i<=MAXNUM;i++)
			num_list[i]=original_list[i];
	mysum.Reset();
	if(MAXNUM > 1200)
	{
		mysum.AddArray(num_list,10);
		mysum.AddArray(num_list+10,100);
		mysum.AddArray(num_list+110,1000);
		mysum.AddArray(num_list+1110,MAXNUM-1110);
	}
	else
		mysum.AddArray(num_list,MAXNUM);
	result_OnlineExactSum_AddArray=mysum.GetSum();
	printf("OnlineExactSum:AddArray is done\n");

	printf("iFastSum                   = %.20e\n",result_iFastSum);
	printf("OnlineExactSum (direct)    = %.20e\n",result_OnlineExactSum);
	printf("OnlineExactSum (AddNumber) = %.20e\n",result_OnlineExactSum_AddNumber);
	printf("OnlineExactSum (AddArray)  = %.20e\n",result_OnlineExactSum_AddArray);

	delete[]num_list;
	delete[]original_list;
	return 0;
} 
