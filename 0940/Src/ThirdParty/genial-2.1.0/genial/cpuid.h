//GENIAL - GENeric Image & Array Library
//Copyright (C) 2007  Patrick LAURENT
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifndef CPUID_H
#define CPUID_H

#include <string.h>


struct cpuid_reg { unsigned int eax,ebx,ecx,edx; };

#if defined(__ICL) || defined(_MSC_VER)
#include <windows.h>
inline cpuid_reg _cpuid(unsigned int a, unsigned int c=0)
{
  cpuid_reg r;
  __try 
  {
    _asm  
    {
      mov eax, a
      mov ecx, c
      cpuid
      mov r.eax, eax
      mov r.ebx, ebx
      mov r.ecx, ecx
      mov r.edx, edx
    }
  } __except (EXCEPTION_EXECUTE_HANDLER) { r.eax=r.ebx=r.ecx=r.edx=0; }
  return r;
}
#elif defined(__GNUC__) && (defined(i386)||defined(__x86_64__))
inline cpuid_reg _cpuid(unsigned int a, unsigned int c=0)
{
  cpuid_reg r;
  __asm__("cpuid":"=a"(r.eax),"=b"(r.ebx),"=c"(r.ecx),"=d"(r.edx):"a"(a),"c"(c));
  return r;
}
#else
inline cpuid_reg _cpuid(unsigned int a, unsigned int c=0)
{
  cpuid_reg r={0,0,0,0};
  return r;
}
#endif

inline cpuid_reg cpuid(unsigned int a, unsigned int c=0)
{
  if (a<=_cpuid(0).eax) 
    return _cpuid(a,c);
  cpuid_reg r={0,0,0,0};
  return r;
}

inline cpuid_reg cpuid_ext(unsigned int a, unsigned int c=0)
{
  a|=0x80000000;
  if (a<=_cpuid(0x80000000).eax) 
    return _cpuid(a,c);
  cpuid_reg r={0,0,0,0};
  return r;
}

static char _cpu_vendor[13];
static char _cpu_name[48];
inline char *cpu_vendor() { unsigned int *p=(unsigned int *)_cpu_vendor; cpuid_reg reg=cpuid(0); p[0]=reg.ebx; p[1]=reg.edx; p[2]=reg.ecx; _cpu_vendor[12]='\0'; return _cpu_vendor; }
inline char *cpu_name  () { unsigned int *p=(unsigned int *)_cpu_name; cpuid_reg reg; reg=cpuid_ext(2); p[0]=reg.eax; p[1]=reg.ebx; p[2]=reg.ecx;  p[3]=reg.edx; reg=cpuid_ext(3); p[4]=reg.eax; p[5]=reg.ebx; p[6]=reg.ecx;  p[7]=reg.edx; reg=cpuid_ext(4); p[8]=reg.eax; p[9]=reg.ebx; p[10]=reg.ecx;  p[11]=reg.edx; return _cpu_name; }
inline bool is_intel_cpu() { return strcmp(cpu_vendor(),"GenuineIntel")==0; }
inline bool is_amd_cpu  () { return strcmp(cpu_vendor(),"AuthenticAMD")==0; }
inline bool is_mmx_cpu  () { return (cpuid(1).edx&0x00800000)!=0; }
inline bool is_sse_cpu  () { return (cpuid(1).edx&0x02000000)!=0; }
inline bool is_sse2_cpu () { return (cpuid(1).edx&0x04000000)!=0; }
inline bool is_sse3_cpu () { return (cpuid(1).ecx&0x00000001)!=0; }

inline bool is_3dnow_cpu() { return (cpuid_ext(1).edx&0x80000000)!=0; }

inline unsigned int num_l1_threads      () { return ((cpuid(4,1).eax>>14)&0xFFF)+1; }
inline unsigned int num_l2_threads      () { return ((cpuid(4,2).eax>>14)&0xFFF)+1; }
inline unsigned int l1_cache_size       () { cpuid_reg reg=cpuid(4,1); return (((reg.ebx>>22)&0x3FF)+1)*(((reg.ebx>>12)&0x3FF)+1)*((reg.ebx&0xFFF)+1)*(reg.ecx+1); }
inline unsigned int l2_cache_size       () { cpuid_reg reg=cpuid(4,2); return (((reg.ebx>>22)&0x3FF)+1)*(((reg.ebx>>12)&0x3FF)+1)*((reg.ebx&0xFFF)+1)*(reg.ecx+1); }
inline unsigned int num_logical_cores   () { cpuid_reg reg=cpuid(1  ); return (reg.edx&0x10000000)?(reg.ebx>>16)&0xFF:1; }
inline unsigned int num_physical_cores  () { return (!is_amd_cpu()) ? ((cpuid(4,0).eax>>26)&0x03F)+1 : num_logical_cores(); }
//inline unsigned int num_physical_cores  () { return ((cpuid(4,0).eax>>26)&0x03F)+1; }

#if (defined(__ICL) || defined(_MSC_VER)) && defined(WIN32)
  #include <windows.h>
  inline int num_processors()
  {
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors;
  }
#elif (defined(__MACOSX__) || defined(__APPLE__))
  #include <Multiprocessing.h>
  inline int num_processors() { return MPProcessorsScheduled(); }
#else
  inline int num_processors() { return num_logical_cores(); }
#endif

inline int num_physical_processors() { return num_processors()*num_physical_cores()/num_logical_cores(); }


#endif


//#include "cpuid.h"
//#include <iostream>
//using namespace std;
//int main()
//{
//  cout << "name  \t" << cpu_name    () << endl;
//  cout << "vendor\t" << cpu_vendor  () << endl;
//  cout << "intel \t" << is_intel_cpu() << endl;
//  cout << "amd   \t" << is_amd_cpu  () << endl;
//  cout << endl;
//  cout << "mmx   \t" << is_mmx_cpu  () << endl;
//  cout << "sse   \t" << is_sse_cpu  () << endl;
//  cout << "sse2  \t" << is_sse2_cpu () << endl;
//  cout << "sse3  \t" << is_sse3_cpu () << endl;
//  cout << "3dnow \t" << is_3dnow_cpu() << endl;
//  cout << endl;
//  cout << "L1 cache\t" << l1_cache_size() << endl;
//  cout << "L2 cache\t" << l2_cache_size() << endl;
//  cout << endl;
//  cout << "processors    \t" << num_processors    () << endl;
//  cout << "logical  cores\t" << num_logical_cores () << endl;
//  cout << "physical cores\t" << num_physical_cores() << endl;
//}

