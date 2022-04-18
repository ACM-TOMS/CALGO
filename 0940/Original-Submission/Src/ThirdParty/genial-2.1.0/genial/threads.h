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


#ifndef THREADS_H
#define THREADS_H


#if defined(OPENMP)
#include <omp.h>
#endif

#if !defined(PTHREADS) && !defined(WINTHREADS) && !defined(MPTASKS)
  #if defined(__ICL) || defined(_MSC_VER)
    #define WINTHREADS
  #elif defined(__MACOSX__) || defined(__APPLE__)
    #define MPTASKS
  #else
    #define PTHREADS
  #endif
#endif

#if defined(PTHREADS)
  #include <pthread.h>
  #include <semaphore.h>
#elif defined(WINTHREADS)
  #ifdef _WIN32_WINNT
    #undef _WIN32_WINNT
  #endif
  #define _WIN32_WINNT 0x0403
  #include <windows.h>
  #include <winbase.h>
  #include <process.h>
  #undef min
  #undef max
#elif defined(MPTASKS)
  #include <Multiprocessing.h>
#endif

#include <functional>
#include <exception>
#include <cassert>
#include <list>
#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>

#include "memory.h"

template<class G> class Vector;
template<class G> class Matrix;

//group=Multithreading Functions

#ifdef HIDE_FROM_DOCJET
#else
namespace gmt
{
#endif

using namespace std;



#if (defined(__ICL) || defined(_MSC_VER)) && defined(WIN32)
template<class T> inline T gmt_min(const T &x, const T &y) { return __min(x,y); }
#else
template<class T> inline T gmt_min(const T &x, const T &y) { return min(x,y); }
#endif

#ifndef MAX_THREADS
#define MAX_THREADS 0
#endif

#if defined(OPENMP)
//{unsecret}
//Summary: Returns the maximum number of threads simultaneously executed in algorihms.
inline int max_threads() { return omp_get_max_threads(); }
//{unsecret}
//Summary: Gets or sets the number of threads simultaneously executed in algorihms.
inline int  num_threads() { return omp_get_num_threads(); }
//{unsecret}
inline void num_threads(int n) { return omp_set_num_threads(n); }
#else
#if MAX_THREADS==0
#include "cpuid.h"
//static int numthreads = num_processors();
//static int maxthreads = num_processors();
static int numthreads = num_physical_cores();
static int maxthreads = num_physical_cores();
inline int max_threads() { return maxthreads; }
#else
static int numthreads = MAX_THREADS;
inline int max_threads() { return MAX_THREADS; }
#endif
inline int num_threads() { return numthreads; }
inline void num_threads(int n) { numthreads=n; }
#endif

class thread_error : public std::exception
{
  private:
    int err;

  public:
    inline thread_error() : err(0) {}
    inline thread_error(int e) : err(e) {}

    inline ~thread_error() throw() {}
    inline int error() const { return err; }
    virtual const char* what() const throw() { return "thread error"; }
};

//{unsecret}
//{group:Multithreading Classes}
//Summary: Scoped lock
//Remarks: Acquires a lock on a mutex and automatically releases it on destruction
//See: ^mutex^, ^spin_mutex^, ^condition_mutex^
template<class M>
class scoped_lock
{
  public:
    typedef M mutex_type;

  private:
    mutex_type &mut;

  public:
    inline scoped_lock(mutex_type &m) : mut(m) { mut.acquire(); }
    inline ~scoped_lock() { mut.release(); }
};


#if defined(WINTHREADS)

inline void assert_thread(HANDLE h) { assert(h!=0); }
inline void assert_thread(BOOL   b) { assert(b   ); }
inline void check_thread (HANDLE h) { if (h!=0) return; throw thread_error(GetLastError()); }
inline void check_thread (BOOL   b) { if (b   ) return; throw thread_error(GetLastError()); }


//{unsecret}
//{group:Multithreading Classes}
//Summary: Protects shared data structures from concurrent modifications
//Remarks: Prefer the use of a scoped_lock rather than explicitly calling the acquire and release methods.
//Example:
//  int var;
//  typedef gmt::mutex mutex_type;
//  mutex_type mut;
//
//  inline set_var(int x)
//  {
//    mutex_type::scoped_lock lock(mut);
//    var=x;
//  };
//See: ^spin_mutex^, ^condition_mutex^, ^scoped_lock^
class mutex
{
  public:
    typedef mutex self;
    typedef scoped_lock<self> scoped_lock;

  private:
    HANDLE h;

  public:
    inline mutex() { check_thread(h=CreateMutex(NULL,FALSE,NULL)); }
    inline ~mutex() { assert_thread(CloseHandle(h)); }

    //Summary: Waits until a lock can be acquired
    inline void acquire()     { assert_thread(WaitForSingleObject(h,INFINITE)==WAIT_OBJECT_0); }
    //Summary: Returns true if it can immediately acquire a lock, else returns false
    inline bool try_acquire() { return WaitForSingleObject(h,0)==WAIT_OBJECT_0; }
    //Summary: Releases the lock
    inline void release()     { assert_thread(ReleaseMutex(h)); }

  private:
    mutex(const self &);
    const self &operator=(const self &);
};

//{unsecret}
//{group:Multithreading Classes}
//Summary: Mutex spinning in user space
//Remarks: For short waits, spining in user space is faster than sleeping.
//See: ^mutex^, ^condition_mutex^, ^scoped_lock^
class spin_mutex
{
  public:
    typedef spin_mutex self;
    typedef scoped_lock<self> scoped_lock;

  private:
    CRITICAL_SECTION mut;

  public:
    inline spin_mutex() { InitializeCriticalSection(&mut); }
    inline ~spin_mutex() { DeleteCriticalSection(&mut); }
    //inline spin_mutex() { check_thread(InitializeCriticalSectionAndSpinCount(&mut,0x80000400)); }

    //Summary: Waits until a lock can be acquired
    inline void acquire() { EnterCriticalSection(&mut); }
    //Summary: Releases the lock
    inline void release() { LeaveCriticalSection(&mut); }

  private:
    spin_mutex(const self &);
    const self &operator=(const self &);
};

//{unsecret}
//{group:Multithreading Classes}
//Summary: Suspends execution and relinquishes the processors until some predicate on shared data is satisfied
//Remarks: You would in most cases prefer to use a condition_mutex
//See: ^condition_mutex^
class condition
{
  public:
    typedef condition self;

  private:
    HANDLE h;

    condition(const self &);
    const self &operator=(const self &);

  public:
    inline condition() { check_thread(h=CreateEvent(NULL,false,false,NULL)); } //auto-reset, nonsignaled
    inline ~condition() {  assert_thread(CloseHandle(h)); }

    //Summary: Restarts one of the threads that are waiting on the condition variable
    inline void signal() { assert_thread(SetEvent  (h)); }
    //Summary: Releases the mutex, and waits for the condition variable to be signaled
    template<typename Mutex> inline void wait(Mutex &mut) { mut.release(); assert_thread(WaitForSingleObject(h,INFINITE)==WAIT_OBJECT_0); mut.acquire(); }
    template<typename Mutex> inline bool wait(Mutex &mut, int dt) { mut.release(); DWORD w=WaitForSingleObject(h,dt); assert_thread(w==WAIT_OBJECT_0 || w==WAIT_TIMEOUT); mut.acquire(); return w!=WAIT_TIMEOUT; }
};

//{unsecret}
//{group:Multithreading Classes}
//Summary: Gathers a spin_mutex and a condition together
//Example:
//  int x,y;
//  typedef gmt::condition_mutex mutex_type;
//  mutex_type cond;
//
//  //First thread
//  {
//    mutex_type::scoped_lock lock(cond);
//    while (x<=y) cond.wait();
//    // operate on x and y
//  }
//
//  // Second thread
//  {
//    mutex_type::scoped_lock lock(cond);
//    // modify x and y
//    if (x>y) cond.signal();
//  }
//See: ^condition^, ^mutex^, ^spin_mutex^
class condition_mutex
{
  public:
    typedef condition_mutex self;
    typedef scoped_lock<self> scoped_lock;

    typedef spin_mutex mutex_type;
    typedef condition condition_type;

  private:
    mutex_type mut;
    condition_type cond;

  public:
    //Summary: Waits until a lock can be acquired
    inline void acquire() { mut.acquire(); }
    //Summary: Releases the lock
    inline void release() { mut.release(); }

    //Summary: Restarts one of the threads that are waiting on the condition variable
    inline void signal() { cond.signal(); }
    //Summary: Waits for the condition variable to be signaled
    inline void wait() { cond.wait(mut); }
    inline bool wait(int dt) { return cond.wait(mut,dt); }
};

//{unsecret}
//{group:Multithreading Classes}
//Summary: Counter for resources shared between threads
class semaphore
{
  public:
    typedef semaphore self;
    typedef scoped_lock<self> scoped_lock;

  private:
    HANDLE h;

    semaphore(const self &);
    const self &operator=(const self &);

  public:
    //Remarks: The count associated with the semaphore is set initially to n, or zero if not specified
    inline semaphore()      { assert_thread(h=CreateSemaphore(NULL,0,std::numeric_limits<long>::max(),NULL)); };
    inline semaphore(int n) { assert_thread(h=CreateSemaphore(NULL,n,std::numeric_limits<long>::max(),NULL)); };
    inline ~semaphore() { assert_thread(CloseHandle(h)); }

    //Summary: Either releases one or multiple threads if there are any waiting, or increments the count if not enough were waiting
    inline void release()      { assert_thread(ReleaseSemaphore(h,1,NULL)); }
    inline void release(int n) { assert_thread(ReleaseSemaphore(h,n,NULL)); }
    //Summary: Waits until it can decrement the count
    inline void acquire()      { assert_thread(WaitForSingleObject(h,INFINITE)==WAIT_OBJECT_0); }
    //Summary: Returns the number with which it could immediately decrement the count
    inline int try_acquire()      { return WaitForSingleObject(h,0)==WAIT_OBJECT_0; }
    inline int try_acquire(int n) { for (int i=0; i<n; ++i) if (!try_acquire()) return i; return n; }
};

template<class F>
unsigned int __stdcall thread_proxy(void *f)
{
  (*static_cast<F *>(f))();
  return 0;
}

//{unsecret}
//{group:Multithreading Classes}
//Summary: Concurrent execution
//Example:
//  struct func0
//  {
//    void operator() const { ... }
//  };
//
//  gmt::thread<func0> th0(func0());
//See: ^condition_mutex^
template<class F>
class thread
{
  public:
    typedef thread self;
    typedef F function_type;

  private:
    function_type func;
    HANDLE th;

  public:
    //Summary: Creates a new thread that executes concurrently with the calling thread.
    //Remarks: The new thread applies the function object f.
    inline thread() : func() { init(); }
    template<class A> inline explicit thread(      A &a) : func(a) { init(); }
    template<class A> inline explicit thread(const A &a) : func(a) { init(); }
    template<class A,class B> inline explicit thread(      A &a, const B &b) : func(a,b) { init(); }
    template<class A,class B> inline explicit thread(const A &a, const B &b) : func(a,b) { init(); }

    inline thread(const self &x) : func(x.func) { init(); }
    inline self &operator=(const self &x) { func=x.func; init(); }

    //Summary: Waits until the executed function returns, and then closes the thread.
    inline ~thread() { join(); CloseHandle(th); th=0; }

    //Summary: Waits until the executed function returns.
    inline void join() { WaitForSingleObject(th,INFINITE); }

  private:
    inline void init() { unsigned int id;  th=reinterpret_cast<HANDLE>(_beginthreadex(0,0,&thread_proxy<F>,(void *)&func,0,&id)); if (!th) throw thread_error(); }
};

#elif defined(MPTASKS)

inline void assert_thread(OSStatus s) { assert(s==noErr); }
inline void check_thread( OSStatus s) { if (s==noErr) return; throw thread_error(s); }

class mutex
{
  public:
    typedef mutex self;
    typedef scoped_lock<self> scoped_lock;

  private:
    MPCriticalRegionID mut;

  public:
    inline mutex() { check_thread(MPCreateCriticalRegion(&mut)); }
    inline ~mutex() { MPDeleteCriticalRegion(mut); }

    inline void acquire() { MPEnterCriticalRegion(mut,kDurationForever); }
    inline void release() { MPExitCriticalRegion(mut); }

  private:
    mutex(const self &);
    const self &operator=(const self &);
};

typedef mutex spin_mutex;


class condition
{
  public:
    typedef condition self;

  private:
    MPEventID event;

    condition(const self &);
    const self &operator=(const self &);

  public:
    inline condition () { assert_thread(MPCreateEvent(&event)); }
    inline ~condition() { assert_thread(MPDeleteEvent( event)); }

    inline void signal() { assert_thread(MPSetEvent(event,1)); }
    template<typename Mutex> inline void wait(Mutex &mut) { mut.release(); assert_thread(MPWaitForEvent(event,NULL,kDurationForever)); mut.acquire(); }
};


class condition_mutex
{
  public:
    typedef condition_mutex self;
    typedef scoped_lock<self> scoped_lock;

    typedef spin_mutex mutex_type;
    typedef condition condition_type;

  private:
    mutex_type mut;
    condition_type cond;

  public:
    inline void acquire() { mut.acquire(); }
    inline void release() { mut.release(); }

    inline void signal() { cond.signal(); }
    inline void wait() { cond.wait(mut); }
};


class semaphore
{
  public:
    typedef semaphore self;
    typedef scoped_lock<self> scoped_lock;

  private:
    MPSemaphoreID sem;

    semaphore(const self &);
    const self &operator=(const self &);

  public:
    inline semaphore(     ) { assert_thread(MPCreateSemaphore(std::numeric_limits<int>::max(),0,&sem)); };
    inline semaphore(int n) { assert_thread(MPCreateSemaphore(std::numeric_limits<int>::max(),n,&sem)); };
    inline ~semaphore() { assert_thread(MPDeleteSemaphore(sem)); }

    inline void release()      { check_thread(MPSignalSemaphore(sem)); }
    inline void release(int n) { for (int i=0; i<n; ++i) release(); }

    inline void acquire()      { assert_thread(MPWaitOnSemaphore(sem,kDurationForever)); }
    inline int try_acquire()      { return MPWaitOnSemaphore(sem,kDurationImmediate)==noErr; }
    inline int try_acquire(int n) { for (int i=0; i<n; ++i) if (!try_acquire()) return i; return n; }
};


template<class F>
OSStatus thread_proxy(void *f)
{
  (*static_cast<F *>(f))();
  return 0;
}

template<class F>
class thread
{
  public:
    typedef thread self;
    typedef F function_type;

  private:
    function_type func;
    MPQueueID queue;
    MPTaskID task;
    
  public:
    inline thread() : func() { init(); }
    template<class A> inline explicit thread(      A &a) : func(a) { init(); }
    template<class A> inline explicit thread(const A &a) : func(a) { init(); }
    template<class A,class B> inline explicit thread(      A &a, const B &b) : func(a,b) { init(); }
    template<class A,class B> inline explicit thread(const A &a, const B &b) : func(a,b) { init(); }

    inline thread(const self &x) : func(x.func) { init(); }
    inline self &operator=(const self &x) { func=x.func; init(); }

    inline ~thread() { join(); MPDeleteQueue(queue); }

    inline void join() { check_thread(MPWaitOnQueue(queue,NULL, NULL,NULL,kDurationForever)); }

  private:
    void init() 
    { 
      static bool b = MPLibraryIsLoaded();
      if (!b) throw thread_error();
      check_thread(MPCreateQueue(&queue));
      OSStatus s=MPCreateTask(&thread_proxy<F>,(void *)&func,0,queue,NULL,NULL,0,&task);
      if (s!=noErr) { MPDeleteQueue(queue); throw thread_error(s); }
    }
};

#elif defined(PTHREADS)

inline void assert_thread(int r) { assert(r==0); }
inline void check_thread(int r) { if (r==0) return; throw thread_error(r); }

class mutex
{
  public:
    typedef mutex self;
    typedef scoped_lock<self> scoped_lock;

  private:
    pthread_mutex_t mut;

  public:
    inline mutex() { pthread_mutex_init(&mut,NULL); }
    inline ~mutex() { assert_thread(pthread_mutex_destroy(&mut)); }

    inline void acquire() { assert_thread(pthread_mutex_lock(&mut)); }
    inline void release() { assert_thread(pthread_mutex_unlock(&mut)); }

    inline pthread_mutex_t       &get()       { return mut; }
    inline const pthread_mutex_t &get() const { return mut; }

  private:
    mutex(const self &);
    const self &operator=(const self &);
};

#if !defined(__ICL) && !defined(_MSC_VER)
typedef mutex spin_mutex;
#else
class spin_mutex
{
  public:
    typedef spin_mutex self;
    typedef scoped_lock<self> scoped_lock;

  private:
    pthread_spinlock_t mut;

  public:
    inline spin_mutex() { pthread_spin_init(&mut,PTHREAD_PROCESS_PRIVATE); }
    inline ~spin_mutex() { assert_thread(pthread_spin_destroy(&mut)); }

    inline void acquire() { assert_thread(pthread_spin_lock(&mut)); }
    inline void release() { assert_thread(pthread_spin_unlock(&mut)); }

    inline pthread_spinlock_t       &get()       { return mut; }
    inline const pthread_spinlock_t &get() const { return mut; }

  private:
    spin_mutex(const self &);
    const self &operator=(const self &);
};
#endif

class condition
{
  public:
    typedef condition self;

  private:
    pthread_cond_t cond;

    condition(const self &);
    const self &operator=(const self &);

  public:
    inline condition () { check_thread(pthread_cond_init(&cond,0)); }
    inline ~condition() { pthread_cond_destroy(&cond); }

    template<typename Mutex> inline void wait(Mutex &mut) { pthread_cond_wait(&cond,&mut.get()); }
    inline void signal() { pthread_cond_signal(&cond); }
};

class condition_mutex
{
  public:
    typedef condition_mutex self;
    typedef scoped_lock<self> scoped_lock;

    typedef mutex mutex_type;
    typedef condition condition_type;

  private:
    mutex_type mut;
    condition_type cond;

  public:
    inline ~condition_mutex() {}

    inline void acquire() { mut.acquire(); }
    inline void release() { mut.release(); }

    inline void signal() { cond.signal(); }
    inline void wait() { cond.wait(mut); }
};

#if (defined(__MACOSX__) || defined(__APPLE__))
class semaphore
{
  public:
    typedef semaphore self;
    typedef scoped_lock<self> scoped_lock;

    typedef condition_mutex mutex_type;

  private:
    int count;
    int waiting;
    mutex_type cond;

    semaphore(const self &);
    const self &operator=(const self &);

  public:
    inline semaphore(     ) : count(0), waiting(0) { }
    inline semaphore(int n) : count(n), waiting(0) { }
    inline ~semaphore() {}

    inline void release()      { mutex_type::scoped_lock lock(cond); ++count; if (waiting>0) cond.signal(); }
    inline void release(int n) { for (int i=0; i<n; ++i) release(); }

    inline void acquire()      { mutex_type::scoped_lock lock(cond);  ++waiting; while (count==0) cond.wait();  --waiting;  --count; }
    inline void acquire(int n) { mutex_type::scoped_lock lock(cond); waiting+=n; while (count< n) cond.wait(); waiting-=n; count-=n; }

    inline int try_acquire()      { mutex_type::scoped_lock lock(cond); if (count>0) { --count; return 1; } else return 0; }
    inline int try_acquire(int n) { mutex_type::scoped_lock lock(cond); if (count>=n) count-=n; else { n=count; count=0; } return n; }
};
#else
class semaphore
{
  public:
    typedef semaphore self;
    typedef scoped_lock<self> scoped_lock;

  private:
    sem_t sem;

    semaphore(const self &);
    const self &operator=(const self &);

  public:
    inline semaphore(     ) { assert_thread(sem_init(&sem,0,0)); };
    inline semaphore(int n) { assert_thread(sem_init(&sem,0,n)); };
    inline ~semaphore() { assert_thread(sem_destroy(&sem)); }

    inline void release()      { assert_thread(sem_post(&sem)); }
#if defined(__GNUC__)
    inline void release(int n) { for (int i=0; i<n; ++i) release(); }
#else
    inline void release(int n) { assert_thread(sem_post_multiple(&sem,n)); }
#endif

    inline void acquire()      { assert_thread(sem_wait(&sem)); }
    inline int try_acquire()      { return sem_trywait(&sem)!=-1; }
    inline int try_acquire(int n) { for (int i=0; i<n; ++i) if (!try_acquire()) return i; return n; }

    inline unsigned int value() { int sval; assert_thread(sem_getvalue(&sem, &sval)); return sval; }
};
#endif



template<class F>
void *thread_proxy(void *f)
{
  (*static_cast<F *>(f))();
  return 0;
}

template<class F>
class thread
{
  public:
    typedef thread self;
    typedef F function_type;

  private:
    function_type func;
    pthread_t th;

  public:
    inline thread() : func() { init(); }
    inline explicit thread(const F &f) : func(f) { init(); }
    template<class A> inline explicit thread(      A &a) : func(a) { init(); }
    template<class A> inline explicit thread(const A &a) : func(a) { init(); }
    template<class A,class B> inline explicit thread(      A &a, const B &b) : func(a,b) { init(); }
    template<class A,class B> inline explicit thread(const A &a, const B &b) : func(a,b) { init(); }

    inline thread(const self &x) : func(x.func) { init(); }
    inline self &operator=(const self &x) { func=x.func; init(); return *this; }

    inline ~thread() { join(); pthread_detach(th); }

    inline function_type       &function()       { return func; }
    inline function_type const &function() const { return func; }

    inline void join() { pthread_join(th,0); }

    inline void init() { check_thread(pthread_create(&th,0,&thread_proxy<function_type>,(void *)&func)); }
};

#endif


template <typename It>
class task_data
{
  public:
    typedef It iterator;

  private:
    iterator it1, it2;

  public:
    inline task_data(iterator begin, iterator end) : it1(begin), it2(end) {}
    inline iterator begin() const { return it1; }
    inline iterator end  () const { return it2; }
};

template<class RanIt>
class task_manager
{
  public:
    typedef task_manager self;
    typedef RanIt iterator;
    typedef task_data<iterator> result_type;
    typedef condition_mutex mutex_type;

  protected:
    iterator it;
    int count;
    int grain;
    int nthreads;
    static mutex_type mut;

  public:
    inline task_manager(iterator begin,int n,int grainsize, int nth) : it(begin),count(n),grain(grainsize),nthreads(nth) { }

    inline int size() const { return count; }
    inline int grain_size() const { return grain; }
    inline int nchuncks() const { return (size()-1)/grain_size()+1; }

    inline result_type operator()()
    {
      typename mutex_type::scoped_lock lock(mut);
      int n=gmt_min(size(),grain_size());
      iterator begin=it;
      it+=n; count-=n;
      return result_type(begin,it);
    };

    inline void release()
    {
      typename mutex_type::scoped_lock lock(mut);
      if (!--nthreads) mut.signal();
    }

    inline void join()
    {
      typename mutex_type::scoped_lock lock(mut);
      while (nthreads!=0) mut.wait();
    }
};

template<class RanIt> typename task_manager<RanIt>::mutex_type task_manager<RanIt>::mut;

template<class F,class RanIt>
class reduction_task_manager : public task_manager<RanIt>
{
  public:
    typedef reduction_task_manager self;
    typedef task_manager<RanIt> base;
    typedef F function_type;

    typedef typename base::iterator iterator;
    typedef typename base::mutex_type mutex_type;

    using base::mut;

  private:
    function_type &func;

  public:
    inline reduction_task_manager(function_type &f,iterator begin,int count,int grainsize, int nth) : base(begin,count,grainsize,nth), func(f) { }

    using base::join;

    inline void join(const function_type &f)
    {
      typename mutex_type::scoped_lock lock(mut);
      func.join(f);
    };
};

template<class M,class F>
class managed_function
{
  public:
    typedef managed_function self;
    typedef M manager_type;
    typedef F function_type;

  private:
    mutable manager_type &manager;
    function_type func;

  public:
    inline managed_function(const self &x) : manager(x.manager), func(x.func) {}
    inline managed_function(manager_type &m, const function_type &f) : manager(m), func(f) {}

    inline void join() { manager.join(); }

    inline void operator()()
    {
      for (typename manager_type::result_type data=manager(); data.begin()!=data.end(); data=manager())
        func(data.begin(),data.end());
      manager.release();
    }
};

template<class M,class F>
class reduced_function
{
  public:
    typedef reduced_function self;
    typedef M manager_type;
    typedef F function_type;

  private:
    mutable manager_type &manager;
    function_type func;

  public:
    inline reduced_function(const self &x) : manager(x.manager), func(x.func) {}
    inline reduced_function(manager_type &m, const function_type &f) : manager(m), func(f) {}

    inline void join() { manager.join(); }

    inline void operator()()
    {
      for (typename manager_type::result_type data=manager(); data.begin()!=data.end(); data=manager())
        func(data.begin(),data.end());
      manager.join(func);
      manager.release();
    }
};


template<class F>
void function_proxy(void *pf)
{
  (*static_cast<F *>(pf))();
}

struct proxy_function
{
  public:
    typedef proxy_function self;
    typedef void (*function_type)(void *);
    typedef double value_type;

  private:
    function_type func;
    int count;
    value_type *param;

  public:
    inline proxy_function() : func(NULL), count(0), param(NULL) {}

    template<class F> inline explicit proxy_function(const F &f) : func(function_proxy<F>), count((sizeof(F)-1)/sizeof(value_type)+1), param(new value_type[count])
    {
      value_type *p = (value_type *)&f;
      copy(p,p+count,param);
    }

    inline proxy_function(const self &x) : func(x.func), count(x.count), param(new value_type[count])
    {
      copy(x.param,x.param+count,param);
    }

    inline self &operator=(const self &x)
    {
      func=x.func;
      if (count!=x.count) { count=x.count; delete param; param=new value_type[count]; }
      copy(x.param,x.param+count,param);
      return *this;
    }

    template<class F> inline self &operator=(const F &f)
    {
      func=function_proxy<F>;
      int n=(sizeof(F)-1)/sizeof(value_type)+1;
      if (count!=n) { count=n; delete param; param=new value_type[count]; }
      value_type *p = (value_type *)&f;
      copy(p,p+count,param);
      return *this;
    }

    inline ~proxy_function() { delete param; }
    inline void operator()() { func(param); }
};

//struct proxy_function
//{
//  public:
//    typedef proxy_function self;
//    typedef void (*function_type)(void *);
//
//  private:
//    function_type func;
//    void *param;
//
//  public:
//    inline proxy_function() : func(NULL), param(NULL) {}
//    template<class F> explicit proxy_function(const F &f) : func(function_proxy<F>), param(&f) {}
//    inline proxy_function(const self &x) : func(x.func), param(x.param) {}
//
//    self &operator=(const self &x) { func=x.func; param=x.param; return *this; }
//    template<class F> self &operator=(const F &f) { func=function_proxy<F>; param=(void *)&f; return *this; }
//
//    inline void operator()() { func(param); }
//};

struct thread_manager
{
  public:
    typedef thread_manager self;
    typedef scoped_lock<self> lock;
    typedef proxy_function function_type;
    typedef condition_mutex condition_type;

  private:
    bool b;
    int count;
    //semaphore threads;
    semaphore reserved;
    condition_type cond;
    function_type func;

  public:
    inline thread_manager() : b(true), count(0) { }

    inline ~thread_manager() { }

    //inline void acquire    (int n) { threads.acquire(n); }
    //inline int  try_acquire(int n) { return threads.try_acquire(n); }

    inline void acquire    (int n) { }
    inline int  try_acquire(int n) {  return gmt_min(n,num_threads()); }

    inline void set(int n)
    {
      condition_type::scoped_lock lock(cond);
      while (count!=0) cond.wait();
      count+=n;
      b=false;
      reserved.release(n);
    }

    template<class F> inline void set(int n,const F &f)
    {
      if (n==0) return;
      condition_type::scoped_lock lock(cond);
      while (count!=0) cond.wait();
      count+=n;
      func=f;
      b=true;
      reserved.release(n);
    }

    inline bool get(function_type &f)
    {
      //threads.release();
      reserved.acquire();
      condition_type::scoped_lock lock(cond);
      if (--count==0) cond.signal();
      f=func;
      return b;
    }
};

template<class M>
class managed_thread_function
{
  public:
    typedef managed_thread_function self;
    typedef M manager_type;
    typedef typename manager_type::function_type function_type;

  private:
    manager_type &manag;

  public:
    inline managed_thread_function(manager_type &m) : manag(m) {}

    inline manager_type       &manager()       { return manag; }
    inline const manager_type &manager() const { return manag; }

    inline void operator()()
    {
      __aligned_stack  // GCC preserves the stack alignment, but some pthread implementations obviously do not.
      function_type func;
      while (manag.get(func))
        func();
    }
};

template<class M>
class managed_thread
{
  public:
    typedef managed_thread self;
    typedef M manager_type;
    typedef thread<managed_thread_function<manager_type> > thread_type;

  private:
    thread_type th;

  public:
    inline managed_thread(manager_type &manager) : th(manager) {}
};

template<class M>
class managed_thread_group
{
  public:
    typedef managed_thread_group self;
    typedef M manager_type;
    typedef managed_thread<manager_type> thread_type;

  private:
    manager_type manag;
    vector<thread_type*> threads;

  public:
    inline managed_thread_group(int n) : manag(), threads(0)
    {
      for (int i=0; i<n; ++i)
        threads.push_back(new thread_type(manag));
    }

    inline ~managed_thread_group()
    {
      manag.set(size());
      while (!threads.empty())
      {
        delete threads.back();
        threads.pop_back();
      }
    }

    inline int size() const { return threads.size(); }

    inline void acquire   (int n) { manag.acquire(n); }
    inline int try_acquire(int n) { return manag.try_acquire(n); }
    template<class F> inline void set(int n,const F &f) { manag.set(n,f); }
};

#if !defined(OPENMP)
static managed_thread_group<thread_manager> thread_pool(max_threads());
#endif


//{unsecret}
//Summary: Computes parallel iteration of a function over a range of values
//Arguments:
//  begin - Integer or random access iterator addressing the position of the first element in the range to be operated on.
//  end   - Integer or random access iterator addressing the position one past the final element in the range operated on.
//  f     - User-defined function object that is applied on subranges.
//  grain - Number of iteration for a reasonable chunk to deal out to a thread.
//          If no value given, the range is subdivided in num_threads() chunks.
//  X     - Range to be operated on.
//Remarks: The function object must model the following requirements:
//
//    @untitled table
//    Copy constructor      F::F(const F &)
//    Destructor            F::~F(const F &)
//    Apply on subranges    template<class RanIt> F::operator()(RanIt,RanIt) const
//Example:
//  template<class V>
//  struct fill_function
//  {
//    V val;
//    fill_function(const V &v) : val(v) {}
//    template<class It> void operator()(It begin, It end) const { fill(begin,end,val); }
//  };
//
//  template<class RanIt,class T>
//  void parallel_fill(RanIt begin, RanIt end, const T &val, int grain=10000)
//  {
//    gmt::parallel_for(X.begin(),X.end(), fill_function<T>(val),grain);
//  }
//See: ^parallel_copy^, ^parallel_reduce^
template<class RanIt,class F> void parallel_for(RanIt begin,RanIt end,const F &f,int grain)
{
#if MAX_THREADS==1
  f(begin,end);
#elif defined(OPENMP)
  int n=end-begin;
  int nchunks=(n-1)/grain+1;
  if (nchunks<=1) return f(begin,end);

  #pragma omp parallel for schedule(dynamic,1) if (nchunks>1)
  for (int i=0; i<nchunks; ++i)
  {
    int d=i*grain;
    RanIt it=begin+d;
    f(it,it+gmt_min(n-d,grain));
  }
#else
  typedef task_manager<RanIt> manager_type;
  typedef managed_function<manager_type,const F> function_type;

  int n=end-begin;
  int nchunks=(n-1)/grain+1;
  int nthreads=thread_pool.try_acquire(nchunks);
  if (nthreads<=1) return f(begin,end);

  manager_type manag(begin,n,grain,nthreads);
  function_type func(manag,f);

  thread_pool.set(nthreads,func);
  func.join();
#endif
}
//{unsecret}
template<class G,class F> void parallel_for(const Vector<G> &X, const F &f,int grain) { parallel_for(X.begin(),X.end(),f,grain); }
//{unsecret}
template<class G,class F> void parallel_for(const Vector<G> &X, const F &f)           { parallel_for(X.begin(),X.end(),f,(X.nelms()-1)/num_threads()+1); }
//{unsecret}
template<class G,class F> void parallel_for(const Matrix<G> &X, const F &f,int grain) { parallel_for(X.begin(),X.end(),f,grain); }
//{unsecret}
template<class G,class F> void parallel_for(const Matrix<G> &X, const F &f)           { parallel_for(X.begin(),X.end(),f,(X.nelms()-1)/num_threads()+1); }

//{unsecret}
//Summary: Computes parallel reduction over a range of values
//Arguments:
//  begin - Integer or random access iterator addressing the position of the first element in the range to be operated on.
//  end   - Integer or random access iterator addressing the position one past the final element in the range operated on.
//  f     - User-defined function object that is applied on subranges.
//  grain - Number of iteration for a reasonable chunk to deal out to a thread.
//          If no value given, the range is subdivided in num_threads() chunks.
//  X     - Range to be operated on.
//Remarks:
//  The function object must model the following requirements:
//
//    @untitled table
//    Copy constructor                F::F(const F &)
//    Destructor                      F::~F(const F &)
//    Merge the results of 2 threads  F::join(const F &) const
//    Apply on subranges              template<class RanIt> F::operator()(RanIt,RanIt) const
//
//  Any modified variable, that is used to store intermediate results, has to be declared as mutable.
//Example:
//  template<class V>
//  struct accumulate_function
//  {
//    mutable V val;
//    accumulate_function() : val(0) {}
//    template<class It> void operator()(It begin, It end) const { val=accumulate(begin,end,val); }
//    void join(const accumulate_function &x) const { val+=x.val; }
//  };
//
//  template<class RanIt,class T>
//  T parallel_accumulate(RanIt begin, RanIt end, const T &val, int grain=10000)
//  {
//    accumulate_function<T> func;
//    gmt::parallel_reduce(begin,end,func,grain);
//    return val+func.val;
//  }
//See: ^parallel_copy^, ^parallel_for^
template<class RanIt,class F> void parallel_reduce(RanIt begin,RanIt end,const F &f,int grain)
{
#if MAX_THREADS==1
  f(begin,end);
#elif defined(OPENMP)
  int n=end-begin;
  int nchunks=(n-1)/grain+1;

  const F func;
  #pragma omp parallel if (nchunks>1) firstprivate(func)
  {
    #pragma omp for schedule(dynamic,1)
    for (int i=0; i<nchunks; ++i)
    {
      int d=i*grain;
      RanIt it=begin+d;
      func(it,it+gmt_min(n-d,grain));
    }
    #pragma omp critical
    f.join(func);
  }
#else
  typedef reduction_task_manager<const F,RanIt> manager_type;
  typedef reduced_function<manager_type,const F> function_type;

  int n=end-begin;
  int nchunks=(n-1)/grain+1;
  int nthreads=thread_pool.try_acquire(nchunks);
  if (nthreads<=1) return f(begin,end);

  manager_type manag(f,begin,n,grain,nthreads);
  function_type func(manag,f);

  thread_pool.set(nthreads,func);
  func.join();
#endif
}
//{unsecret}
template<class G,class F> void parallel_reduce(const Vector<G> &X, const F &f,int grain) { parallel_reduce(X.begin(),X.end(),f,grain); }
//{unsecret}
template<class G,class F> void parallel_reduce(const Vector<G> &X, const F &f)           { parallel_reduce(X.begin(),X.end(),f,(X.nelms()-1)/num_threads()+1); }
//{unsecret}
template<class G,class F> void parallel_reduce(const Matrix<G> &X, const F &f,int grain) { parallel_reduce(X.begin(),X.end(),f,grain); }
//{unsecret}
template<class G,class F> void parallel_reduce(const Matrix<G> &X, const F &f)           { parallel_reduce(X.begin(),X.end(),f,(X.nelms()-1)/num_threads()+1); }




template<class RanIt, class OutIt>
struct copy_function
{
  RanIt it1;
  OutIt it2;
  copy_function(RanIt begin, OutIt out) : it1(begin), it2(out) {}
  template<class It> void operator()(It begin, It end) const { copy(begin,end,it2+(begin-it1)); }
};

//{unsecret}
//Summary: Assigns the values of elements from a source range to a destination range.
//Arguments:
//  begin - Random access iterator addressing the position of the first element in the source range.
//  end   - Random access iterator addressing the position one past the final element in the source range.
//  out   - Random access iterator addressing the position of the first element in the destination range.
//  grain - Number of iteration for a reasonable chunk to deal out to a thread.
//          If no value given, the range is subdivided in num_threads() chunks.
//  X     - Source range.
//  Y     - Destination range.
//Example:
//  DenseMatrix<ucharImage::index_type>::self M;
//  gmt::parallel_copy(motion<8,8>(X,Y,32,32),M,X.ncols()/8);
//Return: An output iterator addressing the position that is one past the final element in the destination range
//See: ^parallel_for^, ^parallel_reduce^
template<class RanIt,class OutIt> OutIt parallel_copy(RanIt begin,RanIt end,OutIt out,int grain)
{
  parallel_for(begin,end,copy_function<RanIt,OutIt>(begin,out),grain);
  return out+(end-begin);
}
//{unsecret}
template<class G1,class G2> void parallel_copy(const Vector<G1> &X, Vector<G2> &Y,int grain) { Y.resize(X.size()); Y.set_lower_bound(X.lower_bound()); parallel_copy(X.begin(),X.end(),Y.begin(),grain); }
//{unsecret}
template<class G1,class G2> void parallel_copy(const Vector<G1> &X, Vector<G2> &Y)           { Y.resize(X.size()); Y.set_lower_bound(X.lower_bound()); parallel_copy(X.begin(),X.end(),Y.begin(),(X.nelms()-1)/num_threads()+1); }
//{unsecret}
template<class G1,class G2> void parallel_copy(const Matrix<G1> &X, Matrix<G2> &Y,int grain) { Y.resize(X.size()); Y.set_lower_bound(X.lower_bound()); parallel_copy(X.begin(),X.end(),Y.begin(),grain); }
//{unsecret}
template<class G1,class G2> void parallel_copy(const Matrix<G1> &X, Matrix<G2> &Y)           { Y.resize(X.size()); Y.set_lower_bound(X.lower_bound()); parallel_copy(X.begin(),X.end(),Y.begin(),(X.nelms()-1)/num_threads()+1); }


#ifdef HIDE_FROM_DOCJET
#else
} //namespace gmt
#endif

#endif
