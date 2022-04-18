with System;

package Unix_Resource_Usage is

  type Process_Times is private;

  type times_enum is (self, children);

  type page_seconds is new integer;

  --   expressed in units of pages * clock ticks (1 tick = 1/50 second).
  --   The value is calculated by summing the number of shared memory
  --   pages in use each time the internal system clock ticks, and then
  --   averaging over 1 second intervals.

  function Get_Process_Times ( who : times_enum := self ) return Process_Times;

  function Total_Time_of ( times: in Process_Times ) return duration;

  function User_CPU_Time_Of ( times: in Process_Times ) return duration;

  function System_CPU_Time_Of ( times: in Process_Times ) return duration;

  function Max_Resident_Set_Size_of ( times: in Process_Times ) return natural;

  function Shared_Pages_Value_of ( times: in Process_Times )
                                 return page_seconds;

  -- DESCRIPTION :
  --   returns the amount of memory used by the text segment which was
  --   also shared among other processes.

  function Unshared_Data_Pages_Value_of ( times: in Process_Times )
                                        return page_seconds;

  -- DESCRIPTION :
  --   returns the amount of unshared memory residing in the data segment
  --   of the process.

  function Stack_Pages_Value_of ( times: in Process_Times ) return page_seconds;

  -- DESCRIPTION :
  --   returns the amount of unshared memory residing in the stack segment
  --   of the process

  function Non_IO_Page_Faults_of ( times: in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of page faults serviced without any I/O activity;
  --   here I/O activity is avoided by "reclaiming" a page frame from the
  --   list of pages awaiting reallocation.

  function IO_Page_Faults_of ( times: in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of page faults serviced which required I/O activity.

  function Swaps_of ( times : in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of times the process was swapped out of main memory.

  function Input_Blocks_of ( times : in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of times the file system had to perform input.

  function Output_Blocks_of ( times : in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of times the file system had to perform output.

  function Socket_Messages_Sent_of ( times : in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of messages sent over sockets.

  function Socket_Messages_Received_of ( times : in Process_Times )
                                       return natural;
  -- DESCRIPTION :
  --   returns the number of messages received over sockets.

  function Signals_Delivered_of ( times : in Process_Times ) return natural;

  -- DESCRIPTION :
  --   returns the number of signals delivered.

  function Voluntary_Context_Switches_of ( times: in Process_Times )
                                         return natural;
  -- DESCRIPTION :
  --   returns the number of times a context switch resulted due to a process
  --   voluntarily giving up the processor before its time slice was completed
  --   (usually to await availability of a resource).

  function Involuntary_Context_Switches_of ( times: in Process_Times )
                                           return natural;
  -- DESCRIPTION :
  --   returns the number of times a context switch resulted due to a
  --   higher priority process becoming runnable or because the current
  --   process exceeded its time slice.

private

  type timeval is record
    tv_sec : integer;        -- Ada integer is C/SunOS long
    tv_usec : integer;       -- Ada integer is C/SunOS long
  end record;

  type rusage is record
    ru_utime : timeval;
    ru_stime : timeval;
    ru_maxrss : integer;
    ru_ixrss : integer;    -- integral shared text memory size
    ru_idrss : integer;    -- integral unshared data size
    ru_isrss : integer;    -- integral unshared stack size
    ru_minflt : integer;   -- page reclaims
    ru_majflt : integer;   -- page faults
    ru_nswap : integer;    -- swaps
    ru_inblock : integer;  -- block input operations
    ru_outblock : integer; -- block output operations
    ru_msgsnd : integer;   -- messages sent
    ru_msgrcv : integer;   -- messages received
    ru_nsignals : integer; -- signals received
    ru_nvcsw : integer;    -- voluntary context switches
    ru_nivcsw : integer;   -- involuntary context switches
  end record;

  type process_times is new rusage;

  pragma inline (get_process_times, total_time_of, user_cpu_time_of,
                 system_cpu_time_of, max_resident_set_size_of,
                 shared_pages_value_of, unshared_data_pages_value_of,
                 stack_pages_value_of, non_io_page_faults_of,
                 io_page_faults_of, swaps_of, input_blocks_of,
                 output_blocks_of, socket_messages_sent_of,
                 socket_messages_received_of, signals_delivered_of,
                 voluntary_context_switches_of,
                 involuntary_context_switches_of);

end Unix_Resource_Usage;
