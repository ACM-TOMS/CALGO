package body Unix_Resource_Usage is

  package C_Interfaces is

    times_map : constant array (times_enum)
    		of integer
	      := (self => 0, children => -1);

    function getrusage (who : integer;
    			rusage : system.address)
	return integer;
    pragma interface (C, getrusage);

    function timeval_to_duration (tv : timeval)
        return duration;

  end C_Interfaces;

    function Get_Process_Times (who : times_enum := self)
	return Process_Times
    is
      answer : Process_Times;
      c_result : integer;
    begin
      c_result := C_Interfaces.getrusage
			(who => C_Interfaces.times_map(who),
			 rusage => answer'address);
      if (c_result = -1) then
	raise program_error;	-- something broke in Unix!
      else
	return answer;
      end if;
    end Get_Process_Times;

    function Total_Time_of (Times: in Process_Times)
        return duration
    is
    begin
      return user_cpu_time_of (times) + system_cpu_time_of (times);
    end;

    function User_CPU_Time_Of (Times: in Process_Times)
        return Duration
    is
    begin
      return C_Interfaces.timeval_to_duration (times.ru_utime);
    end User_CPU_Time_Of;

    function System_CPU_Time_Of (Times: in Process_Times)
        return Duration
    is
    begin
      return C_Interfaces.timeval_to_duration (times.ru_stime);
    end System_CPU_Time_Of;

    function Max_Resident_Set_Size_of (Times: in Process_Times)
        return natural
    is
    begin
      return times.ru_maxrss;
    end max_resident_set_size_of;

    function Shared_Pages_Value_of (Times: in Process_Times)
        return page_seconds
    is
    begin
      return page_seconds(times.ru_ixrss);
    end;

    function Unshared_Data_Pages_Value_of (Times: in Process_Times)
    	return page_seconds
    is
    begin
      return page_seconds(times.ru_idrss);
    end;

    function Stack_Pages_Value_of (Times: in Process_Times)
        return page_seconds
    is
    begin
      return page_seconds(times.ru_isrss);
    end;

    function Non_IO_Page_Faults_of (Times: in Process_Times)
        return natural
    is
    begin
      return times.ru_minflt;
    end;

    function IO_Page_Faults_of (Times: in Process_Times)
    	return natural
    is
    begin
      return times.ru_majflt;
    end;

    function Swaps_of (Times : in Process_Times)
        return natural
    is
    begin
      return times.ru_nswap;
    end;

    function Input_Blocks_of (Times : in Process_Times)
        return natural
    is
    begin
      return times.ru_inblock;
    end;

    function Output_Blocks_of (Times : in Process_Times)
        return natural
    is
    begin
      return times.ru_outblock;
    end;

    function Socket_Messages_Sent_of (Times : in Process_Times)
    	return natural
    is
    begin
      return times.ru_msgsnd;
    end;

    function Socket_Messages_Received_of (Times : in Process_Times)
        return natural
    is
    begin
      return times.ru_msgrcv;
    end;

    function Signals_Delivered_of (Times : in Process_Times)
        return natural
    is
    begin
      return times.ru_nsignals;
    end;

    function Voluntary_Context_Switches_of (Times: in Process_Times)
        return natural
    is
    begin
      return times.ru_nvcsw;
    end;

    function Involuntary_Context_Switches_of (Times: in Process_Times)
        return natural
    is
    begin
      return times.ru_nivcsw;
    end;

  package body C_Interfaces is

    function timeval_to_duration (tv : timeval)
        return duration
    is
      answer : duration;
    begin
      -- on a sun:
        answer := duration(tv.tv_sec) + duration(tv.tv_usec)/1_000_000;
      -- on a dec:
      -- answer := duration(tv.tv_sec);
      -- if float(tv.tv_usec) <= duration'large
      --  then answer := answer + duration(tv.tv_usec)/1_000_000;
      --  else answer := answer + duration(tv.tv_usec/100)/10_000;
      -- end if;
      -- because of the strange fact that on a dec
      -- duration'large is about 2.14E+5 < 1_000_000;
      -- with the following trials, only the seconds were printed:
      --  answer := duration(tv.tv_sec + tv.tv_usec/1_000_000);
      --  answer := duration(tv.tv_sec) + duration(tv.tv_usec/1_000_000);
      return answer;
    end timeval_to_duration;

  end C_Interfaces;

end Unix_Resource_Usage;
