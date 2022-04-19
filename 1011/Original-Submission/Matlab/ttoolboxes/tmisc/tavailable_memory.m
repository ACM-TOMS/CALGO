function freemem = tavailable_memory()
% AVAILABLE_MEMORY returns the memory which can be (pontentially)
% used by Matlab. Works for Windows, Mac and Linux.
% Variable freemem is expressed in BYTE.
    % Taken from: The JSRToolbox by Raphael Jungers
    try
        if(ispc)
            [~,memStats] = memory;
            freemem = memStats.PhysicalMemory.Available;% Built-in function.
        elseif(isunix && not(ismac))
            [~,w] = unix('free -b | grep Mem'); % Excutes free on UNIX system and takes "memory" line (in Byte, -b).
            stats = str2double(regexp(w, '[0-9]*', 'match')); % Splits string into numerical array and ignores char.
            freemem = stats(3) + 0.75*stats(end) ; % Takes the 3 col value (free memory) + 3/4 of cache.
        elseif(ismac)
            [~,w] = unix('vm_stat | grep free'); % Excute vm_stat on MAC system and takes "free" line.
            spaces = strfind(w,' '); % Finds the numeric value.
            freemem = str2double(w(spaces(end):end))*4096; % vm_stat is expressed in pages : we have to multiply it by 4096.
        end
    catch WE %#ok<NASGU>
        freemem = 1e9; % default : 1 GB 
    end
end