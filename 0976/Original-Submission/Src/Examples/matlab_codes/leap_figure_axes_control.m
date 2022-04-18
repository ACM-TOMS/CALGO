function [] = leap_figure_axes_control(varargin)

global config

config.dead_zone = zeros(3,1);


cameratoolbar


t0 = clock;
while 1
    
    num_pointers = establish_number_pointers();
    
    
    
    
    if num_pointers == 0
        
        if etime(clock, t0)>15
%             display('leap timeout');
            break;
        end
        
        
        continue;
    
    else
%         display('resetting clock');
        t0 = clock;
    end
    
    
    
    if (establish_origin(num_pointers)~=1)
        continue;
    end
    
    
    
    switch num_pointers

        case 1
            orbit_camera(1);
            
        case 2
            glide(2);
            
        case 3
            view_pan(3);
            
        case 4
            
        otherwise
            
            
    end
    
    
    
    
end
    

end




function num_pointers = establish_number_pointers()

% display('establishing number of pointers');
prev_num_pointers = 0;
t_start = clock;
while 1
    
    
    
    
    f = matleap(1);
    
    curr_num_pointers = num_extended(f);
    if curr_num_pointers==prev_num_pointers%,curr_num_pointers~=0
        if etime(clock,t_start)>1
            break
        end
    else
        t_start = clock;
    end
        
    prev_num_pointers = curr_num_pointers;
    
end

num_pointers = curr_num_pointers;


end



function glide(required_num_fingers)
global config
% display('gliding');

av_origin = (config.origin(1,:)+ config.origin(2,:))/2;
while 1;
    
%     display('gliding camera');
    current_frame = matleap(1); 
    
	if num_extended(current_frame) ~= required_num_fingers
        break
	end
    
	curr_av_pos = get_average_position(current_frame);
	
	
    a = (curr_av_pos-av_origin)/100;
%     a = min((current_frame.pointables.position-config.origin),100)/100;
    a(abs(a)<0.10) = 0;
    
    
%     pause
    %get dtheta, dphi from the frame
%     dx = 1e-1;
%     dy = 1e-2;
    dz = 1e-1;
    camdolly(0,0,dz*a(3),'fixtarget');
    drawnow;
end

end


function view_pan(required_num_fingers)

global config


% av_origin = average_origin();
av_origin = mean(config.origin);%(config.origin(1,:)+ config.origin(2,:)+ config.origin(3,:))/3;
while 1;
    
   
%     display('orbiting camera');
    current_frame = matleap(1);

	if num_extended(current_frame) ~= required_num_fingers
        break
	end
    
    curr_av_pos = mean(cat(1,current_frame.pointables.position));
    

    a = (curr_av_pos-av_origin)/100;
%     a = min((current_frame.pointables.position-config.origin),100)/100;
    a(abs(a)<0.10) = 0;
    

    dtheta = 1;
    dphi = 1;
    campan(dtheta*a(1),dphi*a(2));
    drawnow;
end

end



function orbit_camera(required_num_fingers)
global config

while 1
    
%     display('orbiting camera');
    current_frame = matleap(1);
    
    if num_extended(current_frame) ~= required_num_fingers
        break
    end
    
    a = (get_average_position(current_frame)-config.origin)/100;

	a(abs(a)<0.10) = 0;
    
    
    %get dtheta, dphi from the frame
    dtheta = 10;
    dphi = 5;
    camorbit(dtheta*a(1),dphi*a(2));
    drawnow;
end


end


function success = establish_origin(num_pointers)
global config

display(sprintf('establishing origin for %i pointers',num_pointers));

config.origin = zeros(num_pointers,3);

t0 = clock;
have_lock = false;

positions = zeros(60,3,num_pointers); % why 60, again?
num_consecutive_frames = 0;
num_nonconsecutive_frames = 0;
prevframeid = 0;
while ~have_lock
    
    current_frame = matleap(1);
    
    if current_frame.id == prevframeid
        continue;
    else
        prevframeid= current_frame.id;
    end
    
    if num_extended(current_frame) ~= num_pointers
        t0 = clock;
        num_consecutive_frames = 1;
        num_nonconsecutive_frames = num_nonconsecutive_frames+1;
%         display('resetting cuz num_pointers');
    else
        num_nonconsecutive_frames = 1;
        num_consecutive_frames = num_consecutive_frames+1;
        
		extended_counter = 0;
		for ii = 1:num_pointers
			if current_frame.pointables(ii).is_extended
				extended_counter = extended_counter+1;
				positions(num_consecutive_frames,:,extended_counter) = current_frame.pointables.position;
				
			end
		end
        
    end
    
    
    
    if and(etime(clock,t0) > 0.5, num_consecutive_frames>=15) %, num_successful_frames
%         positions(:,:,1); % first colon is for all frames.  last
        
        for ii = 1:num_pointers
            config.origin(ii,:) = mean(positions(1:num_consecutive_frames,:,ii));
        end
        

        if (1)
            success = true;
            have_lock = true;
%             display('have origin');
        else
%             t0 = t0/2;
%             t0 = 0;
        end
        
       
    end
    
    if num_nonconsecutive_frames > 20
        success = false;
        break
    end
        
end


%pop up a dialog indicating establishing origin.  

%dialog remains up while doing so.


%when done, close dialog and return the origin
end




function curr_av_pos = get_average_position(current_frame)
curr_av_pos = zeros(1,3);

for ii = 1:num_extended(current_frame)
	if current_frame.pointables(ii).is_extended
		curr_av_pos = curr_av_pos + current_frame.pointables(ii).position;
	end
end
curr_av_pos = curr_av_pos/num_extended(current_frame);

end
	
function num = num_extended(frame)
num = 0;

for ii = 1:length(frame.pointables)
	if frame.pointables(ii).is_extended
		num = num+1;
	end
end


end
