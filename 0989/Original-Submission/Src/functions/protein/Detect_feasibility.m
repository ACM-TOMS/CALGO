function [Detect_fea]=Detect_feasibility(Pos,pos,Mdirect)

% The function is used to detect the feasibility of next candidate direction for backtracking-based repairing.
% The detailed detection method is described in paper:
% -----Benhui Chen, Jinglu Hu, "A Hybrid EDA for Protein Folding Based on HP Model", IEEJ Trans. on Electrical and Electronics Engineering, Vol.5 No.4, pp.459-466, July, 2010.

% Last version 7/31/2010.     


%INPUTS
% Pos: Configuration of the sequence being repaired;
% pos: position in the sequence being inspected;
% Mdirect: the next candidate direction for folding (0--left, 1-- front, 2--right);

%OUTPUT
% Detect_fea ==1 : the next candidate direction is feasible.

global InitConf;
Detect_fea=0;

for i=1:pos
    coor(i,:)=Pos(i,:);
end

if size(coor,1)==length(InitConf) % Folding is complete. 
    Detect_fea=1;
    return;    
end

remnant=length(InitConf)-size(coor,1);

% Define current four boundaries
curr_left=min(coor(:,1))-1;
curr_right=max(coor(:,1))+1;
curr_top=max(coor(:,2))+1;
curr_bottom=max(coor(:,2))-1;

% Calculate the coordinate of the next candidate direction
  i = size(coor,1)+1;

 if(coor(i-1,2)==coor(i-2,2))
  if (Mdirect==0)            % Left
   new_point(1,1) = coor(i-1,1);  
   new_point(1,2) = coor(i-1,2) + (coor(i-1,1)-coor(i-2,1));
  elseif (Mdirect==1)        % Front 
   new_point(1,1) = coor(i-1,1) + (coor(i-1,1)-coor(i-2,1));  
   new_point(1,2) = coor(i-1,2);
  else                         % Right
   new_point(1,1) = coor(i-1,1);  
   new_point(1,2) = coor(i-1,2) - (coor(i-1,1)-coor(i-2,1));
  end
  end
if (coor(i-1,1)==coor(i-2,1))
  if (Mdirect==0)            % Left
    new_point(1,2) = coor(i-1,2);  
    new_point(1,1) = coor(i-1,1) -  (coor(i-1,2)-coor(i-2,2));  
  elseif (Mdirect==1)        % Front 
    new_point(1,2) = coor(i-1,2) +  (coor(i-1,2)-coor(i-2,2));  
    new_point(1,1) = coor(i-1,1);
  else                         % Right
    new_point(1,2) = coor(i-1,2);  
    new_point(1,1) = coor(i-1,1) + (coor(i-1,2)-coor(i-2,2));
  end
 end

%% Floodfill strategy is used to detect empty positions.

if Position_over(coor,new_point)==1  %% Check overlapping
   Detect_fea=0;
   return;
end
%% Initialization
Count_coor=0;
queue=new_point;
head=1;
rear=1;

while (head <= rear)  %% the queue is not empty
   curr_point=queue(head,:); %% Set the first element of queue as the current detected position
   if Position_over(coor,curr_point)==0  %% Check overlapping
        coor=[coor;curr_point];
        Count_coor=Count_coor+1; 
   end
   if (Count_coor >= remnant) || (curr_point(1,1)<=curr_left) || (curr_point(1,1)>=curr_right) || (curr_point(1,2)<=curr_bottom) || (curr_point(1,2)>=curr_top)
       %% The number of feasible positions larger than the remnant sequence, or the detection procedure meet the current boundaries
       Detect_fea=1;
       return;
   end
   head=head+1; % Remove the first element from queue
   west=[curr_point(1,1)-1 curr_point(1,2)];   %% Define the four neighbours
   east=[curr_point(1,1)+1 curr_point(1,2)];
   north=[curr_point(1,1) curr_point(1,2)+1];
   south=[curr_point(1,1) curr_point(1,2)-1];
   
   if Position_over(coor,west)==0  %% Check overlapping
        coor=[coor;west]; 
        Count_coor=Count_coor+1; 
        queue=[queue;west];
        rear=rear+1; % Add the position to the rear of queue
   end
   if Position_over(coor,east)==0  
        coor=[coor;east];  
        Count_coor=Count_coor+1; 
        queue=[queue;east];
        rear=rear+1; 
   end   
   if Position_over(coor,north)==0  
        coor=[coor;north]; 
        Count_coor=Count_coor+1;
        queue=[queue;north];
        rear=rear+1; 
   end   
   if Position_over(coor,south)==0 
        coor=[coor;south]; 
        Count_coor=Count_coor+1; 
        queue=[queue;south];
        rear=rear+1; 
    end
end

function [Over]=Position_over(Pos1,new_point1)
% Check the overlapping 
Over=0;
for i=1:size(Pos1,1)
    if Pos1(i,:)==new_point1
        Over=1;
        return;
    end
end














