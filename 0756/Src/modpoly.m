function [w,beta,indx] = modpoly(w,beta)
%MODPOLY Modify a polygon.
%	[WNEW,BETANEW] = MODPOLY(W,BETA) plots the polygon given by W
%	and BETA and allows the user to change it with the mouse.  At
%	the start, MODPOLY allows you to move vertices.  Move the cursor
%	over a vertex you want to move, hold down the left mouse button,
%	drag it to its new location, and release.  The vertex changes
%	color and affected sides become dashed as you move the vertex.
%	
%	To delete a vertex, press the Delete button.  The pointer
%	will change to a fleur.  After your next click and release, the
%	selected vertex will be deleted and you will return to movement
%	mode.  To cancel a requested deletion, press Delete again.
%	
%	The Add button works similarly.  To add, press the button and
%	then click and release on a polygon side.  A vertex will be
%	added to the middle of the side, the polygon is redrawn, and you
%	return to movement mode.
%	
%	Infinite vertices cannot be moved(!), deleted, or added.  When
%	moving the neighbor of an infinite vertex, the angle at infinity
%	is kept constant.  When you delete a neighbor of infinity, the
%	turn at the deleted vertex is lost and the angle at infinity
%	changes.  You cannot delete a vertex with two infinite
%	neighbors.  When you add a vertex to an infinite side, the new
%	vertex appears at a "reasonable" distance from its finite
%	neighbor.
%	
%       [WNEW,BETANEW,IDX] = MODPOLY(W,BETA) also returns an index
%       vector to help keep track of additions and deletions.  IDX has
%       the same length as WNEW, and if IDX(J) is an integer, it gives
%       the index that WNEW(J) had in the original W.  If WNEW(J) was
%       added, then IDX(J) is NaN.
%	
%	Note: MODPOLY makes no attempt to keep the polygon "legal."  You
%	can easily create things which are not polygons, or change
%	infinite vertices into unrecognized finite ones.  
%	
%	See also DRAWPOLY, PLOTPOLY.
%	
%	Written by Toby Driscoll.  Last updated 6/1/95.

global sc_hs sc_hv sc_k sc_w sc_beta sc_idx
n = length(sc_w);
ptr = get(gcf,'pointer');

if ~isstr(w)				% initial call
  % Draw polygon and initialize global vars
  sc_w = w(:);				% vertices
  sc_beta = beta(:);			% angles
  hold off
  sc_hs = plotpoly(w,beta);		% side handles
  hold on
  sc_hv = zeros(n,1);			% vertex handles
  for j=find(~isinf(w))'
    sc_hv(j) = plot(real(w(j)),imag(w(j)),'.','mark',22);
  end
  sc_idx = (1:length(w))';		% indices
  oldptr = ptr;
  set(gcf,'pointer','circle')
  % Create uicontrols
  pb_done = uicontrol('style','push','string','Done','pos',[5,5,60,22],...
      'call','set(get(gcf,''currentobj''),''user'',1)');
  set(pb_done,'user',0);
  pb_del = uicontrol('style','push','string','Delete','pos',[5,30,60,22],...
      'call','modpoly(''delete'');');
  pb_add = uicontrol('style','push','string','Add','pos',[5,55,60,22],...
      'call','modpoly(''add'');');
  set(gcf,'windowbuttondown','modpoly(''down'');')
  % Run in place until finished
  while ~get(pb_done,'user')
    drawnow
  end
  set(gcf,'windowbuttondown','')
  % Recover new info and clean up
  w = sc_w;
  beta = sc_beta;
  hold off
  delete(pb_done)
  delete(pb_del)
  delete(pb_add)
  set(gcf,'pointer',oldptr)
  plotpoly(w,beta)
  indx = sc_idx;
  clear sc_hs sc_hv sc_k sc_w sc_beta sc_idx

elseif strcmp(w,'down')			% button down
  h = get(gcf,'currentobj');
  % Act only if h is a line object
  if strcmp(get(h,'type'),'line')
    if ~strcmp(ptr,'crosshair')		% move or delete
      if strcmp(get(h,'linesty'),'.')	% vertex?
	sc_k = find(h==sc_hv);
	colr = get(gca,'colororder');
	set(h,'color',colr(2,:))
	set(sc_hs([sc_k,rem(sc_k-2+n,n)+1]),'linesty','--')
	if strcmp(ptr,'circle')
	  % Mouse movement needed only when moving vertices
	  set(gcf,'windowbuttonmotion','modpoly(''move'');')
	end
	set(gcf,'windowbuttonup','modpoly(''up'');')
      end
    else				% insert
      if strcmp(get(h,'linesty'),'-')	% edge?
	sc_k = find(h==sc_hs);
	set(sc_hs(sc_k),'linesty','--')
	set(gcf,'windowbuttonup','modpoly(''up'');')
      end	
    end
  end

elseif strcmp(w,'move')			% mouse move
  z = get(gca,'currentpoint');
  k = sc_k;
  set(sc_hv(k),'xd',z(1,1),'yd',z(1,2))
  % Must handle case of infinite predecessor/successor separately.
  j = rem(k,n)+1;			% successor
  if isinf(sc_w(j))
    xd = get(sc_hs(k),'xd');
    yd = get(sc_hs(k),'yd');
    phi = atan2(diff(yd),diff(xd));
    r = sqrt(diff(xd)^2+diff(yd)^2);
    y = sc_w(k) + [0,r*exp(i*phi)];
  else
    y = [z(1,1)+i*z(1,2),sc_w(j)];
    phi = angle(-diff(y)/diff(sc_w([j,k])));
    sc_beta(k) = sc_beta(k)-phi/pi;
    sc_beta(j) = sc_beta(j)+phi/pi;
  end
  set(sc_hs(k),'xd',real(y),'yd',imag(y))
  j = rem(k-2+n,n)+1;			% predecessor
  if isinf(sc_w(j))
    xd = get(sc_hs(j),'xd');
    yd = get(sc_hs(j),'yd');
    phi = atan2(-diff(yd),-diff(xd));
    r = sqrt(diff(xd)^2+diff(yd)^2);
    y = sc_w(k) + [r*exp(i*phi),0];
  else
    y = [sc_w(j),z(1,1)+i*z(1,2)];
    phi = angle(diff(y)/diff(sc_w([j,k])));
    sc_beta(k) = sc_beta(k)+phi/pi;
    sc_beta(j) = sc_beta(j)-phi/pi;
  end
  % Make change effective
  set(sc_hs(j),'xd',real(y),'yd',imag(y))
  drawnow
  sc_w(k) = z(1,1)+i*z(1,2);

elseif strcmp(w,'up')			% button up
  set(sc_hs([sc_k,rem(sc_k-2+n,n)+1]),'linesty','-')
  set(gcf,'windowbuttonup','')
  if strcmp(ptr,'circle')
    % Moved a vertex.  Just clean up.
    colr = get(gca,'colororder');
    set(sc_hv(sc_k),'color',colr(1,:))
    set(gcf,'windowbuttonmotion','')
  elseif strcmp(ptr,'fleur') 		% Delete...
    colr = get(gca,'colororder');
    set(sc_hv(sc_k),'color',colr(1,:))
    set(gcf,'pointer','circle')
    idx = rem(sc_k+(-2:1)+n-1,n)+1;	% 2 back, here, and 1 forward
    infb = isinf(sc_w(idx(2)));
    infa = isinf(sc_w(idx(4)));
    if n <= 3 | (infa & infb)
      return				% do nothing
    elseif ~infb & ~infa 
      % Finite neighborhs; easy.
      v = get(sc_hs(idx(1)),'xdata')+i*get(sc_hs(idx(1)),'ydata');
      v(3:4) = get(sc_hs(idx(4)),'xdata')+i*get(sc_hs(idx(4)),'ydata');
      b = scangle(v);
      sc_beta(idx([2,4])) = b(2:3);
      v = v(2:3);
    else
      % An infinite neighbor
      axlim = axis;
      r = sqrt(diff(axlim(1:2))^2+diff(axlim(3:4))^2);
      x = get(sc_hs(sc_k*infb + idx(2)*infa),'xdata');
      y = get(sc_hs(sc_k*infb + idx(2)*infa),'ydata');
      ang = atan2(diff(y),diff(x)) + (pi*infb);
      j = (idx(2)*infb) + (idx(4)*infa);
      sc_beta(j) = sc_beta(j) + sc_beta(sc_k);
      if infb
	v = sc_w(idx(4)) + [1.1*r*exp(i*ang),0];
      else
	v = sc_w(idx(2)) + [0,1.1*r*exp(i*ang)];
      end
    end
      
    set(sc_hs(idx(2)),...
	'xdata',real(v),'ydata',imag(v),'linesty','-')
    delete(sc_hv(sc_k))
    sc_hv(sc_k) = [];
    sc_w(sc_k) = [];
    delete(sc_hs(sc_k))
    sc_hs(sc_k) = [];
    sc_beta(sc_k) = [];
    sc_idx(sc_k) = [];
   
  elseif strcmp(ptr,'crosshair') 	% Add...
    [wn,bn] = scaddvtx(sc_w,sc_beta,sc_k);
    sc_w = wn(:);
    sc_beta = bn(:);
    hold off
    sc_hs = plotpoly(wn,bn);
    hold on
    sc_hv = zeros(n+1,1);
    for j=find(~isinf(wn))'
      sc_hv(j) = plot(real(wn(j)),imag(wn(j)),'.','mark',22);
    end
    sc_idx = [sc_idx(1:sc_k);NaN;sc_idx(sc_k+1:n)];
    set(gcf,'pointer','circle')
    
  end

elseif strcmp(w,'delete') 		% toggle delete state
  if ~strcmp(ptr,'fleur')
    set(gcf,'pointer','fleur')
  else
    set(gcf,'pointer','circle')
  end

elseif strcmp(w,'add')			% toggle add state
  if ~strcmp(ptr,'crosshair')
    set(gcf,'pointer','crosshair')
  else
    set(gcf,'pointer','circle')
  end

end
