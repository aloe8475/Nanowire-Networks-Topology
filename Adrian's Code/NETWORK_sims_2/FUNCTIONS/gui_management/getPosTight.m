function pos=getPosTight(outerpos,ti)
%%reduce blank space between actual axes position and space ocuppied
% by the figure, tight inset is the size of axes

left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
pos = [left bottom ax_width ax_height];

end