function plotMovie(x,y,xx,frameRate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compiles and saves a movie in which a marker traces along a 
% plot.
%
% ARGUMENTS: 
% x,y - point in the graph.
% xx  - entries in x that get a frame (values, not indices).
% frameRate - movie's frame rate.
% 
% OUTPUT:
% none.
% 
% REQUIRES:
% none. 
%
% USAGE:
%{
    plotMovie(t,conductance,timestamps,10);
%}
%
% Author:
% Ido Marcus  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    fprintf('\nCompiling movie...\n');
    
    v = VideoWriter('plotMovie.mp4','MPEG-4');
    v.FrameRate = frameRate;
    v.Quality = 100;
        
    open(v);
    for i = 1 : length(xx)
        progressBar(i,length(xx));
    
        % Produce the figure (every time, inefficiently):
        frameFig = figure('visible','off');
        plot(x,y);
        hold on;
        scatter(xx(i),y(x==xx(i)),150,'red','filled','d'); 
        %title('Conductance as function of time');
        %xlabel('Time (sec)');
        %ylabel('Conductance (S)');
        grid on;
        
        % Save it:
        writeVideo(v,getframe(frameFig));
        
        close(frameFig);
    end
    close(v);
    
    fprintf('\nDone.\n');
end
