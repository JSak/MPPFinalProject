clear all; close all; clc; format compact; format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the minimum distance from any point on the film to the airfoil.
% Written by: Jordan Sakakeeny
% Date: September 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Import files
files = dir('snapshot-0.001.dat');

% Get the coordinates for the airfoil
airfoil_thickness = 12;
film_thickness    = 0.003; % m 
y_airfoil = @(x) 0.5*(5*airfoil_thickness/100*(0.2969.*sqrt(2*(x-film_thickness))...
        -0.1260.*(2*(x-film_thickness))...
        -0.3516.*(2*(x-film_thickness)).^2 ...
        +0.2843.*(2*(x-film_thickness)).^3 ...
        -0.1015.*(2*(x-film_thickness)).^4)); 

% Loop through files
for file = files'
    clearvars -except files file airfoil_thickness y_airfoil
    % Import thin film interface data
    filename = file.name;
    timestep_coord = importdata(filename);
    for i = 1:length(timestep_coord)
        if timestep_coord(i,1) < 0.54
            x_interface(i) = timestep_coord(i,1);
            y_interface(i) = timestep_coord(i,2);
        end
    end
%     plot(x_interface,y_interface,'.')
%     pause
    % For every point on the interface, loop through points on the airfoil
    % and get the distance, delta, between them. The airfoil is divided up
    % into as many points as are in the interface, O(5,000). Then, output
    % the minimum delta distance and the coordinates.
    numIter = length(x_interface);
%     x_airfoil = logspace(1E-15,0.3-(1E-15),numIter);

    x_airfoil = zeros(numIter,1);
    x_airfoil(2) = 1E-6;
    for v = 3:2600
        x_airfoil(v) = 1.004*x_airfoil(v-1);
    end
    x_airfoil_back = linspace(x_airfoil(2500),0.253668-1E-12,numIter-2500);
    for v = 2501:numIter
        x_airfoil(v) = x_airfoil_back(v-2500);
    end
    S_0_a = zeros(numIter,1);
    for c = 1:numIter
        if c > 1
                if x_airfoil(c) >= x_airfoil(c-1)
                    S_0_a(c) = S_0_a(c-1)+sqrt((x_airfoil(c)-x_airfoil(c-1))^2+...
                    (y_airfoil(x_airfoil(c))-y_airfoil(x_airfoil(c-1)))^2);
                elseif x_airfoil(c) < x_airfoil(c-1)
                    S_0_a(c) = S_0_a(c-1)-sqrt((x_airfoil(c)-x_airfoil(c-1))^2+...
                    (y_airfoil(x_airfoil(c))-y_airfoil(x_airfoil(c-1)))^2);
                end
        end
    end
    
    for l = 1:numIter %length(x_interface)
        % Divide up the airfoil
%         x_airfoil = linspace(1E-15,0.3-(1E-15),numIter);%length(x_interface));
        min_delta = 1;
%         band = 2700;
%         if l < band
%             lower_bound = 1;
%             if (l + band) > numIter
%                 upper_bound = numIter -1;
%             else
%                 upper_bound = l + band;
%             end
%         elseif l >= (band) 
%             % && l <= (numIter - (band + 1))
%             % lower_bound =  l - band;
%             % upper_bound =  l + band
%             % elseif l > (numIter - (band + 1))
%             lower_bound = numIter - band;
%             upper_bound = numIter - 1;
%         end
        lower_bound = 1;
        upper_bound = numIter;
        
        for m = lower_bound:upper_bound
            
            % Get the distance between one point on the interface and every
            % point on the airfoil
            temp_delta = sqrt((x_airfoil(m)-x_interface(l))^2 + ...
                (y_airfoil(x_airfoil(m))-y_interface(l))^2);
            % Get the minimum distance between the point on the interface
            % and the airfoil
            % The position of the clean airfoil design point is then evaluated as
            % the straight-line distance between the two design points
            if temp_delta <= min_delta
                min_delta = temp_delta;
                s_airfoil_min_temp = S_0_a(m);
            end
%             plot(m,temp_delta,'.r')
%             hold on
%             plot(m,min_delta,'.b')
%             hold on
        end
%         plot(min_delta,'.b')
%         hold on
%         pause
        s_airfoil_min(l) = s_airfoil_min_temp;
        delta(l) = min_delta;
    end
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(s_airfoil_min,delta,'.b')
    xlabel('\bf S (m)')
    ylabel('\bf Thickness, \delta (m)')
    legend('Current Thin Film Interface')
    title(filename)
    daspect([4 0.75 1])
    axis([0 0.26 0 0.005])
    hold on
end
    
    

%% Saves the figures
h=findobj('type','figure') % find the handles of the opened figures
folder='C:\Users\Jordan\Documents\ThinFilmLayout\full_airfoil\Re4E5\Film\run6\SDelta'  % Desination folder
for k=1:numel(h)
  filename=sprintf('image%d.png',k)
  file=fullfile(folder,filename)
  saveas(h(k),file)
end

toc

    