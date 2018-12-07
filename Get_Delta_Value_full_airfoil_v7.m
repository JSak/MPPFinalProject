clear all; close all; clc; format compact; format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the minimum distance from any point on the film to the airfoil.
% Written by: Jordan Sakakeeny
% Date: September 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
files = dir('snapshot-0.001.dat');

% Get the coordinates for the airfoil
airfoil_thickness = 12;
film_thickness    = 0.003; % m 
% y_airfoil = @(x) 0.5*(5*airfoil_thickness/100*(0.2969.*sqrt(2*(x-film_thickness))...
%         -0.1260.*(2*(x-film_thickness))...
%         -0.3516.*(2*(x-film_thickness)).^2 ...
%         +0.2843.*(2*(x-film_thickness)).^3 ...
%         -0.1015.*(2*(x-film_thickness)).^4)); 
y_airfoil = @(x) 0.5*(5*airfoil_thickness/100*(0.2969.*sqrt(2*(x-film_thickness))+...
    (((-0.1015.*(2*(x-film_thickness))+...
    0.2843).*(2*(x-film_thickness)) + ...
    -0.3516).*(2*(x-film_thickness)) + ...
    -0.1260).*(2*(x-film_thickness))));

iter = 1;
% Loop throughclearfiles
num_files = length(files');
file_counter = 1;
for file = files'
    fprintf('Load data: ')
%     tic
    clearvars -except files file airfoil_thickness y_airfoil iter num_files file_counter s_airfoil_min delta S_0_a %S_0_a_total 
    % Import thin film interface data
    filename = file.name;
    timestep_coord = importdata(filename);
    x_interface = timestep_coord(:,1);
    y_interface = timestep_coord(:,2);

    % For every point on the interface, loop through points on the airfoil
    % and get the distance, delta, between them. The airfoil is divided up
    % into as many points as are in the interface, O(5,000). Then, output
    % the minimum delta distance and the coordinates.
    numIter = length(x_interface);
    x_airfoil = gpuArray.linspace(0.003,2*0.253668-1E-12,numIter);
    y_airfoil_calc = y_airfoil(x_airfoil);
%     x_mid_idx = length(x_airfoil)/2;
%     toc
    fprintf('Calculate S_0_a: ')
%     tic
    S_0_a(1,file_counter) = 0.0;
    S_0_a_gpu = gpuArray(S_0_a(:,file_counter));
    s_airfoil_min(:,file_counter) = 0.0;
    s_airfoil_min_gpu = gpuArray(s_airfoil_min(:,file_counter));
    delta(:,file_counter) = 0.0;
    delta_gpu = gpuArray(delta(:,file_counter));
    for c = 1:numIter
        if c > 1
                if x_airfoil(c) >= x_airfoil(c-1)
                    S_0_a_gpu(c) = S_0_a_gpu(c-1)+sqrt((x_airfoil(c)-x_airfoil(c-1))^2+...
                    (y_airfoil(x_airfoil(c))-y_airfoil(x_airfoil(c-1)))^2);
                elseif x_airfoil(c) < x_airfoil(c-1)
                    S_0_a_gpu(c) = S_0_a_gpu(c-1)-sqrt((x_airfoil(c)-x_airfoil(c-1))^2+...
                    (y_airfoil(x_airfoil(c))-y_airfoil(x_airfoil(c-1)))^2);
                end
        end
    end
%     toc
    fprintf('Main loop: \n')
%     tic
%     parfor l = 1:numIter 
        
    % TEST SECTION
%     s_airfoil_min       = zeros(numIter,num_files);
%     delta               = zeros(numIter,num_files);
%     x_airfoil_gpu       = gpuArray(x_airfoil);
%     y_airfoil_calc_gpu  = gpuArray(y_airfoil_calc);
%     x_interface_gpu     = gpuArray(x_interface);
%     y_interface_gpu     = gpuArray(y_interface);
%     S_0_a_gpu           = gpuArray(S_0_a(:,file_counter));
%     s_airfoil_min_gpu   = gpuArray(s_airfoil_min(:,file_counter));
%     delta_gpu           = gpuArray(delta(:,file_counter));
    
    parfor l = 1:numIter 
        min_delta = 1;
        lower_bound = 1;
        upper_bound = numIter; 
        
        % Vectorize the for loop and it's much, much faster!
        temp_delta = sqrt((x_airfoil-x_interface(l)).^2 + ...
            (y_airfoil_calc-y_interface(l)).^2);
        [min_delta,min_idx] = min(real(temp_delta(:)));
        s_airfoil_min_gpu(l) = S_0_a_gpu(min_idx);
        delta_gpu(l) = min_delta;
%         s_airfoil_min(l,file_counter) = S_0_a_gpu(min_idx,file_counter);
%         delta(l,file_counter) = min_delta;
    end
    S_0_a(1:length(S_0_a_gpu),file_counter) = zeros(length(S_0_a_gpu),1);
    s_airfoil_min(1:length(s_airfoil_min_gpu'),file_counter) = zeros(length(s_airfoil_min_gpu'),1);
    delta(1:length(delta_gpu'),file_counter) = zeros(length(delta_gpu'),1);
    
    S_0_a(:,file_counter) = gather(S_0_a_gpu(:));
    s_airfoil_min(:,file_counter) = gather(s_airfoil_min_gpu(:));
    delta(:,file_counter) = gather(delta_gpu(:));
    
    file_counter = file_counter + 1;
end

for gg = 1:num_files
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(s_airfoil_min(:,gg),delta(:,gg),'.b')
    xlabel('\bf S (m)')
    ylabel('\bf Thickness, \delta (m)')
    legend('Current Thin Film Interface')
    title(files(gg).name)
    daspect([4 0.75 1])
    axis([0 0.54 0 0.01])
    hold on
end

toc
fprintf('End of main loop. \n')


%% Saves the figures
% tic
% h=findobj('type','figure') % find the handles of the opened figures
% folder='C:\Users\Jordan\Documents\ThinFilmLayout\full_airfoil\Re4E5\Film\run6\SDelta'  % Desination folder
% for k=numel(h):-1:1
%   filename=sprintf('image%d.png',k)
%   file=fullfile(folder,filename)
%   saveas(h(k),file)
% end
% 
% toc
% 

    