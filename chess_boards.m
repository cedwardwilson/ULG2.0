%% Script for making 2D chess-board style geometries

X = 16   % Note that this will give locations in terms of indices not 
Y = 16;   % spatial coordinates, therefore lambdas will only be satisfied 
           % if resolution is same in all directions 

expnr = '206';
area = X*Y;
n = 2; % This must divide both X and Y
special = true;
ibl = true;

%% Special case 50:50
if special == true
    % Green
    xl = 0;
    xu = X;
    yl = 0;
    yu = Y/2;
    gyl = Y/2;
    gyu = Y;
    gzl = 1;
    gzu = 2;
    if ibl == true
        buildings = [gyl,gyu,xl,xu,gzl,gzu; yl,yu,xl,xu,1,1];
        green_outline = [gyl,gyu,xl,xu];
    else 
        buildings = [xl,xu,gyl,gyu,gzl,gzu; xl,xu,yl,yu,1,1];
        green_outline = [xl,xu,gyl,gyu];
    end
else

    %% General case, square board
    
    divisor = Y/n;
    xpoints = 0:divisor:X;
    ypoints = 0:divisor:Y;
    buildings = [];
    green_outline = [];
    
    for i=1:length(xpoints)-1
        for j = 1:length(ypoints)-1
            if mod(i,2) == mod(j,2)
                buildings = [buildings; xpoints(i), xpoints(i+1), ypoints(j), ypoints(j+1), 1,2];
                green_outline = [green_outline,; xpoints(i), xpoints(i+1), ypoints(j), ypoints(j+1)];
            else 
                buildings = [buildings; xpoints(i), xpoints(i+1), ypoints(j), ypoints(j+1), 1,1];
            end
        end
    end
end

%% Try transforming the coordinates to see if preprocessing params can match
% Here we have to add to the x and y lower coords and take away from the z
% upper coord such that the preprocessing recreates the geom we want.
% Note that these corrections are very different to those needed to convert
% from the geometries made in Python. I expect this is because Python
% indexes from 0 not 1...

buildings(:,1) = buildings(:,1)+1;
buildings(:,3) = buildings(:,3)+1;
green_outline(:,1) = green_outline(:,1)+1;
green_outline(:,3) = green_outline(:,3)+1;
%buildings(:,6) = buildings(:,6)-1;
%% Write out the geometries for uDALES to use
blocks_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/buildings.' expnr '.mat'];
outline_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/green_oultine.' expnr '.mat'];
params_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/generated_params.' expnr '.mat'];
area_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/area.' expnr '.mat'];
save(blocks_file_path, 'buildings');
save(outline_file_path, 'green_outline');
save(area_file_path,'area');

