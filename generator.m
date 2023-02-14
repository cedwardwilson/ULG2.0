%% ULG that makes random vegetated urban areas, this code struggles for lp+lv>0.3
clear all

%% Define the input parameters 
X = 32;                    % Domain length
Y  = 32;                   % Domain width
lp = 0.25;                   % Built fraction
lv = 0.1;                   % Vegetated fraction
lf = 0.2;                  % Frontal aspect ratio
minx = 5;                   % Narrowest a building can be in x
miny = 5;                   % Narrowest a building can be in y
minxg = 2;                  % Narrowest a green can be in x
minyg = 2;                  % Narrowest a green can be in y
minz = 5;                   % Minimum building height
maxz = 10;                  % Max building height
nbuild = 4;                % Number of buildings
ngreen = 2;                 % Number of green spaces
max_chunk = 100;            % Limits how much the frontal surface area of building can increase at any one time
paddingb = 1;               % Minimum gap between buildings
paddingg = 1;               % Minimum gap between green space and anything else
expnr = '253';
tiled = true;
xtiles = 3;
ytiles = 2;

area   = X*Y;
builtarea = lp*area;
greenarea = lv*area;
front = lf*area;

%% Make the buildings
% Make the original building block
builtwidth = floor(sqrt(builtarea));
builtlength = floor(builtarea/builtwidth);
builtblocks = [0,builtwidth,0,builtlength,0];

% Cut into sub-blocks until we have as many as wanted
numbuilt = 1;
while numbuilt < nbuild
    builtblocks = cutfromlist(builtblocks, minx,miny);
    numbuilt = numbuilt + 1;
end    

%% Give the buildings height
%heights = heightadder1(builtblocks,minz,maxz,front,max_chunk)
heights = heightadder2(builtblocks,minz,maxz,front)

%% Check the frontal ratio is satisfied
total_front = sum(heights.*builtblocks(:,4));
ratio = total_front/front;
final_heights = round(heights./ratio);
new_ratio = sum(final_heights.*builtblocks(:,4))/front; % Rescale if need be

%% Add the building heights to the block list that describes them 
builtblocks(:,5) = final_heights+1;   % Add 1 to heights because the buildings start at z = 1 not z = 0

%% Make green patches
gnwidth = floor(sqrt(greenarea));
gnlength = floor(greenarea/gnwidth);
ptchs = [0,gnwidth,0,gnlength];
numptch = 1;
while numptch < ngreen
    ptchs = cutfromlist(ptchs, minxg,minyg);
    numptch = numptch + 1;
end    

ptchs_tall = ones(length(ptchs(:,1)),5);  %Prepare list in which green has height
ptchs_tall(:,1:4) = ptchs(:,1:4); 

%% Placing the buildings and green space randomly
topo = zeros(Y, X);
first_topo = topo;       % Originally no blocks have been added
first_list = [];         % Originally no blocks have been added

% Add the buildings
[second_list, second_topo,~] = placer(builtblocks,paddingb,X,Y,first_topo,first_list);

% Add the green around the buildings
[final_list, final_topo, green_outline] = placer(ptchs_tall,paddingg, X,Y,second_topo,second_list);


%% Make green areas the tallest to be consistent with uDALES preprocessing 
maxh = max(final_list(:,6));
final_list(final_list(:,6) == 1,6) = maxh+1; 
buildings = final_list;

%% Determine the actual lambda values we have achieved
lams = lambcalc(buildings,area);

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
buildings(:,6) = buildings(:,6)-1;

%% Make repeated tiles if desired
if tiled == true
    repeated_buildings = [];
    repeated_green_outline = [];
    for xadd = 1:xtiles
        for yadd = 1:ytiles
            for b = 1:length(buildings(:,1))
                temp = [];
                new_xl = buildings(b,1)+(xadd-1)*X;
                new_xu = buildings(b,2)+(xadd-1)*X;
                new_yl = buildings(b,3)+(yadd-1)*Y;
                new_yu = buildings(b,4)+(yadd-1)*Y;
                temp = [new_xl, new_xu, new_yl,new_yu, buildings(b,5), buildings(b,6)];
                repeated_buildings = [repeated_buildings; temp];
            end 
            for g = 1:length(green_outline(:,1))
                temp = [];
                new_xl = green_outline(g,1)+(xadd-1)*X;
                new_xu = green_outline(g,2)+(xadd-1)*X;
                new_yl = green_outline(g,3)+(yadd-1)*Y;
                new_yu = green_outline(g,4)+(yadd-1)*Y;
                temp = [new_xl, new_xu, new_yl,new_yu];
                repeated_green_outline = [repeated_green_outline; temp];    
            end
        end
    end  
    buildings = repeated_buildings;
    green_outline = repeated_green_outline;
end     

%% Write out the geometries for uDALES to use
blocks_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/buildings.' expnr '.mat'];
outline_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/green_oultine.' expnr '.mat'];
params_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/generated_params.' expnr '.mat'];
area_file_path = ['/media/chris/Project1/uDALES_veg/experiments/' expnr '/area.' expnr '.mat'];
save(blocks_file_path, 'buildings');
save(outline_file_path, 'green_outline');
save(params_file_path,'lams');
save(area_file_path,'area');

function newblocks = cut(block, minx, miny) % Cuts a given block randomly i
% nto 2 rectangles with minimum width and length
    disp('cutting')
    xl = 0;
    xu = block(2);
    yl = 0;
    yu = block(4);

    if xu < 2*minx & yu < 2*miny
        disp('too thin');
    elseif yu < 2*miny
        disp('too thin in y, cut along constant x');
        cutdir = 1;
    elseif xu < 2*minx    
        disp('too thin in x, force cut along constant y');
        cutdir = 2
    else
        cutdir = randi(2);
    end 

    if cutdir == 1
        disp('cut along constant x, cutting in y direction')
        cutline = randi([minx,xu-minx]);
        newblock1 = [0,cutline,0,yu];
        newblock2 = [0,xu-cutline, 0,yu];
    elseif cutdir == 2
        disp('cut along constant y, cutting in x direction')
        cutline = randi([miny,yu-miny]);
        newblock1 = [0,xu,0,cutline];
        newblock2 = [0,xu, 0,yu- cutline];
    else 
        disp('too thing to cut')
    end 
    newblocks = [newblock1; newblock2];
end    

function newlist = cutfromlist(list, minx,miny) % Takes a list of blocks,
% picks one and cuts it, picks new block if current is too small to cut 
    n = size(list,1);
    cutind = randi([1,n]);
    block2cut = list(cutind,:);
    xu = block2cut(2);
    yu = block2cut(4);
    while xu < 2*minx & yu<2*miny
        disp('block too small to cut in x or y')
        cutind = randi([1,n]);
        block2cut = list(cutind,:);
        xu = block2cut(2);
        yu = block2cut(4);
    end    
    list(cutind,:) = [];
    newblocks = cut(block2cut, minx,miny);
    newlist = [list;newblocks];
end    

function [new_list, new_topo, outline] = placer(blocks,padding,X,Y,old_topo,old_list)
% Takes a set of blocks, and places them randomly on a pre-existing
% arrangement, each block being no closer to another than a padding
% distance. It then returns a list consisting of the original blocks and
% the new ones that have been added (the list is in [xl,xu,yl,yu,zl,zu;
% ...] form. It also returns a new topo with all the blocks in it as well
% as an outline, giving the locations of the edges of the new blocks that
% were added.
    blocks(:,1:4) = blocks(:,1:4)+1; %To avoid indexing from zero issues
    empty = zeros(Y,X);
    outline = [];
    i = 1;
    topo = old_topo; 
    while i <= length(blocks(:,1))
        attempts = 0;
        blk = blocks(i,:);
        xlmax = X - (blk(2)+padding);
        ylmax = Y - (blk(4)+padding);
        fit = false;
        while fit == false
            attempts = attempts + 1;
            if attempts > 10000
                i = 1;
                topo = old_topo;
                outline = [];
                break
            end
            disp([attempts]);
            disp(i);
            trial = empty;
            xlmax 
            ylmax
            xadd = randi([padding,xlmax]);
            yadd = randi([padding,ylmax]);
            shft(1) = blk(1)+xadd;
            shft(2) = blk(2)+xadd;
            shft(3) = blk(3)+yadd;
            shft(4) = blk(4)+yadd;
            trial(shft(3):shft(4),shft(1):shft(2)) = blk(5);
            padded = empty;
            padded(shft(3)-padding:shft(4)+padding,shft(1)-padding:shft(2)+padding) = blk(5);
            if (padded.*topo == 0)
                topo = topo + trial;
                old_list = [old_list;shft(1),shft(2),shft(3), shft(4), 1, blk(5)];
                fit = true;
                outline = [outline; shft(1),shft(2),shft(3), shft(4)];
                i = i+1;
            end
        end
    end 
    new_list = old_list;
    new_topo = topo;
end

function lams = lambcalc(blocks,area)
    buildings = blocks(blocks(:,6) ~= max(blocks(:,6)),:); 
    builtareas = (buildings(:,2)-buildings(:,1)).*(buildings(:,4)-buildings(:,3));
    builtarea = sum(builtareas);
    lamp = builtarea/area;

    frontareas = (buildings(:,4)-buildings(:,3)).*(buildings(:,6)-buildings(:,5));
    frontarea = sum(frontareas);
    lamf = frontarea/area;

    greens = blocks(blocks(:,6) == max(blocks(:,6)),:);
    greenareas = (greens(:,2)-greens(:,1)).*(greens(:,4)-greens(:,3));
    greenarea = sum(greenareas);
    lamv = greenarea/area;

    lams = [lamp,lamf,lamv];
end

function new_heights = heightadder1(blocks, minz,maxz, front,max_chunk)
% Here we pick a random area to add then a random building to add it to
    side_lens = blocks(:,4);
    %tot_len = sum(side_lens);
    min_areas = side_lens*minz;                     % Frontal areas based on minimum heights
    remaining = front - sum(min_areas);             % How much area we need to add until lf is satisfied
    heights = minz*ones(length(blocks(:,1)),1);     % Each block is initially its minimum height
    while remaining > 0
        chunk = randi(max_chunk);                   % An amount of area we will take of remaining and add to a building   
        blockind = randi(length(blocks(:,1)));             % The block to add to
        extra_height = floor(chunk/blocks(blockind,4));    % The height we will add to that block
        if heights(blockind)+extra_height <= maxz               % Make sure block doesn't get too tall
            heights(blockind) = heights(blockind)+extra_height;
            remaining = remaining - chunk;
        end    
    end
    new_heights = heights;
end

function new_heights = heightadder2(blocks, minz,maxz, front)
% Here we pick a random building then add a random height to it. 
    side_lens = blocks(:,4);
    %tot_len = sum(side_lens);
    min_areas = side_lens*minz;                     % Frontal areas based on minimum heights
    remaining = front - sum(min_areas);             % How much area we need to add until lf is satisfied
    heights = minz*ones(length(blocks(:,1)),1);     % Each block is initially its minimum height
    blockind_list = 1:length(blocks(:,4));
    while remaining > 0
        if length(blockind_list) == 0
            disp('No more blocks we can add to.')
            break
        end     
        if remaining < min(blocks(:,4))
            disp('Insufficient area remaining to add 1m to smallest building.')
            break
        end    
        blockind_ind = randi(length(blockind_list));
        blockind = blockind_list(blockind_ind);

        max_add = maxz-heights(blockind);
        if max_add <= 0 
            disp('block max height')
            blockind_list(blockind_list == blockind) = [];
            continue
        end
        blockwidth = blocks(blockind,4);
        max_add = min(max_add,floor(remaining/blockwidth));
        if max_add <=0
            disp('Not enough remaining to add to this block')
            blockind_list(blockind_list == blockind) = [];
            continue
        end    
        height_add = randi(max_add);
        chunk = height_add*blockwidth;
        heights(blockind) = heights(blockind)+height_add;
        remaining = remaining - chunk;
    end
    new_heights = heights
end