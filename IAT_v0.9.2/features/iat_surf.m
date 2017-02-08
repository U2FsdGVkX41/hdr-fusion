function [descriptors, locs] = iat_surf(image,sizeD)
% [DESCS, LOCS] = IAT_SURF(IMAGE, SIZED)
% IAT_SURF executes the SURF (Speeded Uo Robust Features) algorithm [1].
% If the Computer Vision Toolbox is installed, the function calls the
% MATLAB's own implementation. Otherwise, the function calls the executable
% from <iat_dir>/SURF folder, created by the original code provided by the 
% authors of [1]. This is executable in Windows and Linux platforms.
% Note that the two implementations are not identical and may produce
% different results.
%
% -->Input
% IMAGE:                The input image
% SIZED:                The descriptor length (default: 64). While 
%                       Mathworks's implementation accepts only 64 and 128, 
%                       original implementations creates descriptors of 
%                       length 36, 64 and 128.
%
% -->Output
% DESCS:                An NxSIZED array with SURF descriptors (one per row).
% LOCS:                 The keypoints described by DESCS. It is an
%                       array of size Nx5, being each row of the form
%                       [x y a b c], where (x,y) is the keypoint and
%                       [a b;b c] is the ellipse matrix that
%                       describes it's surrounding area.
%
% References
% [1] H.Bay, A.Ess, T.Tuytelaars, L.V.Gool, "SURF: Speeded Up
% Robust Features", Computer Vision and Image Understanding (CVIU),
% Vol. 110, No. 3, pp. 346--359, 2008
%
% -------------------
% Authors: Georgios Evangelidis, Panagiotis Anatolitis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios.evangelidis@inria.fr> or
% <anatolitis@ceid.upatras.gr>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(image,3)>1
    image = rgb2gray(image);
end

if nargin==1
    sizeD=64;
end

if exist(fullfile(matlabroot,'toolbox','vision','vision','detectSURFFeatures.m'),'file') % Mathworks implementation
    if nargin>1
        if ~(sizeD==64 || sizeD==128)
            error('iat_surf: SURF size in Mathworks implementation can be 64 or 128');
        end
    end
    
    points = detectSURFFeatures(image);
    [descriptors, valid_points] = extractFeatures(image, points, 'SURFSize', sizeD);
    x = valid_points.Location(:,1);
    y = valid_points.Location(:,2);
    s = valid_points.Scale; 
    s = s.*s;
    locs = [x y 1./(s) zeros(length(x),1) 1./s];
    
    
else % Bay et al's implementation
    
    
    % Path of executable file
    surfPath=fileparts(mfilename('fullpath'));
    surfPath=[surfPath filesep 'SURF' filesep];
    
    
    [rows, cols] = size(image);
    
    % Convert image into PGM imagefile, readable by executable
    f = fopen([surfPath 'tmp.pgm'],'w');
    pgmFile=[surfPath 'tmp.pgm'];
    surfFile=[surfPath 'img.surf'];
    
    if f == -1
        error('iat_surf: File tmp.pgm cannot be created.');
    end
    fprintf(f, 'P5\n%d\n%d\n255\n', cols, rows);
    fwrite(f, image', 'uint8');
    fclose(f);
    
    % Call the executable file
    if ispc || (isunix&&~ismac)
        
        if ispc
        command=['!' surfPath 'surfWINDLLDemo.exe'];
        else
          command=['!' surfPath 'surf.ln'];
        end
            
        
        if sizeD==64
            sizeCommand=' -in 4';
        elseif sizeD==36;
            sizeCommand=' -in 3';
        elseif sizeD==128;
            sizeCommand=' -in 4 -e';
        else
            error('iat_surf: Unknown descriptor size. Type "help iat_surf"');
        end
        
        % command = [command ' -i tmp.pgm -o img.surf ' sizeCommand];
        command = [command ' -i ' pgmFile ' -o ' surfFile sizeCommand];
       
        %command
        eval(command);
        
        
        g=fopen(surfFile,'r');
        if g == -1
            error('iat_surf: SURF file cannot be accessed.');
        end
        [header, count] = fscanf(g, '%d %d', [1 2]);
        if count ~= 2
            error('iat_surf: Invalid keypoint file beginning.');
        end
        num = header(2);
        len = header(1);
        
        if nargin==1
            if (len ~= (sizeD+1)) % we call the executable without -nl flag
                error(['iat_surf: Unknown descriptor size. Type "help iat_surf"']);
            end
        end
        
        % Create the two output matrices (use known size for efficiency)
        locs = double(zeros(num, 5));
        descriptors = double(zeros(num, sizeD));
        
        % Parse tmp.key
        for i = 1:num
            [vector, count] = fscanf(g, '%f %f %f %f %f', [1 5]); %row col scale ori
            if count ~= 5
                error('iat_surf: Unknown SURF file format');
            end
            locs(i, :) = vector(1, :);
            
            [descrip, count] = fscanf(g, '%f', [1 len]);
            
            if (count ~= len)
                error('iat_surf: Unknown SURF file value.');
            end
            descriptors(i, :) = descrip(2:end);
        end
        fclose(g);
        
    else
        error('iat_surf: appropriate version of Mathworks Computer Vision Toolbox is not installed. Execution of original implementation is possible on Windows/Linux platform');
    end
    
end

if ~isempty(descriptors)
descriptors = double(descriptors);
locs = double(locs);
end
