function [err, wI, dI, e, wI_no_pmodel, wI_no_gmodel] = DIRT_Warp(gmodel, g, ...
    pmodel, p, color, sR, sC, tI, mode, unfilled_pixels_warped, sI, unfilled_pixels_diff)
% [err, wI, dI, e, wI_no_pmodel, wI_no_gmodel] = DIRT_Warp(gmodel, g,
%    pmodel, p, color, sR, sC, tI, mode, unfilled_pixels_warped, sI, unfilled_pixels_diff)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~exist('mode', 'var'), mode = false; end;

% transfer the pixels of interest in sR, sC from the source image to the
% target one
[tR, tC] = DIRT_GModel_Apply(gmodel, g, sR, sC);

if nargout > 4
    wI_no_gmodel = DIRT_PModel_Apply(pmodel, color, p, tI);
end;

err = false;

if ~mode
    % mode == false: means that all pixels of interest should be given a value,
    % and that the warped and error images are vectors
    
    % warp the target image tI
    if ~color
        wI = interp2(tI, tC, tR);
    else
        wI(:,1) = interp2(tI(:,:,1), tC, tR);
        wI(:,2) = interp2(tI(:,:,2), tC, tR);
        wI(:,3) = interp2(tI(:,:,3), tC, tR);
    end;
    % check that all pixels are filled in
    if any(isnan(wI))
        err = true;
%        error('the warp went outside the image');
    end;
    if nargout > 3
        wI_no_pmodel = wI;
    end;
    % apply the photometric registration
    wI = DIRT_PModel_Apply(pmodel, color, p, wI);
    if nargout > 1
        % compute the error image dI
        if ~color
            % in the grey-level case, dI is a (Rn x 1) vector
            dI = wI - sI;
        else
            % in the color case, dI is a (3*Rn x 1) vector
            dI = vect( (wI - sI)' );
        end;
    end;
    if nargout > 2
        % compute the RMS e of the error image dI
        if ~color
            e = norm(dI, 'fro') / sqrt(length(wI));
        else
            e = norm(dI, 'fro') / sqrt(3*length(wI));
        end;
    end;
else
    % mode == true: means that some pixels might be unfilled, and are
    % filled in with a default value, unfilled_pixels_warped in the warped
    % image and unfilled_pixels_diff in the difference image
    
    % warp the target image tI
    if ~color
        wI_tmp = interp2(tI, tC, tR);
    else
        wI_tmp(:,:,1) = interp2(tI(:,:,1), tC, tR);
        wI_tmp(:,:,2) = interp2(tI(:,:,2), tC, tR);
        wI_tmp(:,:,3) = interp2(tI(:,:,3), tC, tR);
    end;
    % find the filled and unfilled pixels
    if ~color
        filledPix_id = find(~isnan(wI_tmp));
        unfilledPix_id = find(isnan(wI_tmp));
    else
        filledPix_id = find(~isnan(wI_tmp(:, :, 1)));
        unfilledPix_id = find(isnan(wI_tmp(:, :, 1)));
    end;
    % copy the filled pixels
    if ~color
        % initialize the warped image with the default value
        wI = ones(size(sI)) * unfilled_pixels_warped;
        wI(filledPix_id) = wI_tmp(filledPix_id);
    else
        for c = 1:3
            % initialize the warped image with the default value
            wI_c = ones(size(sI, 1), size(sI, 2)) * unfilled_pixels_warped(c);
            wI_tmp_c = wI_tmp(:, :, c);
            wI_c(filledPix_id) = wI_tmp_c(filledPix_id);
            wI(:, :, c) = wI_c;
        end;
    end;
    if nargout > 3
        wI_no_pmodel = wI;
        % default value
        if ~color
            wI_no_pmodel(unfilledPix_id) = unfilled_pixels_warped;
        else
            for c = 1:3
                wI_c = wI(:, :, c);
                wI_c(unfilledPix_id) = unfilled_pixels_warped(c);
                wI_no_pmodel(:, :, c) = wI_c;
            end;
        end;
    end;
    % apply the photometric registration
    wI = DIRT_PModel_Apply(pmodel, color, p, wI);
    % default value
    if ~color
        wI(unfilledPix_id) = unfilled_pixels_warped;
    else
        for c = 1:3
            wI_c = wI(:, :, c);
            wI_c(unfilledPix_id) = unfilled_pixels_warped(c);
            wI(:, :, c) = wI_c;
        end;
    end;
    if nargout > 1
        % compute the error image dI
        dI = wI - sI;
        if ~color
            dI(unfilledPix_id) = unfilled_pixels_diff;
        else
            for c = 1:3
                dI_c = dI(:, :, c);
                dI_c(unfilledPix_id) = unfilled_pixels_diff(c);
                dI(:, :, c) = dI_c;
            end;
        end;
    end;
    if nargout > 2
        % compute the RMS e of the error image dI
        if ~color
            e = norm(dI(filledPix_id), 'fro') / sqrt(length(filledPix_id));
        else
            e = norm(dI(filledPix_id), 'fro') / sqrt(3*length(filledPix_id));
        end;
    end;
end;
