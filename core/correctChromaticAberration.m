function newimg = correctChromaticAberration(img,chan,meta,absshift,scale)
%CORRECTCHROMATICABERRATION Correct axial chromatic aberration of a channel.
%   NEWIMG = CORRECTCHROMATICABERRATION(IMG, CHAN, META, ABSSHIFT, SCALE)
%   resamples the z-stack IMG of channel CHAN to compensate for axial
%   chromatic shift relative to the nuclear channel, applying an axial shift
%   ABSSHIFT (microns, default 3) and scaling SCALE (default 0.85). Only
%   applied when META.chromAbFlag is set, the nuclear channel is 'DAPI', and
%   CHAN is not the nuclear channel; otherwise IMG is returned unchanged.

if ~exist('absshift','var')
    absshift = 3; %absolute z shift in microns
end

if ~exist('scale','var')
    scale = 0.85; %z scaling factor
end

%handle chromatic aberration if data is from the dragonfly 40x water
%objective
nucChannel = meta.nucChannel; chromAbFlag = meta.chromAbFlag; channelLabel = meta.channelLabel;
if chromAbFlag && strcmp(channelLabel{nucChannel + 1},'DAPI') && (chan ~= nucChannel)
    disp('correcting')
    nz = size(img,3);
    
    zres = meta.zres;
    shift = absshift/zres;
    
    fprintf('shift = %.3g slices, scale = %.3g\n',shift, scale)

    newimg = zeros(size(img),class(img));
    for zi = 1:nz
        zq = (zi + shift)*scale;
        zlow = floor(zq); zhigh = ceil(zq);
        wlow = zhigh - zq; whigh = zq - zlow;
        if zlow == zhigh
            if zq >= 1 && zq <= nz
                newimg(:,:,zi) = img(:,:,zq);
            end
        else
            if zq > 1 && zq < nz
                newimg(:,:,zi) = uint16(wlow*double(img(:,:,zlow)) + whigh*double(img(:,:,zhigh)));
            elseif zq > nz && zq < nz + 1
                newimg(:,:,zi) = uint16(wlow*double(img(:,:,zlow)));
            elseif zq > 0 && zq < 1 %don't think this is actually possible
                newimg(:,:,zi) = uint16(whigh*double(img(:,:,zhigh)));
            end
        end
    end
    
else
    disp('nuclear channel; no correction')
    newimg = img;
end





end