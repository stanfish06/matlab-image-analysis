function saveTiff(im, fname, mode)
%SAVETIFF Write a uint16 image or z-stack to a (multipage) TIFF file.
%   SAVETIFF(IM, FNAME, MODE) writes IM to FNAME via the Tiff interface.
%   IM is a 2-D image or 3-D stack (uint16 only). MODE is the Tiff open mode,
%   e.g. 'w' (overwrite) or 'a' (append). Each z-slice IM(:,:,k) is written
%   as a separate TIFF directory (page).

    if isa(im,'uint16')
        nbits = 16;
    else
        error('add other cases');
    end

    t = Tiff(fname,mode);
    
    numrows = size(im,1);
    numcols = size(im,2);
    
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.ImageLength = numrows;
    tagstruct.BitsPerSample = nbits;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.ImageWidth = numcols;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    
    t.setTag(tagstruct);
    t.write(im(:,:,1));
    for i=2:size(im,3)
        t.writeDirectory();
        t.setTag(tagstruct);
        t.write(im(:,:,i));
    end
    t.close();   
end