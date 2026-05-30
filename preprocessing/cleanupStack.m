function cleanupStack(dataDir, filename, zcorrection, zcorrfactor, mediansize)
    % background subtraction, intensity correction and median filter per z-slice

    % piece for filenames from micropatternPieVis
    [~,barename,ext] = fileparts(filename);

    % self stitched files start with stitched
    if strcmp(filename(1:8), 'stitched') || strcmp(filename(1:10), 'FFstitched')

        s = strsplit(filename,{'_'});
        prefix = s{1}; %'stitched';
        postfix = [s{2} '_w%.4d_t0000.tif'];
        fnameFormat = [prefix '_' postfix];
        cleanfnameFormat = [prefix '_filtered_' postfix];
        cleanMIPfnameFormat = [prefix '_filtered_MIP_' postfix];
        bioformats_flag = false;
        
    elseif strcmp(ext,'.nd2')
        prefix = barename;
        fnameFormat = [prefix,'_MIP_w%.4d.tif'];
        cleanfnameFormat = [prefix '_filtered_w%.4d_t0000.tif'];
        cleanMIPfnameFormat = [prefix '_filtered_MIP_w%.4d.tif'];
        bioformats_flag = true;
        
    elseif strcmp(ext, '.ims')
        
        s = strsplit(filename,{'_FusionStitcher','.ims'});
        prefix = [s{:}];
        fnameFormat = [prefix,'_MIP_w%.4d.tif'];
        cleanfnameFormat = [prefix '_filtered_w%.4d_t0000.tif'];
        cleanMIPfnameFormat = [prefix '_filtered_MIP_w%.4d.tif'];
        bioformats_flag = true;
        
%     % previous case - DOCUMENT WHAT THIS IS EXPECTING
%     else
%         s = strsplit(filename,{'_Stitched','.tif'});
%         prefix = [s{:}];
%         subDir = [prefix '_zslices'];
%         fnameFormat = [prefix '_MIP_w%.4d.tif'];
    end
    
    % make directory for output
    if ~exist(fullfile(dataDir,'filtered'),'dir')
        mkdir(fullfile(dataDir,'filtered'));
    end
    
    if ~exist('zcorrection','var')
        zcorrection = false;
    end
    if ~exist('zcorrfactor','var')
        zcorrfactor = [1 1 1 1]*1.05;
    end
    if ~exist('mediansize','var')
        mediansize = [5 5];
    end
    
    for channelIdx = 0:3 
        
        if ~bioformats_flag
            fname = sprintf(fnameFormat, channelIdx);
            fullfname = fullfile(dataDir,fname);
        else
            fullfname = fullfile(dataDir,filename);
        end
        
        if exist(fullfname, 'file')
            
            if bioformats_flag
                img = readStack(fullfname); %xyzct image stack
                img = squeeze(img(:,:,channelIdx+1,:,:));
            else
                img = readStack(fullfname);
                img=squeeze(img);
            end

            medstack = 0*img;
%             medstack = zeros(size(img),class(img));
            for zidx = 1:size(img,3)

                % bg subtraction
                if zidx == 1
                    mode = 'overwrite';
                else
                    mode = 'append';
                end
                slice = img(:,:,zidx); 
                slicebg = imopen(slice,strel('disk',50));
                slicebgsub = slice-slicebg;

                % correct intensity
                if zcorrection
                    slicebgsub = slicebgsub*(zcorrfactor(channelIdx+1)^(zidx-1));
                end

                %openstack(:,:,zidx) = imopen(slicebgsub,strel('disk',2));
                medstack(:,:,zidx) = medfilt2(slicebgsub, mediansize);
                tmp = medstack(:,:,zidx);

                cleanfname = sprintf(cleanfnameFormat, channelIdx);
                imwrite(tmp,fullfile(dataDir,'filtered',cleanfname),'WriteMode',mode);
            end

            medMIP = max(medstack,[],3);
            cleanMIPfname = sprintf(cleanMIPfnameFormat, channelIdx);
            imwrite(medMIP,fullfile(dataDir,'filtered',cleanMIPfname));
        
        else
            warning([fullfname 'does not exist']);
        end
    end
end
