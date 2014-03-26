function spt()

%% Set Default Behavior
INT = 1;
%  sptint = 0;
ASK = 1;
% sptpara.f_rate = 
% sptpara.pixl =
% sptpara.mtl =
% sptpara.thresh =
% sptpara.dia =
% sptpara.max_disp =
% sptpara.npts =




%% Get Filenames
[filename,filepath]=uigetfile({'*.tif', 'Tiff File (*.tif)';...
                               '*.stk', 'Metamorph Stack (*.stk)'},... 
                               'Choose Movie(s) to Analyze',...
                               'MultiSelect','on');
cd(filepath);

if iscell(filename)
    ntrials = size(filename,2);
    stackn = cell(1,ntrials);
    for i = 1:ntrials
        stackn{i} = filename{i};
    end

else
    ntrials = 1;
    stackn{1} = filename;
end

%% Setup
for i = 1:ntrials

    % Calculation Parameters
    if ASK == 1

        prompt = {'Frame Rate (hz):',...
            'Pixel Length ({\mu}m):',...
            'Minimum Trajectory Length:',...
            'Particle Threshold:',...
            'Integrated Threshold [optional]:',...
            'Particle Diameter(pixel) [must be odd]:',...
            ['Est. Diffusion Constant ({\mu}m^2/sec)',...
               ' [used to calc max particle displacment]'],...
            'Center Finding Method [0 for 2D Gaussian, 1 for centroid]',...
            ['Window Size (pixel) [Radius of sub-pixel determination ',...
               'window size, Particle Diameter + 4 is a good value. ',...
               'Must be odd!'],['Save filtered image? [0 for no, 1 for',...
               ' in .mat file, 2 for .tiff, 3 for both]'],...
               ['Number of frames a particle is allowed to blink off']}; 
        u_name = 'Input Parameters for Single Particle Tracking';
        numlines = 1;
        defaultanswer = {'','.16','4','','0','3','','1','9','2','0'};
        options.Resize = 'on';
        options.WindowStyle = 'normal';
        options.Interpreter = 'tex';
        user_var = inputdlg(prompt,u_name,numlines,defaultanswer,options);

        sptpara.f_rate = str2double(user_var{1});
        sptpara.pixl = str2double(user_var{2});
        sptpara.mtl = str2double(user_var{3});
        sptpara.thresh = str2double(user_var{4});
        sptpara.IntTh = str2double(user_var{5});
        sptpara.dia = str2double(user_var{6});
        sptpara.estD = str2double(user_var{7});
        sptpara.fitmethod = str2double(user_var{8});
        sptpara.boxr = str2double(user_var{9});
        sptpara.saveFilterImgMode = str2double(user_var{10});
        sptpara.trackMem = str2double(user_var{11});
        
        % set maxdisp to 2 standard deviations above the mean disp
        avgdisp = sqrt(4*sptpara.estD./sptpara.f_rate);
        sptpara.max_disp = round(3*avgdisp./sptpara.pixl);
        
        % fix bad entries
        while true
            if ~rem(sptpara.dia,2)
                
                % Warning Box
                warnstr = 'Particle Diameter must be odd! Please Reset!';
                dlgname = 'SPT Warning!';
                createmode.Interpreter = 'tex';
                createmode.WindowStyle = 'model';
                h = warndlg(warnstr,dlgname,createmode);
                uiwait(h);
                
                % new input
                prompt = 'Particle Diameter(pixel, must be odd):';
                u_name = 'Input Parameters for Single Particle Tracking';
                numlines = 1;
                defaultanswer = {num2str(sptpara.dia)};
                options.Resize = 'on';
                options.WindowStyle = 'normal';
                options.Interpreter = 'tex';
                user_var = inputdlg(prompt,u_name,numlines,defaultanswer,options);
                
                sptpara.dia = str2double(user_var{1});
            else
                break
            end
        end
        
        while true
            if ~rem(sptpara.boxr,2)
                
                % Warning Box
                warnstr = 'Windows Size must be odd! Please Reset!';
                dlgname = 'SPT Warning!';
                createmode.Interpreter = 'tex';
                createmode.WindowStyle = 'model';
                h = warndlg(warnstr,dlgname,createmode);
                uiwait(h);
                
                % new input
                prompt = 'Windows Size(pixel, must be odd):';
                u_name = 'Input Parameters for Interactive Feature Location';
                numlines = 1;
                defaultanswer = {num2str(sptpara.boxr)};
                options.Resize = 'on';
                options.WindowStyle = 'normal';
                options.Interpreter = 'tex';
                user_var = inputdlg(prompt,u_name,numlines,defaultanswer,options);
                
                sptpara.boxr = str2double(user_var{1});
            else
                break
            end
        end
        
        % Same Parameters for all Movies?
        if i == 1
            if ntrials > 1
                qstring = 'Use Same Parameters for all Movies?';
                atitle = 'Input Parameters for Single Particle Tracking';
                button = questdlg(qstring, atitle, 'Yes', 'No', 'Yes');
                
                switch (button)
                    
                    case 'Yes'
                        ASK = 0;
                        
                    case 'No'
                        ASK = 1;
                        
                    otherwise
                        ASK = 0;
                end
                
            else
                ASK = 0;
            end
        end
        
        if INT == 1
            % Interactive Mode?
            qstring = 'Run in Interactive Mode?';
            atitle = 'Single Particle Tracking';
            button = questdlg(qstring, atitle, 'Yes', 'No', 'No');
            switch (button)
                
                case 'Yes'
                    sptint = 1;
                    
                otherwise
                    sptint = 0;
                    
            end
        end
        INT = 0;
    end
    
%% Calculate

    sptpara.file = stackn{i};
    
    % Run Interactive Mode
    if sptint == 1
        
        % If using different parameters for all movies 
        if ASK == 1
            sptpara.intfile = stackn(i);
        else
            sptpara.intfile = stackn;
            sptint = 0;
        end
        
        sptpara = spt_interactive(sptpara);
        sptpara = rmfield(sptpara,'intfile');
    end

    % Track
    disp(['Now Tracking ',sptpara.file]);

    % Load images
    % Get information about the tiff file
    info = imfinfo(sptpara.file);
    imcnt = numel(info);
    
    % read the first plane to get general tags
    imageTags = tiffread(sptpara.file,1,1);
    
    % preallocate space for image
    im.rawImg = zeros([info(1).Height,info(1).Width,imcnt],...
        class(imageTags.data));
    
    % remove data field
    imageTags = rmfield(imageTags,'data');
    im.imageAttr = imageTags;
    
    % read in each plane
    for jIFD = 1:imcnt
        im.rawImg(:,:,jIFD) = imread(sptpara.file, jIFD,...
            'Info', info);
    end

    % find positions   
    [im, poslist] = spt_fndpos(sptpara,im);
    [trajlist, traj] = spt_track(sptpara,poslist);
    
    % Break out image size (for use with spt_plottraj)
    if not(isfield(im.imageAttr,'width'))
        im.imageAttr.width = info(1).Width;
    end
    if not(isfield(im.imageAttr,'height'));
        im.imageAttr.height = info(1).Height;
    end
    sptpara.size = [im.imageAttr.width, im.imageAttr.height];
    
    % write out filtered image file
    [fpath,fname,~] = fileparts(sptpara.file);
    filterFileName = fullfile(fpath,[fname,'_filter.tif']);
    % remove the image if it was previously written
    recycle on;
    delete(filterFileName); 
    if (sptpara.saveFilterImgMode==2)||(sptpara.saveFilterImgMode==3)        
        tiffwrite(filterFileName,im.filterImg,[],'false');
    end
    
    % remove the image fields?
    if (sptpara.saveFilterImgMode==0)||(sptpara.saveFilterImgMode==2) 
        im = rmfield(im,'filterImg');
    end
    im = rmfield(im,'rawImg');
    
    % Save Variables
    save([sptpara.file(1:end-4),'.mat'],'im','traj','sptpara');
end


