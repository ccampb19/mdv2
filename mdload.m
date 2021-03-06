function DATAS = mdload(OPTS)
%MDLOAD loads files for multirho fits
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   DATAS contains that dataset.

pirat = pi/180;
chopidxs = 424:1605;

% Handle some overX edge cases
stitchbit = 0;

% Fix basedir if necessary
if ~any(strcmp(OPTS.basedir(end),{'\','/'}))
    OPTS.basedir(end+1) = filesep;
end

% Generate filenames
filenames = cell(length(OPTS.rhorange),1);
for fnidx = 1:length(OPTS.rhorange)
    % Pad prototype to address case of rho at beginning/end of prototype
    string_proto = strrep(['-' OPTS.filenameprototype '-'],...
        ['-',num2str(OPTS.rhorange(1)),'-'],...
        ['-',num2str(OPTS.rhorange(fnidx)),'-']);
    % Then remove the padding
    filenames{fnidx} = [OPTS.basedir string_proto(2:end-1) '-' OPTS.systemname '.asc'];
end

%% Error check
% There are many potential errors from this file naming scheme. Address the
% most common.
for fnidx = 1:size(filenames,1)
    % Get file name. (not sure why I did this--fileparts() would work just as well)
    testidx = regexp(filenames{fnidx,1},'[^/\\]+$');
    % Compare rhos from filenames with rhorange
    findvar = numel(strfind(filenames{fnidx,1}(testidx:end),[num2str(OPTS.rhorange(fnidx)) '-']));
    if findvar == 0
        % Try without the hyphen, though there should be a delimiter present
        findvar = numel(strfind(filenames{fnidx,1}(testidx:end),num2str(OPTS.rhorange(fnidx))));
        if findvar == 0
            error('Could not replace rho in file prototype');
        elseif findvar ~=1
            % Multiple solutions ('a-14-blah','blah-14','14-blah'). Suggest one.
            error('File names should begin with rho (mm), e.g: "14-ndabsc-..."');
        end
        % Not sure what the error was that led me to restrict rhos to
        % 1 occurrence in the filename. I'm sure it will throw errors again 
        % at some point, unless the root cause has been fixed.
    end
end

%% Load FD Data
for fnidx = 1:length(filenames)
    temp=importdata(filenames{fnidx,1},'\t',16);
    DATAS.phases(fnidx,:,:) = temp.data(:,OPTS.usediodes.*2).*pirat;
    DATAS.amps(fnidx,:,:) = temp.data(:,OPTS.usediodes.*2+1);
end
DATAS.complex = DATAS.amps.*cos(DATAS.phases)...
        +1i.*DATAS.amps.*sin(DATAS.phases);
DATAS.freqs = temp.data(:,1);

%% Filter dark FD data (different from broadband dark data filenames above)?
if isfield(OPTS,'darkname')
    if OPTS.darkreps == 0
        temp = importdata([OPTS.basedir OPTS.darkname '-' OPTS.systemname '.asc'],...
            '\t',16);
        darkamps = temp.data(:,3:2:end);
    else
        for i = 1:OPTS.darkreps
            temp = importdata([OPTS.basedir OPTS.darkname '-' sprintf('%04d',i)...
                '-' OPTS.systemname '.asc'],'\t',16);
            darkamps(:,:,i) = temp.data(:,3:2:end);
        end
    end
    
    % Arbitrary point currently now 10 dB above measured noise floor
    dfloor = 10.^((10+20.*log10(max(darkamps,[],3)))/20)';
    % This is getting ugly. Has to be a better way.
    cutoffs = size(DATAS.amps,2).*ones(size(DATAS.amps,1),size(DATAS.amps,3));
    for f_idx = 1:size(DATAS.amps,1)
        for d_idx = 1:size(DATAS.amps,3)
            testdata = squeeze(DATAS.amps(f_idx,:,d_idx))-dfloor(d_idx,:);
            tries = find(testdata<0);
            triess = strfind(diff(tries),[1 1]);
            if ~isempty(triess)
                cutoffs(f_idx,d_idx) = tries(triess(1));
            end
        end
    end
    % To simplify things, use the cutoff frequency where a normalized
    % line first passes through the cutoff function for a given diode
    % SDS-specific cutoffs could be developed with more complicated code
    normline = (1:f_idx)'.*(size(DATAS.amps,2)./f_idx);
    for d_idx = 1:size(DATAS.amps,3)
        roidx = find((cutoffs(:,d_idx)-normline)<0,1,'first')-1;
        if isempty(roidx)
            roidx = size(cutoffs,1);
        end
        foidx = cutoffs(roidx,d_idx);
        DATAS.cutoff(d_idx,1) = OPTS.rhorange(roidx);
        DATAS.cutoff(d_idx,2) = foidx;
    end
    DATAS.cutoffidxs = cutoffs;
end

%% Load Broadband data?
if OPTS.bb == 1
    disp('Loading broadband data...')
    for fnidx = 1:length(OPTS.bbrhorange)
        % See FD filename replacements for explanatory comments
        string_proto = strrep(['-',OPTS.filenameprototype,'-'],...
            ['-',num2str(OPTS.rhorange(1)),'-'],...
            ['-',num2str(OPTS.bbrhorange(fnidx)),'-']);
        string_proto = string_proto(2:end-1);
        bbfilenames{fnidx,1} = [OPTS.basedir string_proto '-tis.asc'];
        if OPTS.subtractbbdark == 1
            bbfilenames{fnidx,2} = [OPTS.basedir string_proto '-tis-dark.asc'];
        end        
        if OPTS.overx == 1
            bbfilenames{fnidx,3} = [OPTS.basedir string_proto '-x-tis.asc'];
            if OPTS.subtractbbdark == 1
                bbfilenames{fnidx,4} = [OPTS.basedir string_proto '-x-tis-dark.asc'];
            end             
        end
    end
    % Error check
    assert(~strcmp(bbfilenames{1,1},bbfilenames{2,1}),'Error: Filenames not populated correctly')
    
    switch OPTS.sphreps
        case -1
        case 0
            temp = importdata([OPTS.basedir OPTS.sphname '-tis.asc'],...
                '\t',13);
            sphamps = temp.data(chopidxs,2:end);
            sphint = str2num(temp.textdata{6}(24:end));
            srefl = mean(sphamps,2).*(1000/sphint); % counts/s
        otherwise
            for i = 1:OPTS.sphreps
                temp = importdata([OPTS.basedir OPTS.sphname '-' sprintf('%04d',i)...
                    '-tis.asc'],'\t',13);
                sphint(i) = str2num(temp.textdata{6}(24:end));            
                sphamps(:,:,i) = temp.data(chopidxs,2:end).*(1000/sphint(i));
            end
            srefl = mean(mean(sphamps,2),3);
    end
    
    % About time to make a function out of this sequence
    % Leaving sphere cal as an option to study water peak, for now
    if OPTS.sphreps == -1       
    elseif OPTS.subtractbbdark == 0
        temp = importdata([OPTS.basedir OPTS.sphname '-tis-dark.asc'],...
            '\t',13);
        sphdint = str2num(temp.textdata{6}(24:end));      
        sphdrefl = temp.data(chopidxs,2:end).*(1000/sphdint);
        srefl = srefl-mean(sphdrefl,2);
    else
        for i = 1:OPTS.sphreps
            temp = importdata([OPTS.basedir OPTS.sphname '-' sprintf('%04d',i)...
                '-tis-dark.asc'],'\t',13);
            sphdint(i) = str2num(temp.textdata{6}(24:end));            
            sphdamps(:,:,i) = temp.data(chopidxs,2:end).*(1000/sphint(i));
        end
        srefl = srefl-mean(mean(sphdamps,2),3);
    end
    
    % Smooth sphere data?
    if (OPTS.smooth == 1 && OPTS.sphreps ~= -1)
        srefl = spectrumSmoother(srefl,1,3);
        srefl = spectrumSmoother(srefl,3);
    end
    
    % Load measurement reflectances
    for r_idx = 1:size(bbfilenames,1)
        temp=importdata(bbfilenames{r_idx,1},'\t',13);
        inttime = str2num(temp.textdata{6}(24:end));
        unc_refl(:,r_idx) = mean(temp.data(chopidxs,2:end),2).*(1000/inttime); % counts/s
        DATAS.wv = temp.data(chopidxs,1);    

        if OPTS.overx == 1 && ~stitchbit
            % Use double-exposure data to find filter indexes
            temp = importdata(bbfilenames{r_idx,3},'\t',13);
            for g = 1:size(temp.data,2)-1

                try
                    tidx(g,:) = [find(temp.data(:,g+1)==65535,1,'first'),find(temp.data(:,g+1)==65535,1,'last')];
                catch                   
                    % Either we hit the integration time limit, so use the
                    % single exposure...
                    if r_idx > 1
                        warning(['Overexposures ignored for SDS >', num2str(OPTS.bbrhorange(r_idx)-1) ' mm'])
                        stitchbit = r_idx;
                        break
                    else
                        % or spec is configured for online dark correction,
                        % which could interfere with the overX stitching.
                        warning('%s\n%s','Dynamic dark correction used. Skipping overX correction.',...
                            '(Disable option in spectrometer EEPROM to fix.)')
                        stitchbit = -1;
                        break
                    end
                    
                end
            end
            
        end   

        % Here, temp data will be overX data if applicable. Otherwise, single
        % exposure. Find counts below threshold.
        filter_refl(:,r_idx) = mean(temp.data(chopidxs,2:end),2); % for filtering   
        
        if OPTS.overx == 1 && ~stitchbit
            
            temp = importdata(bbfilenames{r_idx,3},'\t',13);
            inttimex = str2num(temp.textdata{6}(24:end));

            % move ~15 nm (35 pixels) from region of CCD saturation
            tempcts = temp.data(chopidxs,2:end);     
            % 65535
            [~,a] = find(tempcts'==65535,1,'first');
            [~,b] = find(tempcts'==65535,1,'last');
            tidxs = [a,b];
             
            if OPTS.binning > 0
                if ~isempty(a)
                    DATAS.oxidxs(:,r_idx) = ...
                        [floor((a-35)/OPTS.binning);ceil((b+35)/OPTS.binning)];
                else                
                    DATAS.oxidxs(:,r_idx) = ...
                        [floor(36/OPTS.binning);ceil((size(tempcts,2)-35)/OPTS.binning)];
                end                
            else
                % No binning
                if ~isempty(a)
                    DATAS.oxidxs(:,r_idx) = [a;b];
                else
                    DATAS.oxidxs(:,r_idx) = [a-35; b+35];
                end                   
            end

            % FIX INTTIME RECORDING IN LBS SOFTWARE
            unc_reflx(:,r_idx) = mean(temp.data(chopidxs,2:end),2).*(1000/(2*inttimex));
            if OPTS.subtractbbdark == 1
                temp=importdata(bbfilenames{r_idx,4},'\t',13);
                % FIX INTTIME RECORDING IN LBS SOFTWARE            
                dreflx(:,r_idx) = mean(temp.data(chopidxs,2:end),2).*(1000/(2*inttimex)); % counts/s  
            end
        end
        
        if OPTS.subtractbbdark == 1
            temp=importdata(bbfilenames{r_idx,2},'\t',13);
            drefl(:,r_idx) = mean(temp.data(chopidxs,2:end),2).*(1000/inttime); % counts/s 
        end
    end
    if OPTS.binning > 0
        b = OPTS.binning;
        nbins = floor(size(unc_refl,1)/b);
        t = zeros(nbins,size(unc_refl,2));            

        for i = 1:nbins
            t(i,:) = mean(unc_refl((b*i-b+1):(b*i),:));
            w(i) = mean(DATAS.wv((b*i-b+1):(b*i)));
        end
        unc_refl = t;
        DATAS.wv = w;

        if exist('unc_reflx','var')
            tx = zeros(nbins,size(unc_reflx,2));
            for i = 1:nbins
                tx(i,:) = mean(unc_reflx((b*i-b+1):(b*i),:));
            end 
            unc_reflx = tx;
        end
        if exist('drefl','var')
            for i = 1:nbins
                t(i,:) = mean(drefl((b*i-b+1):(b*i),:));
            end 
            drefl = t;
        end
        if exist('dreflx','var')
            for i = 1:nbins
                tx(i,:) = mean(dreflx((b*i-b+1):(b*i),:));
            end 
            dreflx = tx;
        end 
        if exist('filter_refl','var')
            for i = 1:nbins
                t(i,:) = mean(filter_refl((b*i-b+1):(b*i),:));
            end
            filter_refl = t;
        end                    
    end
    if OPTS.overx == 1 && stitchbit ~=-1       
        for zzyzx = 1:size(unc_reflx,2)
            unc_refl(1:DATAS.oxidxs(1,zzyzx),zzyzx) = unc_reflx(1:DATAS.oxidxs(1,zzyzx),zzyzx);
            unc_refl(DATAS.oxidxs(2,zzyzx):end,zzyzx) = unc_reflx(DATAS.oxidxs(2,zzyzx):end,zzyzx);
            if OPTS.subtractbbdark == 1
                drefl(1:DATAS.oxidxs(1,zzyzx),zzyzx) = dreflx(1:DATAS.oxidxs(1,zzyzx),zzyzx);
                drefl(DATAS.oxidxs(2,zzyzx):end,zzyzx) = dreflx(DATAS.oxidxs(2,zzyzx):end,zzyzx);
            end
        end
        if OPTS.subtractbbdark == 1
            refl = unc_refl-drefl;
        else
            refl = unc_refl;
        end
    else
        if OPTS.subtractbbdark == 1
            refl = unc_refl-drefl;%.*(1000/inttime); % counts/s
        else
            refl = unc_refl;
        end
    end
    
    % Filter broadband data
    % If spec data (dark-corrected or otherwise) below threshold, don't trust it
    % FIX INTTIME RECORDING IN LBS SOFTWARE            
    [DATAS.bbdarkcols,DATAS.bbdarkrows] = find(filter_refl < OPTS.threshold);
    DATAS.bbdarkidxs = find(filter_refl < OPTS.threshold);
    
    % Also, cut off when reflectance data gets too noisy.
    testdata = (refl-smoothdata(refl,'gaussian',21))./smoothdata(refl,'gaussian',21);
    noisefig = abs(mean(testdata));
    if mod(length(noisefig),2)
        DATAS.ncutoff = find(noisefig == median(noisefig),1,'first');
    else
        DATAS.ncutoff = find(noisefig(2:end) == median(noisefig(2:end)),1,'first');
    end
    if DATAS.ncutoff == 0 || length(OPTS.bbrhorange) < 10
        DATAS.ncutoff = length(noisefig);
    end
    
    if OPTS.smooth == 1
        for r_idx = 1:size(refl,2)
            refl(:,r_idx) = spectrumSmoother(refl(:,r_idx),1,3);
            refl(:,r_idx) = spectrumSmoother(refl(:,r_idx),3);
        end
    end
    DATAS.R = refl;
    if OPTS.sphreps ~= -1
        DATAS.R_sph = refl./srefl;
    end
    
end