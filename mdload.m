function DATAS = mdload(OPTS)
%MDLOAD loads files for multirho fits
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   DATAS contains that dataset.

pirat = pi/180;
chopidxs = 424:1605;

% Fix basedir if necessary
if OPTS.basedir(end) ~= '\'
    OPTS.basedir = [OPTS.basedir '/'];
end

% Generate filenames
filenames = cell(length(OPTS.rhorange),1);
for fnidx = 1:length(OPTS.rhorange)
    string_proto = strrep(OPTS.filenameprototype,...
        num2str(OPTS.rhorange(1)),num2str(OPTS.rhorange(fnidx)));    
    filenames{fnidx} = [OPTS.basedir string_proto '-dcswitch.asc'];

end

%% Error check
for fnidx = 1:size(filenames,1)
    findvar = length(strfind(filenames{fnidx,1},[num2str(OPTS.rhorange(fnidx)) '-']));
    if findvar == 0
        findvar = length(strfind(filenames{fnidx,1},num2str(OPTS.rhorange(fnidx))));
        if findvar == 0
            error('Could not replace rho in file prototype');
        elseif findvar ~=1
            error('Rho can only appear once in filenames');
        end
    elseif findvar ~=1
        error('Rho can only appear once in filenames');
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
        temp = importdata([OPTS.basedir OPTS.darkname '-dcswitch.asc'],...
            '\t',16);
        darkamps = temp.data(:,3:2:end);
    else
        for i = 1:OPTS.darkreps
            temp = importdata([OPTS.basedir OPTS.darkname '-' sprintf('%04d',i)...
                '-dcswitch.asc'],'\t',16);
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
    for fnidx = 1:length(OPTS.bbrhorange)
        string_proto = strrep(OPTS.filenameprototype,...
            num2str(OPTS.rhorange(1)),num2str(OPTS.bbrhorange(fnidx)));        
        bbfilenames{fnidx,1} = [OPTS.basedir string_proto '-tis.asc'];
        if OPTS.subtractdark == 1
            bbfilenames{fnidx,2} = [OPTS.basedir string_proto '-tis-dark.asc'];
        end        
        if OPTS.overx == 1
            bbfilenames{fnidx,3} = [OPTS.basedir string_proto '-x-tis.asc'];
            if OPTS.subtractdark == 1
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
    elseif OPTS.subtractdark == 0
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

        if OPTS.overx == 1
            % Use double-exposure data to find filter indexes
            temp = importdata(bbfilenames{r_idx,3},'\t',13);
            for g = 1:size(temp.data,2)-1
                tidx(g,:) = [find(temp.data(:,g+1)==65535,1,'first'),find(temp.data(:,g+1)==65535,1,'last')];
            end
            
        end               
        filter_refl(:,r_idx) = mean(temp.data(chopidxs,2:end),2); % for filtering   
        
        if OPTS.overx == 1
            temp = importdata(bbfilenames{r_idx,3},'\t',13);
            inttimex = str2num(temp.textdata{6}(24:end));

            % move ~15 nm (35 pixels) from region of CCD saturation
            tempcts = temp.data(chopidxs,2:end);          
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
            if OPTS.subtractdark == 1
                temp=importdata(bbfilenames{r_idx,4},'\t',13);
                % FIX INTTIME RECORDING IN LBS SOFTWARE            
                dreflx(:,r_idx) = mean(temp.data(chopidxs,2:end),2).*(1000/(2*inttimex)); % counts/s  
            end
        end
        
        if OPTS.subtractdark == 1
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
            for i = 1:nbins
                t(i,:) = mean(unc_reflx((b*i-b+1):(b*i),:));
            end 
            unc_reflx = t;
        end
        if exist('drefl','var')
            for i = 1:nbins
                t(i,:) = mean(drefl((b*i-b+1):(b*i),:));
            end 
            drefl = t;
        end
        if exist('dreflx','var')
            for i = 1:nbins
                t(i,:) = mean(dreflx((b*i-b+1):(b*i),:));
            end 
            dreflx = t;
        end 
        if exist('filter_refl','var')
            for i = 1:nbins
                t(i,:) = mean(filter_refl((b*i-b+1):(b*i),:));
            end
            filter_refl = t;
        end                    
    end
    if OPTS.overx == 1               
        for zzyzx = 1:size(unc_refl,2)
            unc_refl(1:DATAS.oxidxs(1,zzyzx),zzyzx) = unc_reflx(1:DATAS.oxidxs(1,zzyzx),zzyzx);
            unc_refl(DATAS.oxidxs(2,zzyzx):end,zzyzx) = unc_reflx(DATAS.oxidxs(2,zzyzx):end,zzyzx);
            if OPTS.subtractdark == 1
                drefl(1:DATAS.oxidxs(1,zzyzx),zzyzx) = dreflx(1:DATAS.oxidxs(1,zzyzx),zzyzx);
                drefl(DATAS.oxidxs(2,zzyzx):end,zzyzx) = dreflx(DATAS.oxidxs(2,zzyzx):end,zzyzx);
            end
        end
        if OPTS.subtractdark == 1
            refl = unc_refl-drefl;
        else
            refl = unc_refl;
        end
    else
        if OPTS.subtractdark == 1
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