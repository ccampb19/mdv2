function DATAS = mdload(OPTS)
%MDLOAD loads files for multirho fits
%   OPTS contains all necessary data to load a mutlidistance dataset.
%   DATAS contains that dataset.

% Fix basedir if necessary
% Currently only implemented for Windows
if OPTS.basedir(end) ~= '\'
    OPTS.basedir = [OPTS.basedir '\'];
end

if exist('mdopt.mat','file') && exist('mddata.mat','file')
    test = load('mdopt.mat');
    if isequal(OPTS,test.OPTS)
        load('mddata.mat','DATAS');
        return;
    end
end

% If files not already loaded, do the rest:
% Generate filenames
filenames = cell(length(OPTS.rhorange),1);
for fnidx = 1:length(OPTS.rhorange)
    filenames{fnidx} = strrep(OPTS.filenameprototype,...
        num2str(OPTS.rhorange(1)),num2str(OPTS.rhorange(fnidx)));
end

%% Error check
for fi = 1:length(filenames)
    findvar = length(strfind(filenames{fi},[num2str(OPTS.rhorange(fi)) '-']));
    if findvar == 0
        findvar = length(strfind(filenames{fi},num2str(OPTS.rhorange(fi))));
        if findvar == 0
            error('Could not replace rho in file prototype');
        elseif findvar ~=1
            error('Rho can only appear once in filenames');
        end
    elseif findvar ~=1
        error('Rho can only appear once in filenames');
    end
end
    
%% Load data
disp('Loading Data...');
pirat = pi/180;
for r_idx = 1:length(filenames)
    temp=importdata([OPTS.basedir filenames{r_idx} '-dcswitch.asc'],'\t',16);
    DATAS.phases(r_idx,:,:) = temp.data(:,OPTS.usediodes.*2).*pirat;
    DATAS.amps(r_idx,:,:) = temp.data(:,OPTS.usediodes.*2+1);
end
DATAS.complex = DATAS.amps.*cos(DATAS.phases)...
        +1i.*DATAS.amps.*sin(DATAS.phases);
DATAS.freqs = temp.data(:,1);

%% Load Broadband data?
if OPTS.bb == 1
    bbfilenames = cell(length(OPTS.bbrhorange),1);
    for fnidx = 1:length(OPTS.bbrhorange)
        bbfilenames{fnidx} = strrep(OPTS.filenameprototype,...
            num2str(OPTS.rhorange(1)),num2str(OPTS.bbrhorange(fnidx)));
    end

    chopidxs = 425:1605;
    temp = importdata([OPTS.basedir OPTS.sphname],'\t',13);
    sphint = str2num(temp.textdata{6}(24:end));
    srefl = temp.data(chopidxs,2).*(1000/sphint); % counts/s
    if OPTS.smooth == 1
        srefl = spectrumSmoother(srefl,1,3);
        srefl = spectrumSmoother(srefl,3);
    end
    for r_idx = 1:length(bbfilenames)
        temp=importdata([OPTS.basedir bbfilenames{r_idx} '-tis.asc'],'\t',13);
        inttime = str2num(temp.textdata{6}(24:end));

%         refl(:,r_idx) = temp.data(:,2);
        refl(:,r_idx) = temp.data(chopidxs,2).*(1000/inttime); % counts/s
        
    end
    if OPTS.smooth == 1
        for r_idx = 1:size(refl,2)
            refl(:,r_idx) = spectrumSmoother(refl(:,r_idx),1,3);
            refl(:,r_idx) = spectrumSmoother(refl(:,r_idx),3);
        end
    end
    DATAS.R = refl;
%     DATAS.R = refl./srefl;
    DATAS.wv = temp.data(chopidxs,1);
    
    % Filter broadband data
    % If spec data (dark-corrected or otherwise) below 8000 cts, don't trust it
    [DATAS.bbdarkcols,DATAS.bbdarkrows] = find(refl < 8000);
    DATAS.bbdarkidxs = find(refl < 8000);
end

%% Filter dark FD data?
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
    
    dfloor = 10.^((3+20.*log10(max(darkamps,[],3)))/20)';
    % This is getting ugly. Has to be a better way.
    cutoffs = size(DATAS.amps,2).*ones(size(DATAS.amps,1),size(DATAS.amps,3));
    for r_idx = 1:size(DATAS.amps,1)
        for d_idx = 1:size(DATAS.amps,3)
            testdata = squeeze(DATAS.amps(r_idx,:,d_idx))-dfloor(d_idx,:);
            tries = find(testdata<0);
            triess = strfind(diff(tries),[1 1]);
            if ~isempty(triess)
                cutoffs(r_idx,d_idx) = tries(triess(1));
            end
        end
    end
    % To simplify things, use the cutoff frequency where a normalized
    % line first passes through the cutoff function for a given diode
    normline = (1:r_idx)'.*(size(DATAS.amps,2)./r_idx);
    for d_idx = 1:size(DATAS.amps,3)
        roidx = find((cutoffs(:,d_idx)-normline)<0,1,'first')-1;
        if isempty(roidx)
            roidx = size(cutoffs,1);
        end
        foidx = cutoffs(roidx,d_idx);
        DATAS.cutoff(d_idx,1) = OPTS.rhorange(roidx);
        DATAS.cutoff(d_idx,2) = foidx;
    end
end
DATAS.cutoffidxs = cutoffs;

disp('Done!')
end

