% example script for implementation of the rho_fingerprint technique
% [Please retain the following traceback notice in your final code]
% Version: 20130501.1 (git: http://git.code.sf.net/p/rflex/fingerprint)
%
% Adaptation: -
%
% This code is the implementation of the 'fingerprint' method to derive Rrs from hyperspectral (ir)radiance measurements:
% Simis, S.G.H. and J. Olsson. Unattended processing of shipborne hyperspectral reflectance measurements. Remote Sensing of Environment, in press. DOI: 10.1016/j.rse.2013.04.001
%
% <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB"><img src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png"></a>
% This work is licensed under a <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>
%
% Rflex hardware/software are described at http://sourceforge.net/p/rflex/wiki/Home/
% This code is maintained in a git repository at http://sourceforge.net/p/rflex/fingerprint/

% See also rho_fingerprint_optimize rho_fingerprint_getfingerprint fminbnd


%% load sample data
clear;
load ('sample data\SampleData_matlab.mat')
Ed = EdSample; Ls = LsSample; Lt = LwSample;
ID = [1:1:12]'; clear *Sample;

% Variables Lt, Ls, Ed, ID, and wl are in memory:
% Lt = radiance recorded by water-pointing sensor (format wavelength (row) x sample (column))
% Ls = sky radiance (same format as Lt)
% Ed = downwelling irradiance (same format as Lt)
% wl = wavelength grid (array) for all of the above, as 1 column
% ID = sample identifier


%% Initialization
% 1. Parameters
native_res = mean(diff(wl));               %spectral interval of input data (monotonous)
fingerprint_res = native_res;              %spectral interval used to lookup fingerprint. For narrow band spectrometers set to 2 (nm) or higher. Ignored if equal to native_res.
rho_Lo = 0.024;                            % lower bound for rho_s. Note that higher bound will be determined dynamically from ratio Lt/Ls

% Exclude the following part(s) of the spectrum, recommend excluding the
% oxygen absorption band around 760 nm. Note that spectrum edges (including NaNs in Lt and Ls) will be excluded automatically
rangeexclude = {find(wl>=750,1,'first'):find(wl>=780,1,'first')}; %add cells to the cell array to include additional sections

% spectral domain where negative Rrs must not occur, this will also set the upper limit of rho_s:
nonnegrange = [find(wl>400,1,'first'):1:find(wl<700,1,'last')]; 

% 2 dynamic parameters (for advanced use)
wlFpRes= [wl(1):fingerprint_res:wl(end)]'; %Generates wavelength grid corresponding to fingerprint_res
bandwidth_opt=round(7.5/native_res);       %Interval over which optimization takes place around each fingerprint feature. 
featureseparator = round(10/mean(diff(wlFpRes)));
edge_width = ceil(bandwidth_opt*native_res/fingerprint_res); 
rangeexcludemask = ones(size(wl)); rangeexcludemask(cell2mat(rangeexclude'))=0;

% 3 initialize rho_s and flag arrays (no need to change)
rho_fingerprint = nan(numel(ID),1);             % will hold sky light reflectance (rho_s) results
rho_fit_error = nan(numel(ID),1);               % fit error from optimization routine, can use for debugging / filtering
rho_fingerprint_wlind = cell(numel(ID),1);      % used to store wavebands of fingerprint results from rho_fingerprint_getfingerprint
rho_fingerprint_ind = cell(numel(ID),1);        % spectral band indices of fingerprint, on the same grid as provided to rho_fingerprint_getfingerprint
rho_fingerprint_nind = nan(numel(ID),1);        % stores the number of fingerprint bands identified (can use for debugging/filtering)
rho_fingerprint_SuspectFlag = zeros(numel(ID),1); % Flag raised if any value of rho within upper and lower limit would result in negative Rrs
rho_fingerprint_HiFlag = zeros(numel(ID),1);    % Flag raised if the fingerprint optmization terminates at the upper limit of rho_s
rho_fingerprint_LoFlag = zeros(numel(ID),1);    % Flag raised if the fingerprint optmization terminates at the lower limit of rho_s
rho_fingerprint_EmptyFlag = zeros(numel(ID),1); % Flag raised if the fingerprint returns no band indices. Rho_s will be NaN.

% 4 fminbnd settings. Optimized for (ir)radiance units of [mW/(m2 nm (sr))] (TriOS RAMSES default calibrated output). Changes usually not required.
fitoptions = optimset('fzero'); fitoptions = optimset('Display','off','TolX',1e-12); %TolX must be low to see small differences in Lu/Ed.

%% start processing
t=tic;tic; %time keeping for large data sets
for i = 1:numel(ID);
    tmpLt = Lt(:,i); %subset of spectrum
    tmpLs = Ls(:,i); %subset of spectrum
    tmpEd = Ed(:,i); %subset of spectrum

    % optional interpolation (only for feature identification)
    if native_res ~= fingerprint_res
        tmpLtFpRes = interp1(wl,tmpLt,wlFpRes);
        tmpLsFpRes = interp1(wl,tmpLs,wlFpRes);
        tmpEdFpRes = interp1(wl,tmpEd,wlFpRes);
    else
        tmpLtFpRes = tmpLt;
        tmpLsFpRes = tmpLs;
        tmpEdFpRes = tmpEd;
    end;
    
    rho_Hi = min(tmpLt(nonnegrange)./tmpLs(nonnegrange)); %define upper rho_s limit (no negative Rrs) 
    
    if rho_Hi<=rho_Lo;
        rho_Hi=rho_Lo*1.1;
        rho_fingerprint_SuspectFlag(i) = 1; %this one is bound to result in negative values
    end;
    
    %retrieve 'fingerprint' band indices; for 'featureseparator' see function help:
    rho_fingerprint_ind{i,1} = rho_fingerprint_getfingerprint(tmpLtFpRes,tmpLsFpRes,featureseparator,edge_width);
    
    %Locate the positions of the fingerprint bands on the original wavelength grid
    for r = 1:numel(rho_fingerprint_ind{i,1});
        [d,match] = min((wl-wlFpRes(rho_fingerprint_ind{i,1}(r))).^2);
        rho_fingerprint_ind{i,1}(r) = match;
    end;

    %remove fingerprint bands in unwanted areas (edges, oxygen band, defined in rangeexclude cell array)
    rho_fingerprint_ind{i,1}(rangeexcludemask(rho_fingerprint_ind{i,1})==0)=[];
    rho_fingerprint_wlind{i,1} = wl(rho_fingerprint_ind{i,1}); %store wavebands used for later inspection/debugging
    rho_fingerprint_nind(i,1) = numel(rho_fingerprint_ind{i,1});  %store number of indices found(can use as quality filter)

    if rho_fingerprint_nind(i,1)==0;
        rho_fingerprint_EmptyFlag(i) = 1;
    else    

    %now optimize rho_s with fminbnd
    [rho_fingerprint(i,1),rho_fingerprint_fiterror(i,1)] = fminbnd(@(r) rho_fingerprint_optimize(r,tmpLt,tmpLs,tmpEd,rho_fingerprint_ind{i,1},bandwidth_opt),rho_Lo,rho_Hi,fitoptions); 

    if rho_fingerprint(i) >= rho_Hi-0.0001 % if returned rho is within narrow limit of rho_Hi, raise flag
        rho_fingerprint_HiFlag(i) = 1;
    end;             
    if rho_fingerprint(i) <= rho_Lo+0.0001 % if returned rho is within narrow limit of rho_Lo, raise flag
        rho_fingerprint_LoFlag(i) = 1;
    end;
    end;    
    
    format short g; %use compact display
    disp(['id=',num2str(i),' progress=',num2str(100*i/numel(ID),'%3.1f'),'%, n=',num2str(rho_fingerprint_nind(i)),', lambda=',num2str(rho_fingerprint_wlind{i,1}')]); %visual progress
end;
rho_fingerprint_noFlags = ((1-rho_fingerprint_HiFlag).*(1-rho_fingerprint_LoFlag).*(1-rho_fingerprint_SuspectFlag).*(1-rho_fingerprint_EmptyFlag));

% Some visual output
disp('Flags: Low High Suspect Empty Pass')
disp([rho_fingerprint_LoFlag rho_fingerprint_HiFlag rho_fingerprint_SuspectFlag rho_fingerprint_EmptyFlag rho_fingerprint_noFlags rho_fingerprint rho_fingerprint_nind])

disp(['bandwidth: ',num2str(bandwidth_opt)])
disp(['featureseparator: ',num2str(featureseparator)])
disp('Flags: Low High Suspect Empty Pass')
disp(sum([rho_fingerprint_LoFlag rho_fingerprint_HiFlag rho_fingerprint_SuspectFlag rho_fingerprint_EmptyFlag rho_fingerprint_noFlags]))

% Calculate Rrs using optimized rho_s:
for i=1:numel(rho_fingerprint);
    Rrs(:,i) = (Lt(:,i) - rho_fingerprint(i)* Ls(:,i))./Ed(:,i); 
end

% Plot results
figure(99);clf; %open a figure to summarize the results
if sum(rho_fingerprint_HiFlag>0);p_hi=plot(wl,Rrs(:,rho_fingerprint_HiFlag==1),'r');hold on;end;
if sum(rho_fingerprint_LoFlag>0);p_lo=plot(wl,Rrs(:,rho_fingerprint_LoFlag==1),'b');hold on;end;
if sum(rho_fingerprint_SuspectFlag>0);p_sp=plot(wl,Rrs(:,rho_fingerprint_SuspectFlag==1),':k');hold on;end;
if sum(rho_fingerprint_noFlags>0);p_no=plot(wl,Rrs(:,rho_fingerprint_noFlags==1),'g');hold on;end;
plot([350 800],[0 0],'k');
xlabel('Wavelength (nm)');ylabel('R_{rs} (sr^{-1})');