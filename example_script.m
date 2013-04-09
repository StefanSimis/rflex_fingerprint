% example script for implementation of the rho_fingerprint technique
% [Please retain the following traceback notice in your final code]
% Version: 20130315.1
%
% Adaptation: -
%
% This code is the implementation of the 'fingerprint' method to derive Rrs from hyperspectral (ir)radiance measurements:
% Simis, S.G.H. and Olsson, J. Unattended processing of shipborne hyperspectral reflectance measurements. Submitted (Oct 2012).
%
% <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB"><img src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png"></a>
% This work is licensed under a <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>
%
% Rflex hardware/software are described at http://sourceforge.net/p/rflex/wiki/Home/
% 
% See also rho_fingerprint_optimize rho_fingerprint_getfingerprint fminbnd


%% load sample data

load ('sample data\SampleData_matlab.mat')
Ed = EdSample; Ls = LsSample; Lt = LwSample;
ID = [1:1:12]'; clear *Sample;

% We expect variables Lt, Ls, Ed, ID, and wl to be in memory:
% Lt = radiance recorded by water-pointing sensor (format wavelength (row) x sample (column))
% Ls = sky radiance (same format as Lt)
% Ed = downwelling irradiance (same format as Lt)
% wl = wavelength grid (array) for all of the above (use interp1 to resample if necessary)
% ID = sample identifier
%% initialize
bandwidth=2;  %proper for 3.3 nm spaced bands (Ramses, Hypersas). Higher resolution -> increase value to approximately bandwidth = round(7.5/resolution)

%initialize variables
rho_fingerprint = nan(numel(ID),1); % will hold sky light reflectance (rho_s) results
rho_fit_error = nan(numel(ID),1);   % fit error from optimization routine, can use for debugging / filtering
rho_fingerprint_wlind = cell(numel(ID),1); % used to store wavebands of fingerprint results from rho_fingerprint_getfingerprint
rho_fingerprint_ind = cell(numel(ID),1);   % spectral band indices of fingerprint, on the same grid as provided to rho_fingerprint_getfingerprint
rho_fingerprint_nind = nan(numel(ID),1);   % used to store the number of fingerprint bands identified (can use for debugging/filtering)
rho_fingerprint_SuspectFlag = zeros(numel(ID),1); % raised if any value of rho within upper and lower limit would result in negative Rrs
rho_fingerprint_HiFlag = zeros(numel(ID),1); % raised if the fingerprint optmization terminates at the upper limit of rho_s
rho_fingerprint_LoFlag = zeros(numel(ID),1); % raised if the fingerprint optmization terminates at the lower limit of rho_s

% fminbnd settings. Optimized for (ir)radiance units of [mW/(m2 nm (sr))] (TriOS RAMSES default calibrated output)
fitoptions = optimset('fzero'); fitoptions = optimset('Display','off','TolX',1e-12); %TolX must be low to see small differences in Lu/Ed.

% algorithm definitions
rho_Lo = 0.024;            % lower bound for r, higher bound will be determined from ratio Lt/Ls

% possibility to exclude part(s) of the spectrum, recommend excluding the oxygen absorption band around 760 nm:
rangeinclude = [bandwidth:find(wl<=740,1,'last'),...
                find(wl>=780,1,'first'):numel(wl)-bandwidth]; 
tmpwl = wl(rangeinclude); % wavelength grid corresponding to selection

% spectral domain where negative Rrs must not occur, this will also set the upper limit of rho_s:
nonnegrange = [find(tmpwl>=400,1,'first'):1:find(tmpwl>=700,1,'first')]; % this uses the subset of wl

%% start processing
t=tic;tic; %time keeping for large data sets
for i = 1:numel(ID);
    if toc>2;
        disp(['Processing, completion: ',num2str(i/numel(ID))]);tic; %update on progress
    end;
    tmpLt = Lt(rangeinclude,i); %subset of spectrum
    tmpLs = Ls(rangeinclude,i); %subset of spectrum
    tmpEd = Ed(rangeinclude,i); %subset of spectrum
    
    rho_Hi = min(tmpLt(nonnegrange)./tmpLs(nonnegrange)); %define upper rho_s limit (no negative Rrs) 
    
    if rho_Hi<=rho_Lo;
        rho_Hi=rho_Lo*1.1;
        rho_fingerprint_SuspectFlag(i) = 1; %this one is bound to result in negative values
    end;
    
    rho_fingerprint_ind{i,1} = rho_fingerprint_getfingerprint(tmpLt,tmpLs,3); %retrieve 'fingerprint' band indices; the value '3' is the 'featureseparator', see function help
    rho_fingerprint_wlind{i,1} = tmpwl(rho_fingerprint_ind{i,1}); %store wavebands used for later inspection/debugging
    rho_fingerprint_nind(i,1) = numel(rho_fingerprint_ind{i,1});  %store number of indices found(can use as quality filter)

    %now optimize rho_s with fminbnd
    [rho_fingerprint(i,1),rho_fingerprint_fiterror(i,1)] = fminbnd(@(r) rho_fingerprint_optimize(r,tmpLt,tmpLs,tmpEd,rho_fingerprint_ind{i,1},bandwidth),rho_Lo,rho_Hi,fitoptions); 

    if rho_fingerprint(i) >= rho_Hi-0.000001; % if returned rho is within a narrow limit of rho_Hi, raise flag
        rho_fingerprint_HiFlag(i) = 1;
    end;             
    if rho_fingerprint(i) <= rho_Lo+0.001 % if returned rho is within a narrow limit of rho_Lo, raise flag
        rho_fingerprint_LoFlag(i) = 1;
    end;

    disp([i rho_fingerprint_nind(i) rho_fingerprint_wlind{i,1}']); %visual progress

end;
toc;

% some more visual output:
disp('Flags: Low High Suspect')
disp([rho_fingerprint_LoFlag rho_fingerprint_HiFlag rho_fingerprint_SuspectFlag rho_fingerprint])


%% Calculate Rrs using optimized rho, plot results
rho_fingerprint_noFlags = ((1-rho_fingerprint_HiFlag).*(1-rho_fingerprint_LoFlag).*(1-rho_fingerprint_SuspectFlag));

% Calculate Rrs using optimized rho_s:
for i=1:numel(rho_fingerprint);
    Rrs(:,i) = (Lt(:,i) - rho_fingerprint(i)* Ls(:,i))./Ed(:,i); 
end

figure(99);clf; %open a figure to summarize the results
if sum(rho_fingerprint_HiFlag>0);p_hi=plot(wl,Rrs(:,rho_fingerprint_HiFlag==1),'r');hold on;end;
if sum(rho_fingerprint_LoFlag>0);p_lo=plot(wl,Rrs(:,rho_fingerprint_LoFlag==1),'b');hold on;end;
if sum(rho_fingerprint_SuspectFlag>0);p_sp=plot(wl,Rrs(:,rho_fingerprint_SuspectFlag==1),':k');hold on;end;
if sum(rho_fingerprint_noFlags>0);p_no=plot(wl,Rrs(:,rho_fingerprint_noFlags==1),'g');hold on;end;
plot([350 800],[0 0],'k');
xlabel('Wavelength (nm)');ylabel('R_{rs} (sr^{-1})');