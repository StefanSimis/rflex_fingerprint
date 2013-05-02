function indices = rho_fingerprint_getfingerprint(Lt,Ls,featureseparator,edge_width)
% Locate peaks and valleys occuring in the same bands of the second derivatives of Ls and Lt
% This function is called once for every sample of Lt and Ls. The retrieved
% indices are then used as input in rho_fingerprint_optimize
% 
% note:
% Input needs to be appropriately spaced, i.e. course bandwidth will not reveal atmospheric gas absorption features 
% whereas too fine resolution will not show co-occurrence in Lt and Ls. For high resolution sensors, interpolate to 2-3 nm intervals
% 
% You can optimize this code: 
% * By setting the hard-coded nfeatures ('search depth') variable so that sufficient indices will be generated. 
% Rule of thumb: spectra spanning 320-950 nm should result in 10-15 indices. 
% Setting nfeatures too high will generate meaningless indices. 
% Try approximately 1/5 of the number of spectral channels (not tested)
%
% SYNTAX: 
% [indices] = rho_fingerprint_getfingerprint(Lt,Ls,featureseparator,edge_width)
% 
% * indices = row indices of spectral features co-occurring in Lt and Ls.
% * Lt = Total radiance measured by water-viewing sensor
% * Ls = Sky radiance
% * featureseparator = interval desired between consecutive indices results: 
% Only the dominant feature in band of width featureseparator will be selected. 
% Example: featureseparator = 3 works for RAMSES 3.3nm input, 10 for WISP-3 measurements, 5 for 2-nm data. 
% * edge_width = area to ignore on either end of the spectrum. 
% This width should be set equal to the bandwidth used in
% rho_fingerprint_optimize. Supports NaN results in spectrum calibration: if NaNs are present at the edges of the spectrum, the ignored area will be extended. 
% 
% [Please retain the following traceback notice in your code]
% Version: 20130501.1 (git: http://git.code.sf.net/p/rflex/fingerprint)
% Adaptation: -
%
% This code is the implementation of the 'fingerprint' method to derive Rrs from hyperspectral (ir)radiance measurements:
% Simis, S.G.H. and J. Olsson. Unattended processing of shipborne hyperspectral reflectance measurements. Remote Sensing of Environment, in press. DOI: 10.1016/j.rse.2013.04.001
%
% <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB"><img src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png"></a>
% This code is licensed under a <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>
%
% Rflex hardware/software are described at http://sourceforge.net/p/rflex/wiki/Home/
% This code is maintained in a git repository at http://sourceforge.net/p/rflex/fingerprint/
% 
% See also rho_fingerprint_optimize

nfeatures = 40; % doubled from previous version due to different sorting (drops and rises merged at start)
ew = edge_width;

%This mask covers edges of the spectrum, including NaNs
%NaNs in Ed are not considered.
edgemask = ones(length(Lt),1);
edgemask(1:ew+find(~isnan(Ls),1,'first')-1)=NaN; %1a edge at start of Ls
edgemask(1:ew+find(~isnan(Lt),1,'first')-1)=NaN; %2a edge at start of Lt
edgemask(2+end-find(~isnan(Ls(end:-1:1)),1,'first')-ew:end)=NaN;%edge at end of Ls
edgemask(2+end-find(~isnan(Lt(end:-1:1)),1,'first')-ew:end)=NaN;%edge at end of Lt

Ltdiff = abs(diff(Lt));
[~,Ltind] = sort(Ltdiff,1,'descend');       %indices sorted by value of Ltdiff
Ltind(isnan(edgemask(Ltind)))=[];           %remove indices where Ltdiff is not a number
Lsdiff = abs(diff(Ls));
[~,Lsind] = sort(Lsdiff,1,'descend');       %indices sorted by value of Lsdiff
Lsind(isnan(edgemask(Lsind)))=[];           %remove indices where Ltdiff is not a number

%identify matching Ls&Lt feature locations from sorted first derivative list
indices = intersect(Ltind(1:nfeatures),Lsind(1:nfeatures));   %Ls and Lt signal changes are strong at these indices
indicesLtd = Ltdiff(indices);               % Ltdiff values associated with drops&rises
indicesLtd = indicesLtd/max(indicesLtd);	% normalized
indicesLsd = Lsdiff(indices);               % Lsdiff values associated with drops&rises
indicesLsd = indicesLsd/max(indicesLsd);    % normalized
scores = indicesLtd+indicesLsd;             % scores for combined Ls and Lt feature strength
nout = [];

while ~isempty(indices)
    [~,hiscore] = max(scores);
    nout = [nout;indices(hiscore)];
    r = [indices(hiscore)-featureseparator:1:indices(hiscore)+featureseparator];
    for i=1:numel(r)
        t = find(indices==r(i));
        scores(t) = [];
        indices(t) = [];
    end;
end;
indices = nout;                               % strongest indices, separated