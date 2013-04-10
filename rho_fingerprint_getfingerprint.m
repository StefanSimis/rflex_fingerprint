function indices = rho_fingerprint_getfingerprint(Lt,Ls,featureseparator)
% Locate extremes occuring in the same bands of the second derivatives of Ls and Lt
% This function is called once for every sample of Lt and Ls. The retrieved
% indices are then used as input in rho_fingerprint_optimize
% 
% note:
% Input needs to be appropriately spaced, i.e. course bandwidth will not reveal atmospheric gas absorption features 
% whereas too fine resolution will not show co-occurrence in Lt and Ls. For high resolution sensors, interpolate to 2-3 nm intervals
% 
% You can optimize this code: 
% * Set number of features to include on first selection, set the nfeatures
% variable (15-20 recommended, if 15-20 indices result chances are Rrs can be fitted well. 
% Note that sensors lacking UV and NIR information will generate fewer bands and Rrs will be less reliable.
%
% SYNTAX: 
% [indices] = rho_fingerprint_getfingerprint(Lt,Ls,featureseparator)
% 
% indices = Resulting row indices of spectral features co-occurring (same spectral band) in Lt and Ls.
% Lt = Total radiance measured by water-viewing sensor
% Ls = Sky radiance
% featureseparator = the waveband interval desired between consecutive results. 
% Results that are grouped together should not all be represented in the output, 
% as this would cause broader features to dominate the optimization in rho_fingerprint_optimize. 
% Example: featureseparator = 3 works for RAMSES 3.3nm input, 10 for WISP-3 measurements, 5 for 2-nm data. 
%
% 
% [Please retain the following traceback notice in your code]
% Version: 20130410.1
% Adaptation: -
%
% This code is the implementation of the 'fingerprint' method to derive Rrs from hyperspectral (ir)radiance measurements:
% Simis, S.G.H. and J Olsson. Unattended processing of shipborne hyperspectral reflectance measurements. Remote Sensing of Environment, In press (Apr 2013).
%
% <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB"><img src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png"></a>
% This work is licensed under a <a href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>
%
% Rflex hardware/software are described at http://sourceforge.net/p/rflex/wiki/Home/
% 
% See also rho_fingerprint_optimize


% nfeatures can be optimized if none or too many results are generated (default = 20)
nfeatures = 20; 

Ltdiff = diff(Lt);
[tmp,ind1] = sort(Ltdiff);
Lsdiff = diff(Ls);
[tmp,ind2] = sort(Lsdiff);
droplo10 = intersect(ind1(1:nfeatures),ind2(1:nfeatures));
drophi10 = intersect(ind1(end-nfeatures-1:end),ind2(end-nfeatures-1:end));

%intersect results are sorted so the diff product will reveal whether
%there are adjacent indices in it. Of every adjacent series, the first will
%have a high diff product, subsequent ones will be 1 or very low. Remove
%these latter ones as instructed through featureseparator variable.

ind = find(diff(droplo10)<=featureseparator);
droplo10(ind) = [];
ind = find(diff(drophi10)<=featureseparator);
drophi10(ind) = [];
ind = [droplo10;drophi10];

%exclude edges of spectrum
tmp = find(ind<featureseparator | ind>numel(Lt)-featureseparator-1);
ind(tmp) = [];

%output
indices = ind;