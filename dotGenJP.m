function [stimDim, h] = dotGenJP(varargin)

% dotGenJP creates magnitude dimension parameters for constructing dot
%   arrays systematically across a 3D space defined by number (N), size
%   (Sz), and spacing (Sp), from the framework developed by DeWind, Adams,
%   Platt, & Brannon (2015), http://dx.doi.org/10.1016/j.cognition.2015.05.016.
%
% [stimDim, h] = dotGenJP(PARAMS, VAL) acceps 'params' and 'val' pairs to
%   produce the stimulus dimension parameters 'stimDim' and the figure
%   handle 'h' for graphically illustrating the parameters on a 3D space. 
%
%   Possible params are:
%       nlim     - [min max] values of number; default is [8 32].
%       rdlim    - [min max] values of r_d, radius of each dot; default is [6 12].
%       rflim    - [min max] values of r_f, radius of the field; default is [120 240].
%       level    - number of levels across N, Sz, and Sp; default is 5.
%       showplot - logical (true or false) to show 3D plot; default is false.
%                  The figure shows both theoretical values and the real
%                  values accounting for rounding error (due to the
%                  diameter following the unit of pixels) of the
%                  dimensional parameters, and the "mean squared error"
%                  value between the two sets of values. 
%
% This function constructs dot arrays by defining Sz as the dimension that
%   varies in area (both individual area and total area) while holding
%   number constant, as was done in Exp. 1 of Park, DeWind, Woldorff, &
%   Brannon (2016), https://doi.org/10.1093/cercor/bhv017, as well as many
%   other studies such as Park (2018), https://doi.org/10.1016/j.dcn.2017.02.011.
%
% For a balanced design, the minimum and the maximum values of nlim, rdlim,
%   and rflim should all the equivalent. To achieve this, because area as
%   in individual area (IA) and field area (FA) is proportionate to the
%   square of radius, if the range of n differs in 4 folds, the range of
%   rd (radius of each dot) and rf (radius of the circular field in which
%   the dots are drawn) should differ in 2 folds. 
% 
% If one wishes to construct an unbalanced and/or skewed design, it is
%  recommended to construct a balanced design first and then hand pick the
%  dimensional parameters to achieve the individual design goal. 
%
% Joonkoo Park August 2021
%

if mod(nargin,2) == 1
    error('Properties must be in STRING-VALUE pairs. See help.');
end

for ivar = 1 : 2 : length(varargin)
    if ~ischar(varargin{ivar})
        error('Invalid property names. Properties must be in STRING-VALUE pairs. See help.');
    end
    switch lower(varargin{ivar})
        case 'nlim'
            nlim = varargin{ivar+1};
        case 'rdlim'
            rdlim = varargin{ivar+1};
        case 'rflim'
            rflim = varargin{ivar+1};
        case 'level'
            level = varargin{ivar+1};
        case 'showplot'
            showplot = varargin{ivar+1};
        otherwise
            error(strcat('Invalid property name, ',varargin{ivar},'.'));
    end
end

% default values when no options are given
if ~exist('nlim','var')
    nlim = [8 32];
end
if ~exist('rdlim','var')
    rdlim = [6 12];
end
if ~exist('rflim','var')
    rflim = [120 240];
end
if ~exist('level','var')
    level = 5;
end
if ~exist('showplot','var')
    showplot = false;
end

% function for removing rounding error
fRound = @(x,p) round(x .* p) ./ p;

% Number (theoretical and real values)
N = 2 .^ linspace(log2(nlim(1)), log2(nlim(2)), level);
N_r = fRound(N,1);

% IA
IA = 2 .^ linspace(log2(pi*rdlim(1)^2),log2(pi*rdlim(2)^2),level);
rd = sqrt(IA/pi);
rd_r = fRound(rd,2);

% FA
FA = 2 .^ linspace(log2(pi*rflim(1)^2),log2(pi*rflim(2)^2),level);
rf = sqrt(FA/pi);
rf_r = fRound(rf,2);

% empty matrix to store magnitude values
magval = nan(3*level,5);    % theoretical values; columns: N, IA, FA, TA, Spar
magval_r = nan(3*level,3);  % real values; columns: N, rd, rf
counter = 1;

for iIA = 1 : level
    for iFA = 1 : level
        for iN = 1 : level
            
            % theoretical values
            magval(counter, 1) = N(iN);
            magval(counter, 2) = IA(iIA);
            magval(counter, 3) = FA(iFA);
            magval(counter, 4) = IA(iIA) * N(iN);  % TA
            magval(counter, 5) = FA(iFA) ./ N(iN); % Spar
            
            % real values
            magval_r(counter, 1) = N_r(iN);
            magval_r(counter, 2) = rd_r(iIA);
            magval_r(counter, 3) = rf_r(iFA);
            
            counter = counter + 1;
        end
    end
end

magvalLog = log2(magval);

logN  = magvalLog(:,1);
logSz = magvalLog(:,2) + magvalLog(:,4);
logSp = magvalLog(:,3) + magvalLog(:,5);

% round up to 10e-3 place
logN  = fRound(logN,100);
logSz = fRound(logSz,100);
logSp = fRound(logSp,100);

% min/max Sz
midTA = IA(1) * N(end);
minlogSz = log2(midTA) + log2(IA(1));
maxlogSz = log2(midTA) + log2(IA(end));

% min/max Sp
midSpar = FA(1) ./ N(1);
minlogSp = log2(midSpar) + log2(FA(1));
maxlogSp = log2(midSpar) + log2(FA(end));

% remove rounding error
minlogSz = fRound(minlogSz,100);
maxlogSz = fRound(maxlogSz,100);
minlogSp = fRound(minlogSp,100);
maxlogSp = fRound(maxlogSp,100);

% index within min and max of logSZ and logSp
idx = (logSz >= minlogSz & logSz <= maxlogSz) & (logSp >= minlogSp & logSp <= maxlogSp);

logN  = logN(idx);
logSz = logSz(idx);
logSp = logSp(idx);

magval = magval(idx,:);
magval_r = magval_r(idx,:);

% logN, logSz, and logSP computed from real values
logN_r  = log2(magval_r(:,1));
logSz_r = log2(pi * magval_r(:,2).^2) + log2(pi * magval_r(:,2).^2 .* magval_r(:,1));
logSp_r = log2(pi * magval_r(:,3).^2) + log2(pi * magval_r(:,3).^2 ./ magval_r(:,1));

if showplot
    h = figure('Position',round(get(0,'screensize')/2));

    subplot(1,2,1);
    scatter3(logSz, logSp, logN);
    xlabel('logSz'); ylabel('logSp'); zlabel('logN');
    title('Theoretical values');
    axis equal;

    msqerr = mean(sqrt( sum(([logN, logSz, logSp] - [logN_r, logSz_r, logSp_r]) .^ 2, 2) ));

    subplot(1,2,2);
    scatter3(logSz_r, logSp_r, logN_r);
    xlabel('logSz'); ylabel('logSp'); zlabel('logN');
    title(sprintf('Real values; msqerr = %.4f',msqerr));
    axis equal;
else
    h = [];
end

% output struct
stimDim.desc = [sprintf('Generated from %s.m at %s. \n',mfilename,datestr(now)), ...
    sprintf('logN/Sz/Sp and magval are theoretical values. \n'), ...
    sprintf('The columns for magval represent: N, IA, FA, TA, Spar. \n'), ...
    sprintf('logN_r/Sz_r/Sp_r and magval_r are real values. \n'), ...
    sprintf('The columns for magval_r represent: N, r_d, r_f. \n'), ...
    ];
    
stimDim.idx = 1:length(logN);

stimDim.logN   = logN;
stimDim.logSz  = logSz;
stimDim.logSp  = logSp;
stimDim.magval = magval;

stimDim.logN_r   = logN;
stimDim.logSz_r  = logSz;
stimDim.logSp_r  = logSp;
stimDim.magval_r = magval_r;


end
