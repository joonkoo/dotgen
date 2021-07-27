
% Create stimulus parameters for a systematic construction of dot arrays
stimDim = dotGenJP();

% Magnitude values containing information about N, r_d, and r_f
magval_r = stimDim.magval_r;

% Buffer allows extra space between the dots
buffer = 1.5;

N = 500;   % width/height of the image in pixels
AxisAdj = ceil(N/2);

nDots = 1;    % number of unique dot arrays generated for each stimulus parameter point

% Struct that will contain all the dot-array parameters
dotArrays = struct([]);

for i = 1 : size(magval_r,1)

    for j = 1 : nDots
        
        % dot coordinates
        for ierr = 1 : 25
            [dts, err] = dotField2GKA(magval_r(i,2) * buffer * ones(magval_r(i,1),1), magval_r(i,3));
            if err == 0
                break;
            end
            if ierr == 25
                warning('dotField2GKA failed.');
            end
        end
        dts = dts + AxisAdj;
        
        % create a background image that's "mulfac" times larger (for anti-aliasing)
        mulfac = 8;
        
        % M = zeros( mulfac*(N-1)+1 );
        M = zeros( mulfac * N );
        dts_aa = dts * mulfac;
        
        M( sub2ind(size(M),dts_aa(:,1),dts_aa(:,2)) ) = 1;
        J = double(bwdist(M) <= magval_r(i,2) * mulfac);
        
        % reduce the image by "mulfac" times and trim <0 and >1 values
        J = imresize(J, 1/mulfac);
        J(J(:)<0) = 0;
        J(J(:)>1) = 1;
        
        dotArrays(i,j).logN  = stimDim.logN(i);
        dotArrays(i,j).logSz = stimDim.logSz(i);
        dotArrays(i,j).logSp = stimDim.logSp(i);
        dotArrays(i,j).num   = magval_r(i,1);
        dotArrays(i,j).r_d   = magval_r(i,2);
        dotArrays(i,j).r_f   = magval_r(i,3);
        dotArrays(i,j).coord = dts - AxisAdj;   % original coordinates
        dotArrays(i,j).img   = J;
        
%         % Save the image as a bmp file
%         fn = sprintf('array_idx_%02g_id_%03g.bmp',i,j);
%         imwrite(J,fullfile('bmps', fn));
        
        figure;
        imshow(J)
        axis square;
        
    end
end
