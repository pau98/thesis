%% 0.1. Initialization
offset = 5; %% Mask offset
nFrames = 60;  %% Number of images

%% 0.2. Segment organoids
% Segment organoids using information from the different channels and store
% them as be tif files
pre_masks(nFrames,offset)

%% 0.3. Contrast enhancing and merging
% Enhance contrast and j
% Join the chosen channels for PIV
pre_contrast(nFrames,4)

%% 0.4. RUN PIVLab and export the results to Workspace
% Import contrast enhanced collagen files (h_t*.tif)
% Import masks (m_t*.tif)
% Analyse (obtain velocities)
% Perform a velocity based validation
% Export results

%% 0.5. Inicialize
load("PIVlab.mat","u_filtered","v_filtered");
szv = size(v_filtered);
nfrms = szv(1); % Obtain number of frames
I0 = imread('h_t001.tif');
[HIm, WIm] = size(I0);
[H, W] = size(v_filtered{1});
dt = 1;
nFrames = nfrms + 1;

%% 0.6. Center velocities and drift correction
drift = zeros(nfrms+1,2);
for frame = 1:nfrms
    drift(frame+1,1) = drift(frame,1)+median([u_filtered{frame,1}(1,:) u_filtered{frame,1}(end,:) u_filtered{frame,1}(2:end-1,1)' u_filtered{frame,1}(2:end-1,end)']);
    drift(frame+1,2) = drift(frame,2)+median([v_filtered{frame,1}(1,:) v_filtered{frame,1}(end,:) v_filtered{frame,1}(2:end-1,1)' v_filtered{frame,1}(2:end-1,end)']);
end

% ---- crop to the common (valid-in-every-frame) field of view --------
% After translation, each frame has a blank border where the shift
% revealed pixels outside the original image. Find the worst-case
% margin needed in each direction across the WHOLE stack and apply that
% single fixed crop to every frame, so no frame retains any fill-value
% ("out of frame") pixels.
sx = -drift(:,1);   % actual x-translation applied to each frame
sy = -drift(:,2);   % actual y-translation applied to each frame

cropLeft   = max(0, ceil(max( sx)));
cropRight  = max(0, ceil(max(-sx)));
cropTop    = max(0, ceil(max( sy)));
cropBottom = max(0, ceil(max(-sy)));

rowRange = (cropTop+1)  : (HIm-cropBottom);
colRange = (cropLeft+1) : (WIm-cropRight);

if isempty(rowRange) || isempty(colRange)
    error(['Drift is larger than the frame size -- cropLeft/Right/Top/Bottom ' ...
           'leave no valid region. Check for a bad centroid fit in some frame ' ...
           '(see fitRadii for outliers) before trusting this result.']);
end

alignedMasks  = false(HIm, WIm, nFrames);
alignedImages = zeros(HIm, WIm, nFrames);
for t = 1:nFrames
    if any(isnan(drift(t,:)))
        warning('Frame %d: skipping alignment (No valid drift value for this frame).', t);
        continue;
    end
    I = im2double(imread(sprintf('h_t%03d.tif',t)));

    % imtranslate takes [dx dy] (x=col, y=row); shift by MINUS the
    % measured drift to bring each frame back onto the reference position
    fileFramePath = sprintf('dc_h_t%03d.tif', t);
    IDriftCorrected = imtranslate(I, [-drift(t,1), -drift(t,2)], 'FillValues', 0);
    imwrite(im2uint16(IDriftCorrected(rowRange, colRange)), fileFramePath, 'Compression', 'none');

    M = imread(sprintf('m_t%03d.tif',t));
    if ~islogical(M), M = M > 127; end   % 1-bit masks (written by pre_masks) are read as logical
    maskFramePath = sprintf('dc_m_t%03d.tif', t);
    MDriftCorrected = imtranslate(M, [-drift(t,1), -drift(t,2)], 'FillValues', 0);
    imwrite(uint8(MDriftCorrected(rowRange, colRange))*255, maskFramePath, 'Compression', 'none');
end

%% 0.7. RUN PIVLab and export the results to Workspace
% Import cropped collagen files (dc_h_t*.tif)
% Import cropped masks (dc_m_t*.tif)
% Analyse (obtain velocities)
% Perform a velocity based validation
% Export results


%% 1.0. Analyse PIV results
clear variables;
load("dc_PIVlab.mat","u_filtered","v_filtered","x","y");
szv = size(v_filtered);
nfrms = szv(1); % Length of velocity time steps (should be one less than nFrames)
I0 = imread('dc_h_t001.tif');
[HIm, WIm] = size(I0);
[H, W] = size(v_filtered{1});
dt = 1;
nFrames = nfrms + 1;
%%
analizePIVresults(nFrames, x, y, nfrms, u_filtered, v_filtered, H, W, WIm, HIm, dt);

%% 1.1 Load analised info
load("PIV_results.mat")
dx = 0.2;   % microns/pixel
dt = 20;    % mins/frame (analizePIVresults above is called with dt = 1, i.e. per frame)

%% 1.2 Compute tight mask (for plots)
pnPlots = perimeter_normals_NOd(nFrames,0);

%% 2.5. Compute central point and Normal stretch
r = zeros(1,nFrames); % Preallocate radius vector
NCN = cell(nFrames,1);

for t = 1:nFrames
    Mask = imread(sprintf('dc_m_t%03d.tif',t));
    Mask = Mask > 127;
    szM = size(Mask);
    A = nnz(Mask);   % pixel count
    r(t) = sqrt(A / pi); % Organoid
    B = bwboundaries(Mask);
    [~, biggest] = max(cellfun(@(c) size(c,1), B));
    boundary = B{biggest};
    xy = [boundary(:,2), boundary(:,1)];
    [cx(t), cy(t), rFit(t)] = fitCircleLS(xy);
    NCN{t} = nX{1}.*C{t,1}.*nX{1}+2*nX{1}.*C{t,2}.*nY{1}+nY{1}.*C{t,3}.*nY{1};
end

%% Radial displacement 
iiframes = linspace(2,nfrms,10);
figure('Name', 'Total radial displacement vs initial distance to boundary');
hold on
for ii = 1:10
    frame =  uint16(iiframes(ii));
    [binCenters, meanRadialDisp, semRadialDisp, nPerBin] = radial_displacement_profile_lagrangian(p0n{frame,1}, p0n{frame,2}, X, Y, disp{frame,1}, disp{frame,2});
    [~, ~, r0] = fitCircleLS([p0n{frame,1} p0n{frame,2}]);

    errorbar((binCenters+r0)*dx, meanRadialDisp*dx, semRadialDisp*dx, 'o-', 'LineWidth', 1.3);
    xline(r0*dx, '--', 'Color', [0.5 0.5 0.5]);
    
end
xlabel('INITIAL distance to organoid boundary (microns)');
    ylabel('mean radial displacement (microns)');
    title('Radial displacement profile');
    grid on

%% Test Smoothing
x = X(1,:);
y = Y(:,1);
outSM = cell(nFrames,1);
for t = 1:nFrames
[in, on] = inpolygon(X, Y, p0n{t,1}, p0n{t,2});
hole = in | on;              % also remove nodes lying exactly on the boundary
U = disp{t,1};
V = disp{t,2};
U(hole) = NaN;
V(hole) = NaN;
outSM{t} = defgrad_pls(x, y, U, V, 1);
end

%%
t = 55
CNN = nX{1}.*outSM{t}.C11.*nX{1} + 2*nX{1}.*outSM{t}.C12.*nY{1} + nY{1}.*outSM{t}.C22.*nY{1};
CNN2 = nX{1}.*C{t,1}.*nX{1} + 2*nX{1}.*C{t,2}.*nY{1} + nY{1}.*C{t,3}.*nY{1};

figure()
subplot(1,2,1)
imagesc(sqrt(CNN))
clim([0.5 1.5])
axis equal tight

subplot(1,2,2)
imagesc(sqrt(CNN2))
clim([0.5 1.5])
axis equal tight

%% Compute fth and traction
% dt = 30;
% dx = 0.346;
% dt = 20;
% dx = 0.206;
params.mu = 0.03; % kPa
params.lm = 1.67; % kPa
params.m = 1;     % kPa
params.ms = 15;
params.mg = 0.01;
params.alp = 10;
tnsn = cell(nFrames,1);
fth = cell(nFrames,1);
for t = 1:nFrames
    C1 = outSM{t}.C11;
    C12 = outSM{t}.C12;
    C2 = outSM{t}.C22;
    [s1, s2] = size(C1);
    W = zeros(s1,s2);
    snn = zeros(s1,s2);
    fthM = zeros(s1,s2);

    for ii = 1:s1
        for jj = 1:s2
            if isnan(C1(ii,jj))
                fthM(ii,jj) = NaN;
            else
            CM = [C1(ii,jj) C12(ii,jj); C12(ii,jj) C2(ii,jj)];
            [WNH,SNH] = NH(CM,params);
            [WSS,SSS] = SS2(CM,params);
            [WMC,SMC] = MC(CM,params);
            if (params.mg > 0)
            W(ii,jj) = WNH + WSS + WMC ;
            S = SNH + SSS + SMC;
            else
            W(ii,jj) = WNH + WSS ;
            S = SNH + SSS;
            end
            % W(ii,jj) = WNH;
            % S1(ii,jj) = SNH ;
            FM = [outSM{t}.F11(ii,jj) outSM{t}.F12(ii,jj); outSM{t}.F21(ii,jj) outSM{t}.F22(ii,jj)];
            J = outSM{t}.J(ii,jj);
            s = FM*S*FM'./J;
            fthM(ii,jj) = -W(ii,jj)/J + nX{1}(ii,jj)*s(1,1)*nX{1}(ii,jj)+nX{1}(ii,jj)*s(1,2)*nY{1}(ii,jj)+...
                nY{1}(ii,jj)*s(2,1)*nX{1}(ii,jj)+nY{1}(ii,jj)*s(2,2)*nY{1}(ii,jj);
            snn(ii,jj) = nX{1}(ii,jj)*s(1,1)*nX{1}(ii,jj)+nX{1}(ii,jj)*s(1,2)*nY{1}(ii,jj)+...
                nY{1}(ii,jj)*s(2,1)*nX{1}(ii,jj)+nY{1}(ii,jj)*s(2,2)*nY{1}(ii,jj);
            end
        end
    end
    tnsn{t} = snn;
    fth{t} = fthM;
end

%% Polar plot of degradation velocity against different variables PIV
plotMode = 3;
tStps = nFrames;
vdWind = 5;
for t = tStps
[bandAxis, fthAlongBand, bandXq, bandYq, interfaceX, interfaceY, info] = field_in_band(p0n{t,1}, p0n{t,2}, X, Y, fth(t), 50, 100,'Aggregate','gaussian','Component','radial');
[bandAxis2, NCNAlongBand, bandXq2, bandYq2, interfaceX2, interfaceY2, info2] = field_in_band(p0n{t,1}, p0n{t,2}, X, Y, C(t,:), 50, 100,'Aggregate','gaussian','Component','radial');
[bandAxis3, tnAlongBand] = field_in_band(p0n{t,1}, p0n{t,2}, X, Y, tnsn(t,:), 50, 100,'Aggregate','gaussian','Component','radial');
[vd_avg, theta_g] = angular_time_average(vd, p0n, cx(1), cy(1), t, vdWind, 100);

[xg, yg] = pol2cart(theta_g, 1);
xg = xg + cx(1);  yg = yg + cy(1);
% two fields on different curves (e.g. two contours), same centre
opts.ZeroRadius = [0.6 0.2 0.3];
opts.Fill = true;
opts.FillMode = 'value';
opts.ColorBy = 'angle';
opts.Names = {'fth', 'f_t_h','v_d(\mu m/h)'};
opts.StarTol = 4;
opts.Fold = 'average';
opts.RScale     = { 'linear','linear','linear'};   % inner ring linear, outer ring symlog
opts.LinThresh  = [NaN,5, NaN];              % C for symlog, in units of f (NaN = auto)
opts.Colormap =  {blueWhiteRed(256),blueWhiteRed(256), parula(256)};
if (plotMode == 1)
    opts.CLim = [ -0.2 0.2; -0.2 0.2; -0.5 8.7];
    opts.RLim =  [ -0.2 0.2; -0.2 0.2; -0.5 8.7];    % field 2
else
    % opts.CLim = [ -0.2 0.2; -350 350; -0.5 8.7];
    % opts.RLim =  [ -0.2 0.2; -10 350; -0.5 2];    % field 2
    opts.CLim = [ -0.2 0.2; -30 30; -0.5 8.7];
    opts.RLim =  [ -0.2 0.2; -0.5 2.5; -0.5 2];    % field 2
end
opts.ShowRaw = false;
opts.LineColor = 'value';    % or any colour, e.g. 'k' (default), [0.2 0.2 0.2]

if (plotMode == 1)
    radialFieldPlot({interfaceX, interfaceX, xg}, {interfaceY, interfaceY, yg}, {0.*NCNAlongBand', sqrt(NCNAlongBand')-1, vd_avg/dt*60*dx}, cx(1), cy(1), opts);
    fontname("Cambria")          % R2023a and later
    fontsize(scale=1.2)
    exportgraphics(gcf, sprintf('vd_NCN_pol_t%03d.pdf',t), 'ContentType', 'vector')
    exportgraphics(gcf, sprintf('vd_NCN_img_pol_t%03d.png',t))
    close all
elseif  (plotMode == 2)
    radialFieldPlot({interfaceX, interfaceX, xg}, {interfaceY, interfaceY, yg}, {0.*NCNAlongBand', tnAlongBand', vd_avg/dt*60*dx}, cx(1), cy(1), opts);
    fontname("Cambria")          % R2023a and later
    fontsize(scale=1.2)
    exportgraphics(gcf, sprintf('vd_fth_pol_t%03d.pdf',t), 'ContentType', 'vector')
    exportgraphics(gcf, sprintf('vd_fth_img_pol_t%03d.png',t))
    close all
else
    radialFieldPlot({interfaceX, interfaceX, xg}, {interfaceY, interfaceY, yg}, {0.*NCNAlongBand', fthAlongBand', vd_avg/dt*60*dx}, cx(1), cy(1), opts);
    fontname("Cambria")          % R2023a and later
    fontsize(scale=1.2)
    exportgraphics(gcf, sprintf('vd_fth_pol_t%03d.pdf',t), 'ContentType', 'vector')
    exportgraphics(gcf, sprintf('vd_fth_img_pol_t%03d.png',t))
    close all
end
end

%% Masks
function pre_masks(nframes,offset)
    for frame= 1:nframes
        c3=imread(sprintf('z_t%03d_c002.tif',frame));
        % c4=imread(sprintf('c11f%03d.tif',frame));
        rdisk = 15;
        se = strel('disk',rdisk);
        c3 = imadjust(c3);
        background = imopen(c3,se);
        I2 = c3 - background;
        % Sum channels and enhance contrast 
        I2=imadjust(I2);
        gfv= 10;
        I = imgaussfilt(I2,gfv);
        Imax = max(max(I));
        I(I > 0.2*Imax) = 0.2*Imax;
        Iadj=imadjust(I);
        imshow(Iadj)
        
        T = imbinarize(I);
        H= imfill(T,'holes');
        
        CC = bwconncomp(H);
        S = regionprops(CC, 'Area');
        L = labelmatrix(CC);
        Mask = ismember(L, find([S.Area] >= 3000));

        % DILATE AND FIND PERIMETER        
        % dilate mask
        SEdil=strel('disk',offset,0);
        Mask = imdilate(Mask,SEdil);
        Mask = imfill(Mask,'holes');
        imwrite(Mask,sprintf('m_t%03d.tif',frame));
    end
end

%% Collagen contrast
function pre_contrast(nframes,chnnl)
    nCh = length(chnnl);
    for frame=1:nframes
        Iadj=imread(sprintf('z_t%03d_c%03d.tif',frame,chnnl(nCh))).*0; 
        for ch = 1:nCh
            % Import channel image
            I=imread(sprintf('z_t%03d_c%03d.tif',frame,chnnl(ch)));
            
            % Enhance contrast 
            I = imadjust(I);
            Iadj = Iadj + imflatfield(I,30)/nCh;
        end
        imwrite(imadjust(Iadj),sprintf('h_t%03d.tif',frame));
    end
end

%% Boundary normal velocity
function pn = perimeter_normals(nframes,offset)
    pn = cell(nframes, 4);
    SEdil = strel('disk',offset,0);
    prevPt = [0 0];
    for frame=1:nframes        
        T = imread(sprintf('dc_m_t%03d.tif',frame)); % Read mask       
        
        Mask = T > 127;
        Mask = imfill(Mask, 'holes');
        
        % DILATE AND FIND PERIMETER        
        % Dilate mask according to offset
        D = imdilate(Mask,SEdil);
        
        % Find perimeter (organoid-collagen interface)
        Bn = bwboundaries(D);
        [~, biggest] = max(cellfun(@(c) size(c,1), Bn));
        boundaryPts = Bn{biggest};        % [row col], closed loop (first==last point)
        boundaryPts = boundaryPts(1:end-1, :);   % drop the closing duplicate
        x = boundaryPts(:,2);
        y = boundaryPts(:,1);
        k = (1:numel(x))';               % already an ordered boundary trace
        [~,pMin] = min((x-prevPt(1)).^2+(y-prevPt(2)).^2);
        prevPt = [x(pMin) y(pMin)];
        x = circshift(x,-pMin+1);
        y = circshift(y,-pMin+1);
        % Take a subset of coordinates; calculate tangent t = [tx,ty] and
        % normal unit vector n=[nx,ny]
        
        perx = x(k(1:5:length(k)));
        pery = y(k(1:5:length(k)));
                
        dx = circshift(perx,-1)-circshift(perx,1);
        dy = circshift(pery,-1)-circshift(pery,1);
        nrmD = sqrt(dx.*dx+dy.*dy);
        tx = dx./nrmD;
        ty = dy./nrmD;
        
        nx = ty;
        ny = -tx;
        
        % Generate cell with offset perimeter and normal vectors
        pn(frame,:) = {perx,pery,nx,ny};
    end
end

function pn = perimeter_normals_NOd(nframes,offset)
    pn = cell(nframes, 4);
    SEdil = strel('disk',offset,0);
    prevPt = [0 0];
    for frame=1:nframes        
        T = imread(sprintf('dc_m_t%03d.tif',frame)); % Read mask       
        
        Mask = T > 127;
        Mask = imfill(Mask, 'holes');
        
        % DILATE AND FIND PERIMETER        
        % Dilate mask according to offset
        D = imdilate(Mask,SEdil);
        
        % Find perimeter (organoid-collagen interface)
        Bn = bwboundaries(D);
        [~, biggest] = max(cellfun(@(c) size(c,1), Bn));
        boundaryPts = Bn{biggest};        % [row col], closed loop (first==last point)
        boundaryPts = boundaryPts(1:end-1, :);   % drop the closing duplicate
        x = boundaryPts(:,2);
        y = boundaryPts(:,1);
        k = (1:numel(x))';               % already an ordered boundary trace
        [~,pMin] = min((x-prevPt(1)).^2+(y-prevPt(2)).^2);
        prevPt = [x(pMin) y(pMin)];
        x = circshift(x,-pMin+1);
        y = circshift(y,-pMin+1);
        % Take a subset of coordinates; calculate tangent t = [tx,ty] and
        % normal unit vector n=[nx,ny]
        
        perx = x(k(1:5:length(k)));
        pery = y(k(1:5:length(k)));
                
        dx = circshift(perx,-1)-circshift(perx,1);
        dy = circshift(pery,-1)-circshift(pery,1);
        nrmD = sqrt(dx.*dx+dy.*dy);
        tx = dx./nrmD;
        ty = dy./nrmD;
        
        nx = ty;
        ny = -tx;
        
        % Generate cell with offset perimeter and normal vectors
        pn(frame,:) = {perx,pery,nx,ny};
    end
end


%%
function [a, b, r] = fitCircleLS(xy)
% FITCIRCLELS  Kasa (1976) algebraic least-squares circle fit.
    x = xy(:,1); y = xy(:,2);
    A = [2*x, 2*y, ones(size(x))];
    d = x.^2 + y.^2;
    p = A \ d;
    a = p(1); b = p(2); c = p(3);
    r = sqrt(c + a^2 + b^2);
end
%%
function [binCenters, meanRadialDisp, semRadialDisp, nPerBin] = radial_displacement_profile_lagrangian(perx0, pery0, X0, Y0, dispX, dispY, varargin)
%RADIAL_DISPLACEMENT_PROFILE_LAGRANGIAN  Mean TOTAL radial displacement of
% material points (tracked from a fixed initial Lagrangian grid), binned
% by each point's distance to the organoid BOUNDARY (given directly as
% boundary points, not a mask) at the initial frame.
%
%   [binCenters, meanRadialDisp, semRadialDisp, nPerBin] = ...
%       radial_displacement_profile_lagrangian(perx0, pery0, X0, Y0, dispX, dispY)
%
% INPUTS
%   perx0, pery0 : organoid boundary points at the REFERENCE frame the
%                  Lagrangian grid was defined on (e.g. pn{1,1}, pn{1,2}
%                  from perimeter_normals.m)
%   X0, Y0       : initial grid point positions (fixed, one value per
%                  tracked point)
%   dispX, dispY : TOTAL (cumulative) displacement of each point, for the
%                  frame you want to plot. ADAPT THE INDEXING to your
%                  actual `disp` array -- e.g. if disp is
%                  (nPoints, dim, nFrames): dispX = disp(:,1,frame);
%                  dispY = disp(:,2,frame); adjust to match your layout.
%
% NAME-VALUE OPTIONS
%   'BinWidth'    distance bin width, in the same units as X0,Y0 (default: 20)
%   'MaxDistance' largest distance to include (default: max observed)
%
% OUTPUT: same layout as before (binCenters, meanRadialDisp,
% semRadialDisp, nPerBin).
%
% METHOD: distance-to-boundary and the local outward normal are computed
% directly from the boundary POINTS (exact point-to-polyline distance +
% inpolygon sign test, with the normal available analytically as
% (query - nearest boundary point)/distance) -- no mask/rasterization
% step at all, unlike the earlier mask-based version of this function.
% Validated numerically against the raster (poly2mask+bwdist) approach
% on a synthetic circle: they agree to within ~0.05-1 px, consistent
% with expected discretization differences (see chat).
%
% Computed ONCE from the fixed initial boundary/grid, same as before --
% no need to repeat this per frame.
%
% REFERENCE-CONFIGURATION CHOICE (see chat, unchanged from before): this
% bins by each point's INITIAL distance to the boundary and projects
% total displacement onto the INITIAL local outward normal. If you want
% the CURRENT-frame-relative version instead, that needs the current
% frame's boundary points and current positions (X0+dispX, Y0+dispY) --
% ask if you want that variant built too.
%
% NOTE: if `disp` was built by cumulatively summing PIV increments that
% haven't been cleaned up yet, run them through denoise_piv_displacement.m
% and piv_background_drift.m BEFORE summing -- this profile directly
% inherits whatever noise/drift is already baked into `disp`.

    p = inputParser;
    addParameter(p, 'BinWidth', 20);
    addParameter(p, 'MaxDistance', []);
    parse(p, varargin{:});
    binWidth = p.Results.BinWidth;
    maxDistance = p.Results.MaxDistance;

    % ---- distance-to-boundary and local outward normal, computed ONCE
    % directly from the boundary points (no mask needed) -----------------
    [distToBoundary0, nx0, ny0] = phiFromBoundary(perx0, pery0, X0, Y0);

    % ---- radial component of TOTAL displacement, projected onto the
    % INITIAL local outward normal ----------------------------------------
    radialDisp = dispX(:) .* nx0(:) + dispY(:) .* ny0(:);

    valid = distToBoundary0(:) > 0 & ~isnan(distToBoundary0(:)) & ~isnan(radialDisp);
    d = distToBoundary0(valid);
    r = radialDisp(valid);

    if isempty(maxDistance)
        maxDistance = max(d);
    end
    edges = 0:binWidth:maxDistance;
    binCenters = edges(1:end-1) + binWidth/2;

    meanRadialDisp = nan(size(binCenters));
    semRadialDisp  = nan(size(binCenters));
    nPerBin = zeros(size(binCenters));
    for i = 1:numel(binCenters)
        inBin = d >= edges(i) & d < edges(i+1);
        nPerBin(i) = nnz(inBin);
        if nPerBin(i) > 0
            meanRadialDisp(i) = mean(r(inBin));
            semRadialDisp(i)  = std(r(inBin)) / sqrt(nPerBin(i));
        end
    end
end
%%
function [phi, nx, ny] = phiFromBoundary(perx, pery, X, Y)
%PHIFROMBOUNDARY  Signed distance (phi) from a set of query points (X,Y)
% to a closed polygon defined by boundary points (perx, pery), computed
% EXACTLY (point-to-polyline distance)
%
%   phi = phiFromBoundary(perx, pery, X, Y)
%   [phi, nx, ny] = phiFromBoundary(perx, pery, X, Y)
%
% INPUTS
%   perx, pery : boundary points (e.g. the perx,pery output of
%                perimeter_normals.m / pn{frame,1}, pn{frame,2}),
%                assumed to trace a closed polygon (does NOT need to be
%                explicitly closed -- the first/last point are connected
%                automatically)
%   X, Y       : query points, any shape (scalar, vector, or a full
%                meshgrid if you want a dense field) -- e.g. your
%                Lagrangian X0,Y0 points
%
% OUTPUT
%   phi    : signed distance, same shape as X,Y. Same sign convention as
%            the rest of this pipeline: phi<0 inside, phi>0 outside.
%   nx,ny  : (optional) outward unit normal direction at each query
%            point, computed ANALYTICALLY as (query - nearest boundary
%            point)/distance -- no gradient/interpolation step needed,
%            unlike the raster+bwdist+gradient approach used elsewhere
%            in this pipeline. Only meaningful for points OUTSIDE the
%            polygon (phi>0); returned as NaN for points with phi<=0,
%            since "outward normal" isn't well-defined there in the
%            same simple sense.
%
% NOTE ON EFFICIENCY: this is O(nQueryPoints x nBoundaryPoints) -- fine
% for a scattered set of query points (e.g. a Lagrangian grid with a few
% hundred/thousand points), but MUCH slower than bwdist if you evaluate
% it over every pixel of a full image (bwdist is a highly optimized
% O(nPixels) algorithm). For a full-image phi field, rasterize with
% poly2mask + bwdist instead (see chat) -- this function is for the
% scattered-point case specifically.

    x1 = perx(:); y1 = pery(:);
    n = numel(x1);
    x2 = [x1(2:end); x1(1)];   % each point's successor, wrapping around
    y2 = [y1(2:end); y1(1)];   % (closes the polygon automatically)

    qx = X(:); qy = Y(:);
    nq = numel(qx);

    d = inf(nq, 1);
    nearestX = nan(nq, 1);
    nearestY = nan(nq, 1);

    for e = 1:n
        dx = x2(e) - x1(e);
        dy = y2(e) - y1(e);
        segLenSq = dx^2 + dy^2;
        if segLenSq == 0
            t = zeros(nq, 1);
        else
            t = ((qx - x1(e))*dx + (qy - y1(e))*dy) / segLenSq;
            t = min(max(t, 0), 1);
        end
        projx = x1(e) + t*dx;
        projy = y1(e) + t*dy;
        segDist = sqrt((qx - projx).^2 + (qy - projy).^2);

        closer = segDist < d;
        d(closer) = segDist(closer);
        nearestX(closer) = projx(closer);
        nearestY(closer) = projy(closer);
    end

    inside = inpolygon(qx, qy, perx, pery);
    phi = d;
    phi(inside) = -phi(inside);
    phi = reshape(phi, size(X));

    if nargout > 1
        nxv = (qx - nearestX) ./ d;
        nyv = (qy - nearestY) ./ d;
        nxv(inside) = NaN;
        nyv(inside) = NaN;
        nx = reshape(nxv, size(X));
        ny = reshape(nyv, size(Y));
    end
end
%%
function [bandAxis, dispAlongBand, bandXq, bandYq, interfaceX,interfaceY] = field_along_band(perx0, pery0, X0, Y0, field, targetDistance, varargin)
%DISPLACEMENT_ALONG_BAND  Plot displacement ALONG a band at a fixed
% distance from the organoid boundary, as a function of angular position
% around the organoid
%
%   [bandAxis, dispAlongBand, bandX, bandY] = ...
%       displacement_along_band(perx0, pery0, X0, Y0, dispX, dispY, targetDistance)
%
% INPUTS
%   perx0, pery0   : organoid boundary points at the reference frame
%                    (e.g. p0n{frame,1}, p0n{frame,2})
%   X0, Y0         : Lagrangian grid point positions (fixed)
%   dispX, dispY   : TOTAL displacement of each point, for the frame you
%                    want to plot
%   targetDistance : the single distance (from the boundary) at which to
%                    extract the band, same units as X0,Y0
%
% NAME-VALUE OPTIONS
%   'Component'   'radial' (default), 'tangential', or 'magnitude' --
%                 which displacement component to plot along the band
%   'NBins'       number of points to resample the band onto, evenly
%                 spaced in angle around the organoid (default: 200)
%   'GridPad'     margin (px) added around the boundary points when
%                 building the internal grid for contour extraction
%                 (default: ceil(targetDistance)+10 -- must be large
%                 enough that the band contour doesn't get clipped)
%
% OUTPUT
%   bandAxis      : 1 x NBins, angle (radians, -pi to pi) around the
%                   organoid's circle-fit center
%   dispAlongBand : 1 x NBins, the requested displacement component,
%                   interpolated from the scattered Lagrangian
%                   displacement data onto the band
%   bandX, bandY  : the band's spatial coordinates (for plotting/QC,
%                   e.g. overlaying the band on the raw image)
%
% METHOD: the band is the exact level set phi(x,y) = targetDistance,
% extracted via contour extraction on a gridded phi field (built internally
% from perx0,pery0 via poly2mask + bwdist) -- this is robust for non-convex
% shapes/larger distances, unlike naively offsetting each boundary point
% outward by a fixed amount, which can self-intersect. The displacement
% field (known only at the scattered Lagrangian points X0,Y0) is
% interpolated onto the band with scatteredInterpolant. The band is
% parametrized by ANGLE around a circle-fit center (same convention as the interface-velocity
% kymograph earlier in this pipeline), not raw contour-trace order,
% since angle gives a consistent, well-defined axis.
%
% CAVEATS carried over from earlier in this pipeline: angle
% parametrization assumes reasonably star-convex geometry around the
% fit center; the local outward normal used for 'radial'/'tangential'
% projection comes from the SAME phi field the band was extracted from,
% so it's consistent by construction, not a separate approximation.

    q = inputParser;
    addParameter(q, 'Component', 'radial');
    addParameter(q, 'NBins', 200);
    addParameter(q, 'GridPad', ceil(targetDistance) + 10);
    parse(q, varargin{:});
    component = q.Results.Component;
    nBins = q.Results.NBins;
    gridPad = q.Results.GridPad;

    % ---- build a local grid around the boundary, just big enough -------
    % (avoids needing the full original image size / a mask at all)
    xMin = min(perx0) - gridPad; xMax = max(perx0) + gridPad;
    yMin = min(pery0) - gridPad; yMax = max(pery0) + gridPad;
    W = ceil(xMax - xMin); H = ceil(yMax - yMin);

    perxShifted = perx0 - xMin;
    peryShifted = pery0 - yMin;
    mask = poly2mask(perxShifted, peryShifted, H, W);
    phi = double(bwdist(mask)) - double(bwdist(~mask));   % phi<0 inside, phi>0 outside
                                                             % (same convention as the rest
                                                             % of this pipeline, e.g.
                                                             % compute_interface_velocity.m)

    % ---- extract the band: the level set phi = targetDistance ----------
    C = contourc(double(phi), [targetDistance targetDistance]);
    % contourc output format: pairs of [level; count] header rows
    % followed by that many [x;y] point columns, possibly repeated for
    % multiple disjoint contour pieces -- take the LONGEST piece
    segments = {};
    i = 1;
    while i < size(C,2)
        cnt = C(2,i);
        segments{end+1} = C(:, i+1:i+cnt); %#ok<AGROW>
        i = i + cnt + 1;
    end
    if isempty(segments)
        error(['No points found at distance %.2f from the boundary. ' ...
               'Increase GridPad, or check targetDistance is within a ' ...
               'reasonable range for this boundary.'], targetDistance);
    end
    [~, biggest] = max(cellfun(@(s) size(s,2), segments));
    band = segments{biggest};   % 2 x N, [x;y] in the SHIFTED grid frame

    bandX = band(1,:)' + xMin;   % shift back to original coordinates
    bandY = band(2,:)' + yMin;

    % ---- cumulative arc length along the (already-ordered) band --------
    segLen = hypot(diff(bandX), diff(bandY));
    s = [0; cumsum(segLen)];
    closingLen = hypot(bandX(1)-bandX(end), bandY(1)-bandY(end));
    totalLength = s(end) + closingLen;   % include the closing segment
 
    % ---- resample onto EQUISPACED arc-length points ---------------------
    sQuery = linspace(0, totalLength, nBins+1);
    sQuery(end) = [];
    bandAxis = sQuery;
 
    sExt = [s; totalLength];
    bandXext = [bandX; bandX(1)];
    bandYext = [bandY; bandY(1)];
    [sExt, ia] = unique(sExt);   % guard against any duplicate arc-length
    bandXext = bandXext(ia);     % values (e.g. a zero-length segment)
    bandYext = bandYext(ia);
 
    bandXq = interp1(sExt, bandXext, bandAxis, 'linear');
    bandYq = interp1(sExt, bandYext, bandAxis, 'linear');
 
    % ---- local outward normal at each resampled band point, from the
    % SAME phi field the band came from (consistent by construction) ----
    % [gx, gy] = gradient(phi);
    % gmag = sqrt(gx.^2 + gy.^2); gmag(gmag==0) = 1;
    % nxField = gx ./ gmag; nyField = gy ./ gmag;
    % nxBand = interp2(1:W, 1:H, nxField, bandXq - xMin, bandYq - yMin, 'linear');
    % nyBand = interp2(1:W, 1:H, nyField, bandXq - xMin, bandYq - yMin, 'linear');
    % figure()
    % quiver(bandXq,bandYq,nxBand,nyBand)
    nyBand2 = circshift(bandXq,-1)-circshift(bandXq,1);
    nxBand2 = -circshift(bandYq,-1)+circshift(bandYq,1);
    nMag = sqrt(nxBand2.*nxBand2+nyBand2.*nyBand2);
    % hold on
    % quiver(bandXq,bandYq,nxBand2./nMag,nyBand2./nMag)
    nxBand = nxBand2./nMag;
    nyBand = nyBand2./nMag;

        % ---- EXACT local outward normal + interface correspondence, via
    % point-to-polyline projection against the ORIGINAL boundary points
    % (replaces the earlier gradient(phi)-based normal -- see chat) -----
    % [nxBand, nyBand, interfaceX, interfaceY] = exactNormalAndNearestPoint(perx0, pery0, bandXq, bandYq);
    [interfaceX,interfaceY]  = exactNearestPoint(perx0, pery0, bandXq, bandYq,nxBand,nyBand);
    % figure()
    % quiver(bandXq,bandYq,nxBand,nyBand)

    % ---- interpolate the scattered Lagrangian displacement field onto
    % the band (scatteredInterpolant handles irregular point clouds) ----
    szFld = size(field);
    if szFld(2) == 1
        fieldV = field{1};
        F = scatteredInterpolant(X0(:), Y0(:), fieldV(:), 'linear', 'none');
        vBand = F(bandXq(:), bandYq(:));
        dispAlongBand = vBand(:);
    elseif szFld(2) == 2
        fieldX = field{1};
        fieldY = field{2};
        Fx = scatteredInterpolant(X0(:), Y0(:), fieldX(:), 'linear', 'none');
        Fy = scatteredInterpolant(X0(:), Y0(:), fieldY(:), 'linear', 'none');
        uBand = Fx(bandXq(:), bandYq(:));
        vBand = Fy(bandXq(:), bandYq(:));
        switch lower(component)
        case 'radial'
            dispAlongBand = uBand(:)'.*nxBand(:)' + vBand(:)'.*nyBand(:)';
        case 'tangential'
            % tangential = normal rotated 90 degrees
            txBand = -nyBand; tyBand = nxBand;
            dispAlongBand = uBand(:)'.*txBand(:)' + vBand(:)'.*tyBand(:)';
        case 'magnitude'
            dispAlongBand = sqrt(uBand(:)'.^2 + vBand(:)'.^2);
            otherwise
            error('Component must be ''radial'', ''tangential'', or ''magnitude''.');
        end
    elseif szFld(2) == 3
        fieldXX = field{1};
        fieldXY = field{2};
        fieldYY = field{3};
        FXX = scatteredInterpolant(X0(:), Y0(:), fieldXX(:), 'linear', 'none');
        FXY = scatteredInterpolant(X0(:), Y0(:), fieldXY(:), 'linear', 'none');
        FYY = scatteredInterpolant(X0(:), Y0(:), fieldYY(:), 'linear', 'none');
        FXXBand = FXX(bandXq(:), bandYq(:));
        FXYBand = FXY(bandXq(:), bandYq(:));
        FYYBand = FYY(bandXq(:), bandYq(:));
        switch lower(component)
        case 'radial'
            dispAlongBand = FXXBand(:)'.*nxBand(:)'.*nxBand(:)' + 2.*FXYBand(:)'.*nxBand(:)'.*nyBand(:)'+FYYBand(:)'.*nyBand(:)'.*nyBand(:)';
        case 'tangential'
            % tangential = normal rotated 90 degrees
            txBand = -nyBand; tyBand = nxBand;
            dispAlongBand = FXXBand(:)'.*txBand(:)'.*txBand(:)' + 2.*FXYBand(:)'.*txBand(:)'.*tyBand(:)'+FYYBand(:)'.*tyBand(:)'.*tyBand(:)';
        case 'magnitude'
            dispAlongBand = FXXBand(:)'.*FYYBand(:)'-FXYBand(:)'.*FXYBand(:)';
            otherwise
            error('Component must be ''radial'', ''tangential'', or ''magnitude''.');
        end
    end

    % 
    % figure('Name', sprintf('Displacement along band at distance %.1f from boundary', targetDistance));
    % subplot(1,2,1);
    % plot(bandAxis, dispAlongBand, 'LineWidth', 1.3);
    % if strcmpi(component,'radial'); yline(0,'--','Color',[0.5 0.5 0.5]); end
    % xlabel('angle around organoid (rad)'); ylabel(sprintf('%s displacement', component));
    % title(sprintf('Band at distance = %.1f', targetDistance));
    % 
    % subplot(1,2,2);
    % plot(perx0, pery0, 'k-', 'LineWidth', 1); hold on;
    % scatter(bandXq, bandYq, 15, dispAlongBand, 'filled');
    % axis equal; set(gca,'YDir','reverse'); colorbar;
    % title('Band location, colored by displacement');
end

function [nearestX, nearestY, tParam] = exactNearestPoint(perx, pery, qx, qy, nx, ny)
%EXACTNEARESTPOINT  For each query point (qx,qy), find the point on the
% closed polyline (perx,pery) that lies along the LINE through (qx,qy)
% in the GIVEN direction (nx,ny) -- i.e. ray/line-to-curve intersection,
% NOT perpendicular-distance projection (that's exactNormalAndNearestPoint,
% used in displacement_along_band.m). The intersection is allowed to fall
% strictly between two consecutive points of (perx,pery), not just at a
% vertex -- it's computed exactly via the line-segment intersection
% formula, not by snapping to the nearest existing point.
%
%   [nearestX, nearestY, tParam] = exactNearestPoint(perx, pery, qx, qy, nx, ny)
%
% INPUTS
%   perx, pery : the closed polyline curve to intersect (closes
%                automatically from the last point back to the first)
%   qx, qy     : query points (ray/line origins)
%   nx, ny     : direction at each query point (need not be unit length;
%                only the direction matters, magnitude cancels out)
%
% OUTPUT
%   nearestX, nearestY : the intersection point CLOSEST to (qx,qy) along
%                         its line (smallest |tParam|), same length as
%                         qx,qy. NaN if that point's line does not
%                         intersect the curve at all.
%   tParam              : signed parameter such that
%                         (nearestX,nearestY) = (qx,qy) + tParam.*(nx,ny)
%                         -- positive if the intersection is in the
%                         direction of (nx,ny), negative if behind it.
%                         Useful as a sanity check that the sign matches
%                         what you expect physically (e.g. should be
%                         consistently positive if (nx,ny) is an outward
%                         normal and the curve is further out).
%
% METHOD (validated against a circle with query points NOT aligned to
% any polygon vertex, confirmed to land correctly between vertices, and
% against hand-computed test cases for the underlying line-segment
% intersection formula -- see chat): for each query point and each
% segment of the curve, solves
%     q + t*[nx,ny] = a + s*(b-a)
% via Cramer's rule for [t,s], keeps intersections with s in [0,1]
% (i.e. actually on that segment), and returns the one with the smallest
% |t| across all segments (closest to the query point along its line, in
% EITHER direction).
%
% If multiple query points share the same underlying curve, this
% function is O(nQueryPoints x nCurveSegments) -- fine for a scattered
% set of points (hundreds/thousands), not intended for full-image use.

    x1 = perx(:); y1 = pery(:);
    n = numel(x1);
    x2 = [x1(2:end); x1(1)];
    y2 = [y1(2:end); y1(1)];

    qx = qx(:); qy = qy(:); nx = nx(:); ny = ny(:);
    nq = numel(qx);

    nearestX = nan(nq, 1);
    nearestY = nan(nq, 1);
    tParam   = nan(nq, 1);
    bestAbsT = inf(nq, 1);

    for e = 1:n
        ax = x1(e); ay = y1(e);
        ex = x2(e) - ax; ey = y2(e) - ay;   % segment direction (b-a)

        denom = ex.*ny - ey.*nx;             % determinant; 0 = ray parallel to segment
        valid = abs(denom) > 1e-12;

        t = nan(nq, 1); s = nan(nq, 1);
        t(valid) = (ex.*(ay - qy(valid)) - ey.*(ax - qx(valid))) ./ denom(valid);
        s(valid) = (nx(valid).*(ay - qy(valid)) - ny(valid).*(ax - qx(valid))) ./ denom(valid);

        inSeg = valid & (s >= -1e-9) & (s <= 1+1e-9);   % small tolerance at segment ends
        better = inSeg & (abs(t) < bestAbsT);

        nearestX(better) = qx(better) + t(better).*nx(better);
        nearestY(better) = qy(better) + t(better).*ny(better);
        tParam(better) = t(better);
        bestAbsT(better) = abs(t(better));
    end

    unmatched = isnan(nearestX);
    if any(unmatched)
        warning(['%d of %d query points'' lines did not intersect the curve at all -- ' ...
                 'check that (nx,ny) actually points toward the curve, and that the ' ...
                 'curve fully encloses/surrounds the query points (or is otherwise ' ...
                 'positioned to be reachable along that direction).'], nnz(unmatched), nq);
    end
end

%%
function w = sgDerivativeWeights(halfWindow, order, dt)
% SGDERIVATIVEWEIGHTS  Fixed linear-filter weights w (length
% 2*halfWindow+1) such that, for uniformly time-sampled data y in a
% window centered at the point of interest, w'*y is the least-squares
% polynomial (Savitzky-Golay-style) estimate of the first derivative at
% the window's center. (Savitzky & Golay, 1964, Anal. Chem. 36(8).)
    idx = (-halfWindow:halfWindow)';
    t = idx * dt;
    A = t.^(0:order);              % design matrix: columns t^0 ... t^order
    Ct = (A' * A) \ A';            % least-squares "hat" matrix
    w = Ct(2, :)';                 % row 2 = coefficient of the linear (t^1) term = derivative
end

%%
function divCmap = defCmap()
divCmap = [ ...
       94    79   162
   90    85   165
   85    91   167
   81    96   170
   76   102   173
   72   108   176
   67   114   178
   63   119   181
   58   125   184
   54   131   187
   51   137   189
   56   142   186
   61   148   184
   66   154   181
   72   160   179
   77   166   177
   82   172   174
   87   178   172
   93   183   169
   98   189   167
  103   195   165
  110   197   165
  117   200   165
  124   203   165
  131   205   165
  138   208   164
  145   211   164
  152   214   164
  159   216   164
  166   219   164
  173   222   164
  179   224   162
  185   227   161
  191   229   160
  197   231   159
  203   234   158
  209   236   156
  215   239   155
  220   241   154
  226   244   153
  231   245   154
  234   246   158
  236   247   161
  239   248   165
  241   249   169
  244   250   173
  246   251   177
  249   252   181
  251   253   185
  254   254   189
  255   253   188
  255   250   183
  255   247   178
  255   244   173
  255   241   167
  254   238   162
  254   235   157
  254   232   152
  254   228   146
  254   225   141
  254   221   136
  254   216   132
  254   211   128
  254   206   124
  254   201   119
  253   196   115
  253   191   111
  253   186   107
  253   181   103
  253   176    98
  252   169    95
  251   163    92
  251   156    89
  250   150    86
  249   143    83
  248   137    80
  247   130    77
  246   123    74
  245   117    71
  244   110    68
  241   105    68
  238   100    69
  235    96    70
  232    91    72
  229    86    73
  226    81    74
  223    77    75
  220    72    76
  216    67    78
  213    62    79
  208    56    78
  202    50    77
  197    44    75
  191    38    74
  186    32    73
  180    26    71
  175    19    70
  169    13    69
  164     7    67
  158     1    66 ] / 255;
end

function analizePIVresults(nFrames, x, y, nfrms, u_filtered, v_filtered, H, W, WIm, HIm, dt)
%% 1.1. Compute displacements 
chi = cell(nFrames,2);
X = double(x{1,1});
Y = double(y{1,1});
chi{1,1} = X;
chi{1,2} = Y;
spatialSmoothing = 1;
for tt = 1:nfrms
    if tt == 1
        Uq = smoothdata2(u_filtered{tt}, 'gaussian', spatialSmoothing);
        Vq = smoothdata2(v_filtered{tt}, 'gaussian', spatialSmoothing);
    else
        % Uq = interp2(X, Y, u_filtered{tt}, chi{tt,1}, chi{tt,2});
        % Vq = interp2(X, Y, v_filtered{tt}, chi{tt,1}, chi{tt,2});
        % griddedInterpolant needs ndgrid format -> transpose meshgrid data
        Fu = griddedInterpolant(X', Y', smoothdata2(u_filtered{tt}', 'gaussian', spatialSmoothing), 'linear', 'linear');
        Fv = griddedInterpolant(X', Y', smoothdata2(v_filtered{tt}', 'gaussian', spatialSmoothing), 'linear', 'linear');
        
        Uq = Fu(chi{tt,1}', chi{tt,2}')';   % transpose in, transpose out
        Vq = Fv(chi{tt,1}', chi{tt,2}')';
    end
    chi{tt+1,1} = chi{tt,1}+double(Uq);
    chi{tt+1,2} = chi{tt,2}+double(Vq);
end
disp = cell(nFrames,2);

for tt = 1:nFrames
    disp{tt,1} = chi{tt,1}-chi{1,1};
    disp{tt,2} = chi{tt,2}-chi{1,2};
end

% Variables
% disp(X,dt*n) = x(dt*n) - X = x(dt*n) - x(0)
% chi(X,dt*n) = x(X,dt*n)
% v(X,dt*(n+1/2))

%% 1.2. Compute perimeter and normals from masks 
pn = perimeter_normals(nFrames,0);

%% 1.3. Compute radius and growth
rFit = zeros(1,nFrames); % Preallocate radius vector
cx = zeros(1,nFrames); % Preallocate x center vector
cy = zeros(1,nFrames); % Preallocate y center vector

for frame = 1:nFrames
    Mask = imread(sprintf('dc_m_t%03d.tif',frame));
    Mask = Mask > 127;
    Bn = bwboundaries(Mask);
    [~, biggest] = max(cellfun(@(c) size(c,1), Bn));
    boundary = Bn{biggest};
    xy = [boundary(:,2), boundary(:,1)];
    [cx(frame), cy(frame), rFit(frame)] = fitCircleLS(xy);
end

%% 1.4. Interface on the reference configuration
% (computed before F: the deformation gradient below needs p0n to mask the
% organoid hole; this block only depends on pn, chi and disp)

p0n = cell(nFrames,4); % Same structure as pn, it is the "hole" left by
% the organoid in the reference configuration
for frame = 1:nFrames
    % Reshape matrices to vectors and removing NaNs
    szrshp = size(double(chi{frame,1}));
    interpx = reshape(double(chi{frame,1}),[szrshp(1)*szrshp(2) 1]); 
    interpy = reshape(double(chi{frame,2}),[szrshp(1)*szrshp(2) 1]);
    interpdx = reshape(double(disp{frame,1}),[szrshp(1)*szrshp(2) 1]);
    interpdy = reshape(double(disp{frame,2}),[szrshp(1)*szrshp(2) 1]);
    
    % Create interpolations for the displacement in the deformed
    % configuration
    dispxInterp = scatteredInterpolant(interpx,interpy,interpdx,"linear","linear");
    dispyInterp = scatteredInterpolant(interpx,interpy,interpdy,"linear","linear");
    
    
    xb = pn{frame,1};
    yb = pn{frame,2};
    dxb = dispxInterp(xb,yb);
    dyb = dispyInterp(xb,yb);
    per0x = xb-dxb;
    per0y = yb-dyb;
    dx = circshift(per0x,-1)-circshift(per0x,1);
    dy = circshift(per0y,-1)-circshift(per0y,1);
    nrmD = sqrt(dx.*dx+dy.*dy);
    t0x = dx./nrmD;
    t0y = dy./nrmD;
    
    n0x = t0y;
    n0y = -t0x;
    p0n(frame,:)={per0x,per0y,n0x,n0y};
end

%% 1.5. Compute gradient deformation tensor
% Deformation gradient by pointwise local least squares (defgrad_pls,
% 2-D Savitzky-Golay differentiator), same settings as the "Test Smoothing"
% section of the main script: (2*mPLS+1)x(2*mPLS+1) window, plane fit, and
% the nodes inside (or on) the reference interface p0n{frame} set to NaN so
% that no derivative stencil crosses the organoid.
mPLS     = 1;        % half-width of the fitting window (grid points)
maskHole = true;     % true = exclude nodes inside p0n{frame} (hole -> NaN)
F = cell(nFrames,4);
C = cell(nFrames,3);
B = cell(nFrames,3);
xPLS = X(1,:);       % reference grid vectors (uniform PIV grid)
yPLS = Y(:,1);
for frame = 1:nFrames
    U = disp{frame,1};
    V = disp{frame,2};
    if maskHole
        [in, on] = inpolygon(X, Y, p0n{frame,1}, p0n{frame,2});
        hole = in | on;  % also remove nodes lying exactly on the boundary
        U(hole) = NaN;
        V(hole) = NaN;
    end
    G = defgrad_pls(xPLS, yPLS, U, V, mPLS);
    FxX = G.F11;  FxY = G.F12;      % F_iJ = delta_iJ + du_i/dX_J
    FyX = G.F21;  FyY = G.F22;
    F{frame,1} = FxX;
    F{frame,2} = FxY;
    F{frame,3} = FyX;
    F{frame,4} = FyY;
    
    % Compute C = F^T*F
    C{frame,1} = FxX.*FxX + FyX.*FyX;
    C{frame,2} = FxX.*FxY + FyX.*FyY;
    C{frame,3} = FxY.*FxY + FyY.*FyY;
    
    B{frame,1} = FxX.*FxX + FxY.*FxY;
    B{frame,2} = FxX.*FyX + FxY.*FyY;
    B{frame,3} = FyX.*FyX + FyY.*FyY;
end

pp = cell(nFrames,6);

for frame = 1:nFrames
    for ii = 1:H
        for jj = 1:W
            Ctemp = [C{frame,1}(ii,jj) C{frame,2}(ii,jj); C{frame,2}(ii,jj) C{frame,3}(ii,jj)];
            if any(~isfinite(Ctemp(:)))      % masked node (hole / no data)
                for kk = 1:6, pp{frame,kk}(ii,jj) = NaN; end
                continue
            end
            [V,D] = eig(Ctemp);
            if D(1,1) > D(2,2)
                pp{frame,1}(ii,jj) = D(1,1);
                pp{frame,2}(ii,jj) = D(2,2);
                pp{frame,3}(ii,jj) = V(1,1);
                pp{frame,4}(ii,jj) = V(2,1);
                pp{frame,5}(ii,jj) = V(1,2);
                pp{frame,6}(ii,jj) = V(2,2);
            else
                pp{frame,2}(ii,jj) = D(1,1);
                pp{frame,1}(ii,jj) = D(2,2);
                pp{frame,5}(ii,jj) = V(1,1);
                pp{frame,6}(ii,jj) = V(2,1);
                pp{frame,3}(ii,jj) = V(1,2);
                pp{frame,4}(ii,jj) = V(2,2);
            end
        end
    end
end

%% 1.6. Normal vector pointing away from the organoid center 
nX = cell(nFrames,1);
nY = cell(nFrames,1);
for t = 1:nFrames
    nX{t} = X-cx(t);
    nY{t} = Y-cy(t);
    mN = sqrt(nX{t}.*nX{t}+nY{t}.*nY{t});
    nX{t} = nX{t}./mN;
    nY{t} = nY{t}./mN;
end

%% 1.7. Normal interface velocity
% Compute the normal velocity of the organoid interface across a
% time-evolving mask stack, with temporal smoothing.
%
% METHOD:
%   1. phi(x,y,t) = signed distance transform of each mask (negative
%      inside, positive outside). Since phi satisfies the eikonal
%      equation |grad phi| = 1, the level-set relation
%           V_n = -(1/|grad phi|) * d(phi)/dt
%      collapses to simply V_n = -d(phi)/dt (Machacek & Danuser 2006,
%      Biophys J 90(4):1439-1452, use the same device).
%   2. d(phi)/dt is estimated with a Savitzky-Golay-style local
%      polynomial fit over a small TEMPORAL window of frames (Savitzky &
%      Golay, 1964, Anal. Chem. 36(8):1627-1639) -- this is what
%      suppresses frame-to-frame noise, verified >10x noise reduction
%      vs naive differencing on a synthetic test (see chat).

% -------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------

sgFrameLen  = 7;      % temporal smoothing window (frames); must be odd,
% >= 2*sgOrder+1.
sgOrder     = 2;      % polynomial order for the temporal fit (2)
vDownSmplng = 5;      % Downsampling of the velocity field interpolator

% -------------------------------------------------------------------

[XGDs, YGDs] = meshgrid(1:vDownSmplng:WIm, 1:vDownSmplng:HIm);
half = floor(sgFrameLen/2);
if mod(sgFrameLen,2) == 0 || sgFrameLen < 2*sgOrder+1
    error('sgFrameLen must be odd and >= 2*sgOrder+1.');
end
[HDs, WDs] = size(XGDs);

% ---- 1. Signed distance transform per frame -------------------------
phi = zeros(HDs, WDs, nFrames);
for t = 1:nFrames
    Mask = imread(sprintf('dc_m_t%03d.tif',t));      
    m = Mask > 167;
    phi(:,:,t) = phiFromBoundary(pn{t,1}, pn{t,2}, XGDs, YGDs);
end

% ---- 2. Savitzky-Golay-style temporal derivative --------------------
% Fixed filter weights for a least-squares polynomial derivative
% estimate at a window's center, for uniformly time-sampled data
% (Savitzky & Golay 1964). Computed directly via least squares -- no
% Signal Processing Toolbox needed.
w = sgDerivativeWeights(half, sgOrder, dt);   % length 2*half+1

dphidt = nan(HDs, WDs, nFrames);
for t = (half+1):(nFrames-half)
    acc = zeros(HDs, WDs);
    for k = -half:half
        acc = acc + w(k+half+1) * phi(:,:,t+k);
    end
    dphidt(:,:,t) = acc;
end

vn_field = -dphidt; % Section: V_n = -d(phi)/dt, valid where |grad phi| ~ 1


%% 1.8. Normal interface velocity at reference configuration
% ---- 1. Signed distance transform per frame -------------------------
phi0 = zeros(HDs, WDs, nFrames);
for t = 1:nFrames
    phi0(:,:,t) = phiFromBoundary(p0n{t,1}, p0n{t,2}, XGDs, YGDs);
end

% ---- 2. Savitzky-Golay-style temporal derivative --------------------
dphi0dt = nan(HDs, WDs, nFrames);
for t = (half+1):(nFrames-half)
    acc = zeros(HDs, WDs);
    for k = -half:half
        acc = acc + w(k+half+1) * phi0(:,:,t+k);
    end
    dphi0dt(:,:,t) = acc;
end

V0n_field = -dphi0dt; % Section: V_n = -d(phi)/dt, valid where |grad phi| ~ 1

%% 1.9. Windowing PIV velocity
u_filtered_cat = cat(3, u_filtered{:});   % Ny x Nx x Nt
Us = smoothdata(u_filtered_cat, 3, 'gaussian', 5);  % W = window length, along time dim
v_filtered_cat = cat(3, v_filtered{:});   % Ny x Nx x Nt
Vs = smoothdata(v_filtered_cat, 3, 'gaussian', 5);  % W = window length, along time dim

%% 1.10. Interface velocities
vn = cell(nFrames,1);
Vd = cell(nFrames,1);
vd = cell(nFrames,1);
vCn = cell(nFrames,1);

for t = (half+1):(nFrames-half)
    vn{t} = interp2(XGDs,YGDs,vn_field(:,:,t),pn{t,1},pn{t,2}); % Normal interface vel
    Vd{t} = interp2(XGDs,YGDs,V0n_field(:,:,t),p0n{t,1},p0n{t,2}); % Reference normal interface vel
    vCx = interp2(X,Y,Us(:,:,t),pn{t,1},pn{t,2}); % u interpolated at interface
    vCy = interp2(X,Y,Vs(:,:,t),pn{t,1},pn{t,2}); % v interpolated at interface
    vCn{t} = vCx.*pn{t,3}+vCy.*pn{t,4}; % Normal interface collagen velocity
    vd{t} = vn{t}-vCn{t};
end
save("PIV_results","vn","Vd","vd","vCn","pn","pp","p0n","chi","disp","F","C","nX","nY","X","Y")
end

%%
function [W,S] = NH(C,params)
    mu = params.mu;
    lm = params.lm;
    iC = inv(C);
    J = sqrt(det(C));
    W = mu/2*(trace(C)-2)-mu*log(J)+lm/2*(log(J))^2;
    S = mu*(eye(2)-iC)+lm*log(J)*iC;
end

function [W,S] = SS(C,params)
    m = params.m;
    ms = params.ms;
    W = m/ms*(exp(ms*(trace(C)-2))-ms*(trace(C)-2)-1);
    S = 2*m*(exp(ms*(trace(C)-2))-1)*eye(2);
end

function [W,S] = SS2(C,params)
    m = params.m;
    ms = params.ms;
    [l1, l2, N1, N2] = principalStretches2D(C(1,1), C(2,2), C(1,2));
    l1 = sqrt(l1);
    l2 = sqrt(l2);
    if (l1 > 1)
        f1 = m/(2*ms*(2*ms+1))*(exp(ms*(l1^2-1))-2*ms*(l1-1)-1);
        df1dl1 = m/(2*ms+1)*(l1*exp(ms*(l1^2-1))-1);
    else
        f1 = m/2*(l1-1)^2;
        df1dl1 = m*(l1-1);
    end
    if (l2 > 1)
        f2 = m/(2*ms*(2*ms+1))*(exp(ms*(l2^2-1))-2*ms*(l2-1)-1);
        df2dl2 = m/(2*ms+1)*(l2*exp(ms*(l2^2-1))-1);
    else
        f2 = m/2*(l2-1)^2;
        df2dl2 = m*(l2-1);
    end
    dl1dC = N1*N1'/(2*l1);
    dl2dC = N2*N2'/(2*l2);
    W = f1 + f2;
    S = 2*(df1dl1*dl1dC + df2dl2*dl2dC);
end

function [W,S] = MC(C,params)
    m = params.m;
    mg = params.mg;
    alp = params.alp;
    [l1, l2, N1, N2] = principalStretches2D(C(1,1), C(2,2), C(1,2));
    t = l1^2*l2^(2+alp)+l2^2*l1^(2+alp)-2;
    if (t < 0)
        W = 0;
        S = zeros(2);
    else
        W = m/(mg^2)*(exp(mg*t)-1-mg*t-mg^2*t^2/2);
        dWdt = m/mg*(exp(mg*t)-1-mg*t);
        dtdC1 = (l2^(2+alp)+(1+alp/2)*l2^2*l1^alp).*(N1*N1');
        dtdC2 = (l1^(2+alp)+(1+alp/2)*l1^2*l2^alp).*(N2*N2');
        S = 2*dWdt*(dtdC1+dtdC2);
    end
end

function [lam1, lam2, N1, N2] = principalStretches2D(C11, C22, C12)
% Principal stretches (lam1 >= lam2) of symmetric 2x2 C = [C11 C12; C12 C22].
% Inputs can be arrays of equal size (vectorized over all tensors).
    m    = 0.5*(C11 + C22);
    r    = hypot(0.5*(C11 - C22), C12);   % sqrt(a^2+b^2) without overflow
    detC = C11.*C22 - C12.^2;

    c1 = m + r;           % larger eigenvalue of C
    c2 = detC ./ c1;      % smaller one, cancellation-free

    lam1 = sqrt(c1);
    lam2 = sqrt(c2);

    if nargout > 2
        th = 0.5*atan2(2*C12, C11 - C22);
        N1 = [cos(th(:)), sin(th(:))].';   % 2 x n, columns = N1 per tensor
        N2 = [-sin(th(:)), cos(th(:))].';
    end
end

%%
function out = defgrad_pls(x, y, U, V, m, varargin)
%DEFGRAD_PLS  Deformation gradient from a gridded 2-D displacement field
%   using pointwise local least-squares fitting (2-D Savitzky-Golay
%   differentiator / "PLS" strain estimation).
%
%   out = defgrad_pls(x, y, U, V, m)
%   out = defgrad_pls(x, y, U, V, m, 'Order', 2, 'MinPoints', 9)
%
%   INPUTS
%     x, y   Grid coordinate vectors in the REFERENCE configuration.
%            numel(x) = size(U,2) (columns), numel(y) = size(U,1) (rows).
%            Spacing must be uniform. A negative spacing (e.g. image rows
%            that increase downward while you want Y up) is handled
%            correctly: derivatives are taken w.r.t. the actual x, y.
%     U, V   Displacement components u(X,Y), v(X,Y), size [ny nx].
%            NaN marks missing data (holes, cracks, outside the specimen).
%     m      Half-width of the fitting window: (2m+1) x (2m+1) points.
%            Larger m -> less noise, more bias (virtual strain gauge size).
%
%   OPTIONS (name/value)
%     'Order'      1 (plane, default) or 2 (full quadratic). On complete,
%                  symmetric windows both give IDENTICAL gradients; order 2
%                  only changes truncated windows (edges/holes), where it
%                  reduces bias but increases noise.
%     'MinPoints'  Minimum number of valid points required for a fit
%                  (default: 2*number of unknowns). Fewer -> NaN output.
%
%   OUTPUT  struct with fields (all [ny nx]):
%     dudX dudY dvdX dvdY   displacement gradient components
%     F11 F12 F21 F22       F = I + grad_X(u),  F_iJ = delta_iJ + du_i/dX_J
%     J                     det(F)
%     C11 C12 C22           right Cauchy-Green tensor C = F'F (E = (C-I)/2)
%     nPts                  number of points used in each fit
%
%   NOTE: assumes a Lagrangian description (u given at reference points),
%   which is the case for subset-based DIC. If u is Eulerian (given on the
%   current configuration, e.g. some PIV/optical-flow outputs), then
%   grad_x(u) = I - inv(F), not F - I.
 
% ---------------- parse options ----------------
p = inputParser;
p.addParameter('Order', 1, @(s) any(s == [1 2]));
p.addParameter('MinPoints', [], @(s) isempty(s) || (isscalar(s) && s > 0));
p.parse(varargin{:});
order = p.Results.Order;
nUnk  = 3 + 3*(order == 2);              % unknowns per component
minPts = p.Results.MinPoints;
if isempty(minPts), minPts = 2*nUnk; end
 
x = x(:).'; y = y(:);
[ny, nx] = size(U);
assert(isequal(size(V), [ny nx]), 'U and V must have the same size');
assert(numel(x) == nx && numel(y) == ny, 'x must match columns, y rows');
hx = mean(diff(x));  hy = mean(diff(y));
assert(all(abs(diff(x) - hx) < 1e-8*abs(hx)) && ...
       all(abs(diff(y) - hy) < 1e-8*abs(hy)), 'grid must be uniform');
 
valid = isfinite(U) & isfinite(V);
 
% ---------------- 1) fast path: complete interior windows ----------------
% For a full symmetric window the LS slope decouples (Sx = Sy = Sxy = 0):
%   a1 = sum(x_k u_k)/sum(x_k^2),   a2 = sum(y_k u_k)/sum(y_k^2)
% which is a fixed linear stencil -> a single correlation (filter2).
[I, Jc] = ndgrid(-m:m, -m:m);            % I = row offset, Jc = column offset
dX = Jc*hx;  dY = I*hy;                  % local coordinates of window points
Kx = dX / sum(dX(:).^2);
Ky = dY / sum(dY(:).^2);
 
U0 = U; U0(~valid) = 0;
V0 = V; V0(~valid) = 0;
dudX = filter2(Kx, U0, 'same');  dudY = filter2(Ky, U0, 'same');
dvdX = filter2(Kx, V0, 'same');  dvdY = filter2(Ky, V0, 'same');
 
nPts = round(filter2(ones(2*m+1), double(valid), 'same'));
full = (nPts == (2*m+1)^2);              % false near borders and holes
 
% ---------------- 2) slow path: truncated windows (edges, holes) ----------
% Weighted-free LS fit over whichever points exist in the window.
todo = find(valid & ~full);
for k = todo.'
    [r, c] = ind2sub([ny nx], k);
    rr = max(1, r-m):min(ny, r+m);
    cc = max(1, c-m):min(nx, c+m);
    [R, C] = ndgrid(rr, cc);
    ok = valid(rr, cc);
    xl = (C(ok) - c)*hx;                 % local coords centred on node
    yl = (R(ok) - r)*hy;
    if order == 1
        A = [ones(size(xl)) xl yl];
    else
        A = [ones(size(xl)) xl yl xl.^2 xl.*yl yl.^2];
    end
    if numel(xl) < minPts || rank(A) < nUnk
        dudX(k) = NaN; dudY(k) = NaN; dvdX(k) = NaN; dvdY(k) = NaN;
        continue
    end
    uu = U(rr, cc); vv = V(rr, cc);
    coef = A \ [uu(ok) vv(ok)];          % QR least squares, both components
    dudX(k) = coef(2,1);  dudY(k) = coef(3,1);
    dvdX(k) = coef(2,2);  dvdY(k) = coef(3,2);
end
 
% nodes with no data of their own -> NaN
dudX(~valid) = NaN; dudY(~valid) = NaN; dvdX(~valid) = NaN; dvdY(~valid) = NaN;
 
% ---------------- 3) kinematics (computed AFTER smoothing) ----------------
F11 = 1 + dudX;  F12 = dudY;
F21 = dvdX;      F22 = 1 + dvdY;
C11 = F11.^2 + F21.^2;
C22 = F12.^2 + F22.^2;
C12 = F11.*F12 + F21.*F22;
 
out = struct('dudX', dudX, 'dudY', dudY, 'dvdX', dvdX, 'dvdY', dvdY, ...
             'F11', F11, 'F12', F12, 'F21', F21, 'F22', F22, ...
             'J', F11.*F22 - F12.*F21, ...
             'C11', C11, 'C12', C12, 'C22', C22, ...
             'nPts', nPts);
end
%%
function [thq, Fq, info] = radialFieldPlot(x, y, F, xc, yc, varargin)
%RADIALFIELDPLOT  Concentric polar plots f_k(theta) of one or more fields
%   sampled on closed curves, all with angles measured from the same
%   user-defined centre (xc, yc). Field 1 is drawn innermost, field K
%   outermost, each in its own annular band with its own radial scale.
%
%   ONE FIELD
%     radialFieldPlot(x, y, f, xc, yc, ...)
%   SEVERAL FIELDS ON THE SAME CURVE (columns of F, N-by-K)
%     radialFieldPlot(x, y, [f1 f2], xc, yc, 'ZeroRadius', [0.5 0.3], ...)
%   SEVERAL FIELDS ON DIFFERENT CURVES (cell arrays)
%     radialFieldPlot({x1,x2}, {y1,y2}, {f1,f2}, xc, yc, ...)
%
%   Inputs
%     x, y, F : curve points in curve order and the field values there.
%               The closing point may or may not be repeated at the end.
%     xc, yc  : common centre from which angles are measured.
%   Options (per-field options take one value per field, inner -> outer)
%     'ZeroRadius': fraction a_k in (0,1) OF EACH BAND where f_k = 0 is
%                   drawn (dashed circle). NaN = automatic limits. Default
%                   NaN for all.
%     'Fill'      : true -> shade between f_k = 0 and f_k(theta), one colour
%                   for f > 0 and one for f < 0 (default false).
%     'Bands'     : K-by-2 [s0 s1] radial extents of each band, as fractions
%                   of the outer radius (0 = centre, 1 = edge). Default:
%                   equal bands separated by 'BandGap'.
%     'BandGap'   : gap between default bands (default 0.06).
%     'Names'     : cellstr of field names for the legend.
%     'StarTol'   : degrees of backtracking tolerated (default 0). If the
%                   curve folds back by at most this angle (see
%                   info.backtrackDeg), it is treated as nearly star-shaped
%                   and the folds are resolved with 'Fold'.
%     'FillMode'  : 'sign' (default) - two colours by sign of f;
%                   'value' - colour varies continuously with f.
%     'ColorBy'   : with 'value': 'angle' (default) - each radial strip has
%                   the colour of f(theta); 'radius' - colour goes from
%                   f = 0 on the zero circle to f(theta) on the curve.
%     'Colormap'  : cell of colormaps (M-by-3), one per field. Default:
%                   blue-white-red centred on 0 if f changes sign, else
%                   parula.
%     'CLim'      : K-by-2 colour limits per field (NaN rows = auto).
%     'RLim'      : K-by-2 radial limits [inner outer] per field, in units of
%                   f (NaN = automatic). Fixes the radius at which each f
%                   value is drawn, so several figures can be compared.
%                   Both given: the band runs exactly from RLim(1) to
%                   RLim(2) (ZeroRadius is then ignored). One given plus
%                   ZeroRadius: the other end follows from ZeroRadius.
%                   Values outside the limits are clipped to the band edge.
%     'Colorbar'  : true (default) - one colour bar per field.
%     'Fold'      : 'average' (default) - mean of all samples in each
%                   angular bin of the output grid; 'remove' - drop the
%                   samples that go backwards in angle.
%     'RScale'    : radial (and colour) scale, one per field or one for all:
%                   'linear' (default); 'log' - log10(f), f > 0 only;
%                   'symlog' - sign(f)*log10(1 + |f|/C), for mixed-sign
%                   data: linear for |f| << C, logarithmic for |f| >> C,
%                   and f = 0 stays at 0 (so ZeroRadius still applies).
%     'LinThresh' : C for 'symlog', in units of f, one per field
%                   (NaN = auto: 1e-2 * max|f|, i.e. two decades of log).
%     'LineColor' : colour of the f(theta) curve: any MATLAB colour (default
%                   'k'), or 'value' - coloured continuously with the same
%                   colormap, limits and scale as the 'value' fill.
%     'LineWidth' : curve width (default 1.5).
%     'LineOutline': with 'value', draw a thin dark line underneath so
%                   near-white parts of the colormap stay visible (default
%                   true).
%     'PosColor', 'NegColor': K-by-3 fill colours (rows reused if fewer).
%     'FillAlpha' : fill transparency (default 0.35).
%     'NTheta'    : points on the uniform angle grid (default 360).
%     'Method'    : interp1 method: 'linear' (default), 'pchip', ...
%     'ShowRaw'   : overlay original samples as dots (default false).
%     'Plot'      : true (default) / false.
%   Outputs
%     thq  : uniform angle grid in [0, 2*pi), radians (NTheta-by-1).
%     Fq   : NTheta-by-K resampled fields (NaN column if field k is not
%            single-valued in theta).
%     info : 1-by-K struct array (winding number, star-shape flag, ...).
%
%   Angle convention: theta = atan2(y - yc, x - xc), counter-clockwise
%   from +x. The chart is drawn on Cartesian axes (works in any release).
%
%   Example
%     t  = linspace(0, 2*pi, 200); t(end) = [];
%     x  = 2*cos(t) + 0.3*cos(3*t);  y = 1.5*sin(t);
%     f1 = sin(2*t) + 0.3;   f2 = 5*cos(3*t) - 1;
%     radialFieldPlot(x, y, [f1(:) f2(:)], 0.2, -0.1, ...
%         'ZeroRadius', [0.5 0.4], 'Fill', true, 'Names', {'p','\omega'});

% ---------- options ------------------------------------------------------
p = inputParser;
addParameter(p, 'NTheta', 360);
addParameter(p, 'Method', 'linear');
addParameter(p, 'Plot', true);
addParameter(p, 'ShowRaw', false);
addParameter(p, 'Fill', false);
addParameter(p, 'FillAlpha', 0.35);
addParameter(p, 'PosColor', [0.85 0.33 0.10; 0.75 0.15 0.45; 0.60 0.50 0.05]);
addParameter(p, 'NegColor', [0.00 0.45 0.74; 0.10 0.60 0.50; 0.45 0.30 0.70]);
addParameter(p, 'ZeroRadius', []);
addParameter(p, 'Bands', []);
addParameter(p, 'BandGap', 0.06);
addParameter(p, 'Names', {});
addParameter(p, 'StarTol', 0);          % degrees of backtracking tolerated
addParameter(p, 'Fold', 'average');     % 'average' | 'remove'
addParameter(p, 'FillMode', 'sign');    % 'sign' | 'value'
addParameter(p, 'ColorBy', 'angle');    % 'angle' | 'radius' (FillMode 'value')
addParameter(p, 'Colormap', {});        % cell, one colormap per field
addParameter(p, 'CLim', []);            % K-by-2 colour limits, NaN = auto
addParameter(p, 'RLim', []);            % K-by-2 radial limits (in f), NaN = auto
addParameter(p, 'Colorbar', true);
addParameter(p, 'RScale', 'linear');    % 'linear' | 'log' | 'symlog'
addParameter(p, 'LinThresh', []);       % symlog threshold C, NaN = auto
addParameter(p, 'LineColor', 'k');      % any colour spec, or 'value'
addParameter(p, 'LineWidth', 1.5);
addParameter(p, 'LineOutline', true);   % thin dark outline under 'value' line
parse(p, varargin{:});
opt = p.Results;

% ---------- normalise inputs to cell arrays, one cell per field ----------
if iscell(F)
    K  = numel(F);
    Fc = F;
else
    if isvector(F), F = F(:); end
    K  = size(F, 2);
    Fc = num2cell(F, 1);
end
if ~iscell(x), x = repmat({x}, 1, K); end
if ~iscell(y), y = repmat({y}, 1, K); end
if numel(x) ~= K || numel(y) ~= K
    error('Give one x and one y per field (or a single shared x, y).');
end

a = opt.ZeroRadius;
if isempty(a), a = NaN(1, K); end
if isscalar(a), a = repmat(a, 1, K); end
if numel(a) ~= K, error('ZeroRadius needs one value per field (NaN = auto).'); end
if any(~isnan(a) & (a <= 0 | a >= 1)), error('ZeroRadius values must lie in (0,1).'); end

rs = opt.RScale;
if ischar(rs) || isstring(rs), rs = cellstr(rs); end
if isscalar(rs), rs = repmat(rs, 1, K); end
if numel(rs) ~= K, error('RScale needs one value per field (or one for all).'); end
if ~all(ismember(lower(rs), {'linear', 'log', 'symlog'}))
    error('RScale must be ''linear'', ''log'' or ''symlog''.');
end
Cs = opt.LinThresh;
if isempty(Cs), Cs = NaN(1, K); end
if isscalar(Cs), Cs = repmat(Cs, 1, K); end
if numel(Cs) ~= K, error('LinThresh needs one value per field (NaN = auto).'); end

names = opt.Names;
if isempty(names), names = arrayfun(@(k) sprintf('f_%d', k), 1:K, 'UniformOutput', false); end

bands = opt.Bands;
if isempty(bands)
    g  = opt.BandGap;
    w  = (1 - (K-1)*g) / K;           % K bands + (K-1) gaps fill [0, 1]
    s0 = (0:K-1).' * (w + g);
    bands = [s0, s0 + w];
end
if ~isequal(size(bands), [K 2]), error('Bands must be K-by-2.'); end

% ---------- per-field angular processing ---------------------------------
thq = linspace(0, 2*pi, opt.NTheta + 1).';  thq(end) = [];
Fq  = NaN(opt.NTheta, K);
for k = 1:K
    [Fq(:, k), info(k)] = angularResample(x{k}, y{k}, Fc{k}, xc, yc, thq, opt, k); %#ok<AGROW>
end

if ~opt.Plot, return; end

% ---------- plot -----------------------------------------------------------
figure;
ax = axes; hold(ax, 'on'); axis(ax, 'equal', 'off');
tc = linspace(0, 2*pi, 361).';
grey = [0.85 0.85 0.85];

% spokes and angle labels (shared by all bands: common centre and angle)
for ang = 0:30:330
    plot(ax, [0 cosd(ang)], [0 sind(ang)], '-', 'Color', grey, 'HandleVisibility', 'off');
    text(ax, 1.08*cosd(ang), 1.08*sind(ang), sprintf('%d°', ang), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end

cbars = {};
for k = 1:K
    s0 = bands(k, 1);  s1 = bands(k, 2);
    star = info(k).resampled;      % f(theta) available on the grid
    if star, allf = [Fq(:, k); info(k).fRaw]; else, allf = info(k).fRaw; end

    % --- radial scale: g = T(f); everything below is drawn in g ----------
    [T, Ti, Ck] = scaleTransform(rs{k}, Cs(k), allf, k);
    isLog = strcmpi(rs{k}, 'log');
    gall = T(allf);
    gmin = min(gall);  gmax = max(gall);

    % --- radial limits of this band (in g; T(0) = 0 for linear/symlog) ---
    if ~isnan(a(k)) && ~isLog
        % g = 0 at fraction a of the band: rmin = -a/(1-a)*rmax,
        % with rmax >= gmax and rmax >= -gmin*(1-a)/a
        rmax = max([gmax, -gmin*(1 - a(k))/a(k), eps]);
        rmin = -a(k)/(1 - a(k)) * rmax;
    else
        if ~isnan(a(k)) && isLog
            warning('Field %d: ZeroRadius ignored - f = 0 does not exist on a log scale.', k);
        end
        rmin = gmin;  rmax = gmax;
        if opt.Fill && ~isLog, rmin = min(rmin, 0); rmax = max(rmax, 0); end
        if rmin == rmax, rmax = rmin + 1; end
    end
    % --- fixed radial limits (RLim, in f) override the automatic ones -----
    rl = [NaN NaN];
    if ~isempty(opt.RLim), rl = opt.RLim(min(k, size(opt.RLim, 1)), :); end
    if any(~isnan(rl))
        rlg  = T(rl);                          % NaN stays NaN
        useA = ~isnan(a(k)) && ~isLog;
        if ~isreal(rlg) || any(isinf(rlg))
            error('Field %d: RLim [%g %g] is not valid on a ''%s'' scale (log needs RLim > 0).', ...
                  k, rl(1), rl(2), rs{k});
        end
        if all(~isnan(rl))
            if useA
                warning('Field %d: both RLim limits given, so ZeroRadius is ignored.', k);
            end
            rmin = rlg(1);  rmax = rlg(2);
        elseif ~isnan(rl(2))                   % outer limit fixed
            rmax = rlg(2);
            if useA, rmin = -a(k)/(1 - a(k)) * rmax; end
        else                                   % inner limit fixed
            rmin = rlg(1);
            if useA, rmax = -(1 - a(k))/a(k) * rmin; end
        end
        if ~(rmin < rmax)
            error('Field %d: radial limits [%g %g] (in f) must be increasing.', ...
                  k, Ti(rmin), Ti(rmax));
        end
        tol  = 1e-12 * (rmax - rmin);
        nOut = nnz(gall < rmin - tol | gall > rmax + tol);
        if nOut > 0
            warning('Field %d: %d values outside the radial limits [%g %g] are clipped to the band edges.', ...
                    k, nOut, Ti(rmin), Ti(rmax));
        end
    end

    gBase = 0;                     % fill baseline in g
    if isLog, gBase = rmin; end    % log: fill from the inner edge of the band
    info(k).rLim  = [rmin rmax];   % in g = T(f)
    info(k).rLimF = Ti([rmin rmax]);
    info(k).band  = [s0 s1];
    info(k).rScale = rs{k};  info(k).linThresh = Ck;

    % value g -> plot radius s in [s0, s1], then polar -> Cartesian
    % (clamped to the band so clipped values never spill into other rings)
    sOf = @(r) s0 + (s1 - s0) .* min(max((r - rmin) ./ (rmax - rmin), 0), 1);
    P   = @(t, r) deal(sOf(r).*cos(t), sOf(r).*sin(t));

    % --- band outline and radial tick labels ------------------------------
    for sb = unique([s0 s1])
        if sb > 0
            plot(ax, sb*cos(tc), sb*sin(tc), '-', 'Color', grey, 'HandleVisibility', 'off');
        end
    end
    [tf, tl] = scaleTicks(rs{k}, info(k).rLimF(1), info(k).rLimF(2), Ck);
    for i = 1:numel(tf)
        sr = sOf(T(tf(i)));
        text(ax, sr*cosd(75), sr*sind(75), tl{i}, ...
             'FontSize', 7, 'Color', [0.35 0.35 0.35]);
    end

    pc = opt.PosColor(min(k, size(opt.PosColor, 1)), :);
    nc = opt.NegColor(min(k, size(opt.NegColor, 1)), :);

    % --- colormap for this field (shared by fill and line) -----------------
    lineByValue = strcmpi(opt.LineColor, 'value');
    fillByValue = opt.Fill && star && strcmpi(opt.FillMode, 'value');
    if lineByValue || fillByValue
        [cmap, clim_k] = fieldColormap(opt, k, allf);   % limits in f
        clim_g = T(clim_k);                               % ... and in g
        if any(~isfinite(clim_g)) || ~isreal(clim_g)
            error('Field %d: CLim [%g %g] is not valid on a ''%s'' scale (log needs CLim > 0).', ...
                  k, clim_k(1), clim_k(2), rs{k});
        end
        info(k).colormap = cmap;  info(k).cLim = clim_k;
        cbars{end+1} = {cmap, clim_g, names{k}, rs{k}, T, clim_k, Ck}; %#ok<AGROW>
    end

    % --- fill -------------------------------------------------------------
    showZero = (~isnan(a(k)) || (opt.Fill && star)) && ~isLog ...
               && rmin <= 0 && rmax >= 0;
    if opt.Fill && ~star
        warning(['Field %d: Fill skipped - f(theta) is not single-valued ', ...
                 '(centre outside the curve or curve not star-shaped).'], k);
    end
    if opt.Fill && star
        % angle grid with exact zero crossings inserted
        % (roots found in f, where the resampled f is piecewise linear;
        %  T is monotone with T(0) = 0, so signs are unchanged)
        [tt, vv] = insertZeroCrossings(thq, Fq(:, k));
        m  = numel(tt);
        gv = T(vv);
        % vertices: outer (on the curve) then inner (on the baseline)
        [ox, oy] = P(tt, gv);
        [ix, iy] = P(tt, gBase * ones(m, 1));
        V = [ox oy; ix iy];
        % one quadrilateral per angular interval: each is a simple polygon
        % because f has a single sign between consecutive (crossing) points
        j   = (1:m-1).';
        Fc4 = [j, j+1, m+j+1, m+j];
        sgn = sign(vv(j) + vv(j+1));          % sign of f on the interval
        keep = sgn ~= 0;

        if strcmpi(opt.FillMode, 'value')
            % ---- continuous colour: per-vertex RGB, interpolated -------
            if strcmpi(opt.ColorBy, 'radius')
                cv = [gv; gBase * ones(m, 1)];  % baseline -> curve
            else
                cv = [gv; gv];                  % same value along each radius
            end
            Crgb = valueToRGB(cv, cmap, clim_g);   % colour also on scale T
            patch(ax, 'Faces', Fc4(keep, :), 'Vertices', V, ...
                  'FaceVertexCData', Crgb, 'FaceColor', 'interp', ...
                  'EdgeColor', 'none', 'FaceAlpha', opt.FillAlpha, ...
                  'HandleVisibility', 'off');
        else
            % ---- two colours by sign --------------------------------------
            C = repmat(pc, m-1, 1);
            C(sgn < 0, :) = repmat(nc, nnz(sgn < 0), 1);
            patch(ax, 'Faces', Fc4(keep, :), 'Vertices', V, ...
                  'FaceVertexCData', C(keep, :), 'FaceColor', 'flat', ...
                  'EdgeColor', 'none', 'FaceAlpha', opt.FillAlpha, ...
                  'HandleVisibility', 'off');
            % legend proxies
            patch(ax, NaN, NaN, pc, 'FaceAlpha', opt.FillAlpha, 'EdgeColor', 'none', ...
                  'DisplayName', [names{k} ' > 0']);
            patch(ax, NaN, NaN, nc, 'FaceAlpha', opt.FillAlpha, 'EdgeColor', 'none', ...
                  'DisplayName', [names{k} ' < 0']);
        end
    end

    % --- curve, samples, zero circle --------------------------------------
    if star
        cth = [thq; thq(1)];  cf = [Fq(:, k); Fq(1, k)];
    else
        cth = [info(k).thetaRaw; info(k).thetaRaw(1)];  % raw, curve order
        cf  = [info(k).fRaw;     info(k).fRaw(1)];
    end
    [cx, cy] = P(cth, T(cf));
    if lineByValue
        if opt.LineOutline
            plot(ax, cx, cy, '-', 'Color', [0.25 0.25 0.25], ...
                 'LineWidth', opt.LineWidth + 1, 'HandleVisibility', 'off');
        end
        % a patch with no face and interpolated edge colours draws a line
        % whose colour varies along it; the curve is closed, so the
        % patch's closing edge is just the curve's own last segment
        patch(ax, 'XData', cx, 'YData', cy, 'FaceColor', 'none', ...
              'EdgeColor', 'interp', ...
              'FaceVertexCData', valueToRGB(T(cf), cmap, clim_g), ...
              'LineWidth', opt.LineWidth, 'DisplayName', names{k});
    else
        plot(ax, cx, cy, '-', 'Color', opt.LineColor, ...
             'LineWidth', opt.LineWidth, 'DisplayName', names{k});
    end
    if opt.ShowRaw
        [sx, sy] = P(info(k).thetaRaw, T(info(k).fRaw));
        plot(ax, sx, sy, '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 6, ...
             'HandleVisibility', 'off');
    end
    if showZero
        [zx, zy] = P(tc, zeros(size(tc)));
        plot(ax, zx, zy, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

title(ax, sprintf('About (%.3g, %.3g); inner \\rightarrow outer: %s', ...
                  xc, yc, strjoin(names, ', ')));
legend(ax, 'Location', 'southoutside');
xlim(ax, [-1.15 1.15]); ylim(ax, [-1.15 1.15]);

% one colour bar per field: each on an invisible axes carrying that field's
% colormap and limits (the fill itself is true colour, so each ring keeps
% its own colormap)
if opt.Colorbar && ~isempty(cbars)
    nb = numel(cbars);
    ax.Position = [0.05 0.12 0.70 - 0.08*(nb-1) 0.80];
    for b = 1:nb
        x0 = ax.Position(1) + ax.Position(3) + 0.02 + (b-1)*0.10;
        hax = axes('Position', ax.Position, 'Visible', 'off', ...
                   'HitTest', 'off', 'HandleVisibility', 'off');
        colormap(hax, cbars{b}{1});
        set(hax, 'CLim', cbars{b}{2});
        cb = colorbar(hax, 'Position', [x0 0.20 0.025 0.60]);
        cb.Label.String = cbars{b}{3};
        if ~strcmpi(cbars{b}{4}, 'linear')
            % bar runs in g = T(f); label it in f
            Tb = cbars{b}{5};  cf = cbars{b}{6};
            [tf, tl] = scaleTicks(cbars{b}{4}, cf(1), cf(2), cbars{b}{7});
            cb.Ticks = Tb(tf);  cb.TickLabels = tl;
        end
    end
    axes(ax);                                   %#ok<LAXES> back to main axes
end
end

% ==========================================================================
function [fq, info] = angularResample(x, y, f, xc, yc, thq, opt, k)
%ANGULARRESAMPLE  f on a closed curve -> f(theta) about (xc, yc).
x = x(:); y = y(:); f = f(:);
if ~isequal(numel(x), numel(y), numel(f))
    error('Field %d: x, y and f must have the same number of elements.', k);
end

% 0. drop a repeated closing point
scale = max(max(x) - min(x), max(y) - min(y));
if hypot(x(end) - x(1), y(end) - y(1)) <= 1e-12 * scale
    x(end) = []; y(end) = []; f(end) = [];
end
N = numel(x);

% 1. translate to the centre
dx = x - xc;  dy = y - yc;  r = hypot(dx, dy);
if any(r <= 1e-12 * scale)
    error('Field %d: the centre lies on the curve (angle undefined).', k);
end

% 2. four-quadrant angle
th = atan2(dy, dx);

% 3. signed angle steps, wrapped into [-pi, pi)
dth = diff([th; th(1)]);
dth = mod(dth + pi, 2*pi) - pi;

% 4. winding number and star-shape test
W    = round(sum(dth) / (2*pi));
star = abs(W) == 1 && (all(dth > 0) || all(dth < 0));

% 4b. backtracking: unwrapped angle along the curve, oriented so it
%     increases overall; run two laps so a fold straddling point 1 is seen
sig = sign(W); if sig == 0, sig = 1; end
Th1 = sig * [0; cumsum(dth(1:end-1))];          % lap 1 (point 1 at 0)
Th  = [Th1; Th1 + 2*pi];                         % lap 2 = lap 1 + 2*pi
lag = cummax(Th) - Th;                           % how far behind the max
lag = lag(N+1:end);                              % lap 2 has full history
backtrack = max(lag);
nearStar  = abs(W) == 1 && rad2deg(backtrack) <= opt.StarTol;

info = struct('windingNumber', W, 'starShaped', star, ...
              'backtrackDeg', rad2deg(backtrack), 'resampled', false, ...
              'nFolded', nnz(lag > 0), ...
              'thetaRaw', mod(th, 2*pi), 'fRaw', f, 'r', r, ...
              'rLim', [], 'rLimF', [], 'band', [], 'colormap', [], 'cLim', [], ...
              'rScale', 'linear', 'linThresh', NaN);

if W == 0
    warning('Field %d: centre is OUTSIDE the curve (winding number 0); f(theta) is not single-valued.', k);
elseif ~star && ~nearStar
    warning(['Field %d: curve folds back by up to %.2f deg (%d of %d points). ', ...
             'Pass ''StarTol'', %.0f (or more) to resolve the folds with ''Fold''.'], ...
             k, rad2deg(backtrack), nnz(lag > 0), N, ceil(rad2deg(backtrack)));
elseif ~star
    fprintf('Field %d: folds up to %.2f deg resolved by ''%s'' (%d points involved).\n', ...
            k, rad2deg(backtrack), opt.Fold, nnz(lag > 0));
end

fq = NaN(size(thq));
if ~(star || nearStar), return; end
info.resampled = true;
M  = numel(thq);
dq = 2*pi / M;
t0 = mod(th, 2*pi);

if star || strcmpi(opt.Fold, 'remove')
    % 5a. (optionally) drop samples behind the running maximum, then
    %     sort by angle and interpolate periodically
    keep = true(N, 1);
    if ~star, keep = lag <= 0; end
    [ths, idx] = sort(t0(keep));
    fk = f(keep);  fs = fk(idx);
    [ths, ~, g] = unique(ths);
    fs = accumarray(g, fs, [], @mean);
    fq = interp1([ths - 2*pi; ths; ths + 2*pi], [fs; fs; fs], thq, opt.Method);
else
    % 5b. 'average': mean of all samples in each bin centred on thq(j),
    %     bin j = [thq(j) - dq/2, thq(j) + dq/2); empty bins interpolated
    b   = mod(round(t0 / dq), M) + 1;
    cnt = accumarray(b, 1,  [M 1]);
    sm  = accumarray(b, f,  [M 1]);
    has = cnt > 0;
    fq(has) = sm(has) ./ cnt(has);
    if any(~has)
        tk = thq(has);  fk = fq(has);
        fq(~has) = interp1([tk - 2*pi; tk; tk + 2*pi], [fk; fk; fk], ...
                           thq(~has), opt.Method);
    end
end
end

% ==========================================================================
function [tt, vv] = insertZeroCrossings(thq, fq)
%INSERTZEROCROSSINGS  Close the periodic grid and insert the linear roots
%   of f between grid points, so every interval has a single sign.
t = [thq(:); thq(1) + 2*pi];
v = [fq(:);  fq(1)];
n = numel(t);
tt = zeros(2*n, 1); vv = zeros(2*n, 1); m = 0;
for k = 1:n-1
    m = m + 1; tt(m) = t(k); vv(m) = v(k);
    if v(k) * v(k+1) < 0
        % linear root of f between t(k) and t(k+1)
        tz = t(k) - v(k) * (t(k+1) - t(k)) / (v(k+1) - v(k));
        m = m + 1; tt(m) = tz; vv(m) = 0;
    end
end
m = m + 1; tt(m) = t(n); vv(m) = v(n);
tt = tt(1:m); vv = vv(1:m);
end

% ==========================================================================
function [cmap, cl] = fieldColormap(opt, k, allf)
%FIELDCOLORMAP  Colormap and colour limits for field k.
fmin = min(allf);  fmax = max(allf);
cl = [NaN NaN];
if ~isempty(opt.CLim), cl = opt.CLim(min(k, size(opt.CLim, 1)), :); end
diverging = fmin < 0 && fmax > 0;
if all(~isnan(cl))
    % fixed limits: choose the colormap from the limits, not the data, so
    % the same CLim gives the same colours in every figure
    diverging = cl(1) < 0 && cl(2) > 0;
    if abs(cl(1) + cl(2)) > 1e-12 * max(abs(cl)) && diverging
        warning(['Field %d: CLim [%g %g] is not symmetric about 0, so white ', ...
                 'in the diverging map is not at f = 0.'], k, cl(1), cl(2));
    end
elseif any(isnan(cl))
    if diverging
        mag = max(abs([fmin fmax]));  cl = [-mag mag];   % 0 at the centre
    else
        cl = [fmin fmax];
    end
end
if cl(1) == cl(2), cl(2) = cl(1) + 1; end
if ~isempty(opt.Colormap) && numel(opt.Colormap) >= k && ~isempty(opt.Colormap{k})
    cmap = opt.Colormap{k};
elseif diverging
    cmap = blueWhiteRed(256);
else
    cmap = parula(256);
end
end

% ==========================================================================
function rgb = valueToRGB(v, cmap, cl)
%VALUETORGB  Linear map of values to colormap rows, clamped to cl.
n  = size(cmap, 1);
u  = (min(max(v, cl(1)), cl(2)) - cl(1)) / (cl(2) - cl(1));   % in [0,1]
rgb = interp1(linspace(0, 1, n).', cmap, u(:));
end

% ==========================================================================
function cmap = blueWhiteRed(n)
%BLUEWHITERED  Diverging colormap: blue (low) - white (0) - red (high).
c  = [0.23 0.30 0.75; 1 1 1; 0.71 0.02 0.15];
cmap = interp1([0; 0.5; 1], c, linspace(0, 1, n).');
end
% ==========================================================================
function [T, Ti, C] = scaleTransform(scale, C, allf, k)
%SCALETRANSFORM  Forward/inverse radial transform g = T(f), f = Ti(g).
switch lower(scale)
    case 'linear'
        T = @(f) f;   Ti = @(g) g;   C = NaN;
    case 'log'
        if any(allf <= 0)
            error(['Field %d: ''log'' needs f > 0 (min f = %g). ', ...
                   'Use ''symlog'' for data with zeros or both signs.'], k, min(allf));
        end
        T = @(f) log10(f);   Ti = @(g) 10.^g;   C = NaN;
    case 'symlog'
        if isnan(C)
            C = 1e-2 * max(abs(allf));        % two decades of log range
            if ~(C > 0), C = 1; end
        end
        if ~(C > 0), error('Field %d: LinThresh must be > 0.', k); end
        % bi-symmetric log (Webber 2013): odd, monotone, T(0) = 0
        T  = @(f) sign(f) .* log10(1 + abs(f) ./ C);
        Ti = @(g) sign(g) .* C .* (10.^abs(g) - 1);
end
end

% ==========================================================================
function [tf, tl] = scaleTicks(scale, lo, hi, C)
%SCALETICKS  Tick values (in f) and labels for the radial axis / colour bar.
switch lower(scale)
    case 'linear'
        tf = [lo hi];
        if lo < 0 && hi > 0, tf = [lo 0 hi]; end
    case 'log'
        tf = 10.^(ceil(log10(lo)):floor(log10(hi)));
        tf = thinTicks(tf, 4);
        if numel(tf) < 2, tf = [lo tf hi]; end
    case 'symlog'
        m  = max(abs([lo hi]));
        d  = 10.^(ceil(log10(C)):floor(log10(m)));   % decades beyond C
        d  = thinTicks(d, 3);
        tf = [-fliplr(d) 0 d];
        tf = tf(tf >= lo & tf <= hi);
end
tf = unique(tf);

% drop labels that would overlap: keep 0, then larger |f| first, and accept
% a tick only if it is >= 12% of the axis length from every kept tick
if ~strcmpi(scale, 'linear') && numel(tf) > 2
    switch lower(scale)
        case 'log',    gp = log10(tf);   glo = log10(lo);  ghi = log10(hi);
        case 'symlog', g  = @(f) sign(f) .* log10(1 + abs(f) ./ C);
                       gp = g(tf);       glo = g(lo);      ghi = g(hi);
    end
    gap = 0.12 * (ghi - glo);
    [~, order] = sortrows([tf(:) ~= 0, -abs(tf(:))]);   % 0 first, then big
    keep = false(size(tf));
    for i = order(:).'
        if ~any(keep) || all(abs(gp(i) - gp(keep)) >= gap), keep(i) = true; end
    end
    tf = tf(keep);
end
tl = arrayfun(@(v) sprintf('%.3g', v), tf, 'UniformOutput', false);
end

function t = thinTicks(t, n)
%THINTICKS  Keep at most n ticks, always including the largest.
if numel(t) > n
    t = t(unique(round(linspace(numel(t), 1, n))));
end
end
%%
function [vd_avg, theta_g, V, nvalid] = angular_time_average(vd, p0n, cx, cy, t0, halfWin, nTheta, opts)
%ANGULAR_TIME_AVERAGE  Average a field defined on a moving curve over time,
%   at fixed polar angles measured from a fixed center (cx, cy).
%
%   [vd_avg, theta_g] = angular_time_average(vd, p0n, cx, cy, t0, halfWin, nTheta)
%
%   Inputs
%     vd       : cell {Nt x 1}, vd{t} is [N_t x k] (field at the N_t curve nodes;
%                k = 1 for a scalar field, k > 1 for vector components)
%     p0n      : cell {Nt x 2}, p0n{t,1} = x coords [N_t x 1], p0n{t,2} = y coords
%     cx, cy   : center. Scalars (e.g. cx(1), cy(1)), or vectors of length Nt
%                if you want the center to follow the time step.
%     t0       : central time step of the averaging window
%     halfWin  : window is t0-halfWin : t0+halfWin (clipped to [1, Nt])
%     nTheta   : number of angular grid points on [-pi, pi)
%     opts     : (optional) struct with fields
%                  .method   'interp' (default) or 'bin'
%                  .periodic true (default, closed curve) / false (open arc)
%                  .tIdx     explicit vector of time indices (overrides t0/halfWin)
%
%   Outputs
%     vd_avg   : [nTheta x k] time-averaged field vs angle
%     theta_g  : [nTheta x 1] angular grid (rad), in [-pi, pi)
%     V        : [nTheta x k x nW] the per-time-step resampled fields
%     nvalid   : [nTheta x k] number of time steps that contributed per angle
%
%   Method (see accompanying notes):
%     1) theta_i(t) = atan2(y_i - cy, x_i - cx)            (polar angle per node)
%     2) resample vd{t} onto the fixed grid theta_g:
%          'interp': linear interpolation in theta (periodic wrap if closed)
%          'bin'   : arc-length-weighted mean of nodes in each angular bin
%     3) vd_avg(theta_g) = (1/n) * sum_t V(theta_g, t)   (NaNs ignored)

if nargin < 8 || isempty(opts), opts = struct; end
if ~isfield(opts,'method'),   opts.method   = 'interp'; end
if ~isfield(opts,'periodic'), opts.periodic = true;     end

Nt = numel(vd);
if isfield(opts,'tIdx') && ~isempty(opts.tIdx)
    tIdx = opts.tIdx(:).';
else
    tIdx = max(1, t0-halfWin) : min(Nt, t0+halfWin);
end
nW = numel(tIdx);

% fixed angular grid (bin centres), uniform on [-pi, pi)
dth     = 2*pi/nTheta;
theta_g = (-pi + dth*(0:nTheta-1)).';          % nTheta x 1
edges   = [theta_g - dth/2; theta_g(end) + dth/2]; % bin edges

% vd{t} is empty for frames outside the Savitzky-Golay range
% (t <= half or t > nFrames-half); those frames are skipped (left NaN)
hasData = ~cellfun(@isempty, vd(tIdx));
if ~any(hasData)
    error('angular_time_average: no vd data in frames %d-%d.', tIdx(1), tIdx(end));
end
if ~all(hasData)
    warning('angular_time_average: frames %s have no vd data and are skipped.', mat2str(tIdx(~hasData)));
end
k = size(vd{tIdx(find(hasData, 1))}, 2);
V = nan(nTheta, k, nW);

for w = 1:nW
    t   = tIdx(w);
    if ~hasData(w), continue; end
    x   = p0n{t,1}(:);
    y   = p0n{t,2}(:);
    f   = vd{t};
    if isvector(f) && k == 1, f = f(:); end

    xc = cx(min(numel(cx), t));  yc = cy(min(numel(cy), t));

    % --- step 1: polar angle of every node, in (-pi, pi]
    th = atan2(y - yc, x - xc);

    switch lower(opts.method)
        case 'interp'
            % --- step 2a: sort by angle and merge duplicate angles (mean)
            [th, ord] = sort(th);
            f = f(ord, :);
            [thu, ~, g] = unique(th);
            fu = zeros(numel(thu), k);
            for c = 1:k
                fu(:,c) = accumarray(g, f(:,c), [], @mean);
            end
            if opts.periodic
                % --- step 2b: periodic extension theta +/- 2*pi
                thx = [thu - 2*pi; thu; thu + 2*pi];
                fx  = [fu; fu; fu];
                V(:,:,w) = interp1(thx, fx, theta_g, 'linear');
            else
                % open arc: no wrap, NaN outside the covered angular range
                V(:,:,w) = interp1(thu, fu, theta_g, 'linear', NaN);
            end

        case 'bin'
            % --- step 2c: arc-length weights (trapezoidal share per node)
            ds  = hypot(diff(x), diff(y));
            if opts.periodic
                dsc = hypot(x(1)-x(end), y(1)-y(end));
                wts = 0.5*([dsc; ds] + [ds; dsc]);
            else
                wts = 0.5*([0; ds] + [ds; 0]);
            end
            % map theta to [-pi-dth/2, pi-dth/2) so the bin around -pi wraps
            thw = mod(th + pi + dth/2, 2*pi) - pi - dth/2;
            b   = discretize(thw, edges);
            ok  = ~isnan(b);
            W   = accumarray(b(ok), wts(ok), [nTheta 1]);
            for c = 1:k
                S = accumarray(b(ok), wts(ok).*f(ok,c), [nTheta 1]);
                tmp = S ./ W;
                tmp(W == 0) = NaN;               % empty bin
                V(:,c,w) = tmp;
            end
        otherwise
            error('opts.method must be ''interp'' or ''bin''.');
    end
end

% --- step 3: time average at fixed angle (ignore NaN = no data)
nvalid = sum(~isnan(V), 3);
vd_avg = sum(V, 3, 'omitnan') ./ nvalid;
vd_avg(nvalid == 0) = NaN;
end
