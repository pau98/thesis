%% Principal stretches on the reference image  --  video version
% Loops over the time steps, draws the principal-stretch glyphs for each one
% (same drawing code as the single-frame script) and writes the frames to an
% MP4 file. The title shows the elapsed time.
%
% Needs in the workspace: X, Y, p0n, pp, dt  (and nFrames, optional).
% The helper functions (splitLoops, normaliseDir, stretchScalar, divergingMap,
% smoothStretchTensor, strainOf, ellipsesNaN, segmentsNaN, valueToRGB,
% drawColouredSegments, ifelse, segStrainLabel) must stay available: either
% as separate .m files or as local functions at the END of this file.

%% ===================== VIDEO PARAMETERS =====================
if ~exist('nFrames', 'var')
    nFrames = size(pp, 1);          % number of time steps with data
end
tList      = 1:nFrames;             % time steps to include in the video
tZero      = 1;                     % time step that corresponds to t = 0
% dx       = 0.3460;                % micrometres per pixel (defined beforehand)
% dt       = 30;                    % minutes between time steps (defined beforehand)

videoFile  = 'stretch_glyphs.mp4';
frameRate  = 20;
quality    = 100;                   % lower this to get a smaller file
figSize    = [1000 900];            % figure size in screen pixels [width height]
exportRes  = 150;                   % dpi of each exported frame
%% ===========================================================

%% ===================== USER PARAMETERS =====================
imgFile  = sprintf('dc_h_t%03d.tif', 1);     % image: always the reference, t = 1

sortInterface = false;              % true = order the points by angle around
                                    % their centroid (only valid for a
                                    % star-shaped interface, unordered points)

ppIsEigC = true;                    % true : pp{t,1:2} are eigenvalues of C = F'F
                                    %        (= lambda^2) -> take the square root
                                    % false: pp{t,1:2} are already the stretches

% --- Glyphs (principal directions) ---
segStep       = 5;                  % draw every segStep-th grid point (1 = all)
showFamily    = [true true];        % [draw lambda_1 dirs, draw lambda_2 dirs]
glyphType     = 'ellipse';          % 'cross'  : segments of length 2 r lambda_a^alpha
                                    %            (a reference cross of size 2r, stretched)
                                    % 'ellipse': image of a reference circle of radius r
                                    %            under U^alpha, + its axes (optional)
                                    % 'strain' : old behaviour, length ~ |e(lambda)|
amplify       = 1;                  % alpha: exaggeration of the deformation.
                                    % 1 = true stretches; alpha>1 multiplies
                                    % the log strain: ln(lambda^alpha) = alpha ln(lambda)
glyphRadius   = [];                 % r [px]. [] = automatic (largest glyph = segFrac*spacing),
                                    % recomputed for EVERY frame. Fix it to a number to
                                    % keep the same glyph scale across the whole video.
showRefCircle = false;              % faint circle of radius r = undeformed size
refColor      = [0.85 0.85 0.85];
ellColor      = [1 1 1];            % ellipse outline colour
showEllAxes   = true;               % principal axes inside the ellipses
nEll          = 48;                 % points per ellipse
strainMeasure = 'biot';             % strain e(lambda) used for the 'colormap' colours
                                    % and for the 'strain' glyph length:
                                    %   'log'   e = ln(lambda)          (Hencky)
                                    %   'biot'  e = lambda - 1
                                    %   'green' e = (lambda^2 - 1)/2
                                    %   'none'  equal length (directions only; not
                                    %           valid with colorMode = 'colormap')
colorMode     = 'colormap';         % 'colormap': colour = e(lambda) through segCmap
                                    % 'sign'    : colour by extension / contraction
                                    % 'family'  : colour by lambda_1 / lambda_2
segCmap       = [];                 % ['colormap'] [] = blue-white-red (e = 0 white),
                                    % or any n x 3 map, e.g. turbo(256)
segCLim       = [-0.2 0.2];         % ['colormap'] [] = [-1 1]*max|e| per frame (symmetric);
                                    % keep it fixed for a video
segCbLocation = 'westoutside';      % ['colormap'] where its colour bar goes
colExt        = [0.95 0.35 0.20];   % extension   (e > 0)        ['sign']
colCon        = [0.20 0.80 0.95];   % contraction (e < 0)        ['sign']
colFam1       = [0.95 0.85 0.10];   % lambda_1 directions        ['family']
colFam2       = [0.20 0.80 0.95];   % lambda_2 directions        ['family']
segLineWidth  = [3.0 2.5];          % line width [pt] of the principal axes [family 1, family 2]
ellLineWidth  = 2.5;                % line width [pt] of the ellipse outlines
segFrac       = 0.8;                % longest segment = segFrac * spacing
segScale      = [];                 % ['strain' only] [] = automatic (per frame), or a number s:
                                    % segment length [px] = s * |e|
hideInsideIncl = true;              % true = no glyphs inside the inclusion
showLegend     = false;             % legend below the figure (off)

% --- Smoothing: use ALL grid points around each drawn glyph ---
smoothMode   = 'U';                 % 'none' : plain subsampling
                                    % 'logU' : average ln U, then exp (log-Euclidean)
                                    % 'U'    : average U directly (arithmetic)
smoothKernel = 'gauss';             % 'box'  : flat window of side smoothWidth
                                    % 'gauss': Gaussian with sigma = smoothWidth
smoothWidth  = [];                  % in grid points. [] = segStep (box, made odd)
                                    %                      or segStep/2 (gauss)
excludeInclFromAverage = true;      % do not average across the inclusion

% --- Scalar background field ---
scalarField  = 'none';              % 'anisotropy' |ln(lambda_1/lambda_2)|
                                    % 'logJ'       ln(lambda_1*lambda_2)  (area change)
                                    % 'eps1'       ln(lambda_1)
                                    % 'eps2'       ln(lambda_2)
                                    % 'coherency'  |lambda_1^2-lambda_2^2|/(lambda_1^2+lambda_2^2)
                                    % 'none'
scalarAlpha  = 0.6;
scalarCLim   = [];                  % [] = automatic per frame; fix it for a video
showColorbar = true;

% --- Checks on the input directions ---
tolUnit  = 1e-3;                    % warn if | |n| - 1 | exceeds this
tolOrtho = 1e-2;                    % warn if |n1 . n2| exceeds this

% --- Inclusion ---
inclFaceColor = [0.20 0.45 0.85];
inclFaceAlpha = 1;
inclEdgeColor = [0.05 0.15 0.40];   % 'none' to hide the outline
inclLineWidth = 1.0;

% --- Image ---
imgAlpha = 1;
imgCLim  = [];                      % [] = [min max] of the image

% --- Axes ---
yDirection = 'normal';              % keep the same as in plot_deformed_mesh.m
axLimits   = [];                    % [] = image extent
%% ===========================================================

% Values that the drawing code fills in when left empty. They are restored
% at the start of every frame, so each frame behaves exactly like a separate
% run of the single-frame script.
autoInit = struct('glyphRadius', {glyphRadius}, 'segScale', {segScale}, ...
                  'segCLim', {segCLim}, 'segCmap', {segCmap}, ...
                  'smoothWidth', {smoothWidth}, 'scalarCLim', {scalarCLim}, ...
                  'imgCLim', {imgCLim}, 'axLimits', {axLimits});

I = imread(imgFile);                % reference image, the same for every frame

vw = VideoWriter(videoFile, 'MPEG-4');
vw.FrameRate = frameRate;
vw.Quality   = quality;
open(vw)

for t = tList

    % ---------- Restore the automatic parameters ----------
    glyphRadius = autoInit.glyphRadius;   segScale   = autoInit.segScale;
    segCLim     = autoInit.segCLim;       segCmap    = autoInit.segCmap;
    smoothWidth = autoInit.smoothWidth;   scalarCLim = autoInit.scalarCLim;
    imgCLim     = autoInit.imgCLim;       axLimits   = autoInit.axLimits;

    % ---------- Data of time step t ----------
    % REFERENCE positions of the points where pp is defined. CHECK: if pp was
    % computed on the fixed (Eulerian) PIV grid at time t, these are NOT X, Y but
    % the pulled-back positions X - u(X,t) of those grid points.
    xP = X;
    yP = Y;

    % Interface in the reference configuration (NaN-separated loops allowed)
    xI = p0n{t,1}(:);
    yI = p0n{t,2}(:);

    lam1 = pp{t,1};   lam2 = pp{t,2};
    if ppIsEigC
        if any(lam1(:) < 0 | lam2(:) < 0)
            warning('t = %d: negative eigenvalues of C found; C must be positive definite. Set to NaN.', t);
        end
        lam1(lam1 < 0) = NaN;   lam2(lam2 < 0) = NaN;
        lam1 = sqrt(lam1);      lam2 = sqrt(lam2);   % lambda_a = sqrt(eig(C))
    end
    n1x  = pp{t,3};   n1y  = pp{t,4};
    n2x  = pp{t,5};   n2y  = pp{t,6};

    % Split the interface into closed loops (NaN-separated), optionally sorted
    loops = splitLoops(xI, yI, sortInterface);

    % ---------- (0) Check and normalise the directions ----------
    [n1x, n1y, dev1] = normaliseDir(n1x, n1y);
    [n2x, n2y, dev2] = normaliseDir(n2x, n2y);
    if max([dev1(:); dev2(:)], [], 'omitnan') > tolUnit
        warning('t = %d: principal directions were not unit vectors (max | |n|-1 | = %.3g). Normalised.', ...
                t, max([dev1(:); dev2(:)], [], 'omitnan'));
    end
    dotN = abs(n1x.*n2x + n1y.*n2y);
    if max(dotN(:), [], 'omitnan') > tolOrtho
        warning('t = %d: n1 and n2 are not orthogonal (max |n1.n2| = %.3g). Check the pp column order.', ...
                t, max(dotN(:), [], 'omitnan'));
    end

    % Fixed figure size -> every exported frame has the same size
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 figSize], ...
                 'InvertHardcopy', 'off');
    ax = axes(fig); hold(ax,'on');
    ax.SortMethod = 'childorder';

    % ---------- (1) Image as truecolor ----------
    Id = double(I(:,:,1));
    if isempty(imgCLim), imgCLim = [min(Id(:)) max(Id(:))]; end
    s  = (Id - imgCLim(1)) / max(diff(imgCLim), eps);
    s  = min(max(s,0),1);
    hImg = image(ax, repmat(s,1,1,3));
    hImg.AlphaData = imgAlpha;

    % ---------- (2) Scalar field ----------
    [F, fLabel, isSigned] = stretchScalar(lam1, lam2, scalarField);
    if ~isempty(F)
        surface(ax, xP, yP, zeros(size(xP)), F, ...
            'FaceColor','interp', 'EdgeColor','none', 'FaceAlpha', scalarAlpha);
        fmax = max(abs(F(:)), [], 'omitnan');
        if isSigned
            colormap(ax, divergingMap(256));
            if isempty(scalarCLim), scalarCLim = [-1 1]*max(fmax, eps); end
        else
            colormap(ax, parula(256));
            if isempty(scalarCLim), scalarCLim = [0 max(fmax, eps)]; end
        end
        clim(ax, scalarCLim);                    % R2022a+; use caxis on older releases
        if showColorbar
            cb = colorbar(ax);
            cb.Label.String = fLabel;
        end
    end

    % ---------- (3) Inclusion (reference configuration) ----------
    for kL = 1:numel(loops)
        patch(ax, loops{kL}(:,1), loops{kL}(:,2), inclFaceColor, 'FaceAlpha', inclFaceAlpha, ...
              'EdgeColor', inclEdgeColor, 'LineWidth', inclLineWidth);
    end

    % ---------- (4) Principal-direction glyphs ----------
    [ny, nx] = size(xP);
    ri = 1:segStep:ny;   ci = 1:segStep:nx;
    xs = xP(ri,ci);   ys = yP(ri,ci);

    % Points inside the inclusion (full grid), even-odd rule over the loops
    okXY    = isfinite(xP) & isfinite(yP);
    insideF = false(size(xP));
    for kL = 1:numel(loops)
        insideF(okXY) = xor(insideF(okXY), ...
            inpolygon(xP(okXY), yP(okXY), loops{kL}(:,1), loops{kL}(:,2)));
    end

    if ~strcmpi(smoothMode, 'none') && (segStep > 1 || ~isempty(smoothWidth))
        switch lower(smoothKernel)
            case 'box'
                if isempty(smoothWidth), smoothWidth = segStep; end
                k = 2*floor(smoothWidth/2) + 1;           % odd -> centred on the point
                K = ones(k);
            case 'gauss'
                if isempty(smoothWidth), smoothWidth = segStep/2; end
                hw = ceil(3*smoothWidth);
                [gi, gj] = ndgrid(-hw:hw);
                K = exp(-(gi.^2 + gj.^2) / (2*smoothWidth^2));
            otherwise
                error('Unknown smoothKernel "%s".', smoothKernel);
        end
        valid = okXY & isfinite(lam1) & isfinite(lam2) & lam1 > 0 & lam2 > 0 ...
                & isfinite(n1x) & isfinite(n1y) & isfinite(n2x) & isfinite(n2y);
        if excludeInclFromAverage, valid = valid & ~insideF; end
        [l1s, l2s, n1xs, n1ys, n2xs, n2ys] = smoothStretchTensor( ...
            lam1, lam2, n1x, n1y, n2x, n2y, valid, K, smoothMode);
        L  = {l1s(ri,ci),  l2s(ri,ci)};
        Nx = {n1xs(ri,ci), n2xs(ri,ci)};
        Ny = {n1ys(ri,ci), n2ys(ri,ci)};
    else
        L  = {lam1(ri,ci), lam2(ri,ci)};
        Nx = {n1x(ri,ci),  n2x(ri,ci)};
        Ny = {n1y(ri,ci),  n2y(ri,ci)};
    end

    % Log strain of each principal stretch: its sign = extension / contraction
    Lg = cellfun(@(l) strainOf(l, 'log'), L, 'UniformOutput', false);

    drawPt = isfinite(xs) & isfinite(ys);
    if hideInsideIncl
        drawPt = drawPt & ~insideF(ri,ci);
    end
    for a = 1:2
        drawPt = drawPt & isfinite(Lg{a}) & isfinite(Nx{a}) & isfinite(Ny{a});
    end

    % Spacing between drawn glyphs (median distance between grid neighbours)
    % (named dxGrid/dyGrid so the pixel size dx [um/px] is not overwritten)
    dxGrid = hypot(diff(xP,1,2), diff(yP,1,2));
    dyGrid = hypot(diff(xP,1,1), diff(yP,1,1));
    spacing = segStep * min(median(dxGrid(:),'omitnan'), median(dyGrid(:),'omitnan'));

    % Half-length of each glyph axis, H{a}
    if strcmpi(glyphType, 'strain')
        E = cellfun(@(l) strainOf(l, strainMeasure), L, 'UniformOutput', false);
        if isempty(segScale)
            emax = 0;
            for a = find(showFamily)
                emax = max(emax, max(abs(E{a}(drawPt)), [], 'omitnan'));
            end
            segScale = segFrac * spacing / max(emax, eps);
        end
        H = cellfun(@(e) 0.5*segScale*abs(e), E, 'UniformOutput', false);
        ttl = sprintf('segment length = %.3g px per unit of %s strain', segScale, strainMeasure);
    else
        % Reference circle of radius r -> ellipse with semi-axes r*lambda_a^alpha
        G = cellfun(@(lg) exp(amplify*lg), Lg, 'UniformOutput', false);
        if isempty(glyphRadius)
            gmax = max([G{1}(drawPt); G{2}(drawPt)], [], 'omitnan');
            glyphRadius = segFrac * spacing / (2*max(gmax, eps));
        end
        H = cellfun(@(g) glyphRadius*g, G, 'UniformOutput', false);
        ttl = sprintf('reference radius r = %.3g \\mum,   amplification \\alpha = %g', ...
                      glyphRadius*dx, amplify);    % r [um] = r [px] * dx [um/px]
    end

    % Optional undeformed reference circles (radius r)
    if showRefCircle && ~strcmpi(glyphType, 'strain')
        one = ones(nnz(drawPt),1);
        [px, py] = ellipsesNaN(xs(drawPt), ys(drawPt), one, 0*one, ...
                               glyphRadius*one, glyphRadius*one, nEll);
        line(ax, px, py, 'Color', refColor, 'LineWidth', 0.5, 'LineStyle', ':');
    end

    hLeg = gobjects(0);
    if strcmpi(glyphType, 'ellipse')
        [px, py] = ellipsesNaN(xs(drawPt), ys(drawPt), Nx{1}(drawPt), Ny{1}(drawPt), ...
                               H{1}(drawPt), H{2}(drawPt), nEll);
        hLeg(end+1) = line(ax, px, py, 'Color', ellColor, 'LineWidth', ellLineWidth, ...
                           'DisplayName', sprintf('image of circle r under U^{%g}', amplify)); %#ok<SAGROW>
        familiesToDraw = find(showFamily & showEllAxes);
    else
        familiesToDraw = find(showFamily);
    end

    famName = {'\lambda_1', '\lambda_2'};

    if strcmpi(colorMode, 'colormap')
        if strcmpi(strainMeasure, 'none')
            error('colorMode = ''colormap'' needs a strainMeasure other than ''none''.');
        end
        Ec = cellfun(@(l) strainOf(l, strainMeasure), L, 'UniformOutput', false);
        if isempty(segCmap), segCmap = divergingMap(256); end
        if isempty(segCLim)
            emax = 0;
            for a = familiesToDraw
                emax = max(emax, max(abs(Ec{a}(drawPt)), [], 'omitnan'));
            end
            segCLim = [-1 1]*max(emax, eps);
        end
        for a = familiesToDraw
            sel = drawPt;
            if ~any(sel(:)), continue; end
            rgb = valueToRGB(Ec{a}(sel), segCLim, segCmap);
            drawColouredSegments(ax, xs(sel), ys(sel), Nx{a}(sel), Ny{a}(sel), ...
                                 H{a}(sel), rgb, segLineWidth(a));
        end
        familiesToDraw = [];                     % skip the single-colour loop below
    end

    for a = familiesToDraw
        e  = Lg{a};                              % sign decides the colour
        ok = drawPt;
        h  = H{a};                               % half length
        switch colorMode
            case 'sign'
                groups = {ok & e >= 0, ok & e < 0};
                cols   = {colExt, colCon};
                names  = {[famName{a} ' > 1 (ext.)'], [famName{a} ' < 1 (contr.)']};
            otherwise   % 'family'
                groups = {ok};
                cols   = {ifelse(a==1, colFam1, colFam2)};
                names  = {famName{a}};
        end
        for g = 1:numel(groups)
            sel = groups{g};
            if ~any(sel(:)), continue; end
            [px, py] = segmentsNaN(xs(sel), ys(sel), Nx{a}(sel), Ny{a}(sel), h(sel));
            hl = line(ax, px, py, 'Color', cols{g}, 'LineWidth', segLineWidth(a), ...
                      'DisplayName', names{g});
            hLeg(end+1) = hl; %#ok<SAGROW>
        end
    end
    if showLegend && ~isempty(hLeg)
        legend(ax, hLeg, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
               'TextColor', 'k', 'Box', 'off');
    end

    % ---------- Title with the elapsed time ----------
    tMin = (t - tZero) * dt;                 % elapsed time in minutes
    title(ax, sprintf('t = %02d:%02d h (on reference image),   %s', ...
                      floor(tMin/60), mod(round(tMin), 60), ttl));

    % ---------- Axes ----------
    uistack(hImg,'bottom');
    ax.YDir = yDirection;
    axis(ax, 'equal');
    if isempty(axLimits)
        axLimits = [0.5, size(I,2)+0.5, 0.5, size(I,1)+0.5];
    end
    axis(ax, axLimits);
    ax.XTick = [];  ax.YTick = [];  box(ax,'off');  ax.XColor = 'none';  ax.YColor = 'none';

    % Colour bar for the segments, on an invisible overlay axes
    if strcmpi(colorMode, 'colormap')
        drawnow;                                 % let ax settle its final position
        ax2 = axes('Position', ax.Position, 'Visible', 'off', 'HitTest', 'off', ...
                   'PickableParts', 'none');
        colormap(ax2, segCmap);
        clim(ax2, segCLim);
        cb2 = colorbar(ax2, 'Location', segCbLocation);
        cb2.Label.String = segStrainLabel(strainMeasure);
        ax2.Position = ax.Position;              % colorbar shrank ax2: realign
        axes(ax);                                %#ok<LAXES> keep ax current
    end

    % ---------- Export frame ----------
    % print (not exportgraphics) keeps the whole figure area, so the frame
    % size does not change when the title or legend change length.
    drawnow;
    print(fig, 'tmp.png', '-dpng', sprintf('-r%d', exportRes));
    img = imread('tmp.png');

    if isempty(vw.Height)                                  % first frame: make it even
        img = img(1:2*floor(end/2), 1:2*floor(end/2), :);
    else                                                   % later frames: match the first one
        if size(img,1) ~= vw.Height || size(img,2) ~= vw.Width
            fprintf('t = %d: frame %dx%d resized to %dx%d\n', t, ...
                    size(img,1), size(img,2), vw.Height, vw.Width);
        end
        rowIdx = min(1:vw.Height, size(img,1));  % crops if too big, repeats last row if too small
        colIdx = min(1:vw.Width,  size(img,2));
        img    = img(rowIdx, colIdx, :);
    end

    writeVideo(vw, img);
    close(fig);
    fprintf('Frame t = %d written\n', t);
end

close(vw);
delete('tmp.png');
fprintf('Video written to %s\n', videoFile);

% ---- Paste your local helper functions below this line if they are not separate files ----
%% ===================== LOCAL FUNCTIONS =====================
function [nx, ny, dev] = normaliseDir(nx, ny)
% n <- n/|n|, and report how far |n| was from 1.
r   = hypot(nx, ny);
dev = abs(r - 1);
r(r == 0) = NaN;                  % zero vector: no direction, not drawn
nx  = nx ./ r;   ny = ny ./ r;
end

function [l1, l2, n1x, n1y, n2x, n2y] = smoothStretchTensor( ...
            l1, l2, n1x, n1y, n2x, n2y, valid, K, mode)
% Average the 2D stretch tensor over the kernel K and re-diagonalise.
%   'U'    : A = U    = l1 N1(x)N1 + l2 N2(x)N2
%   'logU' : A = ln U = ln(l1) N1(x)N1 + ln(l2) N2(x)N2
% (the old names 'V' / 'logV' are accepted: the algebra is identical)
% Components are averaged by normalised convolution, which ignores invalid
% points (NaN, inside the mask) and renormalises the weights at the edges:
%   Abar = conv(K, w.*A) ./ conv(K, w),   w = 1 on valid points, 0 elsewhere.
% Output: l1 >= l2 everywhere (family 1 = larger stretch).
switch lower(mode)
    case {'u','v'},       f1 = l1;       f2 = l2;
    case {'logu','logv'}, f1 = log(l1);  f2 = log(l2);
    otherwise,   error('Unknown smoothMode "%s".', mode);
end
A11 = f1.*n1x.^2    + f2.*n2x.^2;
A22 = f1.*n1y.^2    + f2.*n2y.^2;
A12 = f1.*n1x.*n1y  + f2.*n2x.*n2y;

w = double(valid);
A11(~valid) = 0;  A22(~valid) = 0;  A12(~valid) = 0;
den = conv2(w, K, 'same');
den(den < 1e-12*sum(K(:))) = NaN;          % no valid point in the window
A11 = conv2(A11, K, 'same') ./ den;
A22 = conv2(A22, K, 'same') ./ den;
A12 = conv2(A12, K, 'same') ./ den;

% Closed-form eigen-decomposition of [A11 A12; A12 A22]
m  = 0.5*(A11 + A22);
r  = hypot(0.5*(A11 - A22), A12);
mu1 = m + r;   mu2 = m - r;                % mu1 >= mu2
th = 0.5*atan2(2*A12, A11 - A22);          % angle of the mu1 eigenvector
n1x = cos(th);   n1y = sin(th);
n2x = -sin(th);  n2y = cos(th);
if any(strcmpi(mode, {'logu','logv'}))
    l1 = exp(mu1);   l2 = exp(mu2);
else
    l1 = mu1;        l2 = mu2;             % > 0: convex combination of SPD tensors
end
end

function e = strainOf(lam, measure)
% Principal value of a Seth-Hill strain, E^(m) = (U^m - I)/m, m -> 0 : ln U.
switch lower(measure)
    case 'log',   e = log(lam);               % m = 0 (Hencky)
    case 'biot',  e = lam - 1;                % m = 1
    case 'green', e = 0.5*(lam.^2 - 1);       % m = 2 (Green-Lagrange)
    case 'none',  e = ones(size(lam));        % directions only
    otherwise,    error('Unknown strainMeasure "%s".', measure);
end
e(~(lam > 0)) = NaN;                          % stretches must be positive
end

function [F, lbl, isSigned] = stretchScalar(l1, l2, name)
% Scalar fields of the 2D stretch tensor (see the derivation notes).
l1(~(l1 > 0)) = NaN;   l2(~(l2 > 0)) = NaN;
isSigned = false;
switch lower(name)
    case 'anisotropy'
        F = abs(log(l1 ./ l2));   lbl = '|ln(\lambda_1/\lambda_2)|';
    case 'logj'
        F = log(l1 .* l2);        lbl = 'ln J = ln(\lambda_1\lambda_2)';  isSigned = true;
    case 'eps1'
        F = log(l1);              lbl = 'ln \lambda_1';                   isSigned = true;
    case 'eps2'
        F = log(l2);              lbl = 'ln \lambda_2';                   isSigned = true;
    case 'coherency'
        F = abs(l1.^2 - l2.^2) ./ (l1.^2 + l2.^2);
        lbl = '|\lambda_1^2-\lambda_2^2| / (\lambda_1^2+\lambda_2^2)';
    case 'none'
        F = [];  lbl = '';
    otherwise
        error('Unknown scalarField "%s".', name);
end
end

function [px, py] = segmentsNaN(x, y, nx, ny, h)
% Segment a = x - h n  ->  b = x + h n for every point, joined into one
% NaN-separated polyline: [a1 b1 NaN a2 b2 NaN ...]. One graphics object
% for all segments, and centred on x by construction.
x = x(:).';  y = y(:).';  nx = nx(:).';  ny = ny(:).';  h = h(:).';
px = [x - h.*nx; x + h.*nx; nan(size(x))];
py = [y - h.*ny; y + h.*ny; nan(size(y))];
px = px(:);  py = py(:);
end

function [px, py] = ellipsesNaN(x, y, n1x, n1y, a1, a2, nPts)
% Ellipses centred on (x,y) with semi-axis a1 along n1 and a2 along n2 = n1
% rotated by +90 deg:  p(th) = x + a1 cos(th) n1 + a2 sin(th) n2.
% All ellipses in one NaN-separated polyline.
x = x(:).';  y = y(:).';  n1x = n1x(:).';  n1y = n1y(:).';  a1 = a1(:).';  a2 = a2(:).';
th = linspace(0, 2*pi, nPts+1).';           % closed: last = first
c  = cos(th);  s = sin(th);
px = x + c.*(a1.*n1x) + s.*(a2.*(-n1y));
py = y + c.*(a1.*n1y) + s.*(a2.*( n1x));
px = [px; nan(1, numel(x))];   py = [py; nan(1, numel(x))];
px = px(:);  py = py(:);
end

function rgb = valueToRGB(v, lim, cmap)
% Linear map of v onto the rows of cmap, clamped to lim:
%   k = round(1 + (N-1) * (v - lim(1)) / (lim(2) - lim(1))),  1 <= k <= N
N = size(cmap, 1);
s = (v(:) - lim(1)) / max(lim(2) - lim(1), eps);
k = round(1 + (N-1) * min(max(s, 0), 1));
rgb = cmap(k, :);
end

function h = drawColouredSegments(ax, x, y, nx, ny, hl, rgb, lw)
% All segments in ONE patch object: vertices a_i = x_i - hl_i n_i and
% b_i = x_i + hl_i n_i, faces [a_i b_i]. Per-vertex truecolour with
% EdgeColor 'flat' paints each segment with its own colour and does not
% use the axes colormap.
x = x(:); y = y(:); nx = nx(:); ny = ny(:); hl = hl(:);
n  = numel(x);
Vt = zeros(2*n, 2);
Vt(1:2:end,:) = [x - hl.*nx, y - hl.*ny];
Vt(2:2:end,:) = [x + hl.*nx, y + hl.*ny];
Ct = zeros(2*n, 3);
Ct(1:2:end,:) = rgb;   Ct(2:2:end,:) = rgb;
h = patch(ax, 'Faces', reshape(1:2*n, 2, n).', 'Vertices', Vt, ...
          'FaceVertexCData', Ct, 'EdgeColor', 'flat', 'FaceColor', 'none', ...
          'LineWidth', lw);
end

function s = segStrainLabel(measure)
switch lower(measure)
    case 'log',   s = 'ln \lambda_a';
    case 'biot',  s = '\lambda_a - 1';
    case 'green', s = '(\lambda_a^2 - 1)/2';
    otherwise,    s = 'e(\lambda_a)';
end
end

function cmap = divergingMap(n)
% Blue - white - red map, linear in RGB, for signed fields centred on 0.
lo = [0.23 0.30 0.75];  mid = [1 1 1];  hi = [0.71 0.02 0.15];
s  = linspace(0, 1, n).';
w  = min(2*s, 1);            % 0 -> 1 on the first half
v  = max(2*s - 1, 0);        % 0 -> 1 on the second half
cmap = (s <= 0.5).*((1-w).*lo + w.*mid) + (s > 0.5).*((1-v).*mid + v.*hi);
end

function out = ifelse(cond, a, b)
if cond, out = a; else, out = b; end
end

function loops = splitLoops(x, y, doSort)
% Split NaN-separated point lists into closed loops (k x 2 arrays).
% doSort: order each loop by angle around its centroid (star-shaped only).
brk  = isnan(x) | isnan(y);
id   = cumsum(brk);
loops = {};
for g = unique(id(~brk)).'
    sel = (id == g) & ~brk;
    P = [x(sel) y(sel)];
    if size(P,1) < 3, continue; end
    if doSort
        c = mean(P, 1);
        [~, o] = sort(atan2(P(:,2) - c(2), P(:,1) - c(1)));
        P = P(o,:);
    end
    if any(P(1,:) ~= P(end,:)), P(end+1,:) = P(1,:); end   %#ok<AGROW> close loop
    loops{end+1} = P; %#ok<AGROW>
end
end