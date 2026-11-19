%% Colours
cl1 = [255, 155, 66]/255;
cl2 = [193, 255, 131]/255;
clFill = [122, 57, 0]/255;

v = VideoWriter('def_mesh.mp4', 'MPEG-4');
v.FrameRate = 15;
v.Quality   = 100;      % lower this to get a smaller file
open(v)
for t = 1:nFrames
%% 2.4. Deformation regular mesh from PIV  --  configurable version
% Draws, in ONE axes:
%   (1) the PIV mesh: deformed chi{t,:}, undeformed chi{tRef,:}, or both
%   (2) the inclusion from the binary mask          (filled patch)
%   (3) the grey-level image                        (bottom layer by default)
%   (4) a scale bar and a time stamp                (always on top)
% All mesh/inclusion colours are explicit RGB triplets.

%% ===================== USER PARAMETERS =====================
% t        = 40;                              % time step to display
tRef     = 1;                                % time step of the undeformed grid
imgFile  = sprintf('dc_h_t%03d.tif',  t);    % image

% --- Grid ---
gridStep   = 5;                    % draw every gridStep-th PIV line (1 = all lines)
gridToShow = 'deformed';           % 'deformed' | 'undeformed' | 'both'

% deformed grid style
gridColor     = cl1;  % RGB in [0,1]
gridLineWidth = 1.5;              % points
gridLineStyle = '-';
meshFaceColor = 'none';            % 'none' or RGB, e.g. [1 1 0.6]
meshFaceAlpha = 0.3;               % only used if meshFaceColor is RGB

% undeformed grid style
refColor      = cl2;
refLineWidth  = 1.5;
refLineStyle  = '-';

% --- Inclusion ---
inclFaceColor = [0.20 0.45 0.85];
inclFaceAlpha = 1;
inclEdgeColor = [0.05 0.15 0.40];  % 'none' to hide the outline
inclLineWidth = 1.0;
maskLevel     = [];                % [] = midpoint of mask range (0/1 -> 0.5, 0/255 -> 127.5)

% --- Image ---
imgAlpha    = 1;                   % 0 = invisible, 1 = opaque
imgColormap = gray(256);
imgCLim     = [];                  % [] = [min max] of the image
imageOnTop  = false;               % false = image underneath mesh and inclusion
                                   % true  = image over everything (original look)

% --- Axes ---
yDirection  = 'normal';            % 'normal' = as in your original script
                                   % 'reverse' = usual image orientation
axLimits    = [];                  % [] = image extent, or [xmin xmax ymin ymax]

% --- Physical calibration ---
pixelSize   = dx;              % micrometres per pixel
dtMin       = dt;                  % minutes between consecutive time steps
tZero       = 1;                   % time step that corresponds to t = 0

% --- Scale bar ---
showScaleBar   = true;
scaleBarUm     = 50;               % bar length in micrometres
scaleBarColor  = [1 1 1];          % white
scaleBarHeight = 8;                % bar thickness in image pixels
scaleBarLabel  = true;             % write "50 µm" above the bar
scaleBarPos    = 'bottomright';    % 'bottomright' | 'bottomleft'

% --- Time stamp ---
showTime      = true;
timeFormat    = 'hhmm';            % 'hhmm' -> "t = 01:30 h"  | 'min' -> "t = 90 min" | 'h' -> "t = 1.5 h"
timeColor     = [1 1 1];
timePos       = 'topleft';         % 'topleft' | 'topright'

% --- Annotation text ---
annotFontSize = 14;                % points (scales with the export resolution)
annotFontWeight = 'bold';
marginFrac    = 0.03;              % margin from the frame edge, fraction of the image size
%% ===========================================================

I    = imread(imgFile);
Xd   = chi{t,1};     Yd = chi{t,2};       % deformed node positions
X0   = chi{tRef,1};  Y0 = chi{tRef,2};    % undeformed node positions

figure('Color','w');
ax = axes; hold(ax,'on');
ax.SortMethod = 'childorder';   % draw in stacking order, so uistack decides what is on top

% ---------- (1) Mesh ----------
% Keep only nodes that belong to complete cells of size gridStep:
% indices 1, 1+k, ..., 1+m*k with m = floor((n-1)/k). The last
% mod(n-1,k) rows/columns would only form partial cells, so they are cut.
[ny, nx] = size(Xd);
nyUse = 1 + floor((ny-1)/gridStep)*gridStep;
nxUse = 1 + floor((nx-1)/gridStep)*gridStep;
ri = 1:gridStep:nyUse;               % rows drawn
ci = 1:gridStep:nxUse;               % columns drawn

hGrid = gobjects(0);  hRef = gobjects(0);  hFace = gobjects(0);

if any(strcmp(gridToShow, {'undeformed','both'}))
    hRef = drawGrid(ax, X0(1:nyUse,1:nxUse), Y0(1:nyUse,1:nxUse), ri, ci, ...
                    refColor, refLineWidth, refLineStyle);
end

if any(strcmp(gridToShow, {'deformed','both'}))
    if ~ischar(meshFaceColor)        % optional filled mesh (cropped to full cells)
        hFace = surface(ax, Xd(1:nyUse,1:nxUse), Yd(1:nyUse,1:nxUse), zeros(nyUse,nxUse), ...
            'FaceColor', meshFaceColor, 'FaceAlpha', meshFaceAlpha, 'EdgeColor', 'none');
    end
    hGrid = drawGrid(ax, Xd(1:nyUse,1:nxUse), Yd(1:nyUse,1:nxUse), ri, ci, ...
                     gridColor, gridLineWidth, gridLineStyle);
end

% ---------- (2) Inclusion ----------
patch(ax,[pnPlots{t,1};pnPlots{t,1}(1)],[pnPlots{t,2};pnPlots{t,2}(1)], inclFaceColor, ...
    'FaceAlpha', inclFaceAlpha, 'EdgeColor', inclEdgeColor, 'LineWidth', inclLineWidth);
% ---------- (3) Image ----------
hImg = imagesc(ax, I);                         % hold is on -> YDir not changed
hImg.AlphaData = imgAlpha;
colormap(ax, imgColormap);
if isempty(imgCLim), imgCLim = double([min(I(:)) max(I(:))]); end
clim(ax, imgCLim);
if imageOnTop, uistack(hImg,'top'); else, uistack(hImg,'bottom'); end

% ---------- Axes ----------
ax.YDir = yDirection;
axis(ax, 'equal');
if isempty(axLimits)
    axLimits = [0.5, size(I,2)+0.5, 0.5, size(I,1)+0.5];
end
axis(ax, axLimits);
ax.Visible = 'off';

% ---------- (4) Scale bar and time stamp ----------
% Everything is drawn in data units (= image pixels), so the bar length is
% exact whatever the figure size or export resolution.
xMin = axLimits(1);  xMax = axLimits(2);
yMin = axLimits(3);  yMax = axLimits(4);
margin = marginFrac * min(xMax - xMin, yMax - yMin);

% "Up on screen" depends on the y-axis direction:
%   YDir 'normal'  -> screen bottom is yMin, screen up is +y
%   YDir 'reverse' -> screen bottom is yMax, screen up is -y
if strcmp(yDirection, 'reverse')
    yScreenBottom = yMax;  yScreenTop = yMin;  up = -1;
else
    yScreenBottom = yMin;  yScreenTop = yMax;  up = +1;
end

hAnnot = gobjects(0);

if showScaleBar
    barLenPx = scaleBarUm / pixelSize;          % 50 / 0.3460 = 144.51 px
    if strcmp(scaleBarPos, 'bottomleft')
        xL = xMin + margin;  xR = xL + barLenPx;
    else
        xR = xMax - margin;  xL = xR - barLenPx;
    end
    y1 = yScreenBottom + up*margin;             % lower edge of the bar on screen
    y2 = y1 + up*scaleBarHeight;                % upper edge of the bar on screen
    hAnnot(end+1) = patch(ax, [xL xR xR xL], [y1 y1 y2 y2], scaleBarColor, ...
        'EdgeColor', 'none');
    if scaleBarLabel
        hAnnot(end+1) = text(ax, (xL+xR)/2, y2 + up*0.5*scaleBarHeight, ...
            sprintf('%g \\mum', scaleBarUm), ...
            'Color', scaleBarColor, 'FontSize', annotFontSize, ...
            'FontWeight', annotFontWeight, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

if showTime
    tMin = (t - tZero) * dtMin;                 % elapsed time in minutes
    switch timeFormat
        case 'min'
            tStr = sprintf('t = %d min', round(tMin));
        case 'h'
            tStr = sprintf('t = %.1f h', tMin/60);
        otherwise   % 'hhmm'
            tStr = sprintf('t = %02d:%02d h', floor(tMin/60), mod(round(tMin), 60));
    end
    if strcmp(timePos, 'topright')
        xT = xMax - margin;  hAlign = 'right';
    else
        xT = xMin + margin;  hAlign = 'left';
    end
    hAnnot(end+1) = text(ax, xT, yScreenTop - up*margin, tStr, ...
        'Color', timeColor, 'FontSize', annotFontSize, ...
        'FontWeight', annotFontWeight, ...
        'HorizontalAlignment', hAlign, 'VerticalAlignment', 'top');
end

for h = hAnnot, uistack(h, 'top'); end         % annotations above image, mesh and inclusion

%% Export frame
exportgraphics(gcf, 'tmp.png', 'Resolution', 200);   % or print(fig,'tmp.png','-dpng','-r200')
img = imread('tmp.png');
if isempty(v.Height)                                   % first frame: just make it even
    img = img(1:2*floor(end/2), 1:2*floor(end/2), :);
else                                                   % later frames: match the first one
    if size(img,1) ~= v.Height || size(img,2) ~= v.Width
        fprintf('t = %d: frame %dx%d resized to %dx%d\n', t, ...
                size(img,1), size(img,2), v.Height, v.Width);
    end
    rows = min(1:v.Height, size(img,1));   % crops if too big, repeats last row if too small
    cols = min(1:v.Width,  size(img,2));
    img  = img(rows, cols, :);
end
writeVideo(v, img);
close all;
end

close(v)

%% ===================== LOCAL FUNCTIONS =====================
function h = drawGrid(ax, X, Y, ri, ci, col, lw, ls)
% Draw rows ri and columns ci of the node arrays X, Y as one line object.
% Every node along a line is kept, so curvature is preserved.
Xh = X(ri,:).';  Yh = Y(ri,:).';     % one column per selected row
Xv = X(:,ci);    Yv = Y(:,ci);       % one column per selected column
Xh(end+1,:) = NaN;  Yh(end+1,:) = NaN;   % NaN separators
Xv(end+1,:) = NaN;  Yv(end+1,:) = NaN;
h = plot(ax, [Xh(:); Xv(:)], [Yh(:); Yv(:)], ls, 'Color', col, 'LineWidth', lw);
end

function h = drawContourPatches(ax, C, fc, fa, ec, lw)
% Parse a contour matrix: each block starts with [level; nVertices],
% followed by nVertices columns of [x; y]. One patch per closed loop.
h = gobjects(0);
k = 1;
while k < size(C,2)
    n  = C(2,k);
    xy = C(:, k+1:k+n);
    h(end+1) = patch(ax, xy(1,:), xy(2,:), fc, ...  %#ok<AGROW>
        'FaceAlpha', fa, 'EdgeColor', ec, 'LineWidth', lw);
    k = k + n + 1;
end
end