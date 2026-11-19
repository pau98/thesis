%% Inclusion contour video with scale bar and time stamp
% Draws, for each time step t:
%   ax1 : grey-level image (bottom)
%   ax2 : inclusion contours  (reference t1 and current t)
%         + scale bar and time stamp (top)
% and writes every frame to an MP4 file.

%% ===================== USER PARAMETERS =====================
% --- Colours ---
cl1    = [255, 155,  66]/255;      % reference contour (t1)
cl2    = [193, 255, 131]/255;      % current contour (t)
clFill = [122,  57,   0]/255;      % reference fill (t1)

% --- Data / frames ---
% nFrames = ...;                   % number of time steps (defined beforehand)
if ~exist('t1', 'var')
    t1 = 1;                        % reference time step for the orange contour
end
imgIndexFixed = 1;                 % image shown in the background:
                                   %   a number -> always that image (as in your original code)
                                   %   []       -> image of the current time step t
axLimitsUser = [];                 % [] = extent of the image (recommended)
                                   % or fixed [xmin xmax ymin ymax] in image pixels

% --- Physical calibration ---
pixelSize = dx;                % micrometres per pixel
dtMin     = dt;                    % minutes between consecutive time steps
tZero     = 1;                     % time step that corresponds to t = 0

% --- Scale bar ---
showScaleBar   = true;
scaleBarUm     = 50;               % bar length in micrometres
scaleBarColor  = [1 1 1];
scaleBarHeightFrac = 0.008;        % bar thickness as a fraction of the image height
                                   % (0.008 x 1024 px ~ 8 px)
scaleBarLabel  = true;             % write "50 µm" above the bar
scaleBarPos    = 'bottomright';    % 'bottomright' | 'bottomleft'

% --- Time stamp ---
showTime   = true;
timeFormat = 'hhmm';               % 'hhmm' -> "t = 01:30 h" | 'min' -> "t = 90 min" | 'h' -> "t = 1.5 h"
timeColor  = [1 1 1];
timePos    = 'topleft';            % 'topleft' | 'topright'

% --- Annotation text ---
annotFontSize   = 14;              % points
annotFontWeight = 'bold';
marginFrac      = 0.03;            % margin from the frame edge, fraction of the image size

% --- Video ---
videoFile  = 'deg_mesh.mp4';
frameRate  = 15;
quality    = 100;                  % lower this to get a smaller file
exportRes  = 200;                  % dpi used by exportgraphics for each frame
%% ===========================================================

v = VideoWriter(videoFile, 'MPEG-4');
v.FrameRate = frameRate;
v.Quality   = quality;
open(v)

for t = 1:nFrames

    fig = figure('Color', 'w');

    % ---------- ax1: image ----------
    ax1 = axes(fig);
    hold(ax1, 'on')
    if isempty(imgIndexFixed), imgIdx = t; else, imgIdx = imgIndexFixed; end
    I = imread(sprintf('dc_h_t%03d.tif', imgIdx));
    
    h = imagesc(ax1, I);
    clim(ax1, double([min(I(:)) max(I(:))]))
    uistack(h, 'bottom')
    colormap(ax1, 'gray')

    % Axes limits follow the size of the image actually loaded
    [imgH, imgW] = size(I, [1 2]);
    if isempty(axLimitsUser)
        axLimits = [0.5, imgW + 0.5, 0.5, imgH + 0.5];   % exact pixel-edge extent of imagesc
    else
        axLimits = axLimitsUser;
    end

    % ---------- ax2: contours ----------
    ax2 = axes(fig);
    hold(ax2, 'on')
    fill(ax2, [p0n{t,1};  p0n{t,1}(1)],  [p0n{t,2};  p0n{t,2}(1)],  'k',    'LineWidth', 1)
    fill(ax2, [p0n{t1,1}; p0n{t1,1}(1)], [p0n{t1,2}; p0n{t1,2}(1)], clFill, 'LineWidth', 1)
    plot(ax2, [p0n{t,1};  p0n{t,1}(1)],  [p0n{t,2};  p0n{t,2}(1)],  'Color', cl2, 'LineWidth', 3)
    plot(ax2, [p0n{t1,1}; p0n{t1,1}(1)], [p0n{t1,2}; p0n{t1,2}(1)], 'Color', cl1, 'LineWidth', 3)
    colormap(ax2, 'default')
    clim(ax2, [1 2])

    % ---------- Axes (both identical and linked) ----------
    linkaxes([ax1, ax2])
    for ax = [ax1, ax2]
        axis(ax, 'equal')
        axis(ax, axLimits)
        ax.Visible = 'off';
        ax.XTick = [];  ax.YTick = [];
    end

    % ---------- Scale bar and time stamp (on ax2, the top axes) ----------
    % Drawn in data units (= image pixels), so the bar length is exact
    % whatever the figure size or export resolution.
    xMin = axLimits(1);  xMax = axLimits(2);
    yMin = axLimits(3);  yMax = axLimits(4);
    margin = marginFrac * min(xMax - xMin, yMax - yMin);
    scaleBarHeight = scaleBarHeightFrac * (yMax - yMin);   % bar thickness in image pixels

    if strcmp(ax2.YDir, 'reverse')
        yScreenBottom = yMax;  yScreenTop = yMin;  up = -1;
    else
        yScreenBottom = yMin;  yScreenTop = yMax;  up = +1;
    end

    if showScaleBar
        barLenPx = scaleBarUm / pixelSize;          % 50 / 0.3460 = 144.51 px
        if t == 1 && barLenPx > 0.5*(xMax - xMin)
            warning('Scale bar (%g um = %.0f px) covers more than half of the image width (%.0f px). Consider a smaller scaleBarUm.', ...
                    scaleBarUm, barLenPx, xMax - xMin);
        end
        if strcmp(scaleBarPos, 'bottomleft')
            xL = xMin + margin;  xR = xL + barLenPx;
        else
            xR = xMax - margin;  xL = xR - barLenPx;
        end
        y1 = yScreenBottom + up*margin;
        y2 = y1 + up*scaleBarHeight;
        patch(ax2, [xL xR xR xL], [y1 y1 y2 y2], scaleBarColor, 'EdgeColor', 'none');
        if scaleBarLabel
            text(ax2, (xL+xR)/2, y2 + up*0.5*scaleBarHeight, ...
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
        text(ax2, xT, yScreenTop - up*margin, tStr, ...
            'Color', timeColor, 'FontSize', annotFontSize, ...
            'FontWeight', annotFontWeight, ...
            'HorizontalAlignment', hAlign, 'VerticalAlignment', 'top');
    end

    % ---------- Export frame ----------
    exportgraphics(fig, 'tmp.png', 'Resolution', exportRes);
    img = imread('tmp.png');

    if isempty(v.Height)                                   % first frame: make it even
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
    close(fig);
end

close(v);
delete('tmp.png');
fprintf('Video written to %s\n', videoFile);