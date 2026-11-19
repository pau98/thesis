cx0 = cx(1);  cy0 = cy(1);                 
Rmax = min([max(X(:))-cx0, cx0-min(X(:)), max(Y(:))-cy0, cy0-min(Y(:))]);
% Rmax = 300;
nCircles = 6;  nTheta = 200;  nSpokes = 16;
linkMode = 'path';      % 'chord' (straight X -> x) or 'path' (trajectory over frames)
showSpokes = false;      % deformed radial material lines

radii = linspace(Rmax/nCircles, Rmax, nCircles);
theta = linspace(0, 2*pi, nTheta);
thS   = linspace(0, 2*pi, nSpokes+1); thS(end) = [];
cols  = lines(nCircles);                   % one colour per circle
colsRef = 0.45*cols + 0.55;                % lighter version for reference

% --- Reference material points (Step 1) ---
[TH, R]   = meshgrid(theta, radii);        % nCircles x nTheta
CX = cx0 + R.*cos(TH);    CY = cy0 + R.*sin(TH);
[THs, Rs] = meshgrid(thS, radii);          % nCircles x nSpokes (link points)
PX = cx0 + Rs.*cos(THs);  PY = cy0 + Rs.*sin(THs);
rr = linspace(0, Rmax, 60);
[RRl, THl] = meshgrid(rr, thS);            % nSpokes x 60 (radial lines)
LX = cx0 + RRl.*cos(THl); LY = cy0 + RRl.*sin(THl);

% --- Map everything through chi (Step 2) ---
dCX = zeros([size(CX) nFrames]); dCY = dCX;
dPX = zeros([size(PX) nFrames]); dPY = dPX;
dLX = zeros([size(LX) nFrames]); dLY = dLX;
for t = 1:nFrames
    dCX(:,:,t) = interp2(X, Y, chi{t,1}, CX, CY, 'spline');
    dCY(:,:,t) = interp2(X, Y, chi{t,2}, CX, CY, 'spline');
    dPX(:,:,t) = interp2(X, Y, chi{t,1}, PX, PY, 'spline');
    dPY(:,:,t) = interp2(X, Y, chi{t,2}, PX, PY, 'spline');
    dLX(:,:,t) = interp2(X, Y, chi{t,1}, LX, LY, 'spline');
    dLY(:,:,t) = interp2(X, Y, chi{t,2}, LX, LY, 'spline');
end

% Fixed axis limits so the GIF doesn't jump between frames
allX = [CX(:); dCX(:)];  allY = [CY(:); dCY(:)];
spanX = max(allX) - min(allX);
spanY = max(allY) - min(allY);
pad = 0.05*max(spanX, spanY);
lims = [min(allX)-pad, max(allX)+pad, min(allY)-pad, max(allY)+pad];


% --- Video ---
videoFile  = 'rad_disp.mp4';
frameRate  = 15;
quality    = 100;                  % lower this to get a smaller file
exportRes  = 200;                  % dpi used by exportgraphics for each frame
v = VideoWriter(videoFile, 'MPEG-4');
v.FrameRate = frameRate;
v.Quality   = quality;
open(v)
for t = 1:nFrames
    figure()
    nexttile;
    % Background image setup
    bgFile  = sprintf('dc_h_t%03d.tif',t);   % <-- your image (or use an already-loaded PIV frame)
    bgAlpha = 0;             % opacity: 0 = invisible, 1 = opaque
    sx = 1;  sy = 1;           % physical units per pixel (1 if X,Y are in pixels)
    
    bg = imread(bgFile);
    if ismatrix(bg)                        % grayscale -> RGB, so the figure colormap
        bg = repmat(bg, [1 1 3]);          % doesn't interfere with other graphics
    end
    bg = im2double(bg);
    % bg = flipud(bg);                       % make row index increase upward
    [Nr, Nc, ~] = size(bg);
    
    % Pixel (u0, v0) of the ORIGINAL (unflipped) image that should sit at (cx0, cy0)
    u0 = cx0;  v0 = cy0;                   % <-- set these (see option below)
    v0f = Nr - v0 + 1;                     % anchor row after flipud
    
    bgXData = cx0 + sx*([1 Nc] - u0);
    bgYData = cy0 + sy*([1 Nr] - v0f);
    

    cla; hold on
        image(bg, 'CData', bg, ...
              'AlphaData', bgAlpha, 'AlphaDataMapping', 'none');
        set(gca, 'YDir', 'normal')         % y up, consistent with the flipud
    for k = 1:nCircles
        % reference circle (light, dashed, same hue as its deformed copy)
        plot(CX(k,:), CY(k,:), '--', 'Color', colsRef(k,:));
    
        % links (Steps 3/4)
        switch linkMode
            case 'chord'
                plot([PX(k,:); dPX(k,:,t)], [PY(k,:); dPY(k,:,t)], ...
                     '-', 'Color', colsRef(k,:), 'LineWidth', 0.8);
            case 'path'
                trajX = permute(dPX(k,:,1:t), [3 2 1]);   % t x nSpokes
                trajY = permute(dPY(k,:,1:t), [3 2 1]);
                plot([PX(k,:); trajX], [PY(k,:); trajY], ...
                     '-', 'Color', colsRef(k,:), 'LineWidth', 0.8);
        end
        plot(PX(k,:), PY(k,:), '.', 'Color', colsRef(k,:), 'MarkerSize', 8);
    
        % deformed circle + its material points
        plot(dCX(k,:,t), dCY(k,:,t), '-', 'Color', cols(k,:), 'LineWidth', 1.4);
        plot(dPX(k,:,t), dPY(k,:,t), 'o', 'Color', cols(k,:), ...
             'MarkerFaceColor', cols(k,:), 'MarkerSize', 3);
        % % theta = 0 point marked to reveal rotation
        % plot(dPX(k,1,t), dPY(k,1,t), 's', 'Color', 'k', ...
        %      'MarkerFaceColor', cols(k,:), 'MarkerSize', 6);
    end
    
    if showSpokes   % deformed radial material lines (Step 5)
        plot(dLX(:,:,t).', dLY(:,:,t).', '-', 'Color', [0.5 0.5 0.5 0.5]);
    end
    
    % plot(cx0, cy0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
    fill(pn{t,1}, pn{t,2}, 'k', 'FaceAlpha', 0.7, 'EdgeColor', 'k', 'LineWidth', 1.5)
    % fill(pn{1,1}, pn{1,2}, 'k', 'FaceAlpha', 0., 'EdgeColor', 'w', 'LineWidth', 1.5,'EdgeAlpha',0.7)
    
    axis equal tight; axis(lims)
    set(gca, 'XTick', [], 'YTick', [])
    tMin = (t - 1) * dt;                   % elapsed time in minutes (t = 1 -> 0 min)
    title(sprintf('t = %02d:%02d h', floor(tMin/60), mod(round(tMin), 60)))
    % ---------- Export frame ----------
    exportgraphics(gcf, 'tmp.png', 'Resolution', exportRes);
    img = imread('tmp.png');

    if isempty(v.Height)                                   % first frame: make it even
        img = img(1:2*floor(end/2), 1:2*floor(end/2), :);
    else                                                   % later frames: match the first one
        if size(img,1) ~= v.Height || size(img,2) ~= v.Width
            fprintf('t = %d: frame %dx%d resized to %dx%d\n', t, ...
                    size(img,1), size(img,2), v.Height, v.Width);
        end
        rowIdx = min(1:v.Height, size(img,1));   % crops if too big, repeats last row if too small
        colIdx = min(1:v.Width,  size(img,2));
        img    = img(rowIdx, colIdx, :);
    end

    writeVideo(v, img);
    close(gcf);
end

close(v);
delete('tmp.png');
fprintf('Video written to %s\n', videoFile);