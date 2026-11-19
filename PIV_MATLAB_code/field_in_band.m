function [bandAxis, fieldAlongBand, bandXq, bandYq, interfaceX, interfaceY, info] = ...
    field_in_band(perx0, pery0, X0, Y0, field, d1, d2, varargin)
%FIELD_IN_BAND  Profile of a field averaged over the band d1 < d(x) < d2,
%where d is the SIGNED distance to the interface (perx0, pery0)
%(d < 0 inside, d > 0 outside), as a function of arc length along the
%interface.
%
%   [bandAxis, fieldAlongBand, bandXq, bandYq, interfaceX, interfaceY, info] = ...
%       field_in_band(perx0, pery0, X0, Y0, field, d1, d2, ...)
%
% METHOD
%   1. Sample points x_i are taken inside the band: either a regular grid
%      of step h ('Sampling','grid', default) or the Lagrangian points
%      themselves ('Sampling','points').
%   2. For every sample the EXACT closest point c_i on the interface
%      polyline is found by point-to-segment projection. This gives
%        d_i = +/-|x_i - c_i|        (sign from inpolygon)
%        s_i = arc length of c_i along the interface, in [0, L)
%        n_i = (x_i - c_i)/d_i       (outward unit normal = grad d)
%   3. The field is evaluated at x_i (scatteredInterpolant, or directly
%      for 'points'), and projected with each sample's OWN normal n_i.
%   4. Samples are aggregated along s, either in NBins equal-width bins
%      ('Aggregate','bin') or with a periodic Gaussian moving window of
%      std 'SigmaS' ('Aggregate','gaussian'). With the grid sampling each
%      sample represents the same area h^2, so the result is the
%      AREA-WEIGHTED mean of the field over the part of the band that
%      projects onto that stretch of interface.
%
% INPUTS
%   perx0, pery0 : interface polyline (closed; closing segment is implied)
%   X0, Y0       : Lagrangian point positions
%   field        : 1x1 cell {scalar}, 1x2 cell {ux, uy}, or 1x3 cell
%                  {Txx, Txy, Tyy} -- same convention as field_along_band
%   d1, d2       : signed distances bounding the band, d1 < d2
%                  (e.g. 5, 15 = outside;  -15, -5 = inside;  -5, 5 = straddling)
%
% NAME-VALUE OPTIONS
%   'Component'  'radial' (default) | 'tangential' | 'magnitude'
%                (for 3-component tensors 'magnitude' returns the
%                determinant Txx*Tyy - Txy^2, as in field_along_band)
%   'NBins'      number of output positions along the interface (200)
%   'Sampling'   'grid' (default) | 'points'
%   'GridStep'   grid spacing h for 'grid' (default (d2-d1)/10)
%   'Aggregate'  'bin' (default) | 'gaussian'
%   'SigmaS'     Gaussian window std in arc-length units
%                (default 2*L/NBins); only used with 'gaussian'
%   'MinCount'   bins/windows with fewer (effective) samples -> NaN (1)
%
% OUTPUTS (row vectors, 1 x NBins)
%   bandAxis        arc length along the INTERFACE at the output positions
%   fieldAlongBand  band-averaged field component
%   bandXq, bandYq  centroid of the samples contributing to each position
%                   (where in the band that value "lives"; for QC plots)
%   interfaceX/Y    interface point at arc length bandAxis
%   info            struct: count, sd, L, sample data, area check
%
% NOTE: because the axis is interface arc length, profiles computed for
% different (d1, d2) share the same axis and are directly comparable.

    q = inputParser;
    addParameter(q, 'Component', 'radial');
    addParameter(q, 'NBins', 200, @(x) isscalar(x) && x >= 1);
    addParameter(q, 'Sampling', 'grid', @(s) any(strcmpi(s, {'grid','points'})));
    addParameter(q, 'GridStep', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(q, 'Aggregate', 'bin', @(s) any(strcmpi(s, {'bin','gaussian'})));
    addParameter(q, 'SigmaS', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(q, 'MinCount', 1, @(x) isscalar(x) && x >= 0);
    parse(q, varargin{:});
    o = q.Results;
    nBins = round(o.NBins);

    if ~(d2 > d1)
        error('field_in_band:badRange', 'Need d1 < d2 (got %.3g, %.3g).', d1, d2);
    end

    % ---- clean the polygon: columns, no zero-length segments ----------
    px = perx0(:);  py = pery0(:);
    keep = [true; hypot(diff(px), diff(py)) > 0];
    px = px(keep);  py = py(keep);
    if hypot(px(1)-px(end), py(1)-py(end)) == 0
        px(end) = [];  py(end) = [];
    end
    M  = numel(px);
    ax = px;  ay = py;                               % segment starts
    ex = circshift(px, -1) - px;                     % segment vectors
    ey = circshift(py, -1) - py;                     % (closing segment included)
    segLen = hypot(ex, ey);
    S0 = [0; cumsum(segLen(1:end-1))];               % arc length at segment starts
    L  = sum(segLen);

    % ---- candidate sample points ------------------------------------
    usePoints = strcmpi(o.Sampling, 'points');
    if usePoints
        qx = X0(:);  qy = Y0(:);
        h  = NaN;
    else
        h = o.GridStep;
        if isempty(h), h = (d2 - d1) / 10; end
        pad = max(abs([d1 d2])) + h;
        [xg, yg] = meshgrid(min(px)-pad : h : max(px)+pad, ...
                            min(py)-pad : h : max(py)+pad);
        qx = xg(:);  qy = yg(:);
    end
    okq = isfinite(qx) & isfinite(qy);
    inside = false(size(qx));
    inside(okq) = inpolygon(qx(okq), qy(okq), px, py);
    cand = okq;
    if d1 >= 0, cand = cand & ~inside; end           % band entirely outside
    if d2 <= 0, cand = cand &  inside; end           % band entirely inside
    idx = find(cand);
    nq  = numel(idx);

    % ---- exact closest point on the polyline (chunked) -------------------
    dAbs = zeros(nq,1);  sNear = zeros(nq,1);
    cx   = zeros(nq,1);  cy    = zeros(nq,1);
    chunk = max(1, floor(1e6 / M));
    axT = ax.';  ayT = ay.';  exT = ex.';  eyT = ey.';  l2T = (segLen.^2).';
    for k0 = 1:chunk:nq
        kk = k0:min(nq, k0 + chunk - 1);
        X  = qx(idx(kk));  Y = qy(idx(kk));
        t  = ((X - axT).*exT + (Y - ayT).*eyT) ./ l2T;   % unconstrained minimiser
        t  = min(max(t, 0), 1);                          % clamp to the segment
        dd = (X - (axT + t.*exT)).^2 + (Y - (ayT + t.*eyT)).^2;
        [dmin, j] = min(dd, [], 2);
        tj = t(sub2ind(size(t), (1:numel(kk))', j));
        dAbs(kk)  = sqrt(dmin);
        sNear(kk) = S0(j) + tj .* segLen(j);
        cx(kk)    = ax(j) + tj .* ex(j);
        cy(kk)    = ay(j) + tj .* ey(j);
    end
    sNear = mod(sNear, L);

    dS  = dAbs;
    dS(inside(idx)) = -dAbs(inside(idx));            % signed distance
    sel = dS > d1 & dS < d2 & dAbs > 0;

    xs = qx(idx(sel));  ys = qy(idx(sel));
    ds = dS(sel);       ss = sNear(sel);
    nx = (xs - cx(sel)) ./ ds;                       % outward unit normal
    ny = (ys - cy(sel)) ./ ds;

    % ---- field values at the samples -------------------------------------
    nc = numel(field);
    if ~ismember(nc, [1 2 3])
        error('field_in_band:field', 'field must be a cell with 1, 2 or 3 entries.');
    end
    comp = cell(1, nc);
    for c = 1:nc
        fc = field{c};
        if usePoints
            comp{c} = fc(idx(sel));
            comp{c} = comp{c}(:);
        else
            g = isfinite(X0(:)) & isfinite(Y0(:)) & isfinite(fc(:));
            F = scatteredInterpolant(X0(g), Y0(g), fc(g), 'linear', 'none');
            comp{c} = F(xs, ys);
        end
    end
    tx = -ny;  ty = nx;                              % tangent = normal rotated +90 deg
    cmpName = lower(o.Component);
    switch nc
        case 1
            val = comp{1};
        case 2
            switch cmpName
                case 'radial',     val = comp{1}.*nx + comp{2}.*ny;
                case 'tangential', val = comp{1}.*tx + comp{2}.*ty;
                case 'magnitude',  val = hypot(comp{1}, comp{2});
                otherwise, error('Component must be ''radial'', ''tangential'', or ''magnitude''.');
            end
        case 3
            switch cmpName
                case 'radial',     val = comp{1}.*nx.^2 + 2*comp{2}.*nx.*ny + comp{3}.*ny.^2;
                case 'tangential', val = comp{1}.*tx.^2 + 2*comp{2}.*tx.*ty + comp{3}.*ty.^2;
                case 'magnitude',  val = comp{1}.*comp{3} - comp{2}.^2;
                otherwise, error('Component must be ''radial'', ''tangential'', or ''magnitude''.');
            end
    end
    good = isfinite(val);                            % drop samples outside data hull
    fracDropped = 1 - mean(good);
    xs = xs(good); ys = ys(good); ss = ss(good); ds = ds(good); val = val(good);

    % ---- aggregate along interface arc length ----------------------------
    dsBin    = L / nBins;
    bandAxis = ((1:nBins) - 0.5) * dsBin;            % bin centres
    switch lower(o.Aggregate)
        case 'bin'
            k   = min(floor(ss / dsBin) + 1, nBins);
            W   = accumarray(k, 1,        [nBins 1]);
            Sv  = accumarray(k, val,      [nBins 1]);
            Sv2 = accumarray(k, val.^2,   [nBins 1]);
            Sx  = accumarray(k, xs,       [nBins 1]);
            Sy  = accumarray(k, ys,       [nBins 1]);
        case 'gaussian'
            sig = o.SigmaS;
            if isempty(sig), sig = 2 * dsBin; end
            W = zeros(nBins,1); Sv = W; Sv2 = W; Sx = W; Sy = W;
            ch = max(1, floor(2e6 / max(1, numel(ss))));
            for b0 = 1:ch:nBins
                bb = b0:min(nBins, b0 + ch - 1);
                K  = wrappedGaussPeriodic(bandAxis(bb).' - ss.', sig, L);
                W(bb)   = K * ones(numel(ss),1);
                Sv(bb)  = K * val;
                Sv2(bb) = K * val.^2;
                Sx(bb)  = K * xs;
                Sy(bb)  = K * ys;
            end
    end
    mu  = Sv ./ W;
    sd  = sqrt(max(Sv2 ./ W - mu.^2, 0));
    bad = W < max(o.MinCount, eps);
    mu(bad) = NaN;  sd(bad) = NaN;

    fieldAlongBand = mu.';
    bandXq = (Sx ./ W).';  bandYq = (Sy ./ W).';
    bandXq(bad) = NaN;     bandYq(bad) = NaN;

    interfaceX = interp1([S0; L], [px; px(1)], bandAxis, 'linear');
    interfaceY = interp1([S0; L], [py; py(1)], bandAxis, 'linear');

    % ---- QC info ------------------------------------------------------
    info.count        = W.';
    info.sd           = sd.';
    info.L            = L;
    info.fracDropped  = fracDropped;
    info.samples      = struct('x', xs, 'y', ys, 'd', ds, 's', ss, 'value', val);
    info.areaSampled  = sum(sel) * h^2;              % NaN for 'points'
    info.areaSteiner  = L*(d2 - d1) + pi*(d2^2 - d1^2);  % exact only within the reach
end

% =====================================================================
function K = wrappedGaussPeriodic(d, sigma, T)
%WRAPPEDGAUSSPERIODIC Unnormalised Gaussian wrapped with period T.
%   d is wrapped to [-T/2, T/2); terms |m| > Lm are dropped where
%   (2*Lm+1)*T/2 >= 8*sigma, so each omitted term is <= exp(-32).
    d  = mod(d + T/2, T) - T/2;
    Lm = max(1, ceil((16*sigma/T - 1) / 2));
    K  = zeros(size(d));
    for m = -Lm:Lm
        K = K + exp(-(d + m*T).^2 / (2*sigma^2));
    end
end
