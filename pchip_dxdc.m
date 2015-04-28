function [dp, data_interp, t, t0, data] = pchip_dxdc(data, varargin)
% returns an array of derivatives of the spline points to perturbations in
% control point positions, dx / dc
%
% Input
% =====
% - data = [x0, y0, ...] -- spline control points
% - dr        (optional) -- derivative step (default = 1)
%
% Output
% ======
% - dp = [dx/dc_x[1], dy/dc_y[1] |...
%         dx/dc_x[2], dy/dc_y[2]  ... ]
%       -- [T * 2 * N ] array,
%           where 
%


[data_interp, t, t0, data] = interp_implicit_pchip(data);
% figure
% plot(x, y)

[N, data_dim] = size(data);
if nargin > 1
    dr  = varargin{1};
else
    dr = 1;
end

for jj = data_dim:-1:1
    for ii = N:-1:1
        data_perturbed = data;
        data_perturbed(ii, jj) = data(ii,jj) + dr;
        
        pp = pchip(t0, data_perturbed(:, jj));
        dp(:, jj, ii) = (ppval(pp, t) - data_interp(:,jj)) / dr;
    end
end
end