function varargout = curvewise_edge_energy(edge_img, ctrl_points, varargin)
% integrates displacement vector for the snake by integrating over all
% curve points

verbose = false;
[deriv_curve_ctrl, curve_points, t, t0] = pchip_dxdc(ctrl_points, varargin{:});

% curve_img_intensity = integrate_over_edge(edge_img, curve_points(:,2), curve_points(:,1) );
for kk = size(edge_img,3):-1:1
    curve_img_intensity(:,:,kk) = interp2(edge_img(:,:,kk), curve_points(:,2), curve_points(:,1) );
end
dr_dt = interp1((1.5:1:numel(t))', diff(curve_points),t,'linear','extrap');

for ii = size(ctrl_points, 1):-1:1
    for jj = size(ctrl_points, 2):-1:1
        curve_energy = - deriv_curve_ctrl(:, jj, ii) .*  dr_dt(:, jj);
        dE_dc(ii, jj, :) = nansum( bsxfun(@times, curve_img_intensity, curve_energy),1 );% ./ sum(curve_energy);
    end
end

dE_dc = dE_dc * numel(t0) / numel(t);

if verbose
    m = 3;
    figure;
    subplot(m,1,1)
    imagesc(edge_img(:,:,1)); hold all; plot(curve_points(:,2), curve_points(:,1), 'w-'); plot(ctrl_points(:,2), ctrl_points(:,1), 'wx')
    xl= get(gca, 'xlim');
    subplot(m,1,2);
    plot(t, curve_img_intensity(:,:,1), 'k-')
    hold all
    set(gca, 'xlim', xl)
    subplot(m,1,3);
    plot(t0, dE_dc(:,1,1), 'r.')    
    set(gca, 'xlim', xl)
end
% Interp2, can give nan's if contour close to border
dE_dc(isnan(dE_dc))=0;

if  size(edge_img,3) == size(ctrl_points, 2)
    % remove junk if needed
    varargout{1} = zeros(size(ctrl_points));
    for kk = 1:size(ctrl_points, 2)
        varargout{1}(:,kk) = dE_dc(:,kk,kk);
    end
else
    varargout = {dE_dc};
end

end