function dE_dc = curvewise_edge_energy(ctrl_points, edge_img, varargin)
% integrates displacement vector for the snake by integrating over all
% curve points

[deriv_curve_ctrl, curve_points, t] = pchip_dxdc(ctrl_points, varargin{:});

curve_img_intensity = integrate_over_edge(edge_img, curve_points(:,1), curve_points(:,2) );
dr_dt = interp1((1.5:1:numel(t))', diff(curve_points),t,'linear','extrap');

for kk = size(ctrl_points, 1):-1:1
    for jj = 2:-1:1
        dE_dc(kk,jj) = sum( curve_img_intensity .*  deriv_curve_ctrl(:,jj,kk) .*  dr_dt(:, jj) );
    end
end

end