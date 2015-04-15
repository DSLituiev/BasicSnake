function [data_interp, t, t0] = interp_implicit_pchip(data)
seglen = sqrt(sum(diff(data,1,1).^2,2));
t0 = [0;cumsum(seglen)];
t = (1:round(sum(seglen)))';

[N, data_dim] = size(data);
data_interp = zeros(numel(t), data_dim);

for jj = 1:data_dim
    pp_ = pchip(t0, data(:, jj));
    data_interp(:, jj) = ppval(pp_, t);
end
end