function [data_interp, t, t0, data] = interp_implicit_pchip(data)

% %% remove points that are too close 
% seglen = sqrt(sum(diff(data,1,1).^2,2));
% data = data( [seglen>1; true], :);
%%
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