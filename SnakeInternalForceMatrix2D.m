function Ainv = SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma, circ)
%
% B=SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma)
%
% inputs,
%   nPoints : The number of snake contour points
%   alpha : membrame energy  (first order)
%   beta : thin plate energy (second order)
%   gamma : Step Size (Time)
%
% outputs,
%   B : The Snake Smoothness regulation matrix
%
% Function is written by D.Kroon University of Twente (July 2010)

% Penta diagonal matrix, one row:
alphas_ = [ 0,    -alpha, 2*alpha, -alpha, 0    ]';
betas_  = [ beta, 4*beta, 6*beta,  4*beta, beta ]';

b = alphas_ + betas_;

eye_ = eye(nPoints);

for ii = 5:-1:1
    A0(:,:,ii) = circshift(eye_,3-ii);
end

A_alpha = A0;
if ~circ
    A_alpha(1:2,:,1:2) = 0;
    A_alpha(1:2,1,1:2) = 1;
    
    A_alpha(end-1:end,:,end-1:end) = 0;
    A_alpha(end-1:end,end,end-1:end) = 1;
    
%     A_alpha(1,:,:) = 0;
%     A_alpha(end,:,:) = 0;
     A_beta = A_alpha;
%     A_beta(2,:,:) = 0;
%     A_beta(end-1,:,:) = 0;
end

A_ = bsxfun(@times, A_alpha, permute( alphas_(:),[2,3,1]) ) + ...
    bsxfun(@times, A_beta, permute(betas_(:),[2,3,1]) ) ;
% Make the penta matrix (for every contour point)
A = sum(A_,3);

% Calculate the inverse
[L, U] = lu(A + gamma .* eye(nPoints,nPoints));
Ainv = inv(U) * inv(L);
