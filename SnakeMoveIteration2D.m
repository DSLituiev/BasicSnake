function P = SnakeMoveIteration2D(B, P, Fext0, gamma, kappa, delta, ForceOnCurve, varargin)
% This function will calculate one iteration of contour Snake movement
%
% P=SnakeMoveIteration2D(S,P,Fext,gamma,kappa)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   Fext : External vector field (from image)
%   gamma : Time step
%   kappa : External (image) field weight
%   delta : Balloon Force weight
%
% outputs,
%   P : The (moved) contour points N x 2;
%
% Function is written by D.Kroon University of Twente (July 2010)

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext0,1));
P(:,2)=min(max(P(:,2),1),size(Fext0,2));

% Get image force on the contour points
if ForceOnCurve
    % open ends only for now
    Fext1 = kappa * curvewise_edge_energy( Fext0, P);
else
    Fext1 = kappa * pointwise_edge_energy( Fext0, P);
end

Closed = false;
% if Closed
% Calculate the baloonforce on the contour points
N = GetContourNormals2D(P, Closed);
Fext2 = delta*N; % zeros(size(N));%
% else
%     Fext2 = zeros(size(Fext1));%
% end
% Update contour positions
ssx = gamma*P(:,1) + Fext1(:,1) + Fext2(:,1);
ssy = gamma*P(:,2) + Fext1(:,2) + Fext2(:,2);
if ~isempty(varargin) && ~isempty(varargin{1})
    fixed = varargin{1};
    if size(fixed,2) == 1
        P_fix = [P(fixed,1), P(fixed,1)];
    elseif size(fixed,2) == 2
        P_fix = P(fixed);
    else
        error('wrong array of fixed points')
    end
end

P(:,1) = B * ssx;
P(:,2) = B * ssy;

if ~isempty(varargin) && ~isempty(varargin{1})
    if size(fixed,2) == 1
        P(fixed, 1) = P_fix(:,1);
        P(fixed, 2) = P_fix(:,2);
    elseif size(fixed,2) == 2
        P(fixed) = P_fix;
    else
        error('wrong array of fixed points')
    end
end

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext0,1));
P(:,2)=min(max(P(:,2),1),size(Fext0,2));


