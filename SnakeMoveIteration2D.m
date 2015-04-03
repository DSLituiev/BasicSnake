function P=SnakeMoveIteration2D(B, P,Fext,gamma,kappa,delta, varargin)
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
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

% Get image force on the contour points
Fext1(:,1)=kappa*interp2(Fext(:,:,1), P(:,2), P(:,1));
Fext1(:,2)=kappa*interp2(Fext(:,:,2), P(:,2), P(:,1));
% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;

% Calculate the baloonforce on the contour points
N=GetContourNormals2D(P);
Fext2=delta*N;

% Update contour positions
ssx = gamma*P(:,1) + Fext1(:,1) + Fext2(:,1);
ssy = gamma*P(:,2) + Fext1(:,2) + Fext2(:,2);
% if ~isempty(varargin) && ~isempty(varargin{1})
%    fixed = varargin{1};
%    if size(fixed,2) == 1       
%        ssx(~fixed) =  gamma*P(fixed,1);
%        ssy(~fixed) =  gamma*P(fixed,2);
%    elseif size(fixed,2) == 2
%        ssx(fixed(:,1)) =  gamma*P(fixed(:,1),1);
%        ssy(fixed(:,2)) =  gamma*P(fixed(:,2),2);
%    else
%        error('wrong array of fixed points')
%    end
% end

P(:,1) = B * ssx;
P(:,2) = B * ssy;
% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

    
