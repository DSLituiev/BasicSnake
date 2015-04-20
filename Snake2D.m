function [P,J]=Snake2D(I,P,Options)
% This function SNAKE implements the basic snake segmentation. A snake is an
% active (moving) contour, in which the points are attracted by edges and
% other boundaries. To keep the contour smooth, an membrame and thin plate
% energy is used as regularization.
%
% [O,J]=Snake2D(I,P,Options)
%
% inputs,
%   I : An Image of type double preferable ranged [0..1]
%   P : List with coordinates descriping the rough contour N x 2
%   Options : A struct with all snake options
%
% outputs,
%   O : List with coordinates of the final contour M x 2
%   J : Binary image with the segmented region
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.nPoints : Number of contour points, default 100
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
%  Options.Sigma1 : Sigma used to calculate image derivatives, default 10
%  Options.Wline : Attraction to lines, if negative to black lines otherwise white
%                    lines , default 0.04
%  Options.Wedge : Attraction to edges, default 2.0
%  Options.Wterm : Attraction to terminations of lines (end points) and
%                    corners, default 0.01
%  Options.Sigma2 : Sigma used to calculate the gradient of the edge energy
%                    image (which gives the image force), default 20
%
% options (Gradient Vector Flow)
%  Options.Mu : Trade of between real edge vectors, and noise vectors,
%                default 0.2. (Warning setting this to high >0.5 gives
%                an instable Vector Flow)
%  Options.GIterations : Number of GVF iterations, default 0
%  Options.Sigma3 : Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.2
%  Options.Delta : Baloon force, default 0.1
%  Options.Kappa : Weight of external image force, default 2
%
%
% Literature:
%   - Michael Kass, Andrew Witkin and Demetri TerzoPoulos "Snakes : Active
%       Contour Models", 1987
%   - Jim Ivins amd John Porrill, "Everything you always wanted to know
%       about snakes (but wer afraid to ask)
%   - Chenyang Xu and Jerry L. Prince, "Gradient Vector Flow: A New
%       external force for Snakes
%
% Example, Basic:
%
%  % Read an image
%   I = imread('testimage.png');
%  % Convert the image to double data type
%   I = im2double(I);
%  % Show the image and select some points with the mouse (at least 4)
%  %figure, imshow(I); [y,x] = getpts;
%   y=[182 233 251 205 169];
%   x=[163 166 207 248 210];
%  % Make an array with the clicked coordinates
%   P=[x(:) y(:)];
%  % Start Snake Process
%   Options=struct;
%   Options.Verbose=true;
%   Options.Iterations=300;
%   [O,J]=Snake2D(I,P,Options);
%  % Show the result
%   Irgb(:,:,1)=I;
%   Irgb(:,:,2)=I;
%   Irgb(:,:,3)=J;
%   figure, imshow(Irgb,[]);
%   hold on; plot([O(:,2);O(1,2)],[O(:,1);O(1,1)]);
%
% Example, GVF:
%   I=im2double(imread('testimage2.png'));
%   x=[96 51 98 202 272 280 182];
%   y=[63 147 242 262 211 97 59];
%   P=[x(:) y(:)];
%   Options=struct;
%   Options.Verbose=true;
%   Options.Iterations=400;
%   Options.Wedge=2;
%   Options.Wline=0;
%   Options.Wterm=0;
%   Options.Kappa=4;
%   Options.Sigma1=8;
%   Options.Sigma2=8;
%   Options.Alpha=0.1;
%   Options.Beta=0.1;
%   Options.Mu=0.2;
%   Options.Delta=-0.1;
%   Options.GIterations=600;
%   [O,J]=Snake2D(I,P,Options);
%
% Function is written by D.Kroon University of Twente (July 2010)

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',100, 'Wline', 0.04,...
    'Wedge', 2, 'Wterm',0.01,'Sigma1',10,'Sigma2', 20,...
    'Alpha', 0.2,'Beta',0.2,'Gamma',1, 'Delta',0.1, 'Kappa',2,...
    'Iterations', 100, 'GIterations', 0,...
    'Mu', 0.2, 'Sigma3', 1, 'Closed', true, 'AbsTol',1e-6, ...
    'Norm', 2, 'MaxStep', 50, 'useAsEnergy', false, ...
    'Fixed', [], 'forceActsUpon', 'points', 'linewidth', 1.2, 'figure', []);
if(~exist('Options','var')),
    Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
        if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))),
        warning('snake:unknownoption','unknown options found');
    end
end

Options.ForceOnCurve = strcmpi(Options.('forceActsUpon'), 'curve');

% The contour must always be clockwise (because of the balloon force)
P=MakeContourClockwise2D(P);

% Make an uniform sampled contour description
if Options.Closed && Options.nPoints ~= size(P,1)
    P=InterpolateContourPoints2D(P, Options.nPoints);
end

% Convert input to double
I = double(I);
% if ~Options.useAsEnergy
    % If color image convert to grayscale
    if(size(I,3)==3), I=rgb2gray(I); end
    
    [ Eext, Fext ] = gvf_energy_force( I, Options );
% else
%     Eext = I;
%     Fext(:,:,1) = Eext;
%     Fext(:,:,2) = Eext;
% end
% Eext = -(cumsum(Fext(:,:,1),1) + cumsum(Fext(:,:,2),2));
% Eext = (ImageDerivatives2D(I,Options.Sigma1,'xx') + ImageDerivatives2D(I,Options.Sigma1,'yy'));

% Fext(:,:,1) = -Fx;
% Fext(:,:,2) = -Fy;

% Fext(:,:,1) = cumsum(-Fext(:,:,1),1);
% Fext(:,:,2) = cumsum(-Fext(:,:,2),2);

%  Fext(:,:,1) = 2*Options.Sigma1.^2 * Eext;
%  Fext(:,:,2) = 2*Options.Sigma1.^2 * Eext;

% Show the image, contour and force field
if(Options.Verbose)
    if~isempty(Options.figure) && isfigure(Options.figure)
        figh = Options.figure;
    else
        figh = figure;
    end
    set(figh,'render','opengl')
    spl(1) = subplot(2,2,1);
    imagesc(Fext(:,:,1));
    hold on;
    title('x-component')
    spl(2) = subplot(2,2,2);
    imagesc(Fext(:,:,2));
    hold all
    title('y-component')
    Q = 0.025;
    set(spl(1:2),'clim', max(abs(quantile(Fext(:), [Q, 1-Q]))) *[-1,1] )
    spl(3) = subplot(2,2,3);
    spacing = 20;
    [x,y]=ndgrid(1:spacing:size(Fext,1),1:spacing:size(Fext,2));
    imagesc(I), hold on; quiver(y,x,...
        Fext(1:spacing:end,1:spacing:end,2),...
        Fext(1:spacing:end,1:spacing:end,1), 'w');
    title('The external force field ')
    spl(4) = subplot(2,2,4);
    title('Snake movement ')
    imagesc(Eext); hold on;
end

% Make the interal force matrix, which constrains the moving points to a
% smooth contour
if size(Options.Fixed, 2) == 1
    Options.Fixed = [Options.Fixed, Options.Fixed];
end

%% plotting
    function [x, y] = cp_ordered(P, Closed)
        if Closed
            x = [P(:,2);P(1,2)];
            y = [P(:,1);P(1,1)];
        else
            x = P(:,2);
            y = P(:,1);
        end
    end

axes(spl(4))
[x0, y0] = cp_ordered(P, Options.Closed);
[data_interp] = interp_implicit_pchip([x0, y0]);
% li = zeros(Options.Iterations,1);
scatter(x0, y0, 5, [0 0.8, 0],'.'); hold all
li(1) = plot(data_interp(:,1), data_interp(:,2), '-','Color',[0 0.8, 0]);
li(2) = plot(data_interp(:,1), data_interp(:,2), '-','Color', 'k', 'linewidth', Options.linewidth);
hold all;
title('Energy and the snake contour movement')

for ii = 1:4
    axes(spl(ii))
    plot(P(:,2),P(:,1),'.','Color',[0, 0.8, 0] );
    hold all;
    h(ii) = scatter(P(:,2),P(:,1),100*pi,'w','.');
end
axes(spl(4))


A_inv = SnakeInternalForceMatrix2D(Options.nPoints, Options.Alpha, Options.Beta, Options.Gamma, Options.Closed);

if Options.ForceOnCurve
    %     Fext_preint = zeros(size(Fext));
    %      Fext_preint(:,:,2) = cumsum(Fext(:,:,2), 2)/200; %Eext; %
    %      Fext_preint(:,:,1) = cumsum(Fext(:,:,1), 1)/200;
    ext_energy_iter_fun = @(x, y)SnakeMoveIteration2D(A_inv, x, Fext, ...
        Options.Gamma, y, Options.Delta, Options.ForceOnCurve, Options.Fixed);
else
    ext_energy_iter_fun = @(x, y)SnakeMoveIteration2D(A_inv, x, Fext, ...
        Options.Gamma, y, Options.Delta, Options.ForceOnCurve, Options.Fixed);
end


ii = 1;
while ii<Options.Iterations
    P_prev = P;
    Options.Kappa = Options.Kappa * ( 1 - 1e-5 );
    
    P = ext_energy_iter_fun(P_prev, Options.Kappa);
    if norm(P_prev - P, Options.Norm)/Options.nPoints < Options.AbsTol
        fprintf('converged in %u iterations!\n', ii)
        break
    end
    
    StepNorm = norm( sqrt(sum((P_prev - P).^2, 2)), Options.Norm);
    if StepNorm > Options.MaxStep
        warning('too big step: %u', round(norm(P_prev - P, Options.Norm)) )
        P = P_prev;
        Options.Kappa = Options.Kappa * (Options.MaxStep/StepNorm);
        continue
    else
        ii = ii+1;
    end
    
    % Show current contour
    if(Options.Verbose)
        c=i/Options.Iterations;
        [x0, y0] = cp_ordered(P, Options.Closed);
        [data_interp] = interp_implicit_pchip([x0, y0]);
        %         set(li, 'xdata',x_, 'ydata', y_, 'Color',[c 1-c 0])
        %         li(i) = line(data_interp(:,1), data_interp(:,2), 'LineStyle', '-','Color', [c, 0.2, 1-c], 'Parent', spl(4));
        set(li(2),  'xdata', data_interp(:,1), 'ydata', data_interp(:,2),...
            'linewidth', Options.linewidth, 'LineStyle', '-','Color', [c, c, c], 'Parent', spl(4));
        set( h , 'xdata', x0, 'ydata', y0, 'zdata', 2*ones(Options.nPoints,1));
        drawnow;
    end
end
fprintf('abs numeric error (%u norm): %f\n', Options.Norm, norm(P_prev - P, Options.Norm)/Options.nPoints)

if(nargout>1)
    J=DrawSegmentedArea2D(P,size(I));
end

end