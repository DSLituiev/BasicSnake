%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interate.m implements the core of the snakes (active contours) algorithm.
% It is called by the snk.m file which is the GUI frontend. If you do not
% want to deal with GUI, look at this file and the included comments which
% explain how active contours work.
%
% To run this code with GUI
%   1. Type guide on the matlab prompt.
%   2. Click on "Go to Existing GUI"
%   3. Select the snk.fig file in the same directory as this file
%   4. Click the green arrow at the top to launch the GUI
%
%   Once the GUI has been launched, you can use snakes by
%   1. Click on "New Image" and load an input image. Samples image are
%   provided.
%   2. Set the smoothing parameter "sigma" or leave it at its default value
%   and click "Filter". This will smooth the image.
%   3. As soon as you click "Filter", cross hairs would appear and using
%   them and left click of you mouse you can pick initial contour location
%   on the image. A red circle would appead everywhere you click and in
%   most cases you should click all the way around the object you want to
%   segement. The last point must be picked using a right-click in order to
%   stop matlab for asking for more points.
%   4. Set the various snake parameters (relative weights of energy terms
%   in the snake objective function) or leave them with their default value
%   and click "Iterate" button. The snake would appear and move as it
%   converges to its low energy state.
%
% Copyright (c) Ritwik Kumar, Harvard University 2010
%               www.seas.harvard.edu/~rkkumar
%
% This code implements Snakes: Active Contour Models by Kass, Witkin and
% Terzopolous incorporating Eline, Eedge and Eterm energy factors. See the
% included report and the paper to learn more.
%
% If you find this useful, also look at Radon-Like Features based
% segmentation in  the following paper:
% Ritwik Kumar, Amelio V. Reina & Hanspeter Pfister, Radon-Like Features 
% and their Application to Connectomics, IEEE Computer Society Workshop %
% on Mathematical Methods in Biomedical Image Analysis (MMBIA) 2010
% http://seas.harvard.edu/~rkkumar
% Its code is also available on MATLAB Central
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xs, ys] = snake_iterate(image, xs, ys, alpha, beta, gamma, kappa, wl, we, wt, iterations)
% image: This is the image data
% xs, ys: The initial snake coordinates
% alpha: Controls tension
% beta: Controls rigidity
% gamma: Step size
% kappa: Controls enegry term
% wl, we, wt: Weights for line, edge and terminal enegy components
% iterations: No. of iteration for which snake is to be moved


%parameters
N = iterations; 
% Calculating size of image
[row col] = size(image);


%Computing external forces

eline = image; %eline is simply the image intensities

[grady,gradx] = gradient(image);
eedge = -1 * sqrt ((gradx .* gradx + grady .* grady)); %eedge is measured by gradient in the image

%masks for taking various derivatives
m1 = [-1 1];
m2 = [-1;1];
m3 = [1 -2 1];
m4 = [1;-2;1];
m5 = [1 -1;-1 1];

cx = conv2(image,m1,'same');
cy = conv2(image,m2,'same');
cxx = conv2(image,m3,'same');
cyy = conv2(image,m4,'same');
cxy = conv2(image,m5,'same');

% for i = 1:row
%     for j= 1:col
%         % eterm as deined in Kass et al Snakes paper
%         eterm(i,j) = (cyy(i,j)*cx(i,j)*cx(i,j) -2 *cxy(i,j)*cx(i,j)*cy(i,j) + cxx(i,j)*cy(i,j)*cy(i,j))/((1+cx(i,j)*cx(i,j) + cy(i,j)*cy(i,j))^1.5);
%     end
% end

eterm = (cyy.*cx.*cx  -2 *cxy.*cx.*cy + cxx.*cy.*cy)./((cx.*cx + cy.*cy ).^1.5);

% imview(eterm);
% imview(abs(eedge));
eext = (wl*eline + we*eedge -wt * eterm)./(wl+we+wt); %eext as a weighted sum of eline, eedge and eterm

[fx, fy] = gradient(eext); %computing the gradient


%initializing the snake
xs=xs(:);
ys=ys(:);
nPoints = numel(xs);    
closed = false;
Ainv = SnakeInternalForceMatrix2D(nPoints, alpha, beta, gamma, closed);
 
%moving the snake in each iteration
    imshow(image,[]); 
    hold on;
    li = plot([xs; xs(1)], [ys; ys(1)], 'r.-');    
    pause(0.001)
    
for i=1:N;
    
    ssx = gamma*xs - kappa*interp2(fx, xs,ys);
    ssy = gamma*ys - kappa*interp2(fy, xs,ys);
    
    %calculating the new position of snake
    xs = Ainv * ssx;
    ys = Ainv * ssy;    
    
    %Displaying the snake in its new position    
    set(li, 'xdata', [xs; xs(1)], 'ydata', [ys; ys(1)]);
    
    pause(0.001)
    if any(isnan(xs)) || any(isnan(ys))
        error('NaN values, iteration %u \n', i)
    end
end;