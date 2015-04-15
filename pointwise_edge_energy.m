function [ Fext1 ] = pointwise_edge_energy( Fext, P )
% returns displacement vector for the snake control points

Fext1(:,1) = interp2(Fext(:,:,1), P(:,2), P(:,1));
Fext1(:,2) = interp2(Fext(:,:,2), P(:,2), P(:,1));

% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;
end

