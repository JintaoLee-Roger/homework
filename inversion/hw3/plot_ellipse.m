function plot_ellipse(DELTA2,C,m)
% set the number of points on the ellipse to generate and plot
% Plot a two-dimensional covariance 
% ellipse about the model parameters, 
% 
% plot_ellipse(DELTA2,C,m)
%       C: the covariance matrix, 
%       DELTA2: $\Delta^2$
%       m: the 2-vector of model parameters.

n = 100;
% construct a vector of n equally-spaced angles from (0,2*pi)
theta = linspace(0, 2 * pi, n)';
% corresponding unit vector
xhat = [cos(theta), sin(theta)];
Cinv = inv(C);
% preallocate output array
r = zeros(n, 2);
for i = 1:n
    % store each (x,y) pair on the confidence ellipse
    % in the corresponding row of r
    r(i, :) = sqrt(DELTA2 / (xhat(i,:) * Cinv * xhat(i,:)')) * xhat(i,:);
end
plot(m(1) + r(:,1), m(2) + r(:,2));
axis equal

