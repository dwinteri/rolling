function ha = draw_ell_axis(E, varargin)

% draw_ell_axis(E) draws the axis of an ellipse E created using
% SPA.get_ellipse.
%
% See also : SPA.get_ellipse, SPA.draw_ellipse, SPA.get_ell_orientation


%==== Display major axis ================================================
k = linspace(-E.l, E.l, 100);
xa = E.x + k*cos(E.theta);
ya = E.y + k*sin(E.theta);


ha = plot(xa,ya, varargin{:});