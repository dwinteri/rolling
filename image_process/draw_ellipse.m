function h = draw_ellipse(E,color)
% SPA.DRAW_ELLIPSE draws an ellipse whose parameters are given by the
%   struct E.
%
%   See also: SPA.get_ellipse, SPA.draw_ell_axis, SPA.get_ell_orientation

t = linspace(0, 2*pi, 200);
x = E.x + E.w*cos(t)*sin(E.theta)+E.l*sin(t)*cos(E.theta);
y = E.y + E.l*sin(t)*sin(E.theta)-E.w*cos(t)*cos(E.theta);


h = plot(x,y, 'color', color);
