function E = get_ellipse(Img, I)
%SPA.GET_ELLIPSE Get the equivalent ellipse of a region in a image
%   E = GET_ELLIPSE(IMG, I) find the ellipse equivalent to the set of pixels I
%       in the image IMG. IMG is a n-by-m image and I is a vector of pixel 
%       index.
%
%   See also: SPA.draw_ellipse, SPA.draw_ell_axis, SPA.get_ell_orientation

if ~exist('I', 'var')
    I = 1:numel(Img);
end

E = struct();

% Compute the moments
E.m00 = moment(Img,I,0,0);
E.m10 = moment(Img,I,1,0);
E.m01 = moment(Img,I,0,1);
E.m11 = moment(Img,I,1,1);
E.m02 = moment(Img,I,0,2);
E.m20 = moment(Img,I,2,0);

% Compute the ellipse properties
E.x = E.m10/E.m00;
E.y = E.m01/E.m00;

a = E.m20/E.m00 - E.x^2;
b = 2*(E.m11/E.m00 - E.x*E.y);
c = E.m02/E.m00 - E.y^2;

E.theta = 1/2*atan(b/(a-c)) + (a<c)*pi/2; %%ATAN   Inverse tangent, result in radians.
E.w = sqrt(6*(a+c-sqrt(b^2+(a-c)^2)))/2;
E.l = sqrt(6*(a+c+sqrt(b^2+(a-c)^2)))/2;



% -------------------------------------------------------------------------
function m = moment(Img, I, p, q)

[j, i] = ind2sub(size(Img),I);

m = sum((i.^p).*(j.^q).*double(Img(I)));

