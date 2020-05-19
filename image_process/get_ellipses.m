function Out = get_ellipses(Img, head_th)

Out = struct();
I = 1:numel(Img);

% --- Find the "head" of the object
[Yh,Xh] = find(Img > head_th);
Out.xh = mean(Xh); % X-Position of the head's center of mass
Out.yh = mean(Yh); % Y-Position of the head's center of mass

% --- Initial Ellipse

Out.E = get_ellipse(Img, I);
X = Out.E.x;
Y = Out.E.y;

% --- Define the pile
P = struct('I', I, 'X', X, 'Y', Y, 'S', Out.E);         % The pile
Q = struct('I', {}, 'X', {}, 'Y', {}, 'S', {});         % The result

for k = 1:numel(P)
    
    [i,j] = ind2sub(size(Img), P(k).I);
    tmp = cos(P(k).S.theta)*(j-P(k).S.x) + sin(P(k).S.theta)*(i-P(k).S.y);
    
    Q(end+1).I = P(k).I(tmp>=0);
    [Q(end).X, Q(end).Y] = ind2sub(size(Img), Q(end).I);
    Q(end).S = get_ellipse(Img, Q(end).I);
    Out.E1 = Q(end).S;
    Out.I1 = Img & reshape(tmp, size(Img))>=0;
    
    Q(end+1).I = P(k).I(tmp<0);
    [Q(end).X, Q(end).Y] = ind2sub(size(Img), Q(end).I);
    Q(end).S = get_ellipse(Img, Q(end).I);
    Out.E2 = Q(end).S;
    Out.I2 = Img & reshape(tmp, size(Img))<0;
    
end
 
% --- Radius of curvature
Out.x0 = (Q(1).S.y-Q(2).S.y - Q(2).S.x/tan(Q(2).S.theta) + Q(1).S.x/tan(Q(1).S.theta))/(1/tan(Q(1).S.theta) - 1/tan(Q(2).S.theta));
Out.y0 = Q(1).S.y + (Q(1).S.x-Out.x0)/tan(Q(1).S.theta);

R = mean(sqrt((X-Out.x0).^2 + (Y-Out.y0).^2));

% --- Sign of curvature
istrigo = (Out.E.x-Out.x0)*(Out.yh-Out.y0) - (Out.E.y-Out.y0)*(Out.xh-Out.x0)<0;
R = R*(2*istrigo-1);

Out.curv = 1/R;

% --- Define E1 as the head ellipse

if (Out.E1.x-Out.xh)^2 + (Out.E1.y-Out.yh)^2 > (Out.E2.x-Out.xh)^2 + (Out.E2.y-Out.yh)^2
    tmp = Out.E1;
    Out.E1 = Out.E2;
    Out.E2 = tmp;
    
    tmp = Out.I1;
    Out.I1 = Out.I2;
    Out.I2 = tmp;
end






