function [xPoints, yPoints, x, yc] = NACABuild(NACA,num)
% Initializing the chord
c = 1;

% Initializing the NACA values using the NACA code
m = 1/100 * str2double(NACA(1));
p = 1/10 * str2double(NACA(2));
t = 1 / 100 * str2double(NACA(3:4));

% Initializing our theta
theta = linspace(-pi,pi,num);

% Mapping x to theta
x = c/2*(1-cos(theta));

x = fliplr(x);

% Calculating the tangent y coordinate
yt = t * c / 0.2 * (.2969*sqrt(x/c) - 0.1260*(x/c)-0.3516*(x/c).^2+0.2843*(x/c).^3 - 0.1036*(x/c).^4);


% If condition that calculates the camber line if one exists
if (p~=0)
    % Calculating the camber line up to the max camber
    yc1 = m*x / p^2 .* (2*p - x/c) .* (x<p*c);

    % Calculating the camber line after the max camber
    yc2 = m*(c-x) / (1-p)^2 .* (1+x/c-2*p) .* (x>=p*c);
    
    % Combining the camber values
    yc = yc1+yc2;

    % Finding the the derivative of the camber before the max camber
    ydc1 = 2*m/p^2*(p-x/c) .* (x<p*c);

    % Finding the the derivative of the camber after the max camber
    ydc2 = 2*m / (1-p)^2*(p-x/c) .* (x>=p*c);

    % Finding the value for exi using the combination of the derivatives
    exi = atan(ydc1+ydc2);
else
    % Setting the camber to 0
    yc = zeros(1,length(x));

    % Finding the angle when the slope is 0
    exi = atan(zeros(1,length(x)));
end

% Calculating upper x values
xU = (x - yt.*sin(exi)) .* (theta>0);

% Calculating the lower x values
xL = (x+yt.*sin(exi)) .* (theta<=0);

% Calculating the upper y values
yU = (yc+yt.*cos(exi)) .* (theta>0);

% Calculating the lower y values
yL = (yc - yt.*cos(exi)) .* (theta<=0);


% Combining into a sing vector
xPoints = xL+xU;
yPoints = yU+yL;



end