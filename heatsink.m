% Clear 
clc; clear; close all; 

% Define Constants
fins = 3;
thickness = 0.5;
height = 1;
n = 4*fins;
length = ((2*fins)-1)*thickness;
vertices = [length/2-i*thickness/2; -length/2-i*thickness/2; zeros(n,1)];

% Find X Vertices
for k = 1:fins
    for j = 1:4
        if j == 1
            vertices(2+j+(4*(k-1))) = ((-length/2)+(2*thickness*(k-1))) + (1i*(thickness/2));
        elseif j == 2
            vertices(2+j+(4*(k-1))) = ((-length/2)+(2*thickness*(k-1))) + (1i*(height));
        elseif j == 3
            vertices(2+j+(4*(k-1))) = (((-length/2)+(2*thickness*(k-1)))+thickness) + (1i*(height));
        elseif j == 4 
            vertices(2+j+(4*(k-1))) = (((-length/2)+(2*thickness*(k-1)))+thickness) + (1i*(thickness/2));
        end
    end
end

% Flip Vertices
vertices = real(vertices)-(1i*imag(vertices));

% Plot Heat Sink
plot([real(vertices); real(vertices(1))], ([imag(vertices); imag(vertices(1))]), 'k', AffectAutoLimits='off')
xlim([-2,2])
ylim([-2,2])
print("Heat Sink", "-r300", '-dpng')

% Do mapping
figure
p = polygon(vertices);
f = hplmap(p);
plot(f);
phi = lapsolve(f);
print("SC Map", "-r300", '-dpng')
hold on 
