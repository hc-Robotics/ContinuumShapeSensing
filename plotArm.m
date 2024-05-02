function Result = plotArm(th0,th1,th2,L)
% calculate the shape of continuum robot in cartesian coordinate

s=0:2e-2:1;

x = zeros(size(s,2),1)';
z = zeros(size(s,2),1)';

for i = 1:size(s,2)
    x(i) = L*integral(@(x)sin(th0.*x+0.5*th1.*x.^2+(1/3)*th2.*x.^3),0,s(i));
    z(i) = L*integral(@(x)cos(th0.*x+0.5*th1.*x.^2+(1/3)*th2.*x.^3),0,s(i));
end
Result = [x;z];
end
