
function y = uexact(t,x)

% if x >= t
%   y = u0(x-t);
% else
%   y = bc(t-x);
% end
a = 2.0;
y=sin(pi*(x-a*t));
