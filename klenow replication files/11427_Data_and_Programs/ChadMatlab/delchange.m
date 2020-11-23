% delchange.m formerly delta.m    delta(x):  Computes the first difference of x.

function x=delchange(x);
[T N] = size(x);
x=x(2:T,:) - x(1:T-1,:);
