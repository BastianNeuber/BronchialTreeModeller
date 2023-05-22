function [x] = APPFNC_ceil2(x)
x = ceil(x);
if x < 2
    x = 2;
end
end

