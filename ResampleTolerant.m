function newvals = ResampleTolerant(vals,length1,length2)
% Wrapper around the resample function that allows it to work even if
% the product of the lengths are long enough to overwhelm the resample.m 
% limit of 2^31
% Works by finding a rational number approximation of the requested length
% ratio... to within a particular tolerance.
% INPUTS
% vals = vector of values to be resampled
% length1 = desired length
% length2 = initial length (often equals length(vals))
%
% Dan Levenstein code made into a function by Brendon Watson
% August 2016


if length1*length2 < 2^31 % if no need to change factors don't
    newvals = resample(vals,length1,length2);
else
    newvals = [1 1];
    
    resamplefact = length1/length2;
    tol = 0.0001;
    while length(newvals(:,1)) ~= length1
        [P,Q] = rat(resamplefact,tol);
        if P==0
            tol = tol/10;
            continue
        end
        if P*Q >=2^20 || tol<1e-300  %Avoid crashing resample...
            vals([1,end],:) = [];
            length2 = length2-2;
            resamplefact = length1/length2;
            tol = 0.0001;
            continue
        end
        newvals = resample(vals,P,Q);
        tol = tol/10;
    end
end

end
