function [ Bool ] = Isregional( NumPR )
% Check is the NumPR correspond to a regional value

Bool=0;

if floor(NumPR/100)-10*floor(NumPR/1000)==9;
    Bool=1;
end

end

