function [ d ] = estimate_distance(result)
   
	d = power(10, result); %d=10^(PL)
    d = transpose(d);
end

