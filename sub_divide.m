function [v_b, v_n ] = sub_divide(v, B) %given a column and row of indices B, returns the subdivision into those in B and not in B
   v_b = zeros(size(B,2),1);
   v_n = zeros(size(v,1) - size(B,2),1);
   temp = v;
   for i=1:size(B,2)
       v_b(i,1) = v(B(1,i)); 
       temp(B(1,i)) = inf; %mark those entries in b from the input vector for later
   end
   j = 1;
   pos = 1;
   while (pos <= size(v_n,1))
       if(~isinf(temp(j,1)))
           v_n(pos,1) = temp(j,1);
           pos = pos +1;
       end
       j = j +1;
   end
end

