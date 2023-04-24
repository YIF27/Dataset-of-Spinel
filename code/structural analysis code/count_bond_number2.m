function [count] = count_bond_number2(matrix, d1, d2)
%matrix is (n*6) double, matrix(:, 6) are the angles
%count the number of angle in the bound of (minDegree, maxDegree)
%return the count
count = 0;
    for i=1:length(matrix(:, 6))
         if  abs(matrix(i, 6)) > 0
            if abs(matrix(i, 4) - d1) < 0.01
                if abs(matrix(i, 5) - d2) < 0.01
                    count = count+1;             
                 end
            end
        end
    end
end
