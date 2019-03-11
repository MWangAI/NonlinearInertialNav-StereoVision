function S_x = Skew(x)
%% Function for Skew symmetric matrix
    [row,col] = size(x);
    
    if ((row == 3)&&(col == 1))
        S_x = [0 -x(3) x(2);
               x(3) 0 -x(1);
              -x(2) x(1) 0];
    else
        disp('Dimentional error') 
    end
end