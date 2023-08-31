function r = F(x,alpha,a)
    r = sqrt(max( 1-alpha^2*(x-a)^2 , 0 ));
end