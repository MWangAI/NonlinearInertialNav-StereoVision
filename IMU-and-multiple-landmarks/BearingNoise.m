
function Xmes=BearingNoise(X,Xval)

% n = Xval*(rand(2,1) - 0.5);




n(1) = Xval*(rand(1,1) - 0.5);
n(2) = Xval*(rand(1,1) - 0.5);

Xnew = [X(1)/X(3)+ n(1), X(2)/X(3)+n(2) 1]';

Xmes = sign(X(3))*Xnew/norm(Xnew);

end
