function Xnew=Tangent_Flight(X,Best_agent,dim,FES)

if rand <=0.25   % high exploration (large tangent flight)
    X=X+tan(rand(1,dim)*pi); 

else   %low exploration (small tangent flight)
    if isequal (Best_agent,X)
        teta=rand(1,dim)*pi/2.5;
        step=0.5*sign(rand(1,dim)-0.5);
        %step=0.5*sign(rand(1,dim)-0.5)/log(1+FES);
        X=X+step.*tan(teta);
    else
        teta=rand(1,dim)*pi/2.5; % small tangent flight
        step=0.5*sign(rand(1,dim)-0.5)*norm(Best_agent-X);
        X=X+step.*tan(teta);

    end
end
Xnew=X;
% id=randi(dim);
% ind=find(rand(1,dim)<=1.5/dim);
% ind=[ind, id];
% Xnew(ind)=X(ind);
% 
% %-------------Check boundries
% Xnew(Xnew>ub)=rand*(ub - lb) + lb;
% Xnew(Xnew<lb)=rand*(ub - lb) + lb;

