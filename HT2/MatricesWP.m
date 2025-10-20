function [W,P] = MatricesWP(P1,Lista_P)
%Return the weigh matrices (W) and the base elements (P),
%Input, P1, point to be approximated, and
%P: list of points or nodes. P1 = [x;y] y P = [x1 x2 x3;y1 y2 y3];
    B = base(P1);
    %e = 0.005;
    % dis_puntos = norm(P-P1);
    % lambda = 1./dis_puntos + e;

    m = length(Lista_P);
    n = length(B);

    %Creation M and P matrices
    W = sparse(m);
    P = zeros(m,n);

    for iter=1:m
        dis_puntos = norm(Lista_P(:,iter)-P1,2);

        % lambda = 1/(dis_puntos + e);
        lambda = 1/dis_puntos;
        W(iter,iter)=sqrt(lambda);
        P(iter,:) = base(Lista_P(:,iter));
    end 
    % M = W*P;
 end

function B = base(P1)
    x = P1(1);
    y = P1(2);
    B = [1 x y x*y x*x y*y];
end
