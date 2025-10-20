function [w] = Weights_wi(P1,P) 
    %Calculation of w_i weights 
    %P1 : Point the interesting
    %Point list (Nodes)

    %Matrix MP whit size mxn whit m > n.
    [W, P] = MatricesWP(P1,P);

    WP = W*P;
    %size MP 
    [m, n] = size(WP);

    %QR Factorization
    [Q, R] = qr(WP);

    %Size Q1 mxn
    Q1 = Q(1:m,1:n);

    %Size R nxn
    R1 = R(1:n,1:n);

    %Matrix phi nxm
    Phi_1 = base(P1);
    Phi = Phi_1';

    %Solve system d = R^{-t}Phi
    % d = inv(R1')*Phi;
    d = (R1')\Phi;

    %Calculo de los Pesos w
    w = W*(Q1*d);

    
    function B = base(P1)
        %Bases
        x = P1(1);
        y = P1(2);
        B = [1 x y x*y x*x y*y];
    end
end


