function out = dIKS1(x)
    % ------------------------------------------------------
    % Implmentation of Jacobian of Ex. 7 on page 487 of 
    %
    %   Izmailov A.F., Kurennoy A.S., Solodov M.V., Critical 
    %   solutions of nonlinear equations: stability issues, 
    %   Math. Program. Ser. B, 2018.
    % ------------------------------------------------------
    out1 = [x(2) x(1) 0];
    out2 = [x(3) 0 x(1)];
    out3 = [0 x(3) x(2)];
    out = [out1;out2;out3];
end
