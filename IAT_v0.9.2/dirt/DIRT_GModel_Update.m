function g = DIRT_GModel_Update(gmodel, g, dg, gopt)
% g = DIRT_GModel_Update(gmodel, g, dg, gopt)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

switch gmodel
    case 'Homography'
        if 1
            % use full inverse of the increment
            dH = eye(3) + [ dg(1) dg(2) dg(3) ; dg(4) dg(5) dg(6) ; dg(7) dg(8) 0 ];
            g = g * inv(dH);
        else
            % use first order approximation of the inverse for the
            % increment
            g = g *( eye(3) - [ dg(1) dg(2) dg(3) ; dg(4) dg(5) dg(6) ; dg(7) dg(8) 0 ] );
        end;
        g = g / norm(g, 'fro');
    case '3PtHomography'
        dH = [ eye(2) gopt.sq ; 0 0 1 ] * [ 1+dg(1) dg(2) 0 ; dg(3) 1+dg(4) 0 ; dg(5) dg(6) 1 ] * [ eye(2) -gopt.sq ; 0 0 1 ];
        g = g * inv(dH);
        g = g / norm(g, 'fro');
    case 'Affine'
        dA = [ 1+dg(1) dg(2) dg(3) ; dg(4) 1+dg(5) dg(6) ; 0 0 1 ];
        g = g * inv(dA);
    case 'Rt'
        co = cos(dg(1));
        si = sin(dg(1));        
        dG = [ co si dg(2) ; -si co dg(3) ; 0 0 1 ];
        g = g * inv(dG);
    otherwise
        error(['unknown geometric model ' gmodel]);
end;
