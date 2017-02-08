function DIRT_optPre_PrintSummary(optPre)
% DIRT_optPre_PrintSummary(optPre)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

fprintf('-------------------------------------------------\n');
fprintf('Precomputation options (optPre structure) summary\n');
fprintf('\n');
fprintf('color:                  %s\n', ABOnOff(optPre.color));
fprintf('\n');
fprintf('goptim:                 %s\n', ABOnOff(optPre.goptim));
fprintf('gmodel:                 %s\n', optPre.gmodel);
fprintf('ng:                     %02d\n', optPre.ng);
fprintf('gnorm:                  %s\n', ABOnOff(optPre.gnorm));
fprintf('\n');
fprintf('poptim:                 %s\n', ABOnOff(optPre.poptim));
fprintf('pmodel:                 %s\n', optPre.pmodel);
fprintf('np:                     %02d\n', optPre.np);
fprintf('pnorm:                  %s\n', ABOnOff(optPre.pnorm));
fprintf('\n');
fprintf('verbose level:          %01d\n', optPre.verb);
fprintf('filter_image_der_sigma: %04.1f\n', optPre.filter_image_der_sigma);
fprintf('debug:                  %s\n', ABOnOff(optPre.debug));
fprintf('debug2:                 %s\n', ABOnOff(optPre.debug2));
fprintf('-------------------------------------------------\n');
