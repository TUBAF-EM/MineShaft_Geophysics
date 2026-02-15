function DIKE2 = dike2(DI,RHO1,RHO2,THICK,XC,YC,XP,YP)
% DIKE2 = DIKE2(DI,RHO1,RHO2,THICK,XC,YC,XP,YP)
%
% resistivity soundings of a current and a potential electrode
% over a vertical dike structure.
%
% *****  INPUT  *****
% DI.....current
% RHO1...resistivity of left vertical layer
% RHO2...resistivity of the middle vertical layer, the dike.
% THICK..thickness of the middle layer, dike.
% XC.....x-coord. of current electrode, parallel to strike direction.
% YC.....y-coord. of current electrode, perpendicular to strike direct.
% XP.....x-coord. of potential electrode, parallel to strike direction.
% YP.....y-coord. of potential electrode, perpendicular to strike direct.
%
% ***** OUTPUT  *****
% DIKE2....voltage measured at the potential electrode
%
% ***** COORDINATE SYSTEM  *****
%             y=0          y=thick        --->y-axis
% -------------+------------+-----------------z=0
%     rho1     I    rho2    I    rho1         x-axis = strike direction
%
% *****  REMARKS  *****
% the formulas are derived analog to the formulas of Telford et al. (Applied
% Geophysics p.683) and Keller & Frischknecht (Electromagnetic Methods in
% Geophysical Prospecting p.181). However the formulas are not the same,
% the formulas in this program are more general. The summation is done with
% the Horner schema. The survey configuration is a combination of this
% current-potential configuration. For further questions ask T. HANSTEIN.

MAX = 100;
PI2 = 2*pi;
X = XP - XC;
X2= X*X;
THICK2 = THICK + THICK;
AK = (RHO2 - RHO1) / (RHO2 + RHO1);
AK2= AK*AK;
CC = 1. + AK;
DD = 1. - AK2;
AA = -AK*DD;
BB = -AK*CC;
if (XC == XP) & (YC == YP)
   DIKE2 = 0.;
   return
elseif YC <= 0.
   S = abs(YC);
   S2 = S + S;
   if YP <= 0.
      A =  YP - YC;
      Y =  (MAX + 1)*THICK2 + S2 - A;
      SUM = 1./sqrt( X2 + Y*Y);
      for M = MAX:-1:1
         Y = M*THICK2 + S2 - A;
         SUM = 1./sqrt( X2 + Y*Y) + AK2*SUM;
      end
      V = 1./sqrt(X2+A*A) + AK/sqrt(X2 + (S2-A).^2) + AA*SUM;
      DIKE2 = DI*RHO1/PI2*V;
      return
   elseif YP <= THICK
      A = abs(YP) + abs(YC);
      Y = (MAX + 1)*THICK2 + S2 - A;
      SUM1 = 1./sqrt(X2 + Y*Y);
      Y = MAX*THICK2 + A;
      SUM2 = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = M*THICK2 + S2 - A;
         SUM1 = 1./sqrt(X2 + Y*Y) + AK2*SUM1;
         Y = (M-1)*THICK2 + A;
         SUM2 = 1./sqrt(X2 + Y*Y) + AK2*SUM2;
      end
      V = BB*SUM1 + CC*SUM2;
      DIKE2 = DI*RHO1/PI2*V;
      return
   else
      A = abs(YP) + abs(YC);
      Y = MAX*THICK2 + A;
      SUM = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = (M-1)*THICK2 + A;
         SUM = 1./sqrt(X2 + Y*Y) + AK2*SUM;
      end
      V = DD*SUM;
      DIKE2 = DI*RHO1/PI2*V;
      return
   end
elseif (YC > 0.) & (YC <= THICK)
   S = YC;
   S2 = S + S;
   if YP <= 0.
      A = abs(YP) + YC;
      Y = MAX*THICK2 + A;
      SUM1 = 1./sqrt(X2 + Y*Y);
      Y = (MAX + 1)*THICK2 - S2 + A;
      SUM2 = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = (M - 1)*THICK2 + A;
         SUM1 = 1./sqrt(X2 + Y*Y) + AK2*SUM1;
         Y = M*THICK2 - S2 + A;
         SUM2 = 1./sqrt(X2 + Y*Y) + AK2*SUM2;
      end
      V = (1. + AK)*( SUM1 - AK*SUM2);
      DIKE2 = DI*RHO1/PI2*V;
      return
   elseif YP <= THICK
      A = YP - YC;
      Y = (MAX + 1)*THICK2 - A;
      SUM1 = 1./sqrt(X2 + Y*Y);
      Y = (MAX + 1)*THICK2 + A;
      SUM2 = 1./sqrt(X2 + Y*Y);
      Y = MAX*THICK2 + S2 + A;
      SUM3 = 1./sqrt(X2 + Y*Y);
      Y = (MAX + 1)*THICK2 - S2 - A;
      SUM4 = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = M*THICK2 - A;
         SUM1 = 1./sqrt(X2 + Y*Y) + AK2*SUM1;
         Y = M*THICK2 + A;
         SUM2 = 1./sqrt(X2 + Y*Y) + AK2*SUM2;
         Y = (M-1)*THICK2 + S2 + A;
         SUM3 = 1./sqrt(X2 + Y*Y) + AK2*SUM3;
         Y = M*THICK2 - S2 - A;
         SUM4 = 1./sqrt(X2 + Y*Y) + AK2*SUM4;
      end
      V = 1./sqrt(X2+A*A) + AK2*(SUM1+SUM2) - AK*(SUM3+SUM4);
      DIKE2 = DI*RHO2/PI2*V;
      return
   else
      A = YP - YC;
      Y = MAX*THICK2 + A;
      SUM1 = 1./sqrt(X2 + Y*Y);
      Y = MAX*THICK2 + S2 + A;
      SUM2 = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = (M-1)*THICK2 + A;
         SUM1 = 1./sqrt(X2 + Y*Y) + AK2*SUM1;
         Y = (M-1)*THICK2 + S2 + A;
         SUM2 = 1./sqrt(X2 + Y*Y) + AK2*SUM2;
      end
      V = (1. + AK)*( SUM1 - AK*SUM2);
      DIKE2 = DI*RHO1/PI2*V;
      return
   end
elseif YC > THICK
   S = YC - THICK;
   S2 = S + S;
   if YP <= 0.
      A = abs(YP) + YC;
      Y = MAX*THICK2 + A;
      SUM = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = (M-1)*THICK2 + A;
         SUM = 1./sqrt(X2 + Y*Y) + AK2*SUM;
      end
      V = DD*SUM;
      DIKE2 = DI*RHO1/PI2*V;
      return
   elseif YP <= THICK
      A = abs(YP - YC);
      Y = MAX*THICK2 + A;
      SUM1 = 1./sqrt(X2 + Y*Y);
      Y = (MAX + 1)*THICK2 + S2 - A;
      SUM2 = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = (M-1)*THICK2 + A;
         SUM1 = 1./sqrt(X2 + Y*Y) + AK2*SUM1;
         Y = M*THICK2 + S2 - A;
         SUM2 = 1./sqrt(X2 + Y*Y) + AK2*SUM2;
      end
      V = CC*( SUM1 - AK*SUM2);
      DIKE2 = DI*RHO1/PI2*V;
      return
   else
      A = YP - YC;
      Y = (MAX + 1)*THICK2 + S2 + A;
      SUM = 1./sqrt(X2 + Y*Y);
      for M = MAX:-1:1
         Y = M*THICK2 + S2 + A;
         SUM = 1./sqrt(X2 + Y*Y) + AK2*SUM;
      end
      V = 1./sqrt(X2 + A*A) + AK/sqrt(X2 + (S2+A).^2) + AA*SUM;
      DIKE2 = DI*RHO1/PI2*V;
      return
   end
end
