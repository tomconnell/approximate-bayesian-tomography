function g=gpoly(x0,z0,xcorn,zcorn,ncorn,rho)

% g = gpoly(x0,z0,xcorn,zcorn,ncorn,rho)
% This function computes the vertical attraction of a 2-D body with polygonal
% cross section. Axes are right-handed system with x-axis parallel to long
% direction of body and z-axis vertical, positive downward.
% Based on Fortran subroutine B7 by: 

% - Blakely (1995), Potential Theory in Gravity and Magnetic Applications, C.U.P., p 378

% x0 and z0 are the coordinates of the observation point. 
% Arrays xcorn and zcorn (each of length ncorn) contain the coordinates of the
% polygon corners, arranged clockwise when viewed with x-axis to right. 
% rho is the bulk (constant) density of body.
% all distance parameters in units of km, rho in units of kg/m^3

% Output = vertical attraction of gravity g, in mgal

% The formulas are a little bit different from the original Talwani's
% equations, but give the same results.

% This can be easily modified to compute gradients and other components of g
% J.C. Afonso
% Masters Course GEOS701
% Macquarie University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=x0./1000;
z0=z0./1000;
xcorn=xcorn./1000;
zcorn=zcorn./1000;
%% Key parameters
gamma=6.673e-11; % graviatational constant
si2mg=1e5;       % factor to convert SI units to CGM units 
km2m=1e3;        % factor to convert km to m 

%% main loop
sum=0;
for n=1:ncorn % loop over number of corners (+ 1 to close polygon)
    if n==ncorn % I'm closing the polygon... last step
        n2=1;
    else
        n2=n+1;
    end
    
    x1=xcorn(n)-x0;  % x-distance between observation point and corner n
    z1=zcorn(n)-z0;  % z-distance between observation point and corner n
    x2=xcorn(n2)-x0; % x-distance between observation point and corner n+1
    z2=zcorn(n2)-z0; % z-distance between observation point and corner n+1 
    r1sq=x1*x1+z1*z1; % Euclidean distance from observation point to n
    r2sq=x2*x2+z2*z2; % Euclidean distance from observation point to n+1
    
    if r1sq==0       % if corner is = observation point => move to next corner
        break
        disp('Field point on corner')
    elseif r2sq==0
        break
        disp('Field point on corner')
    end
    
    denom=z2-z1;
    
    if denom==0     % denominator in alpha below; don't allow zero value
        denom=1E-15; % arbitrary small number (vertical slope of line)
    end
    alpha=(x2-x1)/denom; 
    beta=x1 - alpha*z1;
%    beta=(x1*z2-x2*z1)/denom;            % alternative form for beta
    factor=beta/(1+alpha*alpha);          % factor in eq. 9.11
    term1=0.5*(log(r2sq)-log(r1sq));      % 0.5 is there due to the missing sqrt in r terms
    term2=atan2(z2,x2)-atan2(z1,x1);      % "arctang" factor page 194
    sum=sum+factor*(term1-alpha*term2);   % Final solution, Eq. 9.11
end

g=2*rho*gamma*sum*si2mg*km2m;             % put the solution into variable g

