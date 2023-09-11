function [Qt,Qfrac,dTime]=GMD_TOPMODEL_ode(D,WxmD,WxmU,WbmD,WbmU,obsR,obsQ,...
                        params,cs,DT0,t,SINa,SINb,COSa,COSb,area,AREA,Nr,Nc,cW)              
warning('off','all');
%-------------------------------------------------------------------------- 
%number of observed data points [-]
Nobs                 = size(obsR,1);
%observed rainfall resolution [s]
dTime                = t(2)-t(1);
%--------------------------------------------------------------------------
%load uncertain/input model parameters
[d,Tmax,ep,Srzmax,manNhs,manNch] = unPack_uncertain_parameters(params);
%initialise system
[Sx0,Su0,wtd0]       = initialiseSYS(Nc);
%assemble manning's n coefficient for hillslope vs channel classes
manN                 = 0*Sx0 + manNhs;           %hillslope
manN(Nc-Nr+1:Nc)     = manNch;                   %channels
%assemble max root zone storage (channel doesn't have rootzone)
Smax                 = 0*Sx0 + Srzmax;    
%evapo-transpiration rate [m/s]
ET                   = seasonal_sinewave_evap(DT0,1,dTime,(1:length(obsR))');
%--------------------------------------------------------------------------
%convert total rainfall (e.g., from tipping bucket) to a rate [m/s]
drdt                 = obsR./dTime;
%method of interpolation
METHOD               = 'pchip';
%find tippin bucket rainfall times (midpoints)
tR                   = [0;t];
tR                   = (tR(1:end-1)+tR(2:end))./2;
%gridded interpolation
FDR                  = griddedInterpolant(tR,drdt,METHOD);
FET                  = griddedInterpolant(t,ET  ,METHOD);
%assemple the vector of state variables
V0                   = [Sx0;Su0;wtd0];
%--------------------------------------------------------------------------
%jacobian pattern (tells the ode-solver to evaluate only where JPAT=1 
%to save runtime - see model paper for details)
M1                   = full(diag(ones(Nc,1)));
M2                   = sign(WxmD + D' + WbmD);
M2(M1==1)            = 1;
JPAT                 = [M2 M1 M1; M1 M1 M1; M1 M1 M2];        
%--------------------------------------------------------------------------
%ode solver tolerance
r0                   = drdt;
r0(r0==0)            = [];
abstol               = min(min(r0),1e-6);
%ode-solver Options
OPS                  = odeset('JPattern',JPAT,'InitialStep',1e-8,'maxstep',10*dTime,'abstol',abstol);
%solve using ode15s, suitable for "stiff" system if equations
tic;
[~,V]                = ode15s(@HydroGEM_ode_fun,t,V0,OPS,area,d,Nc,Smax,FET,FDR,manN,D,WxmD...
                              ,WbmD,WxmU,WbmU,cs,SINa,SINb,COSa,COSb,Tmax,Nr,cW,ep);
%--------------------------------------------------------------------------
%in case ode-solver has crashed midway
if size(V,1)<Nobs 
    V                = nan(Nobs,3*Nc);
end
%--------------------------------------------------------------------------
%calculate outlet discharge
%--------------------------------------------------------------------------
[Qt,Qfrac]           = base_flow_fraction(V,d',Tmax',Smax',AREA,area',SINa',SINb',manN',cs,cW',Nc);
%**************************************************************************
function [Qt,Qfrac]=base_flow_fraction(V,d,Tmax,Smax,AREA,area,SINa,SINb,manN,cs,cW,Nc)
e                    = 1e-16;
%disaggregate variables 
Sx                   = V(:,1:Nc);
wtd                  = V(:,2*Nc+1:3*Nc);
%for numerical stability ensure positivity
wtd(wtd<e)           = e;
T                    = Tmax./(1+d.*wtd).^d;
qb                   = SINa.*T./cs;
qb                   = area.*qb*AREA;
%overland storage
Sx                   = Sx -Smax;
Sx(Sx<e)             = e;
%hydraulic radius for channel class
Ss                   = Sx./cW*cs;
R                    = Ss.*cW./(2*Ss+cW);
%manning's velocity
v                    = R.^(2/3).*SINb.^(1/2)./manN;
qo                   = Sx.*v/cs;
qo(qo>Sx)            = Sx(qo>Sx);
qo                   = area.*qo*AREA;
%calculate total flow (hydrograph)
Qt                   = qo(:,Nc) + qb(:,Nc);
Qfrac                = sum(area.*qb,2)./(sum(area.*(qb+qo),2)+e);
%**************************************************************************
%**************************************************************************
%                               subfunctions
%**************************************************************************
%**************************************************************************
%**************************************************************************
function dVdt=HydroGEM_ode_fun(t,V,area,d,Nc,Smax,FET,FDR,mannN,D,WxmD...
    ,WbmD,WxmU,WbmU,cs,SINa,SINb,COSa,COSb,Tmax,Nr,cW,Ep0)
%--------------------------------------------------------------------------
%in case it doesn't converge, end & return NaNs 
if toc>Nc*30
    dVdt=nan(3*Nc,1);
    return
end
%--------------------------------------------------------------------------
%machine precision
e                 = 1e-16;
%--------------------------------------------------------------------------
%                          disaggragte variables
%--------------------------------------------------------------------------
Sx                = V(1:Nc,1);          
Su                = V(Nc+1:2*Nc,1);     
wtd               = V(2*Nc+1:3*Nc,1); 
%--------------------------------------------------------------------------
%          interpolate rainfall and evapo-transpiration rates
%--------------------------------------------------------------------------
%interpolate rainfall intensity based on solver time
Rn                = FDR(t);   Rn(Rn<e) = e;
%interpolate potential evapotranspiration rate
Ep                = FET(t);   Ep(Ep<e) = e; Ep=Ep.*Ep0;
%--------------------------------------------------------------------------
%                       subsurface storage updates
%--------------------------------------------------------------------------
%fraction from saturated zone to the surface
[w0]              = stepfun(e,wtd,e);
qwt               = -w0.*wtd;
wtd               = wtd - qwt;
%fraction from unsaturated zone to the surface
[w1,Su]           = stepfun(Su,wtd,e);
quz               = w1.*(Su-wtd);
Su                = Su - quz;
%--------------------------------------------------------------------------
%                    water table inflows and outflows
%--------------------------------------------------------------------------
%subsurface power-law transmissivity profile
T              = Tmax./(1+d.*wtd).^d;
qb             = T./cs.*SINa;
[qbi,qb]       = distribute_flow(WbmD,area,1e64,qb,Nc,e);
%--------------------------------------------------------------------------
%                  unsaturated zone inflows and outflows
%--------------------------------------------------------------------------
e0                = 1e-8;
%mean vertical hydraulic conductivity in the unsaturated zone
Kbar              = (Tmax-T)./(wtd+e).*(Su./(wtd+e0));
%total vertical flux
qv                = Kbar.*Su./(wtd+e0);  
%ensure it doesn't exceed available Su
[w2,qv]           = stepfun(qv,Su,e);
qv                = w2.*Su+(1-w2).*qv;
%--------------------------------------------------------------------------
%                         surface storage updates
%--------------------------------------------------------------------------
%rootzone/interception storage
[w3,Sx]           = stepfun(Sx,Smax,e);
Srz               = w3.*Smax + (1-w3).*Sx;
%actual evapo-transpiration
Ea                = Ep.*Srz./(Smax + e);
[w4,Ea]           = stepfun(Ea,Srz,e);
Ea                = w4.*Srz + (1-w4).*Ea;
%update surface excess storage available for routing
Sx                = Sx - Srz ;  
%total available storage in the subsurface (saturated + unsaturated)
SD                = wtd-Su;
%vertical saturated hydraulic conductivity at the surface
Kmax              = d.^2*Tmax;
%green-ampt infiltration
qx                = max(Kmax.*(wtd+Sx-Su)./(wtd+e),e);
%now make sure infiltration doesn't exceed available surface storage
[w5,qx]           = stepfun(qx,Sx,e);
qx                = w5.*Sx + (1-w5).*qx;
%also make sure it doesn't exceed available subsurage storage deficit
[w6,qx]           = stepfun(qx,SD,e);
qx                = w6.*SD + (1-w6).*qx;
%--------------------------------------------------------------------------
%                       surface inflows and outflows
%--------------------------------------------------------------------------
%scale surface excess to acount for variable channel width (only in channel HSUs)
%hilslope cW is set equal to cs (DEM resolution), so no scaling occurs there
Ss                = Sx./cW*cs;
%hydraulic radius for channel HSUs (assuming rectangular channel)
R                 = Ss.*cW./(2*Ss+cW);
%hydraulic radius for hillslope HSUs
R(1:Nc-Nr)        = Ss(1:Nc-Nr);
%overland flow diffusion
dSxi              = repmat(Sx',Nc,1)-Sx;
dSxdx             = sum(D.*dSxi,2)./cs;
%diffusion flux
flxx              = SINb - COSb.*dSxdx;
%determine the sign (+:downslope, -:upslope)
wx                = stepfun(flxx,0,e);
%Manning's n downslope flow flux (v*depth/grid resolution)
qoD               = wx.*R.^(2/3).*flxx.^(1/2)./mannN.*Sx./cs;
%upslope flow flux 
qoU               = (1-wx).*R.^(2/3).*abs(flxx).^(1/2)./mannN.*Sx./cs;
%route downslope and upslope flows
[qoiD,qoD]        = distribute_flow(WxmD,area,Sx,qoD,Nc,e);
[qoiU,qoU]        = distribute_flow(WxmU,area,Sx,qoU,Nc,e);
%combine downslope and upslope
qo                = qoD  + qoU;
qoi               = qoiD + qoiU;
%--------------------------------------------------------------------------
%                       update all time derivatives 
%--------------------------------------------------------------------------
dSxdt             = qoi - qo  - qx + quz + qwt + Rn - Ea; %surface storage
dSudt             = qx  - qv  - quz;                      %unsaturated zone 
dwtddt            = qb  - qbi - qv + qwt;          %saturated zone
%aggregate derivatives for the ODE-solver
dVdt              = [dSxdt ; dSudt ; dwtddt];
%**************************************************************************
function [qoi,qo] = distribute_flow(W,area,S,q,Nc,e)
%ensure it does not exceed the available storage
[w7,q]            = stepfun(q,S,e);
qo                = (1-w7).*q + w7.*S;
%flux into units
qq                = area.*qo;
qq(Nc)            = e; %last HSU (outlet) does not give flow to other HSUs
qoi               = (W*qq)./(area+e);
%**************************************************************************
function [F,x]    = stepfun(x,x0,e)
x(x<e)            = e;
%a step function to handle discontinous in the solution
F                 = (1+tanh((x-x0)./e))./2;