$offOrder
$onUelList

Sets
    t year          /0*4/
    i unit type     /A,B,C/
    j pv types /PV/
    a wt types /WT/
    k hour          /1*8760/
;    

Parameter D(k,t);
$gdxin in.gdx
$load D
$gdxIn
;

Parameter WTout(k,t);
$gdxin in.gdx
$load WTout
$gdxIn
;

Parameter PVout(k,t);
$gdxin in.gdx
$load PVout
$gdxIn
;

Scalars

    disc    discount rate        /0.1/
    To       years period       /5/
    alpha   depreciation rate  /0.25/  
;
    
Parameters

    year(t)         year-t                            /0 0,1 1, 2 2, 3 3, 4 4/
    
*Generator
    PG(i)           Generator Capacity (MW)              /A 150, B 250, C 100/
    cost_inv(i)     Cost of investment (R per MW)       /A 300000, B 350000, C 250000/
    cost_fuel(i)    Cost of fuel (R per MWh )          /A 20.409, B 14, C 25.953/
    cost_fixed(i)   Cost of Fixed O&M (R per MW)      /A 12000, B 36000, C 30000/
    cost_var(i)    Cost of variable O&M (R per MWh)  /A 1, B 3, C 2.5/
     
*Renewables  
    cost_ipv(j)     Cost of investment PV (R per MW)     /PV 1500000/
    cost_iwt(a)     Cost of investment WT (R per MW)    /WT 2000000/
   
    cost_fpv(j)     Cost of Fixed O&M (R per MW year)  /PV 29600/
    cost_fwt(a)     Cost of Fixed O&M (R per MW year) /WT 34700/
    
*Emission   
    Ro(t)           Carbon trading cost (R per tonne)          /0 70,1 70,2 70,3 70,4 70/
    e(i)            Emission rating of unit i (Tonne per MWh) /A 0.533, B 0.533, C 0.533/
    Hmax(t)         Maximum allowed emission in year t       /0 0,1 0,2 0,3 0,4 0/
;
    
Positive Variables

    P(i,k,t)    Power of unit i at hour k in year t
    PPV(j,k,t)  Power of PV unit j at hour k in year t
    PWT(a,k,t)  Power of WT unit j at hour k in year t

    Cinv        Total Investment cost
    Cfuel       Total Fuel cost
    Csalv       Total Salvage value
    Com         Total O&M cost
    Ccem        Carbon emission penalty cost R
   
;

Integer Variable

    X(i,t) Total unit i installed in year t
    Z(i,t) New unit i installed in year t
    N(j,t) New unit pv installed in year t
    M(a,t) New unit wt installed in year t
;

Variables F Total Cost;

Equations OBJ, INV, NU, SALV, FC, GCAP, LBALANCE, CONM, PVCAP, WCAP,INIT, CTD, EMISSION;

*Objective Function
OBJ..       F =E= Cinv+Cfuel+Com+Ccem-Csalv;

*Economical Calculation
INV..       Cinv =E= sum(t, (1/power((1+disc),year(t)) * (sum(i, cost_inv(i)*PG(i)*Z(i,t)) + (sum(j,N(j,t)*cost_ipv(j))) + (sum(a,M(a,t)*cost_iwt(a))))));

SALV..      Csalv =E= sum(t, (1/power((1+disc),To) * (sum(i, cost_inv(i)*PG(i)*X(i,t)*power((1-alpha),(To-year(t)))) + (sum(j,N(j,t)*cost_ipv(j))) + (sum(a,M(a,t)*cost_iwt(a))))));

FC..        CFuel =E= sum(t, (1/power((1+disc),year(t)) * sum((i,k), P(i,k,t)*cost_fuel(i))));

CONM..      Com =E= sum(t, (1/power((1+disc),year(t)) * (sum((i,k), cost_fixed(i)*PG(i)*X(i,t) + P(i,k,t)*cost_var(i)) + (sum(j,N(j,t)*cost_fpv(j))) + (sum(a,M(a,t)*cost_fwt(a))))));

CTD..       Ccem =E= sum(t, (sum((i,k),P(i,k,t)*e(i)*Ro(t))));

*Power Generation Constraint
GCAP(i,k,t)..    P(i,k,t) =L= X(i,t)*PG(i);

PVCAP(j,k,t)..   PPV(j,k,t) =L= N(j,t) * PVout(k,t);

WCAP(a,k,t)..    PWT(a,k,t) =L= M(a,t) * WTout(k,t);

LBALANCE(k,t)..  D(k,t) =L= sum(i, P(i,k,t)) +  sum(j,PPV(j,k,t)) +  sum(a,PWT(a,k,t));


*conditions
INIT(i,t)..  X(i,t) =G= 0;
NU(i,t)..   Z(i,t) =E= X(i,t)-X(i,t-1);
EMISSION(t).. Hmax(t) =L= sum((i,k),P(i,k,t)*e(i));
*we chage this constraint from =G= to =L=, such that the emission will not be limited


Model INTEGRATION /all/;

solve INTEGRATION using MIP minimizing F;
display Cinv.l, Cfuel.l, Csalv.l, Com.l, Z.l, X.l, F.l, N.l, M.l, Ccem.l;