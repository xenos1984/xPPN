(* ::Package:: *)

(* ::Title:: *)
(*Generalized Brans-Dicke type scalar-tensor gravity*)


(* ::Chapter:: *)
(*Preliminaries*)


(* ::Section:: *)
(*Load package*)


(* ::Text:: *)
(*This loads the xPPN package into the Mathematica session.*)


<< xAct`xPPN`


(* ::Section:: *)
(*Nicer printing*)


(* ::Text:: *)
(*Use this to automatically convert "dollar" indices in xAct before output is displayed.*)


$PrePrint = ScreenDollarIndices;


(* ::Section:: *)
(*Make rules from equations*)


(* ::Text:: *)
(*These functions help to create xAct rules, either with automatic raising / lowering of indices or without.*)


mkrg[eq_Equal] := MakeRule[Evaluate[List @@ eq], MetricOn -> All, ContractMetrics -> True];
mkr0[eq_Equal] := MakeRule[Evaluate[List @@ eq], MetricOn -> None, ContractMetrics -> False];


(* ::Chapter:: *)
(*Object definitions*)


(* ::Section:: *)
(*Scalar field*)


(* ::Text:: *)
(*The scalar field is defined as a tensor with no indices.*)


DefTensor[psi[], {MfSpacetime}, PrintAs -> "\[Psi]"];


(* ::Section:: *)
(*Background value of the scalar field*)


(* ::Text:: *)
(*For the background value of the scalar field we define a constant.*)


DefConstantSymbol[psi0, PrintAs -> "\[CapitalPsi]"];


(* ::Section:: *)
(*Gravitational constant*)


(* ::Text:: *)
(*Another constant parameter we need to define is the gravitational constant.*)


DefConstantSymbol[kappa, PrintAs -> "\[Kappa]"];


(* ::Section:: *)
(*Parameter function in the action*)


(* ::Text:: *)
(*The kinetic coupling of the scalar field is defined as a scalar function, which depends on the scalar field.*)


DefScalarFunction[omega, PrintAs -> "\[Omega]"];


(* ::Section:: *)
(*Field Equations*)


(* ::Text:: *)
(*These tensors will hold the gravitational (metric and scalar) field equations.*)


DefTensor[MetEq[-T4\[Alpha], -T4\[Beta]], {MfSpacetime}, Symmetric[{1,2}], PrintAs -> "\[ScriptCapitalE]"];
DefTensor[ScalEq[], {MfSpacetime}, PrintAs -> "\[ScriptCapitalE]"];


(* ::Section:: *)
(*Constant coefficients*)


(* ::Text:: *)
(*For convenience, we define some constants, which will be used in the ansatz for solving the field equations.*)


aa[i_] := Module[{sym = Symbol["a" <> ToString[i]]}, If[!ConstantSymbolQ[sym], DefConstantSymbol[sym, PrintAs -> StringJoin["\!\(a\_", ToString[i], "\)"]]]; Return[sym]]


(* ::Chapter:: *)
(*Field equations*)


(* ::Section:: *)
(*Metric field equation*)


(* ::Text:: *)
(*Here we define the trace-reversed metric field equations. These are more convenient to solve, since the time-time component contains only derivatives of the time-time component of the metric, and so the equations decouple immediately.*)


psi[] * RicciCD[-T4\[Alpha], -T4\[Beta]] - CD[-T4\[Alpha]][CD[-T4\[Beta]][psi[]]] - PD[-T4\[Alpha]][psi[]] * PD[-T4\[Beta]][psi[]] * omega[psi[]] / psi[] + omega'[psi[]] * InvMet[T4\[Gamma], T4\[Delta]] * PD[-T4\[Gamma]][psi[]] * PD[-T4\[Delta]][psi[]] * Met[-T4\[Alpha], -T4\[Beta]] / (4 * omega[psi[]] + 6) - kappa^2 * (EnergyMomentum[-T4\[Alpha], -T4\[Beta]] - InvMet[T4\[Gamma], T4\[Delta]] * EnergyMomentum[-T4\[Gamma], -T4\[Delta]] * Met[-T4\[Alpha], -T4\[Beta]] * (omega[psi[]] + 1) / (2 * omega[psi[]] + 3))


(* ::Text:: *)
(*Save the equations for later use.*)


meteqdef = MetEq[-T4\[Alpha], -T4\[Beta]] == %;
meteqru = mkr0[meteqdef];


(* ::Section:: *)
(*Scalar field equation*)


(* ::Text:: *)
(*In the scalar field equation we have already substituted the Ricci scalar, using the trace of the metric field equations, so that only second order derivatives of the scalar field remain.*)


(2 * omega[psi[]] + 3) * InvMet[T4\[Alpha], T4\[Beta]] * CD[-T4\[Alpha]][CD[-T4\[Beta]][psi[]]] + omega'[psi[]] * InvMet[T4\[Alpha], T4\[Beta]] * PD[-T4\[Alpha]][psi[]] * PD[-T4\[Beta]][psi[]] - kappa^2 * InvMet[T4\[Alpha], T4\[Beta]] * EnergyMomentum[-T4\[Alpha], -T4\[Beta]]


(* ::Text:: *)
(*Save the equations for later use.*)


scaleqdef = ScalEq[] == %;
scaleqru = mkr0[scaleqdef];


(* ::Chapter:: *)
(*Post-Newtonian expansion*)


(* ::Section:: *)
(*Rules for the scalar field*)


(* ::Text:: *)
(*When the scalar field is expanded in velocity orders, set the zeroth order to the constant background value, and let the first and third orders vanish.*)


OrderSet[PPN[psi, 0][], psi0];
OrderSet[PPN[psi, 1][], 0];
OrderSet[PPN[psi, 3][], 0];


(* ::Section:: *)
(*3 + 1 split*)


(* ::Text:: *)
(*Split the metric equations into their time and space components.*)


{#, # /. meteqru}&[MetEq[-T4\[Alpha], -T4\[Beta]]];
ChangeCovD[%, CD, PD];
Expand[%];
SpaceTimeSplits[#, {-T4\[Alpha] -> -T3a, -T4\[Beta] -> -T3b}]& /@ %;
Expand[%];
Map[ToCanonical, %, {3}];
Map[SortPDs, %, {3}];
meteq31list = %;
meteq31def = Union[Flatten[MapThread[Equal, %, 2]]];
meteq31ru = Flatten[mkrg /@ %];


(* ::Text:: *)
(*Split the scalar equation into its time and space components.*)


{#, # /. scaleqru}&[ScalEq[]];
ChangeCovD[%, CD, PD];
Expand[%];
SpaceTimeSplit[#, {}]& /@ %;
Expand[%];
ToCanonical /@ %;
SortPDs /@ %;
scaleq31list = %;
scaleq31def = Equal @@ %;
scaleq31ru = Flatten[mkrg[%]];


(* ::Section:: *)
(*Velocity orders*)


(* ::Text:: *)
(*Expand each component of the metric equations into the respective velocity orders.*)


Outer[VelocityOrder, meteq31list, Range[0, 4]];
Map[NoScalar, %, {4}];
Expand[%];
Map[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]&, %, {4}];
Map[ToCanonical, %, {4}];
Map[SortPDs, %, {4}];
meteqvlist = Simplify[%];
meteqvdef = Union[Flatten[MapThread[Equal, %, 3]]]
meteqvru = Flatten[mkrg /@ %];


(* ::Text:: *)
(*Expand the scalar equation into the respective velocity orders.*)


Outer[VelocityOrder, scaleq31list, Range[0, 4]];
Map[NoScalar, %, {2}];
Expand[%];
Map[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]&, %, {2}];
Map[ToCanonical, %, {2}];
Map[SortPDs, %, {2}];
scaleqvlist = Simplify[%];
scaleqvdef = Flatten[MapThread[Equal, %, 1]]
scaleqvru = Flatten[mkrg /@ %];


(* ::Chapter:: *)
(*Solution*)


(* ::Section:: *)
(*Check vacuum equations*)


(* ::Text:: *)
(*The zeroth order equations correspond to the vacuum. Check that they are solved for the assumed background.*)


eqns0 = {PPN[MetEq,0][-LI[0],-LI[0]], PPN[MetEq,0][-T3a,-T3b], PPN[ScalEq,0][]} /. meteqvru /. scaleqvru


(* ::Section:: *)
(*Second order*)


(* ::Text:: *)
(*Extract the second order field equations.*)


eqns2 = FullSimplify[{PPN[MetEq,2][-LI[0],-LI[0]], PPN[MetEq,2][-T3a,-T3b], PPN[ScalEq,2][]} /. meteqvru /. scaleqvru]


(* ::Text:: *)
(*Define an ansatz for the second order metric perturbations.*)


ans2def = {PPN[Met,2][-LI[0],-LI[0]] == aa[1] * PotentialU[], PPN[Met,2][-T3a,-T3b] == aa[2] * PotentialU[] * BkgMetricS3[-T3a,-T3b] + aa[3] * PotentialUU[-T3a,-T3b], PPN[psi,2][] == aa[4] * PotentialU[]}
ans2ru = Flatten[mkrg /@ ans2def];


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


eqns2 /. ans2ru;
PotentialUToChi /@ %;
PotentialUUToChi /@ %;
Expand[%];
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]& /@ %;
PotentialToSource /@ %;
Expand[%];
ToCanonical /@ %;
SortPDs /@ %;
eqnsa2 = FullSimplify[%]


(* ::Text:: *)
(*The equations to be solved are extracted as coefficients of the matter terms. The last condition is the gauge condition, which mandates that the spatial part of the metric should be diagonal.*)


eqnsc2 = FullSimplify[{Coefficient[eqnsa2[[1]], Density[]], Coefficient[eqnsa2[[2]], Density[] * BkgMetricS3[-T3a,-T3b]], Coefficient[eqnsa2[[3]], Density[]], aa[3]}]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations.*)


sola2 = FullSimplify[First[Solve[# == 0& /@ eqnsc2, aa /@ Range[1, 4]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa2 /. sola2]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol2def = ans2def /. sola2
sol2ru = Flatten[mkrg /@ sol2def];


(* ::Text:: *)
(*Check that the metric components we found indeed solve the second order field equations.*)


eqns2 /. sol2ru;
Expand[%];
PotentialToSource /@ %;
ToCanonical /@ %;
SortPDs /@ %;
Simplify[%]


(* ::Section:: *)
(*Third order*)


(* ::Text:: *)
(*Extract the third order field equations.*)


eqns3 = FullSimplify[PPN[MetEq, 3][-LI[0], -T3a] /. meteqvru]


(* ::Text:: *)
(*Define an ansatz for the third order metric perturbations.*)


ans3def = PPN[Met,3][-LI[0], -T3a] == aa[5] * PotentialV[-T3a] + aa[6] * PotentialW[-T3a]
ans3ru = mkrg[ans3def];


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


eqns3 /. ans3ru /. sol2ru;
PotentialWToChiV[%];
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialChiToU[%];
PotentialVToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
eqnsa3 = FullSimplify[%]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations, keeping a gauge freedom, which is left up to the fourth velocity order.*)


sola3 = FullSimplify[First[Solve[{eqnsa3 == 0, aa[6] - aa[5] == aa[0]}, {aa[5], aa[6]}]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa3 /. sola3]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol3def = ans3def /. sola3
sol3ru = mkrg[sol3def];


(* ::Text:: *)
(*Check that the metric components we found indeed solve the third order field equations.*)


eqns3 /. sol2ru /. sol3ru;
PotentialWToChiV[%];
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialChiToU[%];
PotentialVToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
Simplify[%]


(* ::Section:: *)
(*Fourth order*)


(* ::Text:: *)
(*Extract the fourth order field equations.*)


eqns4 = PPN[MetEq,4][-LI[0], -LI[0]] /. meteqvru


(* ::Text:: *)
(*Define an ansatz for the fourth order metric perturbations.*)


ans4def = PPN[Met,4][-LI[0], -LI[0]] == aa[7] * PotentialPhi1[] + aa[8] * PotentialPhi2[] + aa[9] * PotentialPhi3[] + aa[10] * PotentialPhi4[] + aa[11] * PotentialU[]^2
ans4ru = mkrg[ans4def];


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


eqns4 /. ans4ru /. sol2ru /. sol3ru;
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialVToU[%];
PotentialWToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
Expand[%];
eqnsa4 = Simplify[ScreenDollarIndices[%]]


(* ::Text:: *)
(*Coefficient of the pressure.*)


eq1 = Simplify[Coefficient[eqnsa4, Pressure[]]]


(* ::Text:: *)
(*Coefficient of the internal energy.*)


eq2 = Simplify[Coefficient[eqnsa4, Density[] * InternalEnergy[]]]


(* ::Text:: *)
(*Coefficient of the gravitational potential energy.*)


eq3 = Simplify[Coefficient[eqnsa4, Density[] * PotentialU[]]]


(* ::Text:: *)
(*Coefficient of the second time derivative.*)


eq4 = Simplify[Coefficient[eqnsa4, ParamD[TimePar, TimePar][PotentialU[]]]]


(* ::Text:: *)
(*Coefficient of the kinetic energy.*)


eq5 = Simplify[Coefficient[eqnsa4, Density[] * Velocity[-T3a] * Velocity[T3a]]]


(* ::Text:: *)
(*Coefficient of the potential term.*)


eq6 = Simplify[Coefficient[eqnsa4, PD[-T3a][PotentialU[]] * PD[T3a][PotentialU[]]]]


(* ::Text:: *)
(*Check that we have fully decomposed the equations.*)


Simplify[Pressure[] * eq1 + Density[] * InternalEnergy[] * eq2 + Density[] * PotentialU[] * eq3 + ParamD[TimePar, TimePar][PotentialU[]] * eq4 + Density[] * Velocity[-T3a] * Velocity[T3a] * eq5 + PD[-T3a][PotentialU[]] * PD[T3a][PotentialU[]] * eq6 - eqnsa4]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations.*)


sola4 = Simplify[First[Solve[# == 0& /@ {eq1, eq2, eq3, eq4, eq5, eq6}, aa /@ Prepend[Range[7, 11], 0]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa4 /. sola4]


(* ::Text:: *)
(*Enhance the third order solution, using the gauge fixing condition determined at the fourth order.*)


sol3def = ans3def /. Simplify[sola3 /. sola4]
sol3ru = mkrg[sol3def];


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol4def = ans4def /. sola4
sol4ru = mkrg[sol4def];


(* ::Text:: *)
(*Check that the metric components we found indeed solve the fourth order field equations.*)


eqns4 /. sol2ru /. sol3ru /. sol4ru;
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialVToU[%];
PotentialWToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
Expand[%];
Simplify[%]


(* ::Chapter:: *)
(*PPN metric and parameters*)


(* ::Section:: *)
(*PPN metric*)


(* ::Text:: *)
(*To read off the PPN parameters, we use the following metric components.*)


metcomp = {PPN[Met,2][-LI[0],-LI[0]], PPN[Met,2][-T3a,-T3b], PPN[Met,3][-LI[0],-T3a], PPN[Met,4][-LI[0],-LI[0]]}


(* ::Text:: *)
(*Insert the solution we obtained into the metric components.*)


metcomp /. sol2ru /. sol3ru /. sol4ru;
ToCanonical[%];
Expand[%];
ppnmet = Simplify[%];
metdef = MapThread[Equal, {metcomp, %}, 1]


(* ::Text:: *)
(*For later comparison, convert the selected components also to their standard form in terms of PPN parameters and potentials.*)


stamet = Simplify[MetricToStandard /@ metcomp]


(* ::Section:: *)
(*Newtonian gravitational constant*)


(* ::Text:: *)
(*To solve for the gravitational constant, compare the second order result with the standard normalization.*)


kappaeq = First[ppnmet] == First[stamet]


(* ::Text:: *)
(*When solving for \[Kappa], make sure to catch the positive root.*)


kappadef = kappa == First[Sqrt[FullSimplify[k2 /. Solve[kappaeq /. kappa -> Sqrt[k2], k2]]]]
kapparu = mkrg[kappadef];


(* ::Section:: *)
(*PPN parameters*)


(* ::Text:: *)
(*To generate equations for the PPN parameters, compare the obtained solution with the standard PPN metric.*)


pareqns = Simplify[ToCanonical[stamet - ppnmet /. kapparu]]


(* ::Text:: *)
(*We will consider the coefficients in front of the following potentials.*)


pots = {PotentialU[] * BkgMetricS3[-T3a,-T3b], PotentialV[-T3a], PotentialW[-T3a], PotentialA[], PotentialU[]^2, PotentialPhiW[], PotentialPhi1[], PotentialPhi2[], PotentialPhi3[], PotentialPhi4[]}


(* ::Text:: *)
(*Extract the coefficients from the difference between our result and the standard PPN metric. These terms must vanish.*)


eqns = DeleteCases[Flatten[Simplify[Outer[Coefficient, pareqns, pots]]],0]


(* ::Text:: *)
(*List of PPN parameters we are solving for.*)


pars = {ParameterBeta, ParameterGamma, ParameterXi, ParameterAlpha1, ParameterAlpha2, ParameterAlpha3, ParameterZeta1, ParameterZeta2, ParameterZeta3, ParameterZeta4}


(* ::Text:: *)
(*Finally, solve the equations and determine the PPN parameters.*)


parsol = FullSimplify[Solve[# == 0& /@ eqns, pars][[1]]]
