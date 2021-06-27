(* ::Package:: *)

(* ::Title:: *)
(*PPN calculation for "New General Relativity" class of teleparallel gravity in xAct*)


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
(*Parameters in the action*)


(* ::Text:: *)
(*In its general form, the action for New General Relativity depends on three parameters.*)


Do[DefConstantSymbol[Symbol["c" <> ToString[i]], PrintAs -> StringJoin["\!\(c\_", ToString[i], "\)"]], {i, 3}]


(* ::Section:: *)
(*Gravitational constant*)


(* ::Text:: *)
(*Another constant parameter we need to define is the gravitational constant.*)


DefConstantSymbol[kappa, PrintAs -> "\[Kappa]"];


(* ::Section:: *)
(*Field Equations*)


(* ::Text:: *)
(*This tensor will hold the gravitational (tetrad) field equations.*)


DefTensor[TetEq[-T4\[Alpha], -T4\[Beta]], {MfSpacetime}, PrintAs -> "\[ScriptCapitalE]"];


(* ::Section:: *)
(*Constant coefficients*)


(* ::Text:: *)
(*For convenience, we define some constants, which will be used in the ansatz for solving the field equations.*)


aa[i_] := Module[{sym = Symbol["a" <> ToString[i]]}, If[!ConstantSymbolQ[sym], DefConstantSymbol[sym, PrintAs -> StringJoin["\!\(a\_", ToString[i], "\)"]]]; Return[sym]]


(* ::Chapter:: *)
(*Field equations and perturbative expansion*)


(* ::Section:: *)
(*Field equation*)


(* ::Text:: *)
(*Here we define the gravitational field equations.*)


TorsionFD[T4\[Gamma], -T4\[Mu], -T4\[Nu]] * TorsionFD[T4\[Delta], -T4\[Rho], -T4\[Sigma]] * (Met[-T4\[Gamma], -T4\[Delta]] * InvMet[T4\[Mu], T4\[Rho]] * InvMet[T4\[Nu], T4\[Sigma]] * c1 + delta[-T4\[Gamma], T4\[Sigma]] * delta[-T4\[Delta], T4\[Nu]] * InvMet[T4\[Mu], T4\[Rho]] * c2 + delta[-T4\[Gamma], T4\[Nu]] * delta[-T4\[Delta], T4\[Sigma]] * InvMet[T4\[Mu], T4\[Rho]] * c3) * Met[-T4\[Alpha], -T4\[Beta]] / 2 +
InvMet[T4\[Gamma], T4\[CurlyEpsilon]] * CD[-T4\[CurlyEpsilon]][2 * Met[-T4\[Beta], -T4\[Delta]] * TorsionFD[T4\[Delta], -T4\[Alpha], -T4\[Gamma]] * c1 + Met[-T4\[Gamma], -T4\[Delta]] * TorsionFD[T4\[Delta], -T4\[Alpha], -T4\[Beta]] * c2 - Met[-T4\[Alpha], -T4\[Delta]] * TorsionFD[T4\[Delta], -T4\[Gamma], -T4\[Beta]] * c2 + Met[-T4\[Alpha], -T4\[Beta]] * TorsionFD[T4\[Delta], -T4\[Delta], -T4\[Gamma]] * c3 - Met[-T4\[Gamma], -T4\[Beta]] * TorsionFD[T4\[Delta], -T4\[Delta], -T4\[Alpha]] * c3] +
InvMet[T4\[Delta], T4\[CurlyEpsilon]] * TorsionFD[T4\[Gamma], -T4\[CurlyEpsilon], -T4\[Alpha]] (Met[-T4\[Beta], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Gamma], -T4\[Delta]] - Met[-T4\[Gamma], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Delta], -T4\[Beta]] + Met[-T4\[Delta], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Gamma], -T4\[Beta]]) * c1 +
(Met[-T4\[Alpha], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Gamma], -T4\[Delta]] * InvMet[T4\[Delta], T4\[CurlyEpsilon]] * (2 * TorsionFD[T4\[Gamma], -T4\[CurlyEpsilon], -T4\[Beta]] - Met[-T4\[Beta], -T4\[Xi]] * InvMet[T4\[Gamma], T4\[Eta]] * TorsionFD[T4\[Xi], -T4\[Eta], -T4\[CurlyEpsilon]]) - InvMet[T4\[Delta], T4\[CurlyEpsilon]] * TorsionFD[T4\[Gamma], -T4\[CurlyEpsilon], -T4\[Alpha]] (Met[-T4\[Beta], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Gamma], -T4\[Delta]] - Met[-T4\[Gamma], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Delta], -T4\[Beta]] + Met[-T4\[Delta], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Gamma], -T4\[Beta]])) * c2 / 2 -
TorsionFD[T4\[Delta], -T4\[Delta], -T4\[Gamma]] * (TorsionFD[T4\[Gamma], -T4\[Alpha], -T4\[Beta]] + InvMet[T4\[Gamma], T4\[Zeta]] * (Met[-T4\[CurlyEpsilon], -T4\[Alpha]] * TorsionFD[T4\[CurlyEpsilon], -T4\[Beta], -T4\[Zeta]] + Met[-T4\[CurlyEpsilon], -T4\[Beta]] * TorsionFD[T4\[CurlyEpsilon], -T4\[Alpha], -T4\[Zeta]])) * c3 / 2 -
kappa^2 * EnergyMomentum[-T4\[Alpha], -T4\[Beta]]


(* ::Text:: *)
(*For illustrative purposes, we decompose the equations into the terms for the three coefficients.*)


Expand[%];
ContractMetric[%];
ToCanonical[%, UseMetricOnVBundle -> None];
Collect[%, {c1, c2, c3}]
Expand[%];
Simplify[%]


(* ::Text:: *)
(*Save the equations for later use.*)


teteqdef = TetEq[-T4\[Alpha], -T4\[Beta]] == %;
teteqru = mkr0[teteqdef];


(* ::Section:: *)
(*3 + 1 split*)


(* ::Text:: *)
(*Split the equations into their time and space components.*)


{#, # /. teteqru}&[TetEq[-T4\[Alpha], -T4\[Beta]]];
ChangeCovD[%, CD, PD];
Expand[%];
SpaceTimeSplits[#, {-T4\[Alpha] -> -T3a, -T4\[Beta] -> -T3b}]& /@ %;
Expand[%];
Map[ToCanonical, %, {3}];
Map[SortPDs, %, {3}];
teteq31list = %;
teteq31def = Union[Flatten[MapThread[Equal, %, 2]]];
teteq31ru = Flatten[mkrg /@ %];


(* ::Section:: *)
(*Velocity orders*)


(* ::Text:: *)
(*Expand each component into the respective velocity orders.*)


Outer[VelocityOrder, teteq31list, Range[0, 4]];
Map[NoScalar, %, {4}];
Expand[%];
Map[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]&, %, {4}];
Map[ToCanonical, %, {4}];
Map[SortPDs, %, {4}];
teteqvlist = Simplify[%];
teteqvdef = Union[Flatten[MapThread[Equal, %, 3]]]
teteqvru = Flatten[mkrg /@ %];


(* ::Chapter:: *)
(*Solution*)


(* ::Section:: *)
(*Check vacuum equations*)


(* ::Text:: *)
(*The zeroth order equations correspond to the vacuum. Check that they are solved for the assumed background.*)


eqns0 = {PPN[TetEq,0][-LI[0],-LI[0]], PPN[TetEq,0][-T3a,-T3b]} /. teteqvru


(* ::Section:: *)
(*Second order*)


(* ::Text:: *)
(*Extract the second order field equations.*)


eqns2 = FullSimplify[{PPN[TetEq,2][-LI[0],-LI[0]], PPN[TetEq,2][-T3a,-T3b]} /. teteqvru]


(* ::Text:: *)
(*Define an ansatz for the second order tetrad perturbations.*)


ans2def = {PPN[Met,2][-LI[0],-LI[0]] == aa[1] * PotentialU[], PPN[Met,2][-T3a,-T3b] == aa[2] * PotentialU[] * BkgMetricS3[-T3a,-T3b] + aa[3] * PotentialUU[-T3a,-T3b], PPN[Asym,2][-T3a,-T3b] == 0}
ans2ru = Flatten[mkrg /@ ans2def];


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


eqns2 /. ans2ru;
Expand[%];
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
ToCanonical /@ %;
PotentialUToChi /@ %;
PotentialUUToChi /@ %;
Expand[%];
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
PotentialToSource /@ %;
Expand[%];
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
ToCanonical /@ %;
SortPDs /@ %;
eqnsa2 = FullSimplify[%]


(* ::Text:: *)
(*The equations to be solved are extracted as coefficients of the matter terms. The last condition is the gauge condition, which mandates that the spatial part of the metric should be diagonal.*)


eqnsc2 = FullSimplify[{Coefficient[eqnsa2[[1]], Density[]], Coefficient[eqnsa2[[2]], Density[] * BkgMetricS3[-T3a,-T3b]], aa[3]}]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations.*)


sola2 = FullSimplify[First[Solve[# == 0& /@ eqnsc2, aa /@ Range[1, 3]]]]


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


eqns3 = FullSimplify[{PPN[TetEq, 3][-T3a, -LI[0]], PPN[TetEq, 3][-LI[0], -T3a]} /. teteqvru]


(* ::Text:: *)
(*Define an ansatz for the third order tetrad perturbations.*)


ans3def = {PPN[Met,3][-LI[0], -T3a] == aa[4] * PotentialV[-T3a] + aa[5] * PotentialW[-T3a], PPN[Asym,3][-LI[0], -T3a] == aa[6] * PotentialV[-T3a] + aa[7] * PotentialW[-T3a]}
ans3ru = Flatten[mkrg /@ ans3def];


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


eqns3 /. ans3ru /. sol2ru;
Expand[%];
PotentialWToChiV /@ %;
Expand[%];
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
PotentialChiToU /@ %;
PotentialVToU /@ %;
PotentialToSource /@ %;
ToCanonical /@ %;
SortPDs /@ %;
eqnsa3 = FullSimplify[%]


(* ::Text:: *)
(*The equations to be solved are extracted as coefficients of the matter terms. The last condition is the gauge condition, which mandates that the spatial part of the metric should be diagonal.*)


eqnsc3 = FullSimplify[Flatten[{
	Coefficient[eqnsa3, Density[] * Velocity[-T3a]],
	Coefficient[eqnsa3, ParamD[TimePar][PD[-T3a][PotentialU[]]]],
	aa[7] - aa[6] - aa[0]
}]]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations, keeping a gauge freedom, which is left up to the fourth velocity order.*)


sola3 = FullSimplify[First[Solve[# == 0& /@ eqnsc3, aa /@ Range[4, 7]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa3 /. sola3]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol3def = ans3def /. sola3
sol3ru = Flatten[mkrg /@ sol3def];


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


eqns4 = PPN[TetEq,4][-LI[0], -LI[0]] * (2 * c1 + c2 + 2 * c3) + PPN[TetEq, 4][-T3a, T3a] * c3 /. teteqvru


(* ::Text:: *)
(*Define an ansatz for the fourth order tetrad perturbations.*)


ans4def = PPN[Met,4][-LI[0], -LI[0]] == aa[8] * PotentialPhi1[] + aa[9] * PotentialPhi2[] + aa[10] * PotentialPhi3[] + aa[11] * PotentialPhi4[] + aa[12] * PotentialU[]^2
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


sola4 = Simplify[First[Solve[# == 0& /@ {eq1, eq2, eq3, eq4, eq5, eq6}, aa /@ Prepend[Range[8, 12], 0]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa4 /. sola4]


(* ::Text:: *)
(*Enhance the third order solution, using the gauge fixing condition determined at the fourth order.*)


sol3def = ans3def /. Simplify[sola3 /. sola4]
sol3ru = Flatten[mkrg /@ sol3def];


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


(* ::Text:: *)
(*These parameters show the derivation from general relativity.*)


Factor[{ParameterBeta - 1, ParameterGamma - 1} /. parsol]


(* ::Text:: *)
(*One finds a family of solutions similar to general relativity.*)


Solve[{ParameterBeta == 1, ParameterGamma == 1} /. parsol, {c1, c2, c3}]


(* ::Text:: *)
(*This include in particular the teleparallel equivalent of general relativity.*)


tegr = {c1 -> 1/4, c2 -> 1/2, c3 -> -1}
parsol /. tegr
