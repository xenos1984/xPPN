(* ::Package:: *)

(* ::Title:: *)
(*"Newer General Relativity" class of symmetric teleparallel gravity*)


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
(*In its general form, the action for Newer General Relativity depends on five parameters, for which we introduce two parametrizations.*)


Do[DefConstantSymbol[Symbol["c" <> ToString[i]], PrintAs -> StringJoin["\!\(c\_", ToString[i], "\)"]], {i, 5}]
Do[DefConstantSymbol[Symbol["k" <> ToString[i]], PrintAs -> StringJoin["\!\(k\_", ToString[i], "\)"]], {i, 5}]


(* ::Section:: *)
(*Gravitational constant*)


(* ::Text:: *)
(*Another constant parameter we need to define is the gravitational constant.*)


DefConstantSymbol[kappa, PrintAs -> "\[Kappa]"];


(* ::Section:: *)
(*Field Equations*)


(* ::Text:: *)
(*This tensor will hold the gravitational (metric) field equations.*)


DefTensor[MetEq[-T4\[Alpha], -T4\[Beta]], {MfSpacetime}, Symmetric[{1, 2}], PrintAs -> "\[ScriptCapitalE]"];


(* ::Section:: *)
(*Constant coefficients*)


(* ::Text:: *)
(*For convenience, we define some constants, which will be used in the ansatz for solving the field equations.*)


aa[i_] := Module[{sym = Symbol["a" <> ToString[i]]}, If[!ConstantSymbolQ[sym], DefConstantSymbol[sym, PrintAs -> StringJoin["\!\(a\_", ToString[i], "\)"]]]; Return[sym]]


(* ::Chapter:: *)
(*Field equations and perturbative expansion*)


(* ::Section:: *)
(*Field equations*)


(* ::Text:: *)
(*Here we define the gravitational field equations, starting with the terms without derivatives.*)


feq0n = Met[-T4\[Alpha], -T4\[Beta]] * NonMet[-T4\[Mu], -T4\[Nu], -T4\[Gamma]] * NonMet[-T4\[Delta], -T4\[Epsilon], -T4\[Zeta]] * (
	InvMet[T4\[Mu], T4\[Delta]] * InvMet[T4\[Nu], T4\[Epsilon]] * InvMet[T4\[Gamma], T4\[Zeta]] * c1 +
	InvMet[T4\[Mu], T4\[Epsilon]] * InvMet[T4\[Nu], T4\[Zeta]] * InvMet[T4\[Gamma], T4\[Delta]] * c2 +
	InvMet[T4\[Mu], T4\[Delta]] * InvMet[T4\[Nu], T4\[Gamma]] * InvMet[T4\[Epsilon], T4\[Zeta]] * c3 +
	InvMet[T4\[Mu], T4\[Nu]] * InvMet[T4\[Delta], T4\[Epsilon]] * InvMet[T4\[Gamma], T4\[Zeta]] * c4 +
	InvMet[T4\[Mu], T4\[Nu]] * InvMet[T4\[Zeta], T4\[Epsilon]] * InvMet[T4\[Gamma], T4\[Delta]] * c5
) / 2;
feq1n = (2 * NonMet[-T4\[Gamma], -T4\[Epsilon], -T4\[Alpha]] * NonMet[-T4\[Zeta], -T4\[Delta], -T4\[Beta]] - NonMet[-T4\[Alpha], -T4\[Gamma], -T4\[Epsilon]] * NonMet[-T4\[Beta], -T4\[Delta], -T4\[Zeta]] - NonMet[-T4\[Gamma], -T4\[Epsilon], -T4\[Alpha]] * NonMet[-T4\[Beta], -T4\[Delta], -T4\[Zeta]] - NonMet[-T4\[Gamma], -T4\[Epsilon], -T4\[Beta]] * NonMet[-T4\[Alpha], -T4\[Delta], -T4\[Zeta]]) * InvMet[T4\[Gamma], T4\[Delta]] * InvMet[T4\[Epsilon], T4\[Zeta]] * c1;
feq2n = (NonMet[-T4\[Gamma], -T4\[Epsilon], -T4\[Alpha]] * NonMet[-T4\[Delta], -T4\[Zeta], -T4\[Beta]] - NonMet[-T4\[Alpha], -T4\[Gamma], -T4\[Epsilon]] * NonMet[-T4\[Beta], -T4\[Delta], -T4\[Zeta]] - NonMet[-T4\[Gamma], -T4\[Epsilon], -T4\[Alpha]] * NonMet[-T4\[Beta], -T4\[Delta], -T4\[Zeta]] / 2 - NonMet[-T4\[Gamma], -T4\[Epsilon], -T4\[Beta]] * NonMet[-T4\[Alpha], -T4\[Delta], -T4\[Zeta]] / 2) * InvMet[T4\[Gamma], T4\[Delta]] * InvMet[T4\[Epsilon], T4\[Zeta]] * c2;
feq3n = -NonMet[-T4\[Alpha], -T4\[Gamma], -T4\[Delta]] * NonMet[-T4\[Beta], -T4\[Epsilon], -T4\[Zeta]] * InvMet[T4\[Gamma], T4\[Delta]] * InvMet[T4\[Epsilon], T4\[Zeta]] * c3;
feq4n = ((NonMet[-T4\[Gamma], -T4\[Alpha], -T4\[Beta]] - NonMet[-T4\[Alpha], -T4\[Beta], -T4\[Gamma]] - NonMet[-T4\[Beta], -T4\[Alpha], -T4\[Gamma]]) * NonMet[-T4\[Epsilon], -T4\[Zeta], -T4\[Delta]] + NonMet[-T4\[Gamma], -T4\[Delta], -T4\[Alpha]] * NonMet[-T4\[Epsilon], -T4\[Zeta], -T4\[Beta]] - NonMet[-T4\[Gamma], -T4\[Delta], -T4\[Alpha]] * NonMet[-T4\[Beta], -T4\[Epsilon], -T4\[Zeta]] / 2 - NonMet[-T4\[Gamma], -T4\[Delta], -T4\[Beta]] * NonMet[-T4\[Alpha], -T4\[Epsilon], -T4\[Zeta]] / 2) * InvMet[T4\[Gamma], T4\[Delta]] * InvMet[T4\[Epsilon], T4\[Zeta]] * c4;
feq5n = ((NonMet[-T4\[Gamma], -T4\[Alpha], -T4\[Beta]] - NonMet[-T4\[Alpha], -T4\[Beta], -T4\[Gamma]] - NonMet[-T4\[Beta], -T4\[Alpha], -T4\[Gamma]]) * NonMet[-T4\[Delta], -T4\[Epsilon], -T4\[Zeta]] - NonMet[-T4\[Alpha], -T4\[Gamma], -T4\[Delta]] * NonMet[-T4\[Beta], -T4\[Epsilon], -T4\[Zeta]]) * InvMet[T4\[Gamma], T4\[Delta]] * InvMet[T4\[Epsilon], T4\[Zeta]] * c5 / 2;


(* ::Text:: *)
(*These terms involve derivatives.*)


feq1d = -2 * InvMet[T4\[Gamma], T4\[Delta]] * CD[-T4\[Gamma]][NonMet[-T4\[Delta], -T4\[Alpha], -T4\[Beta]]] * c1;
feq2d = -InvMet[T4\[Gamma], T4\[Delta]] * CD[-T4\[Gamma]][NonMet[-T4\[Alpha], -T4\[Beta], -T4\[Delta]] + NonMet[-T4\[Beta], -T4\[Alpha], -T4\[Delta]]] * c2;
feq3d = -2 * InvMet[T4\[Gamma], T4\[Delta]] * InvMet[T4\[Epsilon], T4\[Zeta]] * Met[-T4\[Alpha], -T4\[Beta]] * CD[-T4\[Gamma]][NonMet[-T4\[Delta], -T4\[Epsilon], -T4\[Zeta]]] * c3;
feq4d = -InvMet[T4\[Epsilon], T4\[Zeta]] * CD[-T4\[Gamma]][NonMet[-T4\[Epsilon], -T4\[Zeta], -T4\[Alpha]] * delta[-T4\[Beta], T4\[Gamma]] + NonMet[-T4\[Epsilon], -T4\[Zeta], -T4\[Beta]] * delta[-T4\[Alpha], T4\[Gamma]]] * c4;
feq5d = -InvMet[T4\[Epsilon], T4\[Zeta]] * CD[-T4\[Gamma]][InvMet[T4\[Gamma], T4\[Delta]] * NonMet[-T4\[Epsilon], -T4\[Zeta], -T4\[Delta]] * Met[-T4\[Alpha], -T4\[Beta]] + NonMet[-T4\[Alpha], -T4\[Epsilon], -T4\[Zeta]] * delta[-T4\[Beta], T4\[Gamma]] / 2 + NonMet[-T4\[Beta], -T4\[Epsilon], -T4\[Zeta]] * delta[-T4\[Alpha], T4\[Gamma]] / 2] * c5;


(* ::Text:: *)
(*We finally combine these terms, together with thr energy-momentum tensor.*)


feq0n + feq1n + feq2n + feq3n + feq4n + feq5n + feq1d + feq2d + feq3d + feq4d + feq5d + EnergyMomentum[-T4\[Alpha], -T4\[Beta]] * kappa^2;
Expand[%];
ContractMetric[%];
ToCanonical[%, UseMetricOnVBundle -> None];
Collect[%, {c1, c2, c3, c4, c5}]
Expand[%];
Simplify[%]


(* ::Text:: *)
(*Save the equations for later use.*)


meteqdef = MetEq[-T4\[Alpha], -T4\[Beta]] == %;
meteqru = mkr0[meteqdef];


(* ::Section:: *)
(*3 + 1 split*)


(* ::Text:: *)
(*Split the equations into their time and space components.*)


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


(* ::Section:: *)
(*Velocity orders*)


(* ::Text:: *)
(*Expand each component into the respective velocity orders.*)


Outer[VelocityOrder, meteq31list, Range[0, 4]];
Map[NoScalar, %, {4}];
Expand[%];
Map[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]&, %, {4}];
Map[ToCanonical, %, {4}];
Map[SortPDs, %, {4}];
meteqvlist = Simplify[%];
meteqvdef = Union[Flatten[MapThread[Equal, %, 3]]];
meteqvru = Flatten[mkrg /@ %];


(* ::Section:: *)
(*Change parametrization*)


(* ::Text:: *)
(*It is more convenient to introduce this change of parameters.*)


parchg = {c1 -> 3 * k5, c2 -> (k1 + k2 + k3 - 2 * k4 - 4 * k5) / 2, c3 -> k2 - k5, c4 -> (-k1 + k2 + k3 - 2 * k4 - 4 * k5) / 2, c5 -> 2 * (k5 + k4 - k2)}


(* ::Text:: *)
(*Apply this to the equations.*)


meteqvlist = Simplify[meteqvlist /. parchg];
meteqvdef = Union[Flatten[MapThread[Equal, %, 3]]];
meteqvru = Flatten[mkrg /@ %];


(* ::Chapter:: *)
(*Solution*)


(* ::Section:: *)
(*Check vacuum equations*)


(* ::Text:: *)
(*The zeroth order equations correspond to the vacuum. Check that they are solved for the assumed background.*)


eqns0 = {PPN[MetEq,0][-LI[0],-LI[0]], PPN[MetEq,0][-T3a,-T3b]} /. meteqvru


(* ::Section:: *)
(*Second order*)


(* ::Text:: *)
(*Extract the second order field equations.*)


eqns2 = FullSimplify[{PPN[MetEq,2][-LI[0],-LI[0]], PPN[MetEq,2][-T3a,-T3b]} /. meteqvru]


(* ::Text:: *)
(*Define an ansatz for the second order perturbations.*)


ans2def = {PPN[Met,2][-LI[0],-LI[0]] == aa[1] * PotentialU[], PPN[Met,2][-T3a,-T3b] == aa[2] * PotentialU[] * BkgMetricS3[-T3a,-T3b] + aa[3] * PotentialUU[-T3a,-T3b], PPN[Xi, 2][T3a] == aa[4] * PD[T3a][PotentialChi[]]}
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


coeff2 = Simplify[{
	Coefficient[eqnsa2[[1]], Density[]] / 8 / Pi,
	Coefficient[eqnsa2[[2]], Density[] * BkgMetricS3[-T3a, -T3b]] / 8 / Pi,
	Times @@ Select[List @@ (eqnsa2[[2]] /. Density[] -> 0), FreeQ[#, PotentialChi] &],
	aa[3]
}]


(* ::Text:: *)
(*We check whether these equations are independent.*)


TableForm[mat = Simplify[Transpose[Table[Coefficient[#, aa[i]] & /@ coeff2, {i, 4}]]]]
Factor[Det[mat]]


(* ::Text:: *)
(*General relativity is included in a degenerate case, which must be studied separately.*)


degen = {k3 -> 0, k4 -> 0};


(* ::Subsection:: *)
(*Non-degenerate case*)


(* ::Text:: *)
(*Solve for the constant coefficients in the equations.*)


sola2 = FullSimplify[First[Solve[# == 0& /@ coeff2, aa /@ Range[1, 4]]]]


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
ToCanonical /@ %;
PotentialUToChi /@ %;
PotentialUUToChi /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
ToCanonical /@ %;
PotentialToSource[%];
SortPDs /@ %;
Simplify[%]


(* ::Subsection:: *)
(*Degenerate case*)


(* ::Text:: *)
(*Solve for the constant coefficients in the equations.*)


sola2d = FullSimplify[First[Solve[# == 0& /@ coeff2 /. degen, aa /@ Range[1, 3]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[Simplify[eqnsa2 /. degen] /. sola2d]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol2defd = ans2def /. sola2d
sol2rud = Flatten[mkrg /@ sol2defd];


(* ::Text:: *)
(*Check that the metric components we found indeed solve the second order field equations.*)


Simplify[eqns2 /. degen] /. sol2rud;
Expand[%];
ToCanonical /@ %;
PotentialUToChi /@ %;
PotentialUUToChi /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
ToCanonical /@ %;
PotentialToSource[%];
SortPDs /@ %;
Simplify[%]


(* ::Section:: *)
(*Third order*)


(* ::Text:: *)
(*Extract the third order field equations.*)


eqns3 = PPN[MetEq, 3][-LI[0], -T3a] /. meteqvru


(* ::Text:: *)
(*Define an ansatz for the third order tetrad perturbations.*)


ans3def = {PPN[Met, 3][-LI[0], -T3a] == aa[5] * PotentialV[-T3a] + aa[6] * PotentialW[-T3a], PPN[Xi, 3][LI[0]] == -aa[7] * ParamD[TimePar][PotentialChi[]]}
ans3ru = Flatten[mkrg /@ ans3def];


(* ::Subsection:: *)
(*Non-degenerate case*)


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


eqns3 /. ans3ru /. sol2ru;
Expand[%];
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
(*The equations to be solved are extracted as coefficients of the matter terms. The last condition is the gauge condition, to be solved later.*)


coeff3 = FullSimplify[Flatten[{
	Coefficient[eqnsa3, Density[] * Velocity[-T3a]],
	Coefficient[eqnsa3, ParamD[TimePar][PD[-T3a][PotentialU[]]]],
	aa[7] - aa[0]
}]]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations, keeping a gauge freedom, which is left up to the fourth velocity order.*)


sola3 = FullSimplify[First[Solve[# == 0& /@ coeff3, aa /@ Range[5, 7]]]]


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
Expand[%];
PotentialWToChiV[%];
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialChiToU[%];
PotentialVToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
Simplify[%]


(* ::Subsection:: *)
(*Degenerate case*)


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


Simplify[eqns3 /. degen] /. ans3ru /. sol2rud;
Expand[%];
PotentialWToChiV[%];
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialChiToU[%];
PotentialVToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
eqnsa3d = FullSimplify[%]


(* ::Text:: *)
(*The equations to be solved are extracted as coefficients of the matter terms. The last condition is the gauge condition, to be solved later.*)


coeff3d = FullSimplify[Flatten[{
	Coefficient[eqnsa3d, Density[] * Velocity[-T3a]],
	Coefficient[eqnsa3d, ParamD[TimePar][PD[-T3a][PotentialU[]]]],
	aa[7] - aa[4] - aa[0]
}]]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations, keeping a gauge freedom, which is left up to the fourth velocity order.*)


sola3d = FullSimplify[First[Solve[# == 0& /@ coeff3d, aa /@ Range[5, 7]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa3d /. sola3d]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol3defd = ans3def /. sola3d
sol3rud = Flatten[mkrg /@ sol3defd];


(* ::Text:: *)
(*Check that the metric components we found indeed solve the third order field equations.*)


Simplify[eqns3 /. degen] /. sol2rud /. sol3rud;
Expand[%];
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
(*Fourth order tetrad*)


(* ::Text:: *)
(*We extract the fourth order field equations.*)


{PD[-T3a][PD[T3a][PPN[MetEq, 4][-LI[0], -LI[0]]]], PD[-T3b][PD[T3b][PPN[MetEq, 4][T3a, -T3a]]], PD[-T3a][PD[-T3b][PPN[MetEq, 4][T3a, T3b]]]};
% /. meteqvru;
Expand /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]& /@ %;
ToCanonical /@ %;
SortPDs /@ %;
eqns4 = %;


(* ::Text:: *)
(*To solve the fourth order equations, we must first eliminate the irrelevant terms.*)


elim = SortPDs /@ ToCanonical /@ {
	PD[-T3a][PD[T3a][PD[-T3b][PD[T3b][PPN[Met, 4][T3c, -T3c]]]]],
	PD[-T3a][PD[T3a][PD[-T3c][PD[-T3b][PPN[Met, 4][T3b, T3c]]]]],
	PD[-T3a][PD[T3a][PD[-T3b][PD[T3b][PD[-T3c][PPN[Xi, 4][T3c]]]]]]
}


(* ::Text:: *)
(*These are the coefficients of the terms to be eliminated.*)


TableForm[cmat = Simplify[Outer[IndexCoefficient, eqns4, elim]]]


(* ::Text:: *)
(*The fourth velocity order of this theory in general exceeds the standard PPN form. First, we list the standard PPN potentials to use in the solution.*)


pot4 = {PotentialU[]^2, PotentialPhi1[], PotentialPhi2[], PotentialPhi3[], PotentialPhi4[], PotentialPhiW[], PotentialA[], PotentialB[]}


(* ::Text:: *)
(*The potentials above are not sufficient for the solution of the general theory. We define the terms that will solve the latter.*)


ddpot4 = {
	PD[-T3a][PD[T3a][Pressure[]]],
	PD[-T3a][PD[T3a][Density[] * InternalEnergy[]]],
	PD[-T3a][PD[T3a][Density[] * Velocity[-T3b] * Velocity[T3b]]],
	PD[-T3a][PD[T3a][PotentialChi[]]] * PD[-T3b][PD[T3b][PD[-T3c][PD[T3c][PD[-T3d][PD[T3d][PotentialChi[]]]]]]],
	PD[-T3a][PD[-T3b][PotentialChi[]]] * PD[T3a][PD[T3b][PD[-T3c][PD[T3c][PD[-T3d][PD[T3d][PotentialChi[]]]]]]],
	PD[-T3a][PD[-T3b][PD[T3b][PotentialChi[]]]] * PD[T3a][PD[-T3c][PD[T3c][PD[-T3d][PD[T3d][PotentialChi[]]]]]],
	PD[-T3a][PD[-T3b][PD[-T3c][PotentialChi[]]]] * PD[T3a][PD[T3b][PD[T3c][PD[-T3d][PD[T3d][PotentialChi[]]]]]],
	PD[-T3a][PD[T3a][PD[-T3b][PD[T3b][PotentialChi[]]]]] * PD[-T3c][PD[T3c][PD[-T3d][PD[T3d][PotentialChi[]]]]],
	PD[-T3a][PD[-T3b][PD[-T3c][PD[T3c][PotentialChi[]]]]] * PD[T3a][PD[T3b][PD[-T3d][PD[T3d][PotentialChi[]]]]],
	PD[-T3a][PD[-T3b][PD[-T3c][PD[-T3d][PotentialChi[]]]]] * PD[T3a][PD[T3b][PD[T3c][PD[T3d][PotentialChi[]]]]],
	PD[-T3a][PD[-T3b][Density[] * Velocity[T3a] * Velocity[T3b]]]
}


(* ::Text:: *)
(*To relate these lists of terms, define coefficients to solve for.*)


zz[i_, j_] := Module[{sym = Symbol["z" <> ToString[i] <> "z" <> ToString[j]]}, If[! ConstantSymbolQ[sym], DefConstantSymbol[sym, PrintAs -> StringJoin["\!\(z\_\(", ToString[i], ",", ToString[j], "\)\)"]]]; Return[sym]]
TableForm[zmat = Array[zz, {8, 11}]]


(* ::Text:: *)
(*Now we can compare the second list to the double Laplacian of the first.*)


PD[-T3a][PD[T3a][PD[-T3b][PD[T3b][pot4]]]] - zmat . ddpot4;
Expand /@ %;
PotentialToSource /@ %;
PotentialChiToU /@ %;
PotentialToSource /@ %;
SeparateMetric[] /@ %;
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
ToCanonical /@ %;
Simplify /@ %;
Expand /@ %;
poteq = ScreenDollarIndices[%];


(* ::Text:: *)
(*These are the non-constant terms whose coefficients must be equal.*)


terms = Union[Times @@ DeleteCases[Select[List @@ #, Not@*ConstantQ], Power[_, a_ /; a < 0]] & /@ Flatten[poteq /. Plus -> List]]


(* ::Text:: *)
(*Then we extract the coefficient equations.*)


peq = Simplify[Outer[Coefficient, poteq, terms]];


(* ::Text:: *)
(*Check that we have included all terms.*)


peq . terms - poteq;
Expand /@ %


(* ::Text:: *)
(*Finally, solve for the coefficients.*)


zsol = Solve[# == 0 & /@ Flatten[peq], Flatten[zmat]][[1]];
TableForm[potrel = zmat /. zsol]


(* ::Text:: *)
(*The ansatz for solving the field equations will be an arbitrary linear combination.*)


a4var = Array[aa, 11, 8];
a4var2 = Array[aa, 8, 19];
ans4def = PD[-T3a][PD[T3a][PD[-T3b][PD[T3b][PPN[Met, 4][-LI[0], -LI[0]]]]]] == a4var . ddpot4;
ans4def2 = PPN[Met, 4][-LI[0], -LI[0]] == a4var2 . pot4;
ans4ru = mkrg[ans4def];
ans4ru2 = mkrg[ans4def2];


(* ::Text:: *)
(*This is how the coefficients of the more general ansatz look like if we restrict to the standard PPN ansatz.*)


varrel = MapThread[Equal, {a4var, a4var2 . potrel}, 1]


(* ::Text:: *)
(*If we eliminate the coefficients of the standard PPN ansatz, we find the conditions that the parameters of the general ansatz must satisfy in order to reduce to the standard PPN case.*)


varcond = Eliminate[varrel, a4var2] /. And -> List


(* ::Subsection:: *)
(*Non-degenerate case*)


(* ::Text:: *)
(*In the non-degenerate case, the coefficient matrix has a one-dimensional kernel.*)


nvec = FullSimplify[NullSpace[Transpose[cmat]][[1]]];
nvec = FullSimplify[Denominator[nvec[[1]]] * nvec]


(* ::Text:: *)
(*This is the equation to solve.*)


eqn4 = Expand[nvec . eqns4];


(* ::Text:: *)
(*We check that the terms we wanted to eliminate are indeed gone.*)


IndexCoefficient[eqn4, #] & /@ elim


(* ::Text:: *)
(*Now we use the general ansatz and convert all potentials and matter terms into a common form.*)


eqn4 /. ans4ru /. sol2ru /. sol3ru;
Expand[%];
PotentialUToChi[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
ToCanonical[%];
PotentialVToU[%];
PotentialWToU[%];
PotentialChiToU[%];
PotentialToSource[%];
PotentialChiToU[%];
PotentialToSource[%];
TimeRhoToEuler[%];
TimeRhoToEuler[%];
TimeVelToEuler[%] /. sol2ru;
PotentialToSource[%];
Expand[%];
ToCanonical[%];
List @@ %;
SeparateMetric[] /@ %;
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
Plus @@ %;
ToCanonical[%];
Simplify[%];
Expand[%];
eqnsa4 = ScreenDollarIndices[%];


(* ::Text:: *)
(*The equations to solve for are the coefficients of the non-constant terms.*)


eqs = Simplify[Coefficient[eqnsa4, #]] & /@ terms;


(* ::Text:: *)
(*Check that we have fully decomposed the equations.*)


Simplify[eqs . terms - eqnsa4]


(* ::Text:: *)
(*To solve the equations, we must supply a gauge condition. The one we use here is satisfied in the standard PPN gauge.*)


sola4 = Simplify[Solve[Append[# == 0 & /@ eqs, aa[11] == aa[15]], Prepend[a4var, aa[0]]][[1]]]


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


eqn4 /. sol2ru /. sol3ru /. sol4ru;
Expand[%];
PotentialUToChi[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
ToCanonical[%];
PotentialVToU[%];
PotentialWToU[%];
PotentialChiToU[%];
PotentialToSource[%];
PotentialChiToU[%];
PotentialToSource[%];
TimeRhoToEuler[%];
TimeRhoToEuler[%];
TimeVelToEuler[%] /. sol2ru;
PotentialToSource[%];
Expand[%];
ToCanonical[%];
List @@ %;
SeparateMetric[] /@ %;
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
Plus @@ %;
ToCanonical[%];
Simplify[%]


(* ::Text:: *)
(*The solution will in general not have the standard PPN form. These are the conditions which must be satisfied.*)


varconda4 = FullSimplify[varcond /. sola4]


(* ::Text:: *)
(*There are three possible solutions.*)


parcond = Simplify[Solve[varconda4, {k5, k1}, Assumptions -> {k3 != 0, k4 != 0}]]


(* ::Text:: *)
(*In these cases, the metric takes the standard PPN form, and we can solve for the coefficients in the corresponding expansion.*)


varrel /. sola4 /. parcond;
Simplify[%];
sola4p = Simplify[Solve[#, a4var2][[1]] & /@ %]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solutions for the metric components.*)


sol4defp = ans4def2 /. sola4p
sol4rup = mkrg /@ sol4defp;


(* ::Text:: *)
(*Finally, we combine all solutions for each of the three branches.*)


soldefp = Transpose[Append[Transpose[MapAt[Simplify, Join[sol2def, sol3def] /. parcond, {All, All, 2}]], sol4defp]]
solrup = Flatten /@ Map[mkrg, soldefp, {2}];


(* ::Subsection:: *)
(*Degenerate case*)


(* ::Text:: *)
(*These are the coefficients of the terms to be eliminated.*)


TableForm[cmatd = Simplify[cmat /. degen]]


(* ::Text:: *)
(*In the degenerate case, the coefficient matrix has a two-dimensional kernel.*)


nvecd = FullSimplify[NullSpace[Transpose[cmatd[[1 ;; 2, 1 ;; 2]]]][[1]]];
nvecd = FullSimplify[Denominator[nvecd[[1]]] * nvecd]


(* ::Text:: *)
(*This is the equation to solve.*)


eqn4d = Expand[nvecd . eqns4[[1 ;; 2]] /. degen];


(* ::Text:: *)
(*We check that the terms we wanted to eliminate are indeed gone.*)


IndexCoefficient[eqn4d, #] & /@ elim


(* ::Text:: *)
(*Now we use the general ansatz and convert all potentials and matter terms into a common form.*)


eqn4d /. ans4ru /. sol2rud /. sol3rud;
Expand[%];
PotentialUToChi[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
ToCanonical[%];
PotentialVToU[%];
PotentialWToU[%];
PotentialChiToU[%];
PotentialToSource[%];
PotentialChiToU[%];
PotentialToSource[%];
TimeRhoToEuler[%];
TimeRhoToEuler[%];
TimeVelToEuler[%] /. sol2rud;
PotentialToSource[%];
Expand[%];
ToCanonical[%];
List @@ %;
SeparateMetric[] /@ %;
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
Plus @@ %;
ToCanonical[%];
Simplify[%];
Expand[%];
eqnsa4d = ScreenDollarIndices[%];


(* ::Text:: *)
(*The equations to solve for are the coefficients of the non-constant terms.*)


eqsd = Simplify[Coefficient[eqnsa4d, #]] & /@ terms;


(* ::Text:: *)
(*Check that we have fully decomposed the equations.*)


Simplify[eqsd . terms - eqnsa4d]


(* ::Text:: *)
(*To solve the equations, we must supply a gauge condition. The one we use here is satisfied in the standard PPN gauge.*)


sola4d = Simplify[Solve[Append[# == 0 & /@ eqsd, aa[11] == aa[15]], Prepend[a4var, aa[0]]][[1]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[eqnsa4d /. sola4d]


(* ::Text:: *)
(*Enhance the third order solution, using the gauge fixing condition determined at the fourth order.*)


sol3defd = ans3def /. Simplify[sola3d /. sola4d]
sol3rud = Flatten[mkrg /@ sol3defd];


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


sol4defd = ans4def /. sola4d
sol4rud = mkrg[sol4defd];


(* ::Text:: *)
(*Check that the metric components we found indeed solve the fourth order field equations.*)


eqn4d /. sol2rud /. sol3rud /. sol4rud;
Expand[%];
PotentialUToChi[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
ToCanonical[%];
PotentialVToU[%];
PotentialWToU[%];
PotentialChiToU[%];
PotentialToSource[%];
PotentialChiToU[%];
PotentialToSource[%];
TimeRhoToEuler[%];
TimeRhoToEuler[%];
TimeVelToEuler[%] /. sol2rud;
PotentialToSource[%];
Expand[%];
ToCanonical[%];
List @@ %;
SeparateMetric[] /@ %;
ToCanonical /@ %;
ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True] & /@ %;
Plus @@ %;
ToCanonical[%];
Simplify[%]


(* ::Text:: *)
(*The solution will in general not have the standard PPN form. These are the conditions which must be satisfied.*)


varconda4d = FullSimplify[varcond /. sola4d]


(* ::Text:: *)
(*There are two possible solutions.*)


parcondd = Simplify[Solve[varconda4d, k2]]


(* ::Text:: *)
(*In these cases, the metric takes the standard PPN form, and we can solve for the coefficients in the corresponding expansion.*)


varrel /. sola4d /. parcondd;
Simplify[%];
sola4pd = Simplify[Solve[#, a4var2][[1]] & /@ %]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solutions for the metric components.*)


sol4defpd = ans4def2 /. sola4pd
sol4rupd = mkrg /@ sol4defpd;


(* ::Text:: *)
(*Finally, we combine all solutions for each of the three branches.*)


soldefpd = Transpose[Append[Transpose[MapAt[Simplify, Join[sol2defd, sol3defd] /. parcondd, {All, All, 2}]], sol4defpd]]
solrupd = Flatten /@ Map[mkrg, soldefpd, {2}];


(* ::Chapter:: *)
(*PPN metric and parameters*)


(* ::Section:: *)
(*PPN metric*)


(* ::Text:: *)
(*To read off the PPN parameters, we use the following metric components.*)


metcomp = {PPN[Met,2][-LI[0],-LI[0]], PPN[Met,2][-T3a,-T3b], PPN[Met,3][-LI[0],-T3a], PPN[Met,4][-LI[0],-LI[0]]}


(* ::Text:: *)
(*Insert the solution we obtained into the metric components.*)


metcomp /. Join[solrup, solrupd];
Map[ToCanonical, %, {2}];
Map[Expand, %, {2}];
Map[Simplify, %, {2}];
ppnmets = %;


(* ::Text:: *)
(*For later comparison, convert the selected components also to their standard form in terms of PPN parameters and potentials.*)


stamet = Simplify[MetricToStandard /@ metcomp]


(* ::Section:: *)
(*Newtonian gravitational constant*)


(* ::Text:: *)
(*To solve for the gravitational constant, compare the second order result with the standard normalization.*)


kappaeq = (First[#] == First[stamet]&) /@ ppnmets


(* ::Text:: *)
(*When solving for \[Kappa], make sure to catch the positive root.*)


kappadef = (kappa == First[Sqrt[FullSimplify[kappa2 /. Solve[# /. kappa -> Sqrt[kappa2], kappa2]]]]&) /@ kappaeq
kapparu = mkrg /@ kappadef;


(* ::Section:: *)
(*PPN parameters*)


(* ::Text:: *)
(*To generate equations for the PPN parameters, compare the obtained solution with the standard PPN metric.*)


pareqns = Simplify[Map[ToCanonical, MapThread[stamet - ReplaceAll[#1, #2]&, {ppnmets, kapparu}, 1], {2}]]


(* ::Text:: *)
(*We will consider the coefficients in front of the following potentials.*)


pots = {PotentialU[] * BkgMetricS3[-T3a,-T3b], PotentialV[-T3a], PotentialW[-T3a], PotentialA[], PotentialU[]^2, PotentialPhiW[], PotentialPhi1[], PotentialPhi2[], PotentialPhi3[], PotentialPhi4[]}


(* ::Text:: *)
(*Extract the coefficients from the difference between our result and the standard PPN metric. These terms must vanish.*)


eqns = DeleteCases[Flatten[Simplify[Outer[Coefficient, #, pots]]], 0]& /@ pareqns;


(* ::Text:: *)
(*List of PPN parameters we are solving for.*)


pars = {ParameterBeta, ParameterGamma, ParameterXi, ParameterAlpha1, ParameterAlpha2, ParameterAlpha3, ParameterZeta1, ParameterZeta2, ParameterZeta3, ParameterZeta4}


(* ::Text:: *)
(*Finally, solve the equations and determine the PPN parameters.*)


TableForm[parsol = FullSimplify[Solve[#, pars][[1]]]& /@ Map[# == 0&, eqns, {2}]]


(* ::Text:: *)
(*These parameters show the derivation from general relativity.*)


TableForm[Factor[{ParameterBeta - 1, ParameterGamma - 1} /. parsol]]
