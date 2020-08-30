(* ::Package:: *)

(* ::Title:: *)
(*PPN calculation for scalar-tensor gravity in xAct*)


(* ::Chapter:: *)
(*Preliminaries*)


(* ::Section:: *)
(*Load package*)


<< xAct`xPPN`


(* ::Section:: *)
(*Nicer printing*)


$PrePrint = ScreenDollarIndices;


(* ::Section:: *)
(*Make rules from equations*)


mkrg[eq_Equal] := MakeRule[Evaluate[List @@ eq], MetricOn -> All, ContractMetrics -> True];
mkr0[eq_Equal] := MakeRule[Evaluate[List @@ eq], MetricOn -> None, ContractMetrics -> False];


(* ::Chapter:: *)
(*Object definitions*)


(* ::Section:: *)
(*Scalar field*)


DefTensor[psi[], {MfSpacetime}, PrintAs -> "\[Psi]"];


(* ::Section:: *)
(*Background value of the scalar field*)


DefConstantSymbol[psi0, PrintAs -> "\[CapitalPsi]"];


(* ::Section:: *)
(*Gravitational constant*)


DefConstantSymbol[kappa, PrintAs -> "\[Kappa]"];


(* ::Section:: *)
(*Parameter function in the action*)


DefScalarFunction[omega, PrintAs -> "\[Omega]"];


(* ::Section:: *)
(*Field Equations*)


DefTensor[MetEq[-T4\[Alpha], -T4\[Beta]], {MfSpacetime}, Symmetric[{1,2}], PrintAs -> "\[ScriptCapitalE]"];
DefTensor[ScalEq[], {MfSpacetime}, PrintAs -> "\[ScriptCapitalE]"];


(* ::Section:: *)
(*Constant coefficients*)


aa[i_] := Module[{sym = Symbol["a" <> ToString[i]]}, If[!ConstantSymbolQ[sym], DefConstantSymbol[sym, PrintAs -> StringJoin["\!\(a\_", ToString[i], "\)"]]]; Return[sym]]


(* ::Chapter:: *)
(*Field equations*)


(* ::Section:: *)
(*Metric field equation*)


psi[] * RicciCD[-T4\[Alpha], -T4\[Beta]] - CD[-T4\[Alpha]][CD[-T4\[Beta]][psi[]]] - PD[-T4\[Alpha]][psi[]] * PD[-T4\[Beta]][psi[]] * omega[psi[]] / psi[] + omega'[psi[]] * InvMet[T4\[Gamma], T4\[Delta]] * PD[-T4\[Gamma]][psi[]] * PD[-T4\[Delta]][psi[]] * Met[-T4\[Alpha], -T4\[Beta]] / (4 * omega[psi[]] + 6) - kappa^2 * (EnergyMomentum[-T4\[Alpha], -T4\[Beta]] - InvMet[T4\[Gamma], T4\[Delta]] * EnergyMomentum[-T4\[Gamma], -T4\[Delta]] * Met[-T4\[Alpha], -T4\[Beta]] * (omega[psi[]] + 1) / (2 * omega[psi[]] + 3))


meteqdef = MetEq[-T4\[Alpha], -T4\[Beta]] == %;
meteqru = mkr0[meteqdef];


(* ::Section:: *)
(*Scalar field equation*)


(2 * omega[psi[]] + 3) * InvMet[T4\[Alpha], T4\[Beta]] * CD[-T4\[Alpha]][CD[-T4\[Beta]][psi[]]] + omega'[psi[]] * InvMet[T4\[Alpha], T4\[Beta]] * PD[-T4\[Alpha]][psi[]] * PD[-T4\[Beta]][psi[]] - kappa^2 * InvMet[T4\[Alpha], T4\[Beta]] * EnergyMomentum[-T4\[Alpha], -T4\[Beta]]


scaleqdef = ScalEq[] == %;
scaleqru = mkr0[scaleqdef];


(* ::Chapter:: *)
(*Post-Newtonian expansion*)


(* ::Section:: *)
(*Rules for the scalar field*)


OrderSet[PPN[psi,0][], psi0];
OrderSet[PPN[psi,1][], 0];
OrderSet[PPN[psi,3][], 0];


(* ::Section:: *)
(*3 + 1 split*)


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


Outer[VelocityOrder, meteq31list, Range[0, 4]];
Map[NoScalar, %, {4}];
Expand[%];
Map[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]&, %, {4}];
Map[ToCanonical, %, {4}];
Map[SortPDs, %, {4}];
meteqvlist = Simplify[%];
meteqvdef = Union[Flatten[MapThread[Equal, %, 3]]]
meteqvru = Flatten[mkrg /@ %];


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


eqns0 = {PPN[MetEq,0][-LI[0],-LI[0]], PPN[MetEq,0][-T3a,-T3b], PPN[ScalEq,0][]} /. meteqvru /. scaleqvru


(* ::Section:: *)
(*Second order*)


eqns2 = FullSimplify[{PPN[MetEq,2][-LI[0],-LI[0]], PPN[MetEq,2][-T3a,-T3b], PPN[ScalEq,2][]} /. meteqvru /. scaleqvru]


ans2def = {PPN[Met,2][-LI[0],-LI[0]] == aa[1] * PotentialU[], PPN[Met,2][-T3a,-T3b] == aa[2] * PotentialU[] * BkgMetricS3[-T3a,-T3b] + aa[3] * PotentialUU[-T3a,-T3b], PPN[psi,2][] == aa[4] * PotentialU[]}
ans2ru = Flatten[mkrg /@ ans2def];


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


eqnsc2 = FullSimplify[{Coefficient[eqnsa2[[1]], Density[]], Coefficient[eqnsa2[[2]], Density[] * BkgMetricS3[-T3a,-T3b]], Coefficient[eqnsa2[[3]], Density[]], aa[3]}]


sola2 = FullSimplify[First[Solve[# == 0& /@ eqnsc2, aa /@ Range[1, 4]]]]


Simplify[eqnsa2 /. sola2]


sol2def = ans2def /. sola2
sol2ru = Flatten[mkrg /@ sol2def];


eqns2/.sol2ru;
Expand[%];
PotentialToSource /@ %;
ToCanonical /@ %;
SortPDs /@ %;
Simplify[%]


(* ::Section:: *)
(*Third order*)


eqns3 = FullSimplify[PPN[MetEq, 3][-LI[0], -T3a] /. meteqvru]


ans3def = PPN[Met,3][-LI[0], -T3a] == aa[5] * PotentialV[-T3a] + aa[6] * PotentialW[-T3a]
ans3ru = mkrg[ans3def];


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


sola3 = FullSimplify[First[Solve[{eqnsa3 == 0, aa[6] - aa[5] == aa[0]}, {aa[5], aa[6]}]]]


Simplify[eqnsa3 /. sola3]


sol3def = ans3def /. sola3
sol3ru = mkrg[sol3def];


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


eqns4 = PPN[MetEq,4][-LI[0], -LI[0]] /. meteqvru


Simplify[Together[eqns4]]


ans4def = PPN[Met,4][-LI[0], -LI[0]] == aa[7] * PotentialPhi1[] + aa[8] * PotentialPhi2[] + aa[9] * PotentialPhi3[] + aa[10] * PotentialPhi4[] + aa[11] * PotentialU[]^2
ans4ru = mkrg[ans4def];


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


eq1 = Simplify[Coefficient[eqnsa4, Pressure[]]]


eq2 = Simplify[Coefficient[eqnsa4, Density[] * InternalEnergy[]]]


eq3 = Simplify[Coefficient[eqnsa4, Density[] * PotentialU[]]]


eq4 = Simplify[Coefficient[eqnsa4, ParamD[TimePar, TimePar][PotentialU[]]]]


eq5 = Simplify[Coefficient[eqnsa4, Density[] * Velocity[-T3a] * Velocity[T3a]]]


eq6 = Simplify[Coefficient[eqnsa4, PD[-T3a][PotentialU[]] * PD[T3a][PotentialU[]]]]


Simplify[Pressure[] * eq1 + Density[] * InternalEnergy[] * eq2 + Density[] * PotentialU[] * eq3 + ParamD[TimePar, TimePar][PotentialU[]] * eq4 + Density[] * Velocity[-T3a] * Velocity[T3a] * eq5 + PD[-T3a][PotentialU[]] * PD[T3a][PotentialU[]] * eq6 - eqnsa4]


sola4 = Simplify[First[Solve[# == 0& /@ {eq1, eq2, eq3, eq4, eq5, eq6}, aa /@ Prepend[Range[7, 11], 0]]]]


Simplify[eqnsa4 /. sola4]


sol3def = ans3def /. Simplify[sola3 /. sola4]
sol3ru = mkrg[sol3def];


sol4def = ans4def /. sola4
sol4ru = mkrg[sol4def];


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


metcomp = {PPN[Met,2][-LI[0],-LI[0]], PPN[Met,2][-T3a,-T3b], PPN[Met,3][-LI[0],-T3a], PPN[Met,4][-LI[0],-LI[0]]}


metcomp /. sol2ru /. sol3ru /. sol4ru;
ToCanonical[%];
Expand[%];
ppnmet = Simplify[%];
metdef = MapThread[Equal, {metcomp, %}, 1]


stamet = Simplify[MetricToStandard /@ metcomp]


(* ::Section:: *)
(*Newtonian gravitational constant*)


kappaeq = First[ppnmet] == First[stamet]


kappadef = kappa == First[Sqrt[FullSimplify[k2 /. Solve[kappaeq /. kappa -> Sqrt[k2], k2]]]]
kapparu = mkrg[kappadef];


(* ::Section:: *)
(*PPN parameters*)


pareqns = Simplify[ToCanonical[stamet - ppnmet /. kapparu]]


pots = {PotentialU[] * BkgMetricS3[-T3a,-T3b], PotentialV[-T3a], PotentialW[-T3a], PotentialA[], PotentialU[]^2, PotentialPhiW[], PotentialPhi1[], PotentialPhi2[], PotentialPhi3[], PotentialPhi4[]}


eqns = DeleteCases[Flatten[Simplify[Outer[Coefficient, pareqns, pots]]],0]


pars = {ParameterBeta, ParameterGamma, ParameterXi, ParameterAlpha1, ParameterAlpha2, ParameterAlpha3, ParameterZeta1, ParameterZeta2, ParameterZeta3, ParameterZeta4}


parsol = FullSimplify[Solve[# == 0& /@ eqns, pars][[1]]]
