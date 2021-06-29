(* ::Package:: *)

(* ::Title:: *)
(*General scalar-torsion gravity*)


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


DefTensor[phi[], {MfSpacetime}, PrintAs -> "\[Phi]"];


(* ::Section:: *)
(*Background value of the scalar field*)


(* ::Text:: *)
(*For the background value of the scalar field we define a constant.*)


DefConstantSymbol[phi0, PrintAs -> "\[CapitalPhi]"];


(* ::Section:: *)
(*Gravitational constant*)


(* ::Text:: *)
(*Another constant parameter we need to define is the gravitational constant.*)


DefConstantSymbol[kappa, PrintAs -> "\[Kappa]"];


(* ::Section:: *)
(*Parameter function in the action*)


(* ::Text:: *)
(*The action is determined by a single function.*)


DefScalarFunction[ll, PrintAs -> "\[ScriptCapitalL]"];


(* ::Section:: *)
(*Field Equations*)


(* ::Text:: *)
(*These tensors will hold the gravitational (metric and scalar) field equations.*)


DefTensor[TetEq[-T4\[Alpha], -T4\[Beta]], {MfSpacetime}, PrintAs -> "\[ScriptCapitalE]"];
DefTensor[ScalEq[], {MfSpacetime}, PrintAs -> "\[ScriptCapitalE]"];


(* ::Section:: *)
(*Constant coefficients*)


(* ::Text:: *)
(*For convenience, we define some constants, which will be used in the ansatz for solving the field equations.*)


aa[i_] := Module[{sym = Symbol["a" <> ToString[i]]}, If[!ConstantSymbolQ[sym], DefConstantSymbol[sym, PrintAs -> StringJoin["\!\(a\_", ToString[i], "\)"]]]; Return[sym]]


(* ::Chapter:: *)
(*Shorthand notations*)


(* ::Section:: *)
(*Formatting of parameter function*)


(* ::Text:: *)
(*Omit the function parameters when printing.*)


Format[ll[_, _, _, phi[]]] := "\[ScriptCapitalL]";


(* ::Text:: *)
(*Derivatives of the parameter function are written with subscript letters.*)


Format[Derivative[ct_, cx_, cy_, c\[Phi]_][ll][_, _, _, phi[]]] := Subscript["\[ScriptCapitalL]", StringJoin[ConstantArray["\[Phi]", c\[Phi]], ConstantArray["T", ct], ConstantArray["X", cx], ConstantArray["Y", cy]]];


(* ::Text:: *)
(*For the zeroth Taylor coefficient, i.e., the cosmological background value, we use an upright letter.*)


Format[ll[0, 0, 0, phi0]] := "L";


(* ::Text:: *)
(*Further Taylor coefficients are again written with subscript letters.*)


Format[Derivative[ct_, cx_, cy_, c\[Phi]_][ll][0, 0, 0, phi0]] :=  Subscript["L", StringJoin[ConstantArray["\[Phi]", c\[Phi]], ConstantArray["T", ct], ConstantArray["X", cx], ConstantArray["Y", cy]]];


(* ::Section:: *)
(*Action terms*)


(* ::Text:: *)
(*Torsion scalar.*)


ts = Scalar[Expand[TorsionFD[T4\[Gamma], -T4\[Mu], -T4\[Nu]] * TorsionFD[T4\[Delta], -T4\[Rho], -T4\[Sigma]] * (
	Met[-T4\[Gamma], -T4\[Delta]] * InvMet[T4\[Mu], T4\[Rho]] * InvMet[T4\[Nu], T4\[Sigma]] / 4 +
	delta[-T4\[Gamma], T4\[Sigma]] * delta[-T4\[Delta], T4\[Nu]] * InvMet[T4\[Mu], T4\[Rho]] / 2 -
	delta[-T4\[Gamma], T4\[Nu]] * delta[-T4\[Delta], T4\[Sigma]] * InvMet[T4\[Mu], T4\[Rho]]
)]]


(* ::Text:: *)
(*Scalar field kinetic term.*)


xx = -Scalar[InvMet[T4\[Gamma], T4\[Delta]] * CD[-T4\[Gamma]][phi[]] * CD[-T4\[Delta]][phi[]]] / 2


(* ::Text:: *)
(*Scalar field coupling term.*)


yy = Scalar[InvMet[T4\[Gamma], T4\[Delta]] * TorsionFD[T4\[Alpha], -T4\[Alpha], -T4\[Gamma]] * CD[-T4\[Delta]][phi[]]]


(* ::Section:: *)
(*Parameter function and its derivatives*)


(* ::Text:: *)
(*For convenience when entering the field equations, we define a few abbreviations.*)


l0 = ll[ts, xx, yy, phi[]]
lt = Derivative[1, 0, 0, 0][ll][ts, xx, yy, phi[]]
lx = Derivative[0, 1, 0, 0][ll][ts, xx, yy, phi[]]
ly = Derivative[0, 0, 1, 0][ll][ts, xx, yy, phi[]]
lphi = Derivative[0, 0, 0, 1][ll][ts, xx, yy, phi[]]


(* ::Section:: *)
(*Taylor coefficients of parameter function*)


(* ::Text:: *)
(*We will use the following abbreviations for Taylor coefficients later when we come to solving the field equations.*)


ll0 = ll[0, 0, 0, phi0]
llphi = Derivative[0, 0, 0, 1][ll][0, 0, 0, phi0]
llphiphi = Derivative[0, 0, 0, 2][ll][0, 0, 0, phi0]
llphiphiphi = Derivative[0, 0, 0, 3][ll][0, 0, 0, phi0]


(* ::Chapter:: *)
(*Field equations*)


(* ::Section:: *)
(*Tetrad field equation*)


(* ::Text:: *)
(*We start by defining the tetrad field equation.*)


-Met[-T4\[Alpha], -T4\[Beta]] * (l0 + InvMet[T4\[Gamma], T4\[Delta]] * CD[-T4\[Gamma]][ CD[-T4\[Delta]][phi[]] * ly]) - 2 * CD[-T4\[Gamma]][lt * (TorsionFD[T4\[Delta], -T4\[Delta], -T4\[Alpha]] * delta[-T4\[Beta], T4\[Gamma]] - TorsionFD[T4\[Delta], -T4\[Delta], -T4\[CurlyEpsilon]] * InvMet[T4\[CurlyEpsilon], T4\[Gamma]] * Met[-T4\[Alpha], -T4\[Beta]] + (ChristoffelFD[T4\[Delta], -T4\[Beta], -T4\[CurlyEpsilon]] - ChristoffelCD[T4\[Delta], -T4\[Beta], -T4\[CurlyEpsilon]]) * Met[-T4\[Alpha], -T4\[Delta]] * InvMet[T4\[CurlyEpsilon], T4\[Gamma]])] - (TorsionFD[T4\[Zeta], -T4\[Zeta], -T4\[Delta]] * (TorsionFD[T4\[Delta], -T4\[Alpha], -T4\[Beta]] + (TorsionFD[T4\[Gamma], -T4\[Beta], -T4\[CurlyEpsilon]] * Met[-T4\[Gamma], -T4\[Alpha]] + TorsionFD[T4\[Gamma], -T4\[Alpha], -T4\[CurlyEpsilon]] * Met[-T4\[Gamma], -T4\[Beta]]) * InvMet[T4\[CurlyEpsilon], T4\[Delta]]) + Met[-T4\[Alpha], -T4\[Zeta]] * TorsionFD[T4\[Zeta], -T4\[Gamma], -T4\[Delta]] * InvMet[T4\[CurlyEpsilon], T4\[Delta]] * (TorsionFD[T4\[Gamma], -T4\[CurlyEpsilon], -T4\[Beta]] - TorsionFD[T4\[Mu], -T4\[Nu], -T4\[CurlyEpsilon]] * Met[-T4\[Mu], -T4\[Beta]] * InvMet[T4\[Nu], T4\[Gamma]] / 2)) * lt - PD[-T4\[Alpha]][phi[]] * PD[-T4\[Beta]][phi[]] * lx + CD[-T4\[Beta]][CD[-T4\[Alpha]][phi[]] * ly] + ((TorsionFD[T4\[Gamma], -T4\[Alpha], -T4\[Beta]] + InvMet[T4\[Gamma], T4\[Delta]] * (Met[-T4\[Alpha], -T4\[CurlyEpsilon]] * TorsionFD[T4\[CurlyEpsilon], -T4\[Beta], -T4\[Delta]] + Met[-T4\[Beta], -T4\[CurlyEpsilon]] * TorsionFD[T4\[CurlyEpsilon], -T4\[Alpha], -T4\[Delta]])) * CD[-T4\[Gamma]][phi[]] + 2 * TorsionFD[T4\[Gamma], -T4\[Gamma], -T4\[Alpha]] * CD[-T4\[Beta]][ phi[]]) * ly / 2 - 2 * kappa^2 * EnergyMomentum[-T4\[Alpha], -T4\[Beta]];


(* ::Text:: *)
(*Next, we trace reverse the equation.*)


% * (delta[-T4\[Tau], T4\[Alpha]] * delta[-T4\[Omega], T4\[Beta]] - InvMet[T4\[Alpha], T4\[Beta]] * Met[-T4\[Tau], -T4\[Omega]] / 2);


(* ::Text:: *)
(*Finally, save the equations for later use.*)


teteqdef = TetEq[-T4\[Tau], -T4\[Omega]] == %;
teteqru = mkr0[teteqdef];


(* ::Section:: *)
(*Scalar field equation*)


(* ::Text:: *)
(*We also need to define the scalar field equation.*)


-lphi + InvMet[T4\[Alpha], T4\[Beta]] * CD[-T4\[Alpha]][TorsionFD[T4\[Gamma], -T4\[Gamma], -T4\[Beta]] * ly - CD[-T4\[Beta]][phi[]] * lx];


(* ::Text:: *)
(*To simplify the equations we solve later, we use a linear combination of the scalar field equation and the trace of the tetrad field equation.*)


TetEq[-T4\[Alpha], -T4\[Beta]] * InvMet[T4\[Alpha], T4\[Beta]] * ly + 4 * % * lt /. teteqru;


(* ::Text:: *)
(*Finally, save the equation for later use.*)


scaleqdef = ScalEq[] == %;
scaleqru = mkr0[scaleqdef];


(* ::Section:: *)
(*Simplify equations*)


(* ::Text:: *)
(*Since using the full equations for the post-Newtonian expansion is a lengthy procedure, we eliminate any terms whose contribution to the field equations would be of higher than fourth velocity order, so that they do not contribute.*)


{TetEq[-T4\[Alpha], -T4\[Beta]], ScalEq[]} /. scaleqru /. teteqru;
ChangeCovD[#, CD, PD] & /@ %;
Expand /@ %;
ContractMetric /@ %;
ToCanonical[#, UseMetricOnVBundle -> None]& /@ %;
NoScalar /@ %;
Select[
	Count[#, TorsionFD[__]] +
	Count[#, ChristoffelFD[__]] +
	Count[#, ChristoffelCD[__]] +
	Count[#, PD[_][TorsionFD[__]]] +
	Count[#, PD[_][ChristoffelFD[__]]] +
	Count[#, PD[_][ChristoffelCD[__]]] +
	Count[#, PD[_][phi[]]] +
	Count[#, PD[_][PD[_][phi[]]]] <= 2&
] /@ %;
PutScalar /@ %;
Simplify /@ %
{teteqdef, scaleqdef} = {TetEq[-T4\[Alpha], -T4\[Beta]] == %[[1]], ScalEq[] == %[[2]]};
teteqru = mkr0[teteqdef];
scaleqru = mkr0[scaleqdef];


(* ::Chapter:: *)
(*Post-Newtonian expansion*)


(* ::Section:: *)
(*Rules for the scalar field*)


(* ::Text:: *)
(*When the scalar field is expanded in velocity orders, set the zeroth order to the constant background value, and let the first and third orders vanish.*)


OrderSet[PPN[phi, 0][], phi0];
OrderSet[PPN[phi, 1][], 0];
OrderSet[PPN[phi, 3][], 0];


(* ::Section:: *)
(*3 + 1 split*)


(* ::Text:: *)
(*Split the tetrad equations into their time and space components.*)


{#, # /. teteqru}&[TetEq[-T4\[Alpha], -T4\[Beta]]];
ChangeCovD[%, CD, PD];
Expand[%];
NoScalar /@ %;
SpaceTimeSplits[#, {-T4\[Alpha] -> -T3a, -T4\[Beta] -> -T3b}]& /@ %;
Expand[%];
Map[ToCanonical, %, {3}];
Map[SortPDs, %, {3}];
teteq31list = %;
teteq31def = Union[Flatten[MapThread[Equal, %, 2]]];
teteq31ru = Flatten[mkrg /@ %];


(* ::Text:: *)
(*Split the scalar equation into its time and space components.*)


{#, # /. scaleqru}&[ScalEq[]];
ChangeCovD[%, CD, PD];
Expand[%];
NoScalar /@ %;
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
(*Expand each component of the tetrad equations into the respective velocity orders.*)


Outer[VelocityOrder, teteq31list, Range[0, 4]];
Map[NoScalar, %, {4}];
Expand[%];
Map[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]&, %, {4}];
Map[ToCanonical, %, {4}];
Map[SortPDs, %, {4}];
teteqvlist = Simplify[%];
teteqvdef = Union[Flatten[MapThread[Equal, %, 3]]]
teteqvru = Flatten[mkrg /@ %];


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
(*Zeroth order tetrad equations*)


(* ::Text:: *)
(*These are the zeroth order tetrad equations.*)


tetppneqs0 = {PPN[TetEq, 0][-LI[0], -LI[0]], PPN[TetEq, 0][-T3a, -T3b]} /. teteqvru


(* ::Text:: *)
(*They are solved by assuming that the parameter function vanishes at the background.*)


tetsol0 = {ll0 -> 0}
tetppneqs0 /. tetsol0


(* ::Section:: *)
(*Zeroth order scalar equation*)


(* ::Text:: *)
(*This is the zeroth order scalar equation.*)


scalppneq0 = PPN[ScalEq, 0][] /. scaleqvru /. tetsol0


(* ::Text:: *)
(*We set the derivative with respect to the scalar field to vanish; otherwise gravity will disappear.*)


scalsol0 = {llphi -> 0}
scalppneq0 /. scalsol0


(* ::Section:: *)
(*Second order scalar field*)


(* ::Text:: *)
(*In the following calculation, we will assume a massless scalar field.*)


massless = {llphiphi -> 0, llphiphiphi -> 0}


(* ::Text:: *)
(*This is the second order scalar equation.*)


scalppneq2 = FullSimplify[ToCanonical[PPN[ScalEq, 2][] /. scaleqvru /. tetsol0 /. scalsol0 /. massless]]


(* ::Text:: *)
(*Define an ansatz for the second order scalar perturbations.*)


scala2def = PPN[phi, 2][] == aa[1] * PotentialU[]
scala2ru = mkrg[scala2def];


(* ::Text:: *)
(*Insert the ansatz into the scalar field equation and convert derivatives of the PPN potentials into matter source terms.*)


scalppneq2 /. scala2ru;
PotentialToSource[%];
scaleq2 = %


(* ::Text:: *)
(*The equation to be solved is extracted as coefficient of the density.*)


scaleqa2 = Coefficient[scaleq2, Density[]]


(* ::Text:: *)
(*Solve for the constant coefficient in the equation.*)


scalasol2 = Solve[scaleqa2 == 0, aa[1]][[1]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equation.*)


Simplify[scaleqa2 /. scalasol2]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the scalar field.*)


scalsol2 = mkrg[scala2def /. scalasol2]


(* ::Text:: *)
(*Check that the scalar field we found indeed solves the second order field equations.*)


scalppneq2 /. scalsol2 // PotentialToSource // Simplify


(* ::Section:: *)
(*Second order tetrad - time part*)


(* ::Text:: *)
(*This is the second order tetrad equation.*)


tetppneq200 = PPN[TetEq, 2][-LI[0], -LI[0]] /. teteqvru /. tetsol0 /. scalsol0 /. scalsol2 // PotentialToSource // ToCanonical // Simplify


(* ::Text:: *)
(*Define an ansatz for the second order tetrad perturbations.*)


teta200def = PPN[Met, 2][-LI[0], -LI[0]] == aa[2] * PotentialU[]
teta200ru = mkrg[teta200def];


(* ::Text:: *)
(*Insert the ansatz into the field equation and convert derivatives of the PPN potentials into matter source terms.*)


tetppneq200 /. teta200ru;
PotentialToSource[%];
teteq200 = %


(* ::Text:: *)
(*The equation to be solved is extracted as coefficient of the density.*)


teteqa200 = Coefficient[teteq200, Density[]]


(* ::Text:: *)
(*Solve for the constant coefficient in the equation.*)


tetasol200 = Solve[teteqa200 == 0, aa[2]][[1]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equation.*)


Simplify[teteqa200 /. tetasol200]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the tetrad component.*)


tetsol200 = mkrg[teta200def /. tetasol200]


(* ::Text:: *)
(*Check that the tetrad we found indeed solves the second order field equations.*)


tetppneq200 /. tetsol200 // PotentialToSource // Simplify


(* ::Section:: *)
(*Second order tetrad - space part*)


(* ::Text:: *)
(*This is the second order tetrad equation.*)


tetppneq2ij = PPN[TetEq, 2][-T3a, -T3b] /. teteqvru /. tetsol0 /. scalsol0 /. tetsol200 /. scalsol2 // PotentialToSource // Expand // ToCanonical // SortPDs // Simplify


(* ::Text:: *)
(*Define an ansatz for the second order tetrad perturbations.*)


teta2ijdef = PPN[Met, 2][-T3a, -T3b] == aa[3] * PotentialU[] * BkgMetricS3[-T3a, -T3b]
teta2ijru = mkrg[teta2ijdef];


(* ::Text:: *)
(*Insert the ansatz into the field equation and convert derivatives of the PPN potentials into matter source terms.*)


tetppneq2ij /. teta2ijru;
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
teteq2ij = %


(* ::Text:: *)
(*The equation to be solved is extracted as coefficient of the density.*)


teteqa2ij = Coefficient[teteq2ij, Density[] * BkgMetricS3[-T3a, -T3b]]


(* ::Text:: *)
(*Solve for the constant coefficient in the equation.*)


tetasol2ij = Solve[teteqa2ij == 0, aa[3]][[1]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equation.*)


Simplify[teteqa2ij /. tetasol2ij]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the tetrad component.*)


tetsol2ij = mkrg[teta2ijdef /. tetasol2ij]


(* ::Text:: *)
(*Check that the tetrad we found indeed solves the second order field equations.*)


tetppneq2ij /. tetsol2ij // PotentialToSource // ToCanonical // Simplify


(* ::Section:: *)
(*Third order*)


(* ::Text:: *)
(*The antisymmetric part of the third order equations is already solved.*)


PPN[TetEq, 3][-T3a, -LI[0]] - PPN[TetEq, 3][-LI[0], -T3a] /. teteqvru // Expand // ToCanonical // SortPDs // Simplify


(* ::Text:: *)
(*Extract the symmetric third order field equations.*)


tetppneq3sym =  PPN[TetEq, 3][-T3a, -LI[0]] + PPN[TetEq, 3][-LI[0], -T3a] /. teteqvru /. tetsol0 /. scalsol2 /. tetsol200 /. tetsol2ij // Expand // ToCanonical // SortPDs // Simplify


(* ::Text:: *)
(*Define an ansatz for the third order tetrad perturbations.*)


teta3def = {PPN[Met, 3][-LI[0], -T3a] == aa[4] * PotentialV[-T3a] + aa[5] * PotentialW[-T3a], PPN[Asym, 3][-LI[0], -T3a] == 0}
teta3ru = Join @@ (mkrg /@ teta3def);


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


tetppneq3sym /. teta3ru;
PotentialWToChiV[%];
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialChiToU[%];
PotentialVToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
teteq3 = Simplify[%]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations, keeping a gauge freedom, which is left up to the fourth velocity order.*)


tetasol3 = FullSimplify[First[Solve[{teteq3 == 0, aa[5] - aa[4] == aa[0]}, {aa[4], aa[5]}]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[teteq3 /. tetasol3]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the metric components.*)


tetsol3 = Join @@ (mkrg /@ teta3def /. tetasol3)


(* ::Text:: *)
(*Check that the tetrad components we found indeed solve the third order field equations.*)


tetppneq3sym /. tetsol3 /. tetsol200;
PotentialWToChiV[%];
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
(*We eliminate the spatial tetrad component by using a suitable linear combination of the tetrad and scalar field equations.*)


tetc = Coefficient[ScreenDollarIndices[PPN[TetEq, 4][-LI[0], -LI[0]] /. teteqvru], PD[-T3a][PD[T3a][PPN[phi, 4][]]]]
scalc = Coefficient[ScreenDollarIndices[PPN[ScalEq, 4][] /. scaleqvru], PD[-T3a][PD[T3a][PPN[phi, 4][]]]]


(* ::Text:: *)
(*Extract the fourth order field equations.*)


PPN[TetEq, 4][-LI[0], -LI[0]] * scalc - PPN[ScalEq, 4][] * tetc;
% /. teteqvru /. scaleqvru /. massless /. tetsol0 /. scalsol0 /. scalsol2 /. tetsol200 /. tetsol2ij /. tetsol3;
Expand[%];
PotentialWToU[%];
PotentialVToU[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
tetppneq400 = %


(* ::Text:: *)
(*Define an ansatz for the fourth order metric perturbations.*)


teta400def = {PPN[Met, 4][-LI[0], -LI[0]] == aa[6] * PotentialU[]^2 + aa[7] * PotentialPhi1[] + aa[8] * PotentialPhi2[] + aa[9] * PotentialPhi3[] + aa[10] * PotentialPhi4[], PPN[Asym, 2][-T3a, -T3b] == 0}
teta400ru = Flatten[mkrg /@ teta400def];


(* ::Text:: *)
(*Insert the ansatz into the field equations and convert derivatives of the PPN potentials into matter source terms.*)


tetppneq400 /. teta400ru;
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialVToU[%];
PotentialWToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%];
ScreenDollarIndices[%];
teteq4 = Simplify[%]


(* ::Text:: *)
(*Coefficient of the pressure.*)


eq1 = Simplify[Coefficient[teteq4, Pressure[]]]


(* ::Text:: *)
(*Coefficient of the internal energy.*)


eq2 = Simplify[Coefficient[teteq4, Density[] * InternalEnergy[]]]


(* ::Text:: *)
(*Coefficient of the gravitational potential energy.*)


eq3 = Simplify[Coefficient[teteq4, Density[] * PotentialU[]]]


(* ::Text:: *)
(*Coefficient of the second time derivative.*)


eq4 = Simplify[Coefficient[teteq4, ParamD[TimePar, TimePar][PotentialU[]]]]


(* ::Text:: *)
(*Coefficient of the kinetic energy.*)


eq5 = Simplify[Coefficient[teteq4, Density[] * Velocity[-T3a] * Velocity[T3a]]]


(* ::Text:: *)
(*Coefficient of the potential term.*)


eq6 = Simplify[Coefficient[teteq4, PD[-T3a][PotentialU[]] * PD[T3a][PotentialU[]]]]


(* ::Text:: *)
(*Check that we have fully decomposed the equations.*)


Simplify[Pressure[] * eq1 + Density[] * InternalEnergy[] * eq2 + Density[] * PotentialU[] * eq3 + ParamD[TimePar, TimePar][PotentialU[]] * eq4 + Density[] * Velocity[-T3a] * Velocity[T3a] * eq5 + PD[-T3a][PotentialU[]] * PD[T3a][PotentialU[]] * eq6 - teteq4]


(* ::Text:: *)
(*Solve for the constant coefficients in the equations.*)


tetasol400 = Simplify[First[Solve[# == 0& /@ {eq1, eq2, eq3, eq4, eq5, eq6}, aa /@ Prepend[Range[6, 10], 0]]]]


(* ::Text:: *)
(*Check that the solution indeed solves the component equations.*)


Simplify[teteq4 /. tetasol400]


(* ::Text:: *)
(*Enhance the third order solution, using the gauge fixing condition determined at the fourth order.*)


tetsol3 =  Simplify[Join @@ (mkrg /@ teta3def /. tetasol3 /. tetasol400)]


(* ::Text:: *)
(*Insert the solution for the constant coefficients into the ansatz, to obtain the solution for the tetrad components.*)


tetsol400 = Flatten[mkrg /@ (teta400def /. tetasol400)]


(* ::Text:: *)
(*Check that the tetrad components we found indeed solve the fourth order field equations.*)


tetppneq400 /. tetsol400 /. tetasol400;
Expand[%];
ContractMetric[%, OverDerivatives -> True, AllowUpperDerivatives -> True];
PotentialVToU[%];
PotentialWToU[%];
PotentialToSource[%];
ToCanonical[%];
SortPDs[%]; Simplify[%]


(* ::Chapter:: *)
(*PPN metric and parameters*)


(* ::Section:: *)
(*PPN metric*)


(* ::Text:: *)
(*To read off the PPN parameters, we use the following metric components.*)


metcomp = {PPN[Met,2][-LI[0],-LI[0]], PPN[Met,2][-T3a,-T3b], PPN[Met,3][-LI[0],-T3a], PPN[Met,4][-LI[0],-LI[0]]}


(* ::Text:: *)
(*Insert the solution we obtained into the metric components.*)


metcomp /. tetsol200 /. tetsol2ij /. tetsol3 /. tetsol400;
ToCanonical[%];
Expand[%];
Simplify[%];
Expand[%];
metlist = Simplify[%];
metdef = MapThread[Equal, {metcomp, metlist}, 1]
metru = Join @@ (mkrg /@ metdef);


(* ::Text:: *)
(*For later comparison, convert the selected components also to their standard form in terms of PPN parameters and potentials.*)


stamet = Simplify[MetricToStandard /@ metcomp]


(* ::Section:: *)
(*Newtonian gravitational constant*)


(* ::Text:: *)
(*To solve for the gravitational constant, compare the second order result with the standard normalization.*)


kappaeq = First[metlist] == First[stamet]


(* ::Text:: *)
(*When solving for \[Kappa], make sure to catch the positive root.*)


kappadef = kappa == First[Sqrt[FullSimplify[k2 /. Solve[kappaeq /. kappa -> Sqrt[k2], k2]]]]
kapparu = mkrg[kappadef];


(* ::Text:: *)
(*Normalize metric components.*)


metlistn = Simplify[Expand[metlist /. kapparu]];
metdefn = MapThread[Equal, {metcomp, metlistn}, 1]
metrun = Join @@ (mkrg /@ metdefn);


(* ::Section:: *)
(*PPN parameters*)


(* ::Text:: *)
(*To generate equations for the PPN parameters, compare the obtained solution with the standard PPN metric.*)


pareqns = Simplify[ToCanonical[stamet - metlistn]]


(* ::Text:: *)
(*We will consider the coefficients in front of the following potentials.*)


pots = {PotentialU[] * BkgMetricS3[-T3a,-T3b], PotentialV[-T3a], PotentialW[-T3a], PotentialA[], PotentialU[]^2, PotentialPhiW[], PotentialPhi1[], PotentialPhi2[], PotentialPhi3[], PotentialPhi4[]}


(* ::Text:: *)
(*Extract the coefficients from the difference between our result and the standard PPN metric. These terms must vanish.*)


eqns = DeleteCases[Flatten[Simplify[Outer[Coefficient, pareqns, pots]]], 0]


(* ::Text:: *)
(*List of PPN parameters we are solving for.*)


pars = {ParameterBeta, ParameterGamma, ParameterXi, ParameterAlpha1, ParameterAlpha2, ParameterAlpha3, ParameterZeta1, ParameterZeta2, ParameterZeta3, ParameterZeta4}


(* ::Text:: *)
(*Finally, solve the equations and determine the PPN parameters.*)


parsol = FullSimplify[Solve[# == 0& /@ eqns, pars][[1]]]


(* ::Text:: *)
(*These parameters show the derivation from general relativity.*)


grdev = Factor[{ParameterBeta - 1, ParameterGamma - 1} /. parsol]


(* ::Text:: *)
(*Note in particular that the deviation vanishes for minimal coupling.*)


grdev /. Derivative[0, 0, 1, 0][ll][0, 0, 0, phi0] -> 0


(* ::Chapter:: *)
(*Simple scalar-torsion gravity*)


(* ::Section:: *)
(*Parameter functions*)


(* ::Text:: *)
(*The simple action depends on the following functions.*)


DefScalarFunction[pfa, PrintAs -> "\[ScriptCapitalA]"]
DefScalarFunction[pfb, PrintAs -> "\[ScriptCapitalB]"]
DefScalarFunction[pfc, PrintAs -> "\[ScriptCapitalC]"]


(* ::Text:: *)
(*We omit the argument in display.*)


Format[pfa[phi[]]] := "\[ScriptCapitalA]"
Format[pfb[phi[]]] := "\[ScriptCapitalB]"
Format[pfc[phi[]]] := "\[ScriptCapitalC]"


(* ::Text:: *)
(*Background values are denotes with upright letters.*)


Format[pfa[phi0]] := "A"
Format[pfb[phi0]] := "B"
Format[pfc[phi0]] := "C"


(* ::Text:: *)
(*Derivatives are denoted with primes.*)


Format[Derivative[c\[Phi]_][pfa][phi0]] := Superscript["A", StringJoin @@ ConstantArray["'", c\[Phi]]]
Format[Derivative[c\[Phi]_][pfb][phi0]] := Superscript["B", StringJoin @@ ConstantArray["'", c\[Phi]]]
Format[Derivative[c\[Phi]_][pfc][phi0]] := Superscript["C", StringJoin @@ ConstantArray["'", c\[Phi]]]


(* ::Text:: *)
(*This replaces the general action with the specific choice we consider now.*)


abcrep = {ll -> ((-pfa[#4] * #1 + 2 * pfb[#4] * #2 + 2 * pfc[#4] * #3) / 2 /\[Kappa]^2 &)}


(* ::Section:: *)
(*PPN parameters*)


(* ::Text:: *)
(*These are the PPN parameters.*)


parsolabc = parsol /. abcrep // Expand // Simplify


(* ::Text:: *)
(*These parameters show the derivation from general relativity.*)


grdevabc = Factor[{ParameterBeta - 1, ParameterGamma - 1} /. parsolabc]


(* ::Text:: *)
(*Note in particular that the deviation vanishes for minimal coupling.*)


grdevabc /. pfc[phi0] -> 0
