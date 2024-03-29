BeginPackage["xAct`xPPN`", {"xAct`xTensor`", "xAct`xPerm`", "xAct`xCore`"}]

MfSpacetime::usage = "MfSpacetime is the spacetime manifold \!\(M\_4\).";
MfSpace::usage = "MfSpace is the space manifold \!\(S\_3\).";
MfTime::usage = "MfTime is the time manifold \!\(T\_1\).";

LorentzMfSpacetime::usage = "LorentzMfSpacetime is the Lorentz bundle \!\(\[DoubleStruckCapitalL]M\_4\) over the spacetime manifold \!\(M\_4\).";
LorentzMfSpace::usage = "LorentzMfSpace is the Lorentz bundle \!\(\[DoubleStruckCapitalL]S\_3\) over the space manifold \!\(S\_3\).";
LorentzMfTime::usage = "LorentzMfTime is the Lorentz bundle \!\(\[DoubleStruckCapitalL]T\_1\) over the time manifold \!\(T\_1\).";

\[ScriptT]::usage = "\[ScriptT] is the generic index on the tangent bundle of the time manifold.";
\[ScriptCapitalT]::usage = "\[ScriptCapitalT] is the generic index on the Lorentz bundle of the time manifold.";

IndexNamesT3 = FromCharacterCode /@ Range[97, 122];
IndexNamesT4 = FromCharacterCode /@ Join[Range[945, 969], {977, 981, 982, 1008, 1009, 1013}];
IndexNamesL3 = FromCharacterCode /@ Range[65, 90];
IndexNamesL4 = FromCharacterCode /@ Join[Range[913, 929], Range[931, 937]];

IndicesT3 = (Symbol[StringJoin["T3", #]]&) /@ IndexNamesT3;
IndicesT4 = (Symbol[StringJoin["T4", #]]&) /@ IndexNamesT4;
IndicesL3 = (Symbol[StringJoin["L3", #]]&) /@ IndexNamesL3;
IndicesL4 = (Symbol[StringJoin["L4", #]]&) /@ IndexNamesL4;

(MessageName[#, "usage"] = StringJoin[ToString[#], " is an index on the tangent bundle of the space manifold."])& /@ IndicesT3;
(MessageName[#, "usage"] = StringJoin[ToString[#], " is an index on the tangent bundle of the spacetime manifold."])& /@ IndicesT4;
(MessageName[#, "usage"] = StringJoin[ToString[#], " is an index on the Lorentz bundle of the space manifold."])& /@ IndicesL3;
(MessageName[#, "usage"] = StringJoin[ToString[#], " is an index on the Lorentz bundle of the spacetime manifold."])& /@ IndicesL4;

TimePar::usage = "TimePar is the generic time parameter on which split objects depend.";

BkgMetricM4::usage = "BkgMetricM4[-\[Mu], -\[Nu]] is the background metric \!\(\[Eta]\_\(\[Mu]\[Nu]\)\) on the spacetime manifold \!\(M\_4\).";
BkgMetricS3::usage = "BkgMetricS3[-i, -j] is the background metric \!\(\[Delta]\_\(ij\)\) on the space manifold \!\(S\_3\).";
BkgMetricT1::usage = "BkgMetricT1[-\[ScriptT]1, -\[ScriptT]2] is the background metric \!\(\[Eta]\_\(00\) = -1\) on the time manifold \!\(T\_1\).";

BkgTetradM4::usage = "BkgTetradM4[\[CapitalGamma], -\[Mu]] is the background tetrad \!\(\(\[CapitalDelta]\^\[CapitalGamma]\)\_\[Mu]\) on the spacetime manifold \!\(M\_4\).";
BkgTetradS3::usage = "BkgTetradS3[A, -i] is the background tetrad \!\(\(\[CapitalDelta]\^A\)\_i\) on the space manifold \!\(S\_3\).";
BkgTetradT1::usage = "BkgTetradT1[\[ScriptCapitalT], -\[ScriptT]] is the background tetrad \!\(\(\[CapitalDelta]\^0\)\_0 = 1\) on the time manifold \!\(T\_1\).";

BkgInvTetradM4::usage = "BkgInvTetradM4[-\[CapitalGamma], \[Mu]] is the inverse background tetrad \!\(\(\[CapitalDelta]\_\[CapitalGamma]\)\^\[Mu]\) on the spacetime manifold \!\(M\_4\).";
BkgInvTetradS3::usage = "BkgInvTetradS3[-A, i] is the inverse background tetrad \!\(\(\[CapitalDelta]\_A\)\^i\) on the space manifold \!\(S\_3\).";
BkgInvTetradT1::usage = "BkgInvTetradT1[-\[ScriptCapitalT], \[ScriptT]] is the inverse background tetrad \!\(\(\[CapitalDelta]\_0\)\^0 = 1\) on the time manifold \!\(T\_1\).";

EnergyMomentum::usage = "EnergyMomentum[-\[Mu], -\[Nu]] is the energy-momentum tensor \!\(\[CapitalTheta]\_\(\[Mu]\[Nu]\)\).";
TREnergyMomentum::usage = "TREnergyMomentum[-\[Mu], -\[Nu]] is the trace-reversed energy-momentum tensor \!\(\(\[CapitalTheta]\&_\)\_\(\[Mu]\[Nu]\)\).";

Density::usage = "Density[] is the rest energy density \[Rho].";
Pressure::usage = "Pressure[] is the pressure p.";
InternalEnergy::usage = "InternalEnergy[] is the specific internal energy \[CapitalPi].";
Velocity::usage = "Velocity[i] is the source matter velocity \!\(v\^i\).";

PotentialChi::usage = "PotentialChi[] is the post-Newtonian potential \[Chi].";
PotentialU::usage = "PotentialU[] is the post-Newtonian potential U.";
PotentialUU::usage = "PotentialUU[-i, -j] is the post-Newtonian potential \!\(U\_\(ij\)\).";
PotentialV::usage = "PotentialV[-i] is the post-Newtonian potential \!\(V\_i\).";
PotentialW::usage = "PotentialW[-i] is the post-Newtonian potential \!\(W\_i\).";
PotentialPhiW::usage = "PotentialPhiW[] is the post-Newtonian potential \!\(\[CapitalPhi]\_W\).";
PotentialPhi1::usage = "PotentialPhi1[] is the post-Newtonian potential \!\(\[CapitalPhi]\_1\).";
PotentialPhi2::usage = "PotentialPhi2[] is the post-Newtonian potential \!\(\[CapitalPhi]\_2\).";
PotentialPhi3::usage = "PotentialPhi3[] is the post-Newtonian potential \!\(\[CapitalPhi]\_3\).";
PotentialPhi4::usage = "PotentialPhi4[] is the post-Newtonian potential \!\(\[CapitalPhi]\_4\).";
PotentialA::usage = "PotentialA[] is the post-Newtonian potential \[ScriptCapitalA].";
PotentialB::usage = "PotentialB[] is the post-Newtonian potential \[ScriptCapitalB].";

ParameterBeta::usage = "ParameterBeta is the post-Newtonian parameter \[Beta].";
ParameterGamma::usage = "ParameterGamma is the post-Newtonian parameter \[Gamma].";
ParameterAlpha1::usage = "ParameterAlpha1 is the post-Newtonian parameter \!\(\[Alpha]\_1\).";
ParameterAlpha2::usage = "ParameterAlpha2 is the post-Newtonian parameter \!\(\[Alpha]\_2\).";
ParameterAlpha3::usage = "ParameterAlpha3 is the post-Newtonian parameter \!\(\[Alpha]\_3\).";
ParameterZeta1::usage = "ParameterZeta1 is the post-Newtonian parameter \!\(\[Zeta]\_1\).";
ParameterZeta2::usage = "ParameterZeta2 is the post-Newtonian parameter \!\(\[Zeta]\_2\).";
ParameterZeta3::usage = "ParameterZeta3 is the post-Newtonian parameter \!\(\[Zeta]\_3\).";
ParameterZeta4::usage = "ParameterZeta4 is the post-Newtonian parameter \!\(\[Zeta]\_4\).";
ParameterXi::usage = "ParameterXi is the post-Newtonian parameter \[Xi].";

SpaceTimeSplit::usage = "SpaceTimeSplit[expr, rules] splits the tensorial expression expr into time and space components, by expanding all dummy spacetime indices into sums of time and space indices and replacing all free spacetime indices with either time or space indices following the given rules.";
SpaceTimeSplits::usage = "SpaceTimeSplits[expr, rules] splits the tensorial expression expr into an array of time and space components using SpaceTimeSplit, where the free indices are converted to time indices or space indices following the given rules.";

VelocityOrder::usage = "VelocityOrder[expr, order] extracts the term of the given velocity order from expr.";
ApplyPPNRules::usage = "ApplyPPNRules[expr, head] applies any rules defined for the post-Newtonian expansion of the tensor with the given head to expr.\nApplyPPNRules[expr] applies the PPN rules for any tensors which appear in expr.";
PPN::usage = "PPN[head][indices] yields the 3+1 split of a tensor with given head and indices.\nPPN[head, order][indices] yields a single term in the perturbative expansion of the aforementioned tensor.";

UsePPNRules::usage = "UsePPNRules is an option to VelocityOrder which specifies whether PPN rules for tensors at particular velocity orders should be applied or not. Possible values are True and False, with True being the default.";

OrderSet::usage = "OrderSet[expr, value] defines a subsitution rule for an expression expr of the form PPN[head, order][indices], where value must be have the same index structure, to be used by ApplyPPNRules.";
OrderUnset::usage = "OrderSet[expr] removes a post-Newtonian substitution rule for expr previously defined by OrderSet.";
OrderClear::usage = "OrderClear[head] removes all post-Newtonian substitution rules defined for the tensor with the given head.";

SortPDs::usage = "SortPDs[expr] sorts derivatives appearing in expr such that they appear in canonical order: spatial derivatives are applied before time derivatives and are sorted lexicographically.";
SortPDsToBox::usage = "SortPDsToBox[expr] sorts derivatives appearing in expr such that pairs of spatial derivatives which combine to d'Alembert or Laplace operators are grouped and applied first.";
SortPDsToDiv::usage = "SortPDsToDiv[expr] sorts derivatives appearing in expr such that spatial derivatives which are contracted with indices of the tensor on which they act are applied first.";
SortPDsToTime::usage = "SortPDsToTime[expr] sorts derivatives appearing in expr such that time derivatives are applied before spatial derivatives.";

PotentialToSource::usage = "PotentialToSource[expr] replaces derivatives of PPN potentials by source terms.";

PotentialChiToU::usage = "PotentialChiToU[expr] applies the substitution \!\(\[Chi]\_\(,ii\) \[RightArrow] -2U\).";
PotentialUToChi::usage = "PotentialUToChi[expr] applies the substitution \!\(U \[RightArrow] - \(1\/2\) \[Chi]\_\(,ii\)\).";
PotentialUToUU::usage = "PotentialUToChi[expr] applies the substitution \!\(U \[RightArrow] U\_\(ii\)\).";
PotentialUUToU::usage = "PotentialUUToU[expr] applies the substitution \!\(U\_\(ii\) \[RightArrow] U\).";
PotentialUUToChi::usage = "PotentialUUToChi[expr] applies the substitution \!\(U\_\(ij\) \[RightArrow] \[Chi]\_\(,ij\) - \(1\/2\) \[Chi]\_\(,kk\) \[Delta]\_\(ij\)\).";

PotentialUToV::usage = "PotentialUToV[expr] applies the substitution \!\(U\_\(,0\) \[RightArrow] -V\_\(i,i\)\).";
PotentialUToW::usage = "PotentialUToW[expr] applies the substitution \!\(U\_\(,0\) \[RightArrow] W\_\(i,i\)\).";
PotentialVToU::usage = "PotentialVToU[expr] applies the substitution \!\(V\_\(i,i\) \[RightArrow] -U\_\(,0\)\).";
PotentialVToW::usage = "PotentialVToW[expr] applies the substitution \!\(V\_\(i,i\) \[RightArrow] -W\_\(i,i\)\).";
PotentialWToU::usage = "PotentialWToU[expr] applies the substitution \!\(W\_\(i,i\) \[RightArrow] U\_\(,0\)\).";
PotentialWToV::usage = "PotentialWToV[expr] applies the substitution \!\(W\_\(i,i\) \[RightArrow] -V\_\(i,i\)\).";

PotentialVToChiW::usage = "PotentialVToChiW[expr] applies the substitution \!\(V\_i \[RightArrow] W\_i + \[Chi]\_\(,0i\)\).";
PotentialWToChiV::usage = "PotentialWToChiV[expr] applies the substitution \!\(W\_i \[RightArrow] V\_i - \[Chi]\_\(,0i\)\).";
PotentialChiToPhiAB::usage = "PotentialChiToPhiAB[expr] applies the substitution \!\(\[Chi]\_\(,00\) \[RightArrow] \[ScriptCapitalA] + \[ScriptCapitalB] - \[CapitalPhi]\_1\).";
PotentialUToPhiAB::usage = "PotentialUToPhiAB[expr] applies the substitution \!\(U\_\(,00\) \[RightArrow] -\(1\/2\)(\[ScriptCapitalA] + \[ScriptCapitalB] - \[CapitalPhi]\_1)\_\(,ii\)\).";

TimeRhoToEuler::usage = "TimeRhoToEuler[expr] applies the substitution \!\(\[Rho]\_\(,0\) \[RightArrow] -(\[Rho] v\_i)\_i\).";
TimeVelToEuler::usage = "TimeVelToEuler[expr] applies the substitution \!\(v\_\(i,0\) \[RightArrow] \(1\/2\) g\&2\_\(00,i\) - v\_\(i,j\)v\_j - \(p\_\(,i\)\/\[Rho]\)\).";
TimePiToEuler::usage = "TimePiToEuler[expr] applies the substitution \!\(\[CapitalPi]\_\(,0\) \[RightArrow] v\_i (\(p\_\(,i\)\/\[Rho]\) - \[CapitalPi]\_\(,i\) - \(1\/2\) g\&2\_\(00,i\) - \(1\/2\) g\&2\_\(jj,i\)) - \(p v\_\(i,i\)\/\[Rho]\) - \(1\/2\) g\&2\_\(ii,0\)\).";

MetricToStandard::usage = "MetricToStandard[expr] converts all metric perturbations into PPN potentials and PPN parameters in the standard gauge.";

Met::usage = "Met[-\[Mu], -\[Nu]] is the physical metric \!\(g\_\(\[Mu]\[Nu]\)\)";
Tet::usage = "Tet[\[CapitalGamma], -\[Mu]] is the physical tetrad \!\(\(\[Theta]\^\[CapitalGamma]\)\_\[Mu]\)";
InvTet::usage = "InvTet[-\[CapitalGamma], \[Mu]] is the inverse physical tetrad \!\(\(e\_\[CapitalGamma]\)\^\[Mu]\)";
Asym::usage = "Asym[-\[Mu], -\[Nu]] is the antisymmetric part \!\(a\_\(\[Mu]\[Nu]\)\) of the tetrad perturbation.";
Tgen::usage = "Tgen[\[Mu], -\[Nu]] is the perturbation \!\(\[Lambda]\^\[Mu]\_\[Nu]\) of the generator of the general teleparallel connection.";
Xi::usage = "Xi[i] is the generating vector field \!\(\[Xi]\^i\) of the symmetric teleparallel connection.";
NonMetND::usage = "NonMetND[-\[Rho], -\[Mu], -\[Nu]] is the nonmetricity tensor \!\(Q\&\[Times]\_\(\[Rho]\[Mu]\[Nu]\)\) of the symmetric teleparallel connection."
NonMetTD::usage = "NonMetTD[-\[Rho], -\[Mu], -\[Nu]] is the nonmetricity tensor \!\(Q\&\[DoubleVerticalBar]\_\(\[Rho]\[Mu]\[Nu]\)\) of the general teleparallel connection."
CD::usage = "CD[-\[Mu]] is the Levi-Civita covariant derivative \!\(\[EmptyDownTriangle]\&\[EmptyCircle]\).";
ND::usage = "ND[-\[Mu]] is the symmetric teleparallel covariant derivative \!\(\[EmptyDownTriangle]\&\[Times]\).";
FD::usage = "FD[-\[Mu]] is the metric teleparallel covariant derivative \!\(\[EmptyDownTriangle]\&\[FilledCircle]\).";
TD::usage = "TD[-\[Mu]] is the general teleparallel covariant derivative \!\(\[EmptyDownTriangle]\&\[DoubleVerticalBar]\).";

Begin["xAct`xPPN`Private`"]

$MaxPPNOrder = 4;

SpacetimeVBundleQ[vb_] := And[VBundleQ[vb], BaseOfVBundle[vb] === MfSpacetime, BaseOfVBundle /@ (List @@ First[SplittingsOfVBundle[vb]]) === {MfTime, MfSpace}];

(* A PPNTensor object represents a tensor, if the following is satisfied: *)
PPNTensor /: xTensorQ[PPNTensor[head_, slots_List, o_ : 0]] := And[
	xTensorQ[head],
	MemberQ[DependenciesOfTensor[head], MfSpacetime],
	Length[slots] == Length[SlotsOfTensor[head]],
	And @@ MapThread[MemberQ[Switch[#1,
		Tangent[MfSpacetime], {Labels, Tangent[MfSpace]},
		-Tangent[MfSpacetime], {-Labels, -Tangent[MfSpace]},
		LorentzMfSpacetime, {Labels, LorentzMfSpace},
		-LorentzMfSpacetime, {-Labels, -LorentzMfSpace},
		_, {#1}
	], #2]&, {SlotsOfTensor[head], slots}, 1],
	IntegerQ[o],
	NonNegative[o]
];

PPNTensor /: SlotsOfTensor[PPNTensor[head_, slots_List, o_ : 0] ? xTensorQ] := slots;
PPNTensor /: DependenciesOfTensor[PPNTensor[head_, slots_List, o_ : 0] ? xTensorQ] := DependenciesOfTensor[head] /. {MfSpacetime -> Sequence[MfSpace, TimePar]};
PPNTensor /: SymmetryGroupOfTensor[PPNTensor[head_, slots_List, o_ : 0] ? xTensorQ] := PPNTensorSymmetry[head, slots];
PPNTensor /: PrintAs[PPNTensor[head_, slots_List] ? xTensorQ] := PrintAs[head];
PPNTensor /: PrintAs[PPNTensor[head_, slots_List, o_] ? xTensorQ] := OverscriptBox[PrintAs[head], o];

PPNTensorSymmetry[_, _] := StrongGenSet[{}, GenSet[]];

ToPPNTensor[head_[inds___]] := PPNTensor[head, MapThread[ReplaceAll[#1, _ ? SpacetimeVBundleQ -> #2]&, {SlotsOfTensor[head], VBundleOfIndex /@ {inds}}, 1]][inds];
ToPPNTensor[head_[inds___], n_] := PPNTensor[head, MapThread[ReplaceAll[#1, _ ? SpacetimeVBundleQ -> #2]&, {SlotsOfTensor[head], VBundleOfIndex /@ {inds}}, 1], n][inds];

PPN[head_ ? xTensorQ][inds___] := ToPPNTensor[head[inds]];
PPN[head_ ? xTensorQ, n_Integer ? NonNegative][inds___] := ToPPNTensor[head[inds], n];

DefTensorBeginning[head_[inds___], deps_, sym_, opts___] := Null;

DefTensorEnd[head_[inds___], deps_, sym_, opts___] := Module[{ntot, nst, pst, pin, sgs, gst, gin, grp, tup, dec, decs, rep, ss, rem, ind, i0, i1, slots, sltf, p, prm, n, vb},
	(* Do we have any indices to split here? *)
	If[Not[Or[deps === MfSpacetime, And[Head[deps] === List, MemberQ[deps, MfSpacetime]]]], Return[]];

	(* How many indices should we split and which? *)
	ntot = Length[SlotsOfTensor[head]];
	pst = Flatten[Position[SlotsOfTensor[head], _ ? SpacetimeVBundleQ | -_ ? SpacetimeVBundleQ, 1]];
	nst = Length[pst];
	pin = Complement[Range[ntot], pst];

	slots = {SlotsOfTensor[head], SlotsOfTensor[head] /. {_ ? SpacetimeVBundleQ -> Labels}, SlotsOfTensor[head] /. {vb_ ? SpacetimeVBundleQ :> SplittingsOfVBundle[vb][[1, 2]]}};
	sltf = (MapThread[Part, {Transpose[slots], # + 2}, 1]&);

	(* Print[{ntot, nst, pst, pin}]; *)

	(* Decompose the symmetry group into that of spacetime indices and some rest. *)
	sgs = SymmetryGroupOfTensor[head];
	gst = Stabilizer[pin, sgs[[2]]];
	gin = Stabilizer[pst, sgs[[2]]];
	grp = Dimino[gst];

	(* Print[{sgs, gin, gst, grp}]; *)

	(* Look at all possible splittings of indices. *)
	tup = Select[Tuples[{-1, 0, 1}, ntot], Flatten[Position[#, -1]] == pin&];
	decs = Sort[Union[Union /@ Outer[PermuteList, tup, List @@ grp, 1]], Order @@ (First /@ {##}) &];
	Do[
		rep = First[dec];
		If[MemberQ[Dimino[Stabilizer[Flatten[Position[rep, 1]], gst]], Times[-1, _?PermQ]],
			(* This index combination vanishes *)
			Do[
				head /: PPNTensor[head, sltf[ind]] := Zero;
				head /: PPNTensor[head, sltf[ind], _] := Zero,
			{ind, dec}],
			(* This index combination does not vanish *)
			i1 = Range[ntot] * rep * (rep + 1) / 2;
			i0 = Range[ntot] - i1;
			ss = SetStabilizer[Flatten[Position[rep, 1]], gst];
			rem = DeleteCases[ss /. {p_Cycles :> PermProduct[TranslatePerm[PermutationFromTo[PermuteList[i0, p] + i1, Range[ntot]], Cycles], p]}, Cycles[]];
			head /: PPNTensorSymmetry[head, sltf[rep]] := Evaluate[xAct`xTensor`Private`SGSofsym[rem]];
			Do[
				prm = First[Sort[Select[List @@ grp, PermuteList[ind, #] == rep &], Order[DeleteCases[PermuteList[Range[ntot] ind, #1], 0], DeleteCases[PermuteList[Range[ntot] ind, #2], 0]] &]];
				head /: PPNTensor[head, sltf[ind]] := Evaluate[With[{body = prm /. {p_Cycles :> PPNTensor[head, sltf[rep]] @@ (Slot /@ PermuteList[Range[ntot], p])}}, body&]];
				head /: PPNTensor[head, sltf[ind], n_] := Evaluate[With[{body = prm /. {p_Cycles :> PPNTensor[head, sltf[rep], n] @@ (Slot /@ PermuteList[Range[ntot], p])}}, body&]],
			{ind, Rest[dec]}]
		],
	{dec, decs}]
];

xTension["xPPN`", DefTensor, "Beginning"] = DefTensorBeginning;
xTension["xPPN`", DefTensor, "End"] = DefTensorEnd;

PPNRules[_, _List, _Integer] := {};

PPNRules[head_, o_Integer] := Flatten[If[SlotsOfTensor[head] == {}, PPNRules[head, {}, o], Outer[PPNRules[head, {##}, o]&, Sequence @@ (List /@ SlotsOfTensor[head] /. {{Tangent[MfSpacetime]} -> {Labels, Tangent[MfSpace]}, {-Tangent[MfSpacetime]} -> {-Labels, -Tangent[MfSpace]}, {LorentzMfSpacetime} -> {Labels, LorentzMfSpace}, {-LorentzMfSpacetime} -> {-Labels, -LorentzMfSpace}})]]];

PPNRules[head_, slots_List] := Flatten[PPNRules[head, slots, #]& /@ Range[0, $MaxPPNOrder]];

PPNRules[head_] := Flatten[PPNRules[head, #]& /@ Range[0, $MaxPPNOrder]];

OrderSet[lhs : PPNTensor[head_, slots_List, o_] ? xTensorQ[inds___], rhs_] := TagSet[head, PPNRules[head, slots, o], MakeRule[{lhs, rhs}, MetricOn -> All, ContractMetrics -> True]];

OrderUnset[PPNTensor[head_, slots_List, o_] ? xTensorQ[inds___]] := TagUnset[head, PPNRules[head, slots, o]];

OrderClear[head_] := Set[UpValues[head], Select[UpValues[head], FreeQ[#, PPNRules]&]];

CreateXiRules[xi_] := (
	OrderSet[PPN[xi, 0][LI[0]], 0];
	OrderSet[PPN[xi, 0][T3a]  , 0];
	OrderSet[PPN[xi, 1][LI[0]], 0];
	OrderSet[PPN[xi, 1][T3a]  , 0];
	OrderSet[PPN[xi, 2][LI[0]], 0];
	OrderSet[PPN[xi, 3][T3a]  , 0];
	OrderSet[PPN[xi, 4][LI[0]], 0];
	Return[xi];
);
(*
CreateCoincRules[nd_, xi_] := Module[{expr, n},
	expr = Christoffel[nd][T4\[Rho], -T4\[Mu], -T4\[Nu]];
	expr = SpaceTimeSplits[expr, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}];
	expr = Union[Flatten[expr], SameTest -> (SameQ @@ (Head /@ {##})&)];
	expr = VelocityOrder[#, 0]& /@ expr;
	OrderSet[#, 0]& /@ expr;

	Do[
		expr = {Christoffel[nd][T4\[Rho], -T4\[Mu], -T4\[Nu]], PD[-T4\[Mu]][PD[-T4\[Nu]][xi[T4\[Rho]]]] - PD[-T4\[Sigma]][xi[T4\[Rho]]] * Christoffel[nd][T4\[Sigma], -T4\[Mu], -T4\[Nu]]};
		expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
		expr = Union[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
		expr = Map[VelocityOrder[#, n]&, expr, {2}];
		expr = Simplify[ToCanonical[expr]];
		MapThread[OrderSet, Transpose[expr], 1],
	{n, $MaxPPNOrder}];

	Return[GiveSymbol[Christoffel, nd]];
];
*)
CreateCoincRules[nd_, xi_] := Module[{expr, n},
	expr = PD[-T4\[Mu]][PD[-T4\[Nu]][xi[T4\[Rho]]]];
	expr = expr + LieD[xi[T4\[Sigma]]][expr] / 2;
	expr = LieDToCovD[expr, PD];
	expr = {Christoffel[nd][T4\[Rho], -T4\[Mu], -T4\[Nu]], expr};
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {4}];
	expr = Union[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Simplify[ToCanonical[expr]];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[GiveSymbol[Christoffel, nd]];
];

CreateTeleRules[td_, tg_] := Module[{expr, n},
	expr = {Christoffel[td][T4\[Rho], -T4\[Nu], -T4\[Mu]], (2 * delta[-T4\[Alpha], T4\[Rho]] - tg[T4\[Rho], -T4\[Alpha]]) * PD[-T4\[Nu]][tg[T4\[Alpha], -T4\[Mu]]]};
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Union[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = expr /. {delta[-LI[0], LI[0]] -> 1} /. {delta[-_, LI[0]] :> 0, delta[-LI[0], _] :> 0};
	MapThread[OrderSet, Transpose[expr], 1];
	Return[GiveSymbol[Christoffel, td]];
];

CreateNonMetRules[nm_, nd_, met_] := Module[{expr},
	expr = {nm[-T4\[Rho], -T4\[Mu], -T4\[Nu]], ChangeCovD[nd[-T4\[Rho]][met[-T4\[Mu], -T4\[Nu]]], nd, PD]};
	expr = SpaceTimeSplits[#, {-T4\[Rho] -> -T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {4}];
	expr = Union[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Simplify[ToCanonical[expr]];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[nm];
];

CreateTetradRules[tet_, met_, asy_, bkt_, bki_, bkm_] := Module[{n, k},
	(* Zeroth order is background tetrad. *)
	OrderSet[PPN[tet, 0][LI[0], -LI[0]], 1];
	OrderSet[PPN[tet, 0][LI[0], -T3a]  , 0];
	OrderSet[PPN[tet, 0][L3A, -LI[0]]  , 0];
	OrderSet[PPN[tet, 0][L3A, -T3a]    , bkt[L3A, -T3a]];

	(* General formula for higher orders. *)
	Do[
		MapThread[OrderSet, {
			{
				PPN[tet, n][LI[0], -LI[0]],
				PPN[tet, n][LI[0], -T3b],
				PPN[tet, n][L3A, -LI[0]],
				PPN[tet, n][L3A, -T3b]
			},
			Simplify[ToCanonical[#]]& /@ {
				-((PPN[met, n][-LI[0], -LI[0]] /. PPNRules[met, {-Labels, -Labels}, n]) + Sum[(PPN[tet, k][LI[0], -LI[0]] /. PPNRules[tet, {Labels, -Labels}, k]) * (PPN[tet, n - k][LI[0], -LI[0]] /. PPNRules[tet, {Labels, -Labels}, n - k]) - bkm[-T3c, -T3d] * bki[-L3C, T3c] * bki[-L3D, T3d] * (PPN[tet, k][L3C, -LI[0]] /. PPNRules[tet, {LorentzMfSpace, -Labels}, k]) * (PPN[tet, n - k][L3D, -LI[0]] /. PPNRules[tet, {LorentzMfSpace, -Labels}, n - k]), {k, n - 1}]) / 2,
				-((PPN[asy, n][-LI[0], -T3b] /. PPNRules[asy, {-Labels, -Tangent[MfSpace]}, n]) + (PPN[met, n][-LI[0], -T3b] /. PPNRules[met, {-Labels, -Tangent[MfSpace]}, n]) + Sum[(PPN[tet, k][LI[0], -LI[0]] /. PPNRules[tet, {Labels, -Labels}, k]) * (PPN[tet, n - k][LI[0], -T3b] /. PPNRules[tet, {Labels, -Tangent[MfSpace]}, n - k]) - bkm[-T3c, -T3d] * bki[-L3C, T3c] * bki[-L3D, T3d] * (PPN[tet, k][L3C, -LI[0]] /. PPNRules[tet, {LorentzMfSpace, -Labels}, k]) * (PPN[tet, n - k][L3D, -T3b] /. PPNRules[tet, {LorentzMfSpace, -Tangent[MfSpace]}, n - k]), {k, n - 1}]) / 2,
				bkt[L3A, -T3e] * bkm[T3e, T3a] * (-(PPN[asy, n][-LI[0], -T3a] /. PPNRules[asy, {-Labels, -Tangent[MfSpace]}, n]) + (PPN[met, n][-LI[0], -T3a] /. PPNRules[met, {-Labels, -Tangent[MfSpace]}, n]) + Sum[(PPN[tet, k][LI[0], -T3a] /. PPNRules[tet, {Labels, -Tangent[MfSpace]}, k]) * (PPN[tet, n - k][LI[0], -LI[0]] /. PPNRules[tet, {Labels, -Labels}, n - k]) - bkm[-T3c, -T3d] * bki[-L3C, T3c] * bki[-L3D, T3d] * (PPN[tet, k][L3C, -T3a] /. PPNRules[tet, {LorentzMfSpace, -Tangent[MfSpace]}, k]) * (PPN[tet, n - k][L3D, -LI[0]] /. PPNRules[tet, {LorentzMfSpace, -Labels}, n - k]), {k, n - 1}]) / 2,
				bkt[L3A, -T3e] * bkm[T3e, T3a] * ((PPN[asy, n][-T3a, -T3b] /. PPNRules[asy, {-Tangent[MfSpace], -Tangent[MfSpace]}, n]) + (PPN[met, n][-T3a, -T3b] /. PPNRules[met, {-Tangent[MfSpace], -Tangent[MfSpace]}, n]) + Sum[(PPN[tet, k][LI[0], -T3a] /. PPNRules[tet, {Labels, -Tangent[MfSpace]}, k]) * (PPN[tet, n - k][LI[0], -T3b] /. PPNRules[tet, {Labels, -Tangent[MfSpace]}, n - k]) - bkm[-T3c, -T3d] * bki[-L3C, T3c] * bki[-L3D, T3d] * (PPN[tet, k][L3C, -T3a] /. PPNRules[tet, {LorentzMfSpace, -Tangent[MfSpace]}, k]) * (PPN[tet, n - k][L3D, -T3b] /. PPNRules[tet, {LorentzMfSpace, -Tangent[MfSpace]}, n - k]), {k, n - 1}]) / 2
			}
		}, 1],
	{n, $MaxPPNOrder}];

	Return[tet];
];

CreateInvTetradRules[itet_, tet_, bkt_] := Module[{n, m},
	(* Zeroth order is background tetrad. *)
	OrderSet[PPN[itet, 0][-LI[0], LI[0]], 1];
	OrderSet[PPN[itet, 0][-LI[0], T3a]  , 0];
	OrderSet[PPN[itet, 0][-L3A, LI[0]]  , 0];
	OrderSet[PPN[itet, 0][-L3A, T3a]    , bkt[-L3A, T3a]];

	(* Recursive formula for higher orders. *)
	Do[
		MapThread[OrderSet, {
			{
				PPN[itet, n][-LI[0], LI[0]],
				PPN[itet, n][-LI[0], T3a],
				PPN[itet, n][-L3A, LI[0]],
				PPN[itet, n][-L3A, T3a]
			},
			Simplify[ToCanonical[#]]& /@ ({
				-Sum[PPN[itet, n - m][-LI[0], LI[0]] * PPN[tet, m][LI[0], -LI[0]] + PPN[itet, n - m][-L3B, LI[0]] * PPN[tet, m][L3B, -LI[0]], {m, n}],
				-Sum[PPN[itet, n - m][-LI[0], T3a] * PPN[tet, m][LI[0], -LI[0]] + PPN[itet, n - m][-L3B, T3a] * PPN[tet, m][L3B, -LI[0]], {m, n}],
				-Sum[PPN[itet, n - m][-LI[0], LI[0]] * PPN[tet, m][LI[0], -T3b] * bkt[-L3A, T3b] + PPN[itet, n - m][-L3B, LI[0]] * PPN[tet, m][L3B, -T3b] * bkt[-L3A, T3b], {m, n}],
				-Sum[PPN[itet, n - m][-LI[0], T3a] * PPN[tet, m][LI[0], -T3b] * bkt[-L3A, T3b] + PPN[itet, n - m][-L3B, T3a] * PPN[tet, m][L3B, -T3b] * bkt[-L3A, T3b], {m, n}]
			} /. PPNRules[tet] /. PPNRules[itet])
		}, 1],
	{n, $MaxPPNOrder}];

	Return[itet];
];

CreateWeitzRules[fd_, tet_, itet_] := Module[{expr},
	expr = {Christoffel[fd][T4\[Rho], -T4\[Mu], -T4\[Nu]], itet[-L4\[CapitalAlpha], T4\[Rho]] * PD[-T4\[Mu]][tet[L4\[CapitalAlpha], -T4\[Nu]]]};
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {4}];
	expr = Union[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Simplify[ToCanonical[expr]];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[GiveSymbol[Christoffel, fd]];
];

CreateTorsionRules[fd_] := Module[{expr},
	expr = {#, TorsionToChristoffel[#]}&[Torsion[fd][T4\[Rho], -T4\[Mu], -T4\[Nu]]];
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {4}];
	expr = Union[DeleteCases[DeleteCases[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], {0, _}], {-_, _}], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Simplify[ToCanonical[expr]];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[Torsion[fd]];
];

CreateConnDiffRules[cd_, xd_] := Module[{expr},
	expr = {#, BreakChristoffel[#]}&[Christoffel[cd, xd][T4\[Rho], -T4\[Mu], -T4\[Nu]]];
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {4}];
	expr = Union[DeleteCases[DeleteCases[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], {0, _}], {-_, _}], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Simplify[ToCanonical[expr]];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[GiveSymbol[Christoffel, cd, xd]];
];

CreateMetricRules[met_, bkg_] := (
	(* Zeroth order is background metric. *)
	OrderSet[PPN[met, 0][-LI[0], -LI[0]], -1];
	OrderSet[PPN[met, 0][-LI[0], -T3a]  , 0];
	OrderSet[PPN[met, 0][-T3a, -T3b]    , bkg[-T3a, -T3b]];

	(* Vanishing components *)
	OrderSet[PPN[met, 1][-LI[0], -LI[0]], 0];
	OrderSet[PPN[met, 1][-LI[0], -T3a]  , 0];
	OrderSet[PPN[met, 1][-T3a, -T3b]    , 0];
	OrderSet[PPN[met, 2][-LI[0], -T3a]  , 0];
	OrderSet[PPN[met, 3][-LI[0], -LI[0]], 0];
	OrderSet[PPN[met, 3][-T3a, -T3b]    , 0];
	OrderSet[PPN[met, 4][-LI[0], -T3a]  , 0];

	Return[met];
);

CreateInvMetricRules[met_, bkg_] := Module[{imet, n, m},
	imet = Inv[met];

	(* Zeroth order is background metric. *)
	OrderSet[PPN[imet, 0][LI[0], LI[0]], -1];
	OrderSet[PPN[imet, 0][LI[0], T3a]  , 0];
	OrderSet[PPN[imet, 0][T3a, T3b]    , bkg[T3a, T3b]];

	(* Recursive formula for higher orders. *)
	Do[
		MapThread[OrderSet, {
			{
				PPN[imet, n][LI[0], LI[0]],
				PPN[imet, n][LI[0], T3a],
				PPN[imet, n][T3a, T3b]
			},
			Simplify[ToCanonical[#]]& /@ ({
				Sum[PPN[imet, n - m][LI[0], LI[0]] * PPN[met, m][-LI[0], -LI[0]] + PPN[imet, n - m][LI[0], T3c] * PPN[met, m][-LI[0], -T3c], {m, n}],
				Sum[PPN[imet, n - m][LI[0], T3a] * PPN[met, m][-LI[0], -LI[0]] + PPN[imet, n - m][T3a, T3c] * PPN[met, m][-LI[0], -T3c], {m, n}],
				-bkg[T3b, T3d] * Sum[PPN[imet, n - m][LI[0], T3a] * PPN[met, m][-LI[0], -T3d] + PPN[imet, n - m][T3a, T3c] * PPN[met, m][-T3c, -T3d], {m, n}]
			} /. PPNRules[met] /. PPNRules[imet])
		}, 1],
	{n, $MaxPPNOrder}];

	Return[imet];
];

CreateDetMetricRules[met_, bkg_] := Module[{dmet},
	dmet = Determinant[met, AIndex];

	OrderSet[PPN[dmet, 0][], -1];
	OrderSet[PPN[dmet, 1][], 0];
	OrderSet[PPN[dmet, 2][], PPN[met, 2][-LI[0], -LI[0]] - bkg[T3a, T3b] * PPN[met, 2][-T3a, -T3b]];
	OrderSet[PPN[dmet, 3][], 0];
	OrderSet[PPN[dmet, 4][], PPN[met, 4][-LI[0], -LI[0]] - bkg[T3a, T3b] * PPN[met, 4][-T3a, -T3b] + PPN[met, 2][-LI[0], -LI[0]] * bkg[T3a, T3b] * PPN[met, 2][-T3a, -T3b] + bkg[T3a, T3b] * bkg[T3c, T3d] * (PPN[met, 2][-T3a, -T3c] * PPN[met, 2][-T3b, -T3d] - PPN[met, 2][-T3a, -T3b] * PPN[met, 2][-T3c, -T3d]) / 2];

	Return[dmet];
];

CreateEpsilonMetricRules[met_, bkg_] := Module[{dmet, emet, ebkg, n},
	dmet = Determinant[met, AIndex];
	{emet, ebkg} = epsilon /@ {met, bkg};

	Do[
		OrderSet[PPN[emet, n][-LI[0], -T3a, -T3b, -T3c], VelocityOrder[SpaceTimeSplit[Sqrt[-dmet[]], {}], n] * ebkg[-T3a, -T3b, -T3c]];
		OrderSet[PPN[emet, n][-T3a, -T3b, -T3c, -T3d], 0],
	{n, 0, $MaxPPNOrder}];

	Return[emet];
];

CreateTetraMetricRules[met_] := Module[{expr},
	expr = {
		{
			Tetra[met][-T4\[Mu], -T4\[Nu], -T4\[Rho], -T4\[Sigma]],
			(met[-T4\[Mu], -T4\[Rho]] * met[-T4\[Nu], -T4\[Sigma]] + met[-T4\[Mu], -T4\[Sigma]] * met[-T4\[Nu], -T4\[Rho]] - met[-T4\[Mu], -T4\[Nu]] * met[-T4\[Rho], -T4\[Sigma]]) / 2 + I * epsilon[Met][-T4\[Mu], -T4\[Nu], -T4\[Rho], -T4\[Sigma]]
		},
		{
			Dagger[Tetra[met]][-T4\[Mu], -T4\[Nu], -T4\[Rho], -T4\[Sigma]],
			(met[-T4\[Mu], -T4\[Rho]] * met[-T4\[Nu], -T4\[Sigma]] + met[-T4\[Mu], -T4\[Sigma]] * met[-T4\[Nu], -T4\[Rho]] - met[-T4\[Mu], -T4\[Nu]] * met[-T4\[Rho], -T4\[Sigma]]) / 2 - I * epsilon[Met][-T4\[Mu], -T4\[Nu], -T4\[Rho], -T4\[Sigma]]
		}
	};
	expr = Map[SpaceTimeSplits[#, {-T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b, -T4\[Rho] -> -T3c, -T4\[Sigma] -> -T3d}]&, expr, {2}];
	expr = Map[Simplify[ToCanonical[#]]&, expr, {6}];
	expr = Union[DeleteCases[DeleteCases[Flatten[Transpose[expr, {5, 6, 1, 2, 3, 4}], 4], {0, _}], {-_, _}], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[{Tetra[met], Dagger[Tetra[met]]}];
];

CreateAsymRules[asy_] := (
	OrderSet[PPN[asy, 0][-LI[0], -T3a], 0];
	OrderSet[PPN[asy, 0][-T3a, -T3b]  , 0];
	OrderSet[PPN[asy, 1][-LI[0], -T3a], 0];
	OrderSet[PPN[asy, 1][-T3a, -T3b]  , 0];
	OrderSet[PPN[asy, 2][-LI[0], -T3a], 0];
	OrderSet[PPN[asy, 3][-T3a, -T3b]  , 0];
	OrderSet[PPN[asy, 4][-LI[0], -T3a], 0];
	Return[asy];
);

CreateTgenRules[tg_] := (
	OrderSet[PPN[tg, 0][LI[0], -LI[0]], 1];
	OrderSet[PPN[tg, 0][T3a, -LI[0]]  , 0];
	OrderSet[PPN[tg, 0][LI[0], -T3a]  , 0];
	OrderSet[PPN[tg, 0][T3a, -T3b]    , delta[-T3b, T3a]];
	OrderSet[PPN[tg, 1][LI[0], -LI[0]], 0];
	OrderSet[PPN[tg, 1][T3a, -LI[0]]  , 0];
	OrderSet[PPN[tg, 1][LI[0], -T3a]  , 0];
	OrderSet[PPN[tg, 1][T3a, -T3b]    , 0];
	OrderSet[PPN[tg, 2][T3a, -LI[0]]  , 0];
	OrderSet[PPN[tg, 2][LI[0], -T3a]  , 0];
	OrderSet[PPN[tg, 3][LI[0], -LI[0]], 0];
	OrderSet[PPN[tg, 3][T3a, -T3b]    , 0];
	OrderSet[PPN[tg, 4][T3a, -LI[0]]  , 0];
	OrderSet[PPN[tg, 4][LI[0], -T3a]  , 0];
	Return[tg];
);

CreateEnMomRules[em_, met_, dens_, pres_, int_, vel_, bkg_] := (
	OrderSet[PPN[em, 0][-LI[0], -LI[0]], 0];
	OrderSet[PPN[em, 0][-LI[0], -T3a]  , 0];
	OrderSet[PPN[em, 0][-T3a, -T3b]    , 0];
	OrderSet[PPN[em, 1][-LI[0], -LI[0]], 0];
	OrderSet[PPN[em, 1][-LI[0], -T3a]  , 0];
	OrderSet[PPN[em, 1][-T3a, -T3b]    , 0];
	OrderSet[PPN[em, 2][-LI[0], -LI[0]], dens[]];
	OrderSet[PPN[em, 2][-LI[0], -T3a]  , 0];
	OrderSet[PPN[em, 2][-T3a, -T3b]    , 0];
	OrderSet[PPN[em, 3][-LI[0], -LI[0]], 0];
	OrderSet[PPN[em, 3][-LI[0], -T3a]  , -dens[] * vel[-T3a]];
	OrderSet[PPN[em, 3][-T3a, -T3b]    , 0];
	OrderSet[PPN[em, 4][-LI[0], -LI[0]], dens[] * (int[] + vel[T3a] * vel[-T3a] - PPN[met, 2][-LI[0], -LI[0]])];
	OrderSet[PPN[em, 4][-LI[0], -T3a]  , 0];
	OrderSet[PPN[em, 4][-T3a, -T3b]    , dens[] * vel[-T3a] * vel[-T3b] + pres[] * bkg[-T3a, -T3b]];
	OrderSet[PPN[em, 5][-LI[0], -LI[0]], 0];
	OrderSet[PPN[em, 5][-LI[0], -T3a]  , -(dens[] * int[] + dens[] * vel[T3b] * vel[-T3b] + pres[]) * vel[-T3a] - dens[] * PPN[met, 3][-LI[0], -T3a] - dens[] * PPN[met, 2][-T3a, -T3b] * vel[T3b]];
	OrderSet[PPN[em, 5][-T3a, -T3b]    , 0];
	Return[em];
);

CreateLeviCivitaRules[cd_, met_] := Module[{expr},
	expr = {#, ChristoffelToGradMetric[#]}&[Christoffel[cd][T4\[Rho], -T4\[Mu], -T4\[Nu]]];
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Union[Flatten[Transpose[expr, {4, 1, 2, 3}], 2], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[GiveSymbol[Christoffel, cd]];
];

CreateRiemannRules[cd_] := Module[{expr},
	expr = {#, ChangeCurvature[#, cd, PD]}&[Riemann[cd][-T4\[Nu], -T4\[Mu], -T4\[Sigma], T4\[Rho]]];
	expr = SpaceTimeSplits[#, {T4\[Rho] -> T3c, -T4\[Sigma] -> -T3d, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {5}];
	expr = Union[DeleteCases[DeleteCases[Flatten[Transpose[expr, {5, 1, 2, 3, 4}], 3], {0, _}], {-_, _}], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[Riemann[cd]];
];

CreateRiemannDownRules[cd_, met_] := Module[{expr},
	expr = {#, ChangeCurvature[RiemannDownToRiemann[#], cd, PD]}&[RiemannDown[cd][-T4\[Nu], -T4\[Mu], -T4\[Sigma], -T4\[Rho]]];
	expr = SpaceTimeSplits[#, {-T4\[Rho] -> -T3c, -T4\[Sigma] -> -T3d, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {5}];
	expr = Union[DeleteCases[DeleteCases[Flatten[Transpose[expr, {5, 1, 2, 3, 4}], 3], {0, _}], {-_, _}], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[RiemannDown[cd]];
];

CreateRicciRules[cd_, met_] := Module[{expr},
	expr = {Ricci[cd][-T4\[Mu], -T4\[Nu]], RiemannDown[cd][-T4\[Mu], -T4\[Rho], -T4\[Nu], -T4\[Sigma]] * Inv[met][T4\[Rho], T4\[Sigma]]};
	expr = SpaceTimeSplits[#, {-T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Union[Flatten[Transpose[expr, {3, 1, 2}], 1], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[Ricci[cd]];
];

CreateRicciScalarRules[cd_, met_] := Module[{expr},
	expr = {SpaceTimeSplits[RicciScalar[cd][], {}], Tr[SpaceTimeSplits[Ricci[cd][-T4\[Mu], -T4\[Nu]], {-T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}] . Transpose[SpaceTimeSplits[Inv[met][T4\[Mu], T4\[Nu]], {T4\[Mu] -> T3a, T4\[Nu] -> T3b}]]]};
	expr = Transpose[Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]]];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[RicciScalar[cd]];
];

CreateEinsteinRules[cd_, met_] := Module[{expr},
	expr = {#, EinsteinToRicci[#]}&[Einstein[cd][-T4\[Mu], -T4\[Nu]]];
	expr = SpaceTimeSplits[#, {-T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Union[Flatten[Transpose[expr, {3, 1, 2}], 1], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[Einstein[cd]];
];

CreateTFRicciRules[cd_, met_] := Module[{expr},
	expr = {TFRicci[cd][-T4\[Mu], -T4\[Nu]], Ricci[cd][-T4\[Mu], -T4\[Nu]] - RicciScalar[cd][] * met[-T4\[Mu], -T4\[Nu]] / 4};
	expr = SpaceTimeSplits[#, {-T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Union[Flatten[Transpose[expr, {3, 1, 2}], 1], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[TFRicci[cd]];
];

CreateWeylRules[cd_, met_] := Module[{expr},
	expr = {#, WeylToRiemann[#]}&[Weyl[cd][-T4\[Nu], -T4\[Mu], -T4\[Sigma], -T4\[Rho]]];
	expr = SpaceTimeSplits[#, {-T4\[Rho] -> -T3c, -T4\[Sigma] -> -T3d, -T4\[Mu] -> -T3a, -T4\[Nu] -> -T3b}]& /@ expr;
	expr = Map[Simplify[ToCanonical[#]]&, expr, {5}];
	expr = Union[DeleteCases[DeleteCases[Flatten[Transpose[expr, {5, 1, 2, 3, 4}], 3], {0, _}], {-_, _}], SameTest -> (SameQ @@ (Head /@ First /@ {##})&)];
	expr = Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]];
	expr = Flatten[Transpose[expr, {2, 3, 1}], 1];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[Weyl[cd]];
];

CreateKretschmannRules[cd_, met_] := Module[{expr},
	expr = {Kretschmann[cd][], RiemannDown[cd][-T4\[Alpha], -T4\[Beta], -T4\[Gamma], -T4\[Delta]] * RiemannDown[cd][-T4\[Mu], -T4\[Nu], -T4\[Rho], -T4\[Sigma]] * Inv[met][T4\[Alpha], T4\[Mu]] * Inv[met][T4\[Beta], T4\[Nu]] * Inv[met][T4\[Gamma], T4\[Rho]] * Inv[met][T4\[Delta], T4\[Sigma]]};
	expr = SpaceTimeSplits[#, {}]& /@ expr;
	expr = Transpose[Outer[VelocityOrder, expr, Range[0, $MaxPPNOrder]]];
	expr = Map[Simplify[ToCanonical[ContractMetric[#, OverDerivatives -> True, AllowUpperDerivatives -> True]]]&, expr, {2}];
	MapThread[OrderSet, Transpose[expr], 1];
	Return[Kretschmann[cd]];
];

SpaceTimeSplit[expr_, reps_] := Module[{fi, res, h, i, x},
	fi = List @@ IndicesOf[Free, Select[$SumVBundles, BaseOfVBundle[#] === MfSpacetime&]][expr];
	If[Not[Sort[fi] === Sort[First /@ reps]], Throw[Message[SpaceTimeSplit::error, "Replacement list must contain replacements for all free spacetime indices."]]];

	(* Trace dummy indices into space and time. *)
	res = ReplaceIndex[expr, reps];
	res = Fold[TraceProductDummy, res, Select[$SumVBundles, BaseOfVBundle[#] === MfSpacetime&]];
	res = res /. Scalar[x_] :> Scalar[Fold[TraceProductDummy, x, Select[$SumVBundles, BaseOfVBundle[#] === MfSpacetime&]]];
	res = res /. (Rule[#, LI[0]]& /@ (List @@ IndicesOf[Dummy, Up, Select[$VBundles, BaseOfVBundle[#] === MfTime&]][res]));
	res = res /. Scalar[x_] :> Scalar[x /. (Rule[#, LI[0]]& /@ (List @@ IndicesOf[Dummy, Up, Select[$VBundles, BaseOfVBundle[#] === MfTime&]][x]))];
	res = res /. {(h_ ? (And[xTensorQ[#], MemberQ[DependenciesOfTensor[#], MfSpacetime]]&))[i___] :> ToPPNTensor[h[i]]};
	res = res /. {PD[-LI[0]] -> ParamD[TimePar]};
	Return[res];
];

SpaceTimeSplits[expr_, reps_] := Module[{fi, res, ind, rep},
	fi = List @@ IndicesOf[Free, Select[$SumVBundles, BaseOfVBundle[#] === MfSpacetime&]][expr];
	If[Not[Sort[fi] === Sort[First /@ reps]], Throw[Message[SpaceTimeSplits::error, "Replacement list must contain replacements for all free spacetime indices."]]];

	inds = First /@ reps;
	rep = ({If[DownIndexQ[#], -1, 1] * LI[0], #}&) /@ Last /@ reps;
	res = Array[SpaceTimeSplit[expr, MapThread[Rule, {inds, MapThread[Part, {rep, {##}}, 1]}, 1]]&, Table[2, {Length[fi]}]];
	Return[res];
];

DefInertHead[VelocityOrder, PrintAs -> "\[ScriptCapitalO]", DefInfo -> False, ContractThrough -> {delta}, LinearQ -> True];

MakeBoxes[VelocityOrder ? xAct`xTensor`Private`HeldInertHeadQ[expr_, n_], StandardForm] := xAct`xTensor`Private`interpretbox[VelocityOrder[expr, n], RowBox[{SubscriptBox[PrintAs[VelocityOrder], ToString[n]], "[", MakeBoxes[expr, StandardForm], "]"}]];

Options[VelocityOrder] = {UsePPNRules -> True};

VelocityOrder[expr_, n_ ? Negative, OptionsPattern[]] := 0;

VelocityOrder[expr_Plus, n_, opt : OptionsPattern[]] := (VelocityOrder[#, n, opt]&) /@ expr;

VelocityOrder[Times[ex0_, ex1__], n_, opt : OptionsPattern[]] := Module[{k}, Sum[VelocityOrder[ex0, k, opt] * VelocityOrder[Times[ex1], n - k, opt], {k, 0, n}]];

VelocityOrder[Power[expr_, p_Integer ? Positive], n_, opt : OptionsPattern[]] := Module[{k}, Sum[VelocityOrder[expr, k, opt] * VelocityOrder[Power[expr, p - 1], n - k, opt], {k, 0, n}]];

VelocityOrder[ParamD[t : TimePar ..][expr_], n_, opt : OptionsPattern[]] := ParamD[t][VelocityOrder[expr, n - Length[{t}], opt]];

VelocityOrder[Scalar[expr_], n_, opt : OptionsPattern[]] := Scalar[VelocityOrder[expr, n, opt]];

VelocityOrder[fkt_ ? ScalarFunctionQ[args__], n_, opt : OptionsPattern[]] := Module[{t, i}, SeriesCoefficient[fkt @@ (Sum[t^i * VelocityOrder[#, i, opt], {i, 0, n}]& /@ {args}), {t, 0, n}]];

VelocityOrder[PD[i_][expr_], n_, opt : OptionsPattern[]] := PD[i][VelocityOrder[expr, n, opt]];

VelocityOrder[PPNTensor[head_, slots_][inds___], n_, OptionsPattern[]] := If[OptionValue[UsePPNRules], ApplyPPNRulesTo, Identity][PPNTensor[head, slots, n][inds]];

VelocityOrder[expr_ ? ConstantQ, n_, OptionsPattern[]] := KroneckerDelta[n, 0] * expr;

VelocityOrder[expr : delta[_, _], n_, OptionsPattern[]] := KroneckerDelta[n, 0] * expr;

VelocityOrder[expr : BkgMetricS3[_, _], n_, OptionsPattern[]] := KroneckerDelta[n, 0] * expr;

VelocityOrder[expr : BkgTetradS3[_, _], n_, OptionsPattern[]] := KroneckerDelta[n, 0] * expr;

VelocityOrder[expr : BkgInvTetradS3[_, _], n_, OptionsPattern[]] := KroneckerDelta[n, 0] * expr;

ApplyPPNRulesTo[expr : PPNTensor[head_, slots_, n_][inds___]] := If[And[Head[PPNRules[head, slots, n]] === List, Or @@ (MatchQ[expr, #]& /@ First /@ PPNRules[head, slots, n])], expr /. PPNRules[head, slots, n] /. pt : PPNTensor[_, {___}, _][___] :> ApplyPPNRulesTo[pt], expr];

ApplyPPNRulesTo[expr : PPNTensor[head_, slots_, n_][inds___], h_] := If[And[head === h, Head[PPNRules[head, slots, n]] === List, Or @@ (MatchQ[expr, #]& /@ First /@ PPNRules[head, slots, n])], expr /. PPNRules[head, slots, n] /. pt : PPNTensor[h, {___}, _][___] :> ApplyPPNRulesTo[pt, h], expr];

ApplyPPNRules[expr_] := expr /. pt : PPNTensor[_, {___}, _][___] :> ApplyPPNRulesTo[pt];

ApplyPPNRules[expr_, h_] := expr /. pt : PPNTensor[h, {___}, _][___] :> ApplyPPNRulesTo[pt, h];

SortPDsToTime[expr_, head_] := Module[{a, x, t},
	Return[expr //. ParamD[t : TimePar ..][PD[a_][x_]] /; Not[FreeQ[x, head]] :> PD[a][ParamD[t][x]]];
];

SortPDsToDiv[expr_, head_] := Module[{ru, a, b, x, t},
	ru = {
		PD[a_][ParamD[t : TimePar ..][x_]] /; Not[FreeQ[x, head[___, ChangeIndex[a], ___]]] :> ParamD[t][PD[a][x]],
		PD[a_][PD[b_][x_]] /; And[FreeQ[x, head[___, ChangeIndex[b], ___]], Not[FreeQ[x, head[___, ChangeIndex[a], ___]]]] :> PD[b][PD[a][x]],
		PD[a_][ParamD[t : TimePar ..][x_]] /; Not[FreeQ[x, PPNTensor[head, __][___, ChangeIndex[a], ___]]] :> ParamD[t][PD[a][x]],
		PD[a_][PD[b_][x_]] /; And[FreeQ[x, PPNTensor[head, __][___, ChangeIndex[b], ___]], Not[FreeQ[x, PPNTensor[head, __][___, ChangeIndex[a], ___]]]] :> PD[b][PD[a][x]]
	};
	Return[expr //. ru];
];

SortPDsToBox[expr_, head_] := Module[{ru, a, b, x, t},
	ru = {
		PD[-a_][PD[a_][ParamD[t : TimePar ..][x_]]] /; Not[FreeQ[x, head]] :> ParamD[t][PD[-a][PD[a][x]]],
		PD[a_][PD[-a_][ParamD[t : TimePar ..][x_]]] /; Not[FreeQ[x, head]] :> ParamD[t][PD[-a][PD[a][x]]],
		PD[a_][ParamD[t : TimePar ..][x_]] /; And[Not[FreeQ[x, head]], Not[FreeQ[x, PD[ChangeIndex[a]][_]]]] :> ParamD[t][PD[a][x]],
		PD[a_][PD[-a_][x_]] /; Not[FreeQ[x, head]] :> PD[-a][PD[a][x]],
		PD[-a_][PD[a_][PD[b_][x_]]] /; And[Not[FreeQ[x, head]], FreeQ[x, PD[ChangeIndex[b]][_]]] :> PD[b][PD[-a][PD[a][x]]],
		PD[a_][PD[-a_][PD[b_][x_]]] /; And[Not[FreeQ[x, head]], FreeQ[x, PD[ChangeIndex[b]][_]]] :> PD[b][PD[-a][PD[a][x]]],
		PD[a_][PD[b_][x_]] /; And[Not[FreeQ[x, head]], FreeQ[x, PD[ChangeIndex[b]][_]], Not[FreeQ[x, PD[ChangeIndex[a]][_]]]] :> PD[b][PD[a][x]]
	};
	Return[expr //. ru];
];

SortPDs[expr_] := Module[{ru, a, b, x, t},
	ru = {
		PD[a_][ParamD[t : TimePar ..][x_]] :> ParamD[t][PD[a][x]],
		PD[b_][PD[a_][x_]] /; DisorderedPairQ[a, b] :> PD[a][PD[b][x]]
	};
	Return[expr //. ru];
];

End[]

$DefInfoQ = False;
$UndefInfoQ = False;

DefManifold[MfTime, 1, {\[ScriptT]}, PrintAs -> "\!\(T\_1\)"];
DefManifold[MfSpace, 3, IndicesT3, PrintAs -> "\!\(S\_3\)"];
DefManifold[MfSpacetime, 4, IndicesT4, PrintAs -> "\!\(M\_4\)"];
SplitManifold[MfSpacetime, {MfTime, MfSpace}];

DefVBundle[LorentzMfTime, MfTime, 1, {\[ScriptCapitalT]}, PrintAs -> StringJoin["\[DoubleStruckCapitalL]", PrintAs[MfTime]]];
DefVBundle[LorentzMfSpace, MfSpace, 3, IndicesL3, PrintAs -> StringJoin["\[DoubleStruckCapitalL]", PrintAs[MfSpace]]];
DefVBundle[LorentzMfSpacetime, MfSpacetime, 4, IndicesL4, PrintAs -> StringJoin["\[DoubleStruckCapitalL]", PrintAs[MfSpacetime]]];
SplitVBundle[LorentzMfSpacetime, {LorentzMfTime, LorentzMfSpace}];

MapThread[UpSet[PrintAs[#1], #2]&, {IndicesT3, IndexNamesT3}, 1];
MapThread[UpSet[PrintAs[#1], #2]&, {IndicesT4, IndexNamesT4}, 1];
MapThread[UpSet[PrintAs[#1], #2]&, {IndicesL3, IndexNamesL3}, 1];
MapThread[UpSet[PrintAs[#1], #2]&, {IndicesL4, IndexNamesL4}, 1];

DefParameter[TimePar, PrintAs -> "0"];

DefTensor[BkgTetradT1[NewIndexIn[LorentzMfTime], -NewIndexIn[Tangent[MfTime]]], {MfTime}, PrintAs -> "\[CapitalDelta]"];
DefTensor[BkgTetradS3[L3A, -T3a], {MfSpace}, PrintAs -> "\[CapitalDelta]"];
DefTensor[BkgTetradM4[L4\[CapitalAlpha], -T4\[Mu]], {MfSpacetime}, PrintAs -> "\[CapitalDelta]"];

DefTensor[BkgInvTetradT1[-NewIndexIn[LorentzMfTime], NewIndexIn[Tangent[MfTime]]], {MfTime}, PrintAs -> "\[CapitalDelta]"];
DefTensor[BkgInvTetradS3[-L3A, T3a], {MfSpace}, PrintAs -> "\[CapitalDelta]"];
DefTensor[BkgInvTetradM4[-L4\[CapitalAlpha], T4\[Mu]], {MfSpacetime}, PrintAs -> "\[CapitalDelta]"];

DefMetric[-1, BkgMetricT1[-NewIndexIn[Tangent[MfTime]], -NewIndexIn[Tangent[MfTime]]], PD, SymbolOfCovD -> {",", "\[PartialD]"},  PrintAs -> "\[Eth]", FlatMetric -> True];
DefMetric[1, BkgMetricS3[-T3a, -T3b], PD, SymbolOfCovD -> {",", "\[PartialD]"}, PrintAs -> "\[Delta]", FlatMetric -> True];
DefProductMetric[BkgMetricM4[-T4\[Mu], -T4\[Nu]], {{Tangent[MfTime], 1}, {Tangent[MfSpace], 1}}, PD, SymbolOfCovD -> {",", "\[PartialD]"}, PrintAs -> "\[Eta]", FlatMetric -> True];

AutomaticRules[BkgInvTetradM4, MakeRule[{BkgInvTetradM4[-L4\[CapitalAlpha], T4\[Mu]] * BkgTetradM4[L4\[CapitalAlpha], -T4\[Nu]], delta[-T4\[Nu], T4\[Mu]]}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[BkgInvTetradM4, MakeRule[{BkgInvTetradM4[-L4\[CapitalAlpha], T4\[Mu]] * BkgTetradM4[L4\[CapitalBeta], -T4\[Mu]], delta[-L4\[CapitalAlpha], L4\[CapitalBeta]]}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[BkgTetradM4, MakeRule[{PD[-T4\[Nu]][BkgTetradM4[L4\[CapitalAlpha], -T4\[Mu]]], 0}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[BkgInvTetradM4, MakeRule[{PD[-T4\[Nu]][BkgInvTetradM4[-L4\[CapitalAlpha], T4\[Mu]]], 0}, MetricOn -> All, ContractMetrics -> True]];

AutomaticRules[BkgInvTetradS3, MakeRule[{BkgInvTetradS3[-L3A, T3a] * BkgTetradS3[L3A, -T3b], delta[-T3b, T3a]}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[BkgInvTetradS3, MakeRule[{BkgInvTetradS3[-L3A, T3a] * BkgTetradS3[L3B, -T3a], delta[-L3A, L3B]}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[BkgTetradS3, MakeRule[{PD[-T3b][BkgTetradS3[L3A, -T3a]], 0}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[BkgInvTetradS3, MakeRule[{PD[-T3b][BkgInvTetradS3[-L3A, T3a]], 0}, MetricOn -> All, ContractMetrics -> True]];

DefMetric[-1, Met[-T4\[Mu], -T4\[Nu]], CD, SymbolOfCovD -> {";", "\!\(\[EmptyDownTriangle]\&\[EmptyCircle]\)"},  PrintAs -> "g"];
DefTensor[Tet[L4\[CapitalAlpha], -T4\[Mu]], {MfSpacetime}, PrintAs -> "\[Theta]"];
DefTensor[InvTet[-L4\[CapitalAlpha], T4\[Mu]], {MfSpacetime}, PrintAs -> "e"];

AutomaticRules[InvTet, MakeRule[{InvTet[-L4\[CapitalAlpha], T4\[Mu]] * Tet[L4\[CapitalAlpha], -T4\[Nu]], delta[-T4\[Nu], T4\[Mu]]}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[InvTet, MakeRule[{InvTet[-L4\[CapitalAlpha], T4\[Mu]] * Tet[L4\[CapitalBeta], -T4\[Mu]], delta[-L4\[CapitalAlpha], L4\[CapitalBeta]]}, MetricOn -> All, ContractMetrics -> True]];
AutomaticRules[InvTet, MakeRule[{PD[-T4\[Nu]][InvTet[-L4\[CapitalAlpha], T4\[Mu]]], -InvTet[-L4\[CapitalAlpha], T4\[Rho]] * InvTet[-L4\[CapitalBeta], T4\[Mu]] * PD[-T4\[Nu]][Tet[L4\[CapitalBeta], -T4\[Rho]]]}, MetricOn -> All, ContractMetrics -> True]];

DefCovD[FD[-T4\[Mu]], LorentzMfSpacetime, SymbolOfCovD -> {"|", "\!\(\[EmptyDownTriangle]\&\[FilledCircle]\)"}, FromMetric -> Met, Torsion -> True, Curvature -> False];
DefCovD[ND[-T4\[Mu]], SymbolOfCovD -> {"#", "\!\(\[EmptyDownTriangle]\&\[Times]\)"}, Torsion -> False, Curvature -> False];
DefCovD[TD[-T4\[Mu]], SymbolOfCovD -> {":", "\!\(\[EmptyDownTriangle]\&\[DoubleVerticalBar]\)"}, Torsion -> True, Curvature -> False];

DefTensor[Asym[-T4\[Mu], -T4\[Nu]], {MfSpacetime}, Antisymmetric[{1, 2}], PrintAs -> "a"];
DefTensor[Tgen[T4\[Mu], -T4\[Nu]], {MfSpacetime}, PrintAs -> "\[Lambda]"];
DefTensor[Xi[T4\[Mu]], {MfSpacetime}, PrintAs -> "\[Xi]"];
DefTensor[NonMetND[-T4\[Rho], -T4\[Mu], -T4\[Nu]], {MfSpacetime}, Symmetric[{2, 3}], PrintAs -> "\!\(Q\&\[Times]\)"];
DefTensor[NonMetTD[-T4\[Rho], -T4\[Mu], -T4\[Nu]], {MfSpacetime}, Symmetric[{2, 3}], PrintAs -> "\!\(Q\&\[DoubleVerticalBar]\)"];

GiveSymbol[Christoffel, CD, ND];
GiveSymbol[Christoffel, CD, FD];
GiveSymbol[Christoffel, FD, ND];
GiveSymbol[Christoffel, CD, TD];
GiveSymbol[Christoffel, FD, TD];
GiveSymbol[Christoffel, ND, TD];

DefTensor[EnergyMomentum[-T4\[Mu], -T4\[Nu]], {MfSpacetime}, Symmetric[{1, 2}], PrintAs -> "\[CapitalTheta]"];
DefTensor[TREnergyMomentum[-T4\[Mu], -T4\[Nu]], {MfSpacetime}, Symmetric[{1, 2}], PrintAs -> "\!\(\[CapitalTheta]\&_\)"];
AutomaticRules[TREnergyMomentum, MakeRule[{TREnergyMomentum[-T4\[Mu], -T4\[Nu]], EnergyMomentum[-T4\[Mu], -T4\[Nu]] - EnergyMomentum[-T4\[Rho], -T4\[Sigma]] * GiveSymbol[Inv, Met][T4\[Rho], T4\[Sigma]] * Met[-T4\[Mu], -T4\[Nu]] / 2}, MetricOn -> All, ContractMetrics -> True]];

DefTensor[Density[], {MfSpace, TimePar}, PrintAs -> "\[Rho]"];
DefTensor[Pressure[], {MfSpace, TimePar}, PrintAs -> "p"];
DefTensor[InternalEnergy[], {MfSpace, TimePar}, PrintAs -> "\[CapitalPi]"];
DefTensor[Velocity[T3a], {MfSpace, TimePar}, PrintAs -> "v"];

DefTensor[PotentialChi[], {MfSpace, TimePar}, PrintAs -> "\[Chi]"];
DefTensor[PotentialU[], {MfSpace, TimePar}, PrintAs -> "U"];
DefTensor[PotentialUU[-T3a, -T3b], {MfSpace, TimePar}, Symmetric[{1, 2}], PrintAs -> "U"];
DefTensor[PotentialV[-T3a], {MfSpace, TimePar}, PrintAs -> "V"];
DefTensor[PotentialW[-T3a], {MfSpace, TimePar}, PrintAs -> "W"];
DefTensor[PotentialPhiW[], {MfSpace, TimePar}, PrintAs -> "\!\(\[CapitalPhi]\_W\)"];
DefTensor[PotentialPhi1[], {MfSpace, TimePar}, PrintAs -> "\!\(\[CapitalPhi]\_1\)"];
DefTensor[PotentialPhi2[], {MfSpace, TimePar}, PrintAs -> "\!\(\[CapitalPhi]\_2\)"];
DefTensor[PotentialPhi3[], {MfSpace, TimePar}, PrintAs -> "\!\(\[CapitalPhi]\_3\)"];
DefTensor[PotentialPhi4[], {MfSpace, TimePar}, PrintAs -> "\!\(\[CapitalPhi]\_4\)"];
DefTensor[PotentialA[], {MfSpace, TimePar}, PrintAs -> "\[ScriptCapitalA]"];
DefTensor[PotentialB[], {MfSpace, TimePar}, PrintAs -> "\[ScriptCapitalB]"];

DefConstantSymbol[ParameterBeta, PrintAs -> "\[Beta]"];
DefConstantSymbol[ParameterGamma, PrintAs -> "\[Gamma]"];
DefConstantSymbol[ParameterAlpha1, PrintAs -> "\!\(\[Alpha]\_1\)"];
DefConstantSymbol[ParameterAlpha2, PrintAs -> "\!\(\[Alpha]\_2\)"];
DefConstantSymbol[ParameterAlpha3, PrintAs -> "\!\(\[Alpha]\_3\)"];
DefConstantSymbol[ParameterZeta1, PrintAs -> "\!\(\[Zeta]\_1\)"];
DefConstantSymbol[ParameterZeta2, PrintAs -> "\!\(\[Zeta]\_2\)"];
DefConstantSymbol[ParameterZeta3, PrintAs -> "\!\(\[Zeta]\_3\)"];
DefConstantSymbol[ParameterZeta4, PrintAs -> "\!\(\[Zeta]\_4\)"];
DefConstantSymbol[ParameterXi, PrintAs -> "\[Xi]"];

Begin["xAct`xPPN`Private`"]

BkgMetricM4 /: PPNTensor[BkgMetricM4, {-Labels, -Labels}] := (-1&);
BkgMetricM4 /: PPNTensor[BkgMetricM4, {-Labels, -TangentMfSpace}] := Zero;
BkgMetricM4 /: PPNTensor[BkgMetricM4, {-TangentMfSpace, -TangentMfSpace}] := BkgMetricS3;

BkgTetradM4 /: PPNTensor[BkgTetradM4, {Labels, -Labels}] := (1&);
BkgTetradM4 /: PPNTensor[BkgTetradM4, {Labels, -TangentMfSpace}] := Zero;
BkgTetradM4 /: PPNTensor[BkgTetradM4, {LorentzMfSpace, -Labels}] := Zero;
BkgTetradM4 /: PPNTensor[BkgTetradM4, {LorentzMfSpace, -TangentMfSpace}] := BkgTetradS3;

BkgInvTetradM4 /: PPNTensor[BkgInvTetradM4, {-Labels, Labels}] := (1&);
BkgInvTetradM4 /: PPNTensor[BkgInvTetradM4, {-Labels, TangentMfSpace}] := Zero;
BkgInvTetradM4 /: PPNTensor[BkgInvTetradM4, {-LorentzMfSpace, Labels}] := Zero;
BkgInvTetradM4 /: PPNTensor[BkgInvTetradM4, {-LorentzMfSpace, TangentMfSpace}] := BkgInvTetradS3;

$PotentialToSourceRules = Flatten[MakeRule[#, MetricOn -> All, ContractMetrics -> True]& /@ {
	{PD[-T3b][PD[T3b][PD[-T3a][PD[T3a][PotentialA[]]]]], 8 * Pi * PD[-T3a][PD[-T3b][Density[] * Velocity[T3a] * Velocity[T3b]]] - 4 * Pi * PD[-T3a][PD[T3a][Density[] * Velocity[-T3b] * Velocity[T3b]]]},
	{PD[-T3b][PD[T3b][PD[-T3a][PD[T3a][PotentialB[]]]]], 8 * Pi * PD[-T3a][PD[T3a][Pressure[]] - Density[] * PD[T3a][PotentialU[]]]},
	{PD[-T3b][PD[T3b][PD[-T3a][PD[T3a][PotentialChi[]]]]], 8 * Pi * Density[]},
	{PD[-T3a][PD[T3a][PotentialU[]]], -4 * Pi * Density[]},
	{PD[-T3b][PD[T3b][PotentialV[-T3a]]], -4 * Pi * Density[] * Velocity[-T3a]},
	{PD[-T3a][PD[T3a][PotentialPhi1[]]], -4 * Pi * Density[] * Velocity[-T3a] * Velocity[T3a]},
	{PD[-T3a][PD[T3a][PotentialPhi2[]]], -4 * Pi * Density[] * PotentialU[]},
	{PD[-T3a][PD[T3a][PotentialPhi3[]]], -4 * Pi * Density[] * InternalEnergy[]},
	{PD[-T3a][PD[T3a][PotentialPhi4[]]], -4 * Pi * Pressure[]},
	{PD[-T3a][PD[T3a][PotentialPhiW[]]], 2 * PD[-T3a][PD[-T3b][PotentialChi[]]] * PD[T3a][PD[T3b][PotentialU[]]] + PD[-T3a][PD[T3a][3 * PotentialPhi2[] - 2 * PotentialU[]^2]]}
}];

$PotentialChiToURules = MakeRule[{PD[-T3a][PD[T3a][PotentialChi[]]], -2 * PotentialU[]}, MetricOn -> All, ContractMetrics -> True];
$PotentialUToChiRules = MakeRule[{PotentialU[], -PD[-T3a][PD[T3a][PotentialChi[]]] / 2}, MetricOn -> All, ContractMetrics -> True];
$PotentialUToUURules = MakeRule[{PotentialU[], PotentialUU[-T3a, T3a]}, MetricOn -> All, ContractMetrics -> True];
$PotentialUUToURules = MakeRule[{PotentialUU[-T3a, T3a], PotentialU[]}, MetricOn -> All, ContractMetrics -> True];
$PotentialUUToChiRules = MakeRule[{PotentialUU[-T3a, -T3b], PD[-T3b][PD[-T3a][PotentialChi[]]] - PD[-T3c][PD[T3c][PotentialChi[]]] * BkgMetricS3[-T3a, -T3b] / 2}, MetricOn -> All, ContractMetrics -> True];
$PotentialUToVRules = MakeRule[{ParamD[TimePar][PotentialU[]], -PD[-T3a][PotentialV[T3a]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialUToWRules = MakeRule[{ParamD[TimePar][PotentialU[]], PD[-T3a][PotentialW[T3a]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialVToURules = MakeRule[{PD[-T3a][PotentialV[T3a]], -ParamD[TimePar][PotentialU[]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialVToWRules = MakeRule[{PD[-T3a][PotentialV[T3a]], -PD[-T3a][PotentialW[T3a]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialWToURules = MakeRule[{PD[-T3a][PotentialW[T3a]], ParamD[TimePar][PotentialU[]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialWToVRules = MakeRule[{PD[-T3a][PotentialW[T3a]], -PD[-T3a][PotentialV[T3a]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialVToChiWRules = MakeRule[{PotentialV[-T3a], PotentialW[-T3a] + ParamD[TimePar][PD[-T3a][PotentialChi[]]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialWToChiVRules = MakeRule[{PotentialW[-T3a], PotentialV[-T3a] - ParamD[TimePar][PD[-T3a][PotentialChi[]]]}, MetricOn -> All, ContractMetrics -> True];
$PotentialChiToPhiABRules = MakeRule[{ParamD[TimePar, TimePar][PotentialChi[]], PotentialA[] + PotentialB[] - PotentialPhi1[]}, MetricOn -> All, ContractMetrics -> True];
$PotentialUToPhiABRules = MakeRule[{ParamD[TimePar, TimePar][PotentialU[]], -PD[-T3a][PD[T3a][PotentialA[] + PotentialB[] - PotentialPhi1[]]] / 2}, MetricOn -> All, ContractMetrics -> True];

With[{rule = Symbol[StringJoin["$", ToString[#], "Rules"]]}, #[expr_] := expr //. rule]& /@ {
	PotentialUToChi,
	PotentialUToUU,
	PotentialUUToU,
	PotentialUUToChi,
	PotentialVToChiW,
	PotentialWToChiV
};

PotentialToSource[expr_] := Fold[SortPDsToBox, expr, {PotentialA, PotentialB, PotentialChi, PotentialU, PotentialV, PotentialPhiW, PotentialPhi1, PotentialPhi2, PotentialPhi3, PotentialPhi4}] //. $PotentialToSourceRules;

PotentialChiToU[expr_] := SortPDsToBox[expr, PotentialChi] //. $PotentialChiToURules;

PotentialVToU[expr_] := SortPDsToDiv[expr, PotentialV] //. $PotentialVToURules;
PotentialVToW[expr_] := SortPDsToDiv[expr, PotentialV] //. $PotentialVToWRules;
PotentialWToU[expr_] := SortPDsToDiv[expr, PotentialW] //. $PotentialWToURules;
PotentialWToV[expr_] := SortPDsToDiv[expr, PotentialW] //. $PotentialWToVRules;

PotentialUToV[expr_] := SortPDsToTime[expr, PotentialU] //. $PotentialUToVRules;
PotentialUToW[expr_] := SortPDsToTime[expr, PotentialU] //. $PotentialUToWRules;
PotentialUToPhiAB[expr_] := SortPDsToTime[expr, PotentialU] //. $PotentialUToPhiABRules;
PotentialChiToPhiAB[expr_] := SortPDsToTime[expr, PotentialChi] //. $PotentialChiToPhiABRules;

StandardMetricRules[met_, bkg_] := Flatten[MakeRule[#, MetricOn -> All, ContractMetrics -> True]& /@ {
	(* Zeroth order is background metric. *)
	{PPNTensor[met, {-Labels, -Labels}, 2][-LI[0], -LI[0]], 2 * PotentialU[]},
	{PPNTensor[met, {-Tangent[MfSpace], -Tangent[MfSpace]}, 2][-T3a, -T3b], 2 * ParameterGamma * PotentialU[] * bkg[-T3a, -T3b]},
	{PPNTensor[met, {-Labels, -Tangent[MfSpace]}, 3][-LI[0], -T3a], -(4 * ParameterGamma + 3 + ParameterAlpha1 - ParameterAlpha2 + ParameterZeta1 - 2 * ParameterXi) * PotentialV[-T3a] / 2 - (1 + ParameterAlpha2 - ParameterZeta1 + 2 * ParameterXi) * PotentialW[-T3a] / 2},
	{PPNTensor[met, {-Labels, -Labels}, 4][-LI[0], -LI[0]], -2 * ParameterBeta * PotentialU[]^2 - 2 * ParameterXi * PotentialPhiW[] + (2 * ParameterGamma + 2 + ParameterAlpha3 + ParameterZeta1 - 2 * ParameterXi) * PotentialPhi1[] + 2 * (3 * ParameterGamma - 2 * ParameterBeta + 1 + ParameterZeta2 + ParameterXi) * PotentialPhi2[] + 2 * (1 + ParameterZeta3) * PotentialPhi3[] + 2 * (3 * ParameterGamma + 3 * ParameterZeta4 - 2 ParameterXi) * PotentialPhi4[] - (ParameterZeta1 - 2 * ParameterXi) * PotentialA[]}
}];

TimeRhoToEuler[expr_] := SortPDsToTime[expr, Density] //. {ParamD[t___, TimePar][Density[]] :> Module[{T3a}, -ParamD[t][PD[-T3a][Density[] * Velocity[T3a]]]]};
TimeVelToEuler[expr_] := SortPDsToTime[expr, Velocity] //. {ParamD[t___, TimePar][Velocity[a_]] :> Module[{T3b}, ParamD[t][PD[a][PPN[Met, 2][-LI[0], -LI[0]]] / 2 - PD[a][Pressure[]] / Density[] - Velocity[T3b] * PD[-T3b][Velocity[a]]]]};
TimePiToEuler[expr_] := SortPDsToTime[expr, InternalEnergy] //. {ParamD[t___, TimePar][InternalEnergy[]] :> Module[{T3a, T3b}, ParamD[t][Velocity[T3a] * (PD[-T3a][Pressure[]] / Density[] - PD[-T3a][InternalEnergy[] + PPN[Met, 2][-LI[0], -LI[0]] / 2 + PPN[Met, 2][-T3b, T3b] / 2]) - Pressure[] * PD[-T3a][Velocity[T3a]] / Density[] - ParamD[TimePar][PPN[Met, 2][-T3a, T3a] / 2]]]};

$DumpSaveFile = FileNameJoin[{$xActDirectory, "xPPN", "PPNRules.mx"}];
$DumpSaveFile2 = FileNameJoin[{$TemporaryDirectory, "PPNRules.mx"}];

If[
	FileExistsQ[$DumpSaveFile],
	Get[$DumpSaveFile],

	$DumpSaveSymbols = Flatten[{
		CreateMetricRules[Met, BkgMetricS3],
		CreateInvMetricRules[Met, BkgMetricS3],
		CreateDetMetricRules[Met, BkgMetricS3],
		CreateEpsilonMetricRules[Met, BkgMetricS3],
		CreateTetraMetricRules[Met],

		CreateLeviCivitaRules[CD, Met],
		CreateRiemannRules[CD],
		CreateRiemannDownRules[CD, Met],
		CreateRicciRules[CD, Met],
		CreateRicciScalarRules[CD, Met],
		CreateEinsteinRules[CD, Met],
		CreateTFRicciRules[CD, Met],
		CreateWeylRules[CD, Met],
		CreateKretschmannRules[CD, Met],

		CreateAsymRules[Asym],
		CreateTetradRules[Tet, Met, Asym, BkgTetradS3, BkgInvTetradS3, BkgMetricS3],
		CreateInvTetradRules[InvTet, Tet, BkgInvTetradS3],
		CreateWeitzRules[FD, Tet, InvTet],
		CreateTorsionRules[FD],

		CreateXiRules[Xi],
		CreateCoincRules[ND, Xi],
		CreateNonMetRules[NonMetND, ND, Met],

		CreateTgenRules[Tgen],
		CreateTeleRules[TD, Tgen],
		CreateTorsionRules[TD],
		CreateNonMetRules[NonMetTD, TD, Met],

		CreateConnDiffRules[CD, FD],
		CreateConnDiffRules[CD, ND],
		CreateConnDiffRules[FD, ND],
		CreateConnDiffRules[CD, TD],
		CreateConnDiffRules[FD, TD],
		CreateConnDiffRules[ND, TD],

		CreateEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3]
	}];

	Check[DumpSave[$DumpSaveFile, Evaluate[$DumpSaveSymbols]],
		Print["Warning: Could not write expression cache ", $DumpSaveFile, ". Either manually copy the file ", $DumpSaveFile2, " to ", $DumpSaveFile, " or proceed without cache and have the PPN expansion calculated whenever the package is loaded."];
		DumpSave[$DumpSaveFile2, Evaluate[$DumpSaveSymbols]],
		DumpSave::noopen
	];
];

MetricToStandard[expr_] := expr //. StandardMetricRules[Met, BkgMetricS3] //. PPNRules[Met];

End[]

EndPackage[]
