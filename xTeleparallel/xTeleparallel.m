BeginPackage["xAct`xPPN`xTeleparallel`", {"xAct`xPPN`xSpacetime`", "xAct`xTensor`", "xAct`xPerm`", "xAct`xCore`"}]

Met::usage = "";
Tet::usage = "";
InvTet::usage = "";
Tau::usage = "";
CD::usage = "";
FD::usage = "";

$DefInfoQ = False;
$UndefInfoQ = False;

DefMetric[-1, Met[-T4\[Mu], -T4\[Nu]], CD, SymbolOfCovD -> {";", "\!\(\[EmptyDownTriangle]\&\[EmptyCircle]\)"},  PrintAs -> "g"];
DefTensor[Tet[L4\[CapitalAlpha], -T4\[Mu]], {MfSpacetime}, PrintAs -> "\[Theta]"];
DefTensor[InvTet[-L4\[CapitalAlpha], T4\[Mu]], {MfSpacetime}, PrintAs -> "e"];
DefTensor[Tau[-T4\[Mu], -T4\[Nu]], {MfSpacetime}, PrintAs -> "\[Tau]"];
DefCovD[FD[-T4\[Mu]], LorentzMfSpacetime, SymbolOfCovD -> {"|", "\!\(\[EmptyDownTriangle]\&\[FilledCircle]\)"}, FromMetric -> Met, Torsion -> True, Curvature -> False];

GiveSymbol[Christoffel, CD, FD];

Begin["xAct`xPPN`Private`"]

AutomaticRules[InvTet, MakeRule[{InvTet[-L4\[CapitalAlpha], T4\[Mu]] * Tet[L4\[CapitalAlpha], -T4\[Nu]], delta[-T4\[Nu], T4\[Mu]]}, MetricOn -> All, ContractMetrics -> True]];

CreateTauRules[Tau, BkgMetricS3];
CreateTetradTauRules[Tet, Tau, BkgTetradS3, BkgMetricS3];
CreateInvTetradRules[InvTet, Tet, BkgInvTetradS3];

CreateMetricTauRules[Met, Tau, BkgMetricS3];
CreateInvMetricRules[Met, BkgMetricS3];

CreateLeviCivitaRules[CD, Met];
CreateRiemannRules[CD];
CreateRiemannDownRules[CD, Met];
CreateRicciRules[CD, Met];
CreateRicciScalarRules[CD, Met];
CreateEinsteinRules[CD, Met];
(*
CreateWeitzRules[FD, Tet, InvTet];
CreateTorsionRules[FD];
CreateContortionRules[CD, FD];
*)
CreateEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3];

MetricToStandard[expr_] := expr //. StandardMetricRules[Met, BkgMetricS3] //. PPNRules[Met];

End[]

EndPackage[]
