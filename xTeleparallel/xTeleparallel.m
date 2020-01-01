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

PPNRules[Tau] ^= PPNTauRules[Tau, BkgMetricS3];
PPNRules[Tet] ^= PPNTetradTauRules[Tet, Tau, BkgTetradS3, BkgMetricS3];
PPNRules[InvTet] ^= PPNInvTetradRules[InvTet, Tet, BkgInvTetradS3];
PPNRules[GiveSymbol[Christoffel, FD]] ^= PPNWeitzRules[FD, Tet, InvTet];
PPNRules[GiveSymbol[Torsion, FD]] ^= PPNTorsionRules[FD];
PPNRules[Met] ^= PPNMetricTauRules[Met, Tau, BkgMetricS3];
PPNRules[GiveSymbol[Inv, Met]] ^= PPNInvMetricRules[Met, BkgMetricS3];
PPNRules[GiveSymbol[Christoffel, CD]] ^= PPNLeviCivitaRules[CD, Met];
PPNRules[GiveSymbol[Riemann, CD]] ^= PPNRiemannRules[CD];
PPNRules[GiveSymbol[RiemannDown, CD]] ^= PPNRiemannDownRules[CD, Met];
PPNRules[GiveSymbol[Ricci, CD]] ^= PPNRicciRules[CD, Met];
PPNRules[GiveSymbol[RicciScalar, CD]] ^= PPNRicciScalarRules[CD, Met];
PPNRules[GiveSymbol[Einstein, CD]] ^= PPNEinsteinRules[CD, Met];
PPNRules[GiveSymbol[Christoffel, CD, FD]] ^= PPNContortionRules[CD, FD];

PPNRules[EnergyMomentum] ^= PPNEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3];

MetricToStandard[expr_] := expr //. Join[StandardMetricRules[Met, BkgMetricS3], PPNMetricRules[Met, BkgMetricS3]];

End[]

EndPackage[]
