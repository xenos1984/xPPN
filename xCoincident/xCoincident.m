BeginPackage["xAct`xPPN`xMetric`", {"xAct`xPPN`xSpacetime`", "xAct`xTensor`", "xAct`xPerm`", "xAct`xCore`"}]

Met::usage = "";
CD::usage = "";
ND::usage = "";

$DefInfoQ = False;
$UndefInfoQ = False;

DefMetric[-1, Met[-T4\[Mu], -T4\[Nu]], CD, SymbolOfCovD -> {";", "\!\(\[EmptyDownTriangle]\&\[EmptyCircle]\)"},  PrintAs -> "g"];
DefCovD[ND[-T4\[Mu]], SymbolOfCovD -> {"#", "\!\(\[EmptyDownTriangle]\&\[Times]\)"}, Torsion -> False, Curvature -> False];

AutomaticRules[TREnergyMomentum, MakeRule[{TREnergyMomentum[-T4\[Mu], -T4\[Nu]], EnergyMomentum[-T4\[Mu], -T4\[Nu]] - EnergyMomentum[-T4\[Rho], -T4\[Sigma]] * GiveSymbol[Inv, Met][T4\[Rho], T4\[Sigma]] * Met[-T4\[Mu], -T4\[Nu]] / 2}, MetricOn -> All, ContractMetrics -> True]];

GiveSymbol[Christoffel, CD, ND];

Begin["xAct`xPPN`Private`"]

PPNRules[Met] ^= PPNMetricRules[Met, BkgMetricS3];
PPNRules[GiveSymbol[Inv, Met]] ^= PPNInvMetricRules[Met, BkgMetricS3];
PPNRules[GiveSymbol[Christoffel, CD]] ^= PPNLeviCivitaRules[CD, Met];
PPNRules[GiveSymbol[Riemann, CD]] ^= PPNRiemannRules[CD];
PPNRules[GiveSymbol[RiemannDown, CD]] ^= PPNRiemannDownRules[CD, Met];
PPNRules[GiveSymbol[Ricci, CD]] ^= PPNRicciRules[CD, Met];
PPNRules[GiveSymbol[RicciScalar, CD]] ^= PPNRicciScalarRules[CD, Met];
PPNRules[GiveSymbol[Einstein, CD]] ^= PPNEinsteinRules[CD, Met];

PPNRules[EnergyMomentum] ^= PPNEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3];

MetricToStandard[expr_] := expr //. Join[StandardMetricRules[Met, BkgMetricS3], PPNMetricRules[Met, BkgMetricS3]];

End[]

EndPackage[]
