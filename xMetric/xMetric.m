BeginPackage["xAct`xPPN`xMetric`", {"xAct`xPPN`xSpacetime`", "xAct`xTensor`", "xAct`xPerm`", "xAct`xCore`"}]

Met::usage = "";
CD::usage = "";

$DefInfoQ = False;
$UndefInfoQ = False;

DefMetric[-1, Met[-T4\[Mu], -T4\[Nu]], CD, SymbolOfCovD -> {";", "\[EmptyDownTriangle]"},  PrintAs -> "g"];

AutomaticRules[TREnergyMomentum, MakeRule[{TREnergyMomentum[-T4\[Mu], -T4\[Nu]], EnergyMomentum[-T4\[Mu], -T4\[Nu]] - EnergyMomentum[-T4\[Rho], -T4\[Sigma]] * GiveSymbol[Inv, Met][T4\[Rho], T4\[Sigma]] * Met[-T4\[Mu], -T4\[Nu]] / 2}, MetricOn -> All, ContractMetrics -> True]];

Begin["xAct`xPPN`Private`"]

CreateMetricRules[Met, BkgMetricS3];
CreateInvMetricRules[Met, BkgMetricS3];
(*PPNRules[GiveSymbol[Christoffel, CD]] ^= PPNLeviCivitaRules[CD, Met];
PPNRules[GiveSymbol[Riemann, CD]] ^= PPNRiemannRules[CD];
PPNRules[GiveSymbol[RiemannDown, CD]] ^= PPNRiemannDownRules[CD, Met];
PPNRules[GiveSymbol[Ricci, CD]] ^= PPNRicciRules[CD, Met];
PPNRules[GiveSymbol[RicciScalar, CD]] ^= PPNRicciScalarRules[CD, Met];
PPNRules[GiveSymbol[Einstein, CD]] ^= PPNEinsteinRules[CD, Met];
*)
CreateEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3];

MetricToStandard[expr_] := expr //. Join[StandardMetricRules[Met, BkgMetricS3], PPNRules[Met]];

End[]

EndPackage[]
