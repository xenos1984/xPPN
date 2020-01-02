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

CreateLeviCivitaRules[CD, Met];
CreateRiemannRules[CD];
CreateRiemannDownRules[CD, Met];
CreateRicciRules[CD, Met];
CreateRicciScalarRules[CD, Met];
CreateEinsteinRules[CD, Met];

CreateEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3];

MetricToStandard[expr_] := expr //. Join[StandardMetricRules[Met, BkgMetricS3], PPNRules[Met]];

End[]

EndPackage[]
