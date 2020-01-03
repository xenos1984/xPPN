BeginPackage["xAct`xPPN`xCoincident`", {"xAct`xPPN`xSpacetime`", "xAct`xTensor`", "xAct`xPerm`", "xAct`xCore`"}]

Met::usage = "";
Xi::usage = "";
CD::usage = "";
ND::usage = "";
NonMet::usage = ""

$DefInfoQ = False;
$UndefInfoQ = False;

DefMetric[-1, Met[-T4\[Mu], -T4\[Nu]], CD, SymbolOfCovD -> {";", "\!\(\[EmptyDownTriangle]\&\[EmptyCircle]\)"},  PrintAs -> "g"];
DefTensor[Xi[T4\[Mu]], {MfSpacetime}, PrintAs -> "\[Xi]"];
DefCovD[ND[-T4\[Mu]], SymbolOfCovD -> {"#", "\!\(\[EmptyDownTriangle]\&\[Times]\)"}, Torsion -> False, Curvature -> False];
DefTensor[NonMet[-T4\[Rho], -T4\[Mu], -T4\[Nu]], {MfSpacetime}, Symmetric[{2, 3}], PrintAs -> "Q"];

AutomaticRules[TREnergyMomentum, MakeRule[{TREnergyMomentum[-T4\[Mu], -T4\[Nu]], EnergyMomentum[-T4\[Mu], -T4\[Nu]] - EnergyMomentum[-T4\[Rho], -T4\[Sigma]] * GiveSymbol[Inv, Met][T4\[Rho], T4\[Sigma]] * Met[-T4\[Mu], -T4\[Nu]] / 2}, MetricOn -> All, ContractMetrics -> True]];

GiveSymbol[Christoffel, CD, ND];

Begin["xAct`xPPN`Private`"]

CreateMetricRules[Met, BkgMetricS3];
CreateInvMetricRules[Met, BkgMetricS3];

CreateXiRules[Xi];
CreateCoincRules[ND, Xi];
CreateConnDiffRules[CD, ND];
CreateNonMetRules[NonMet, ND, Met];

CreateLeviCivitaRules[CD, Met];
CreateRiemannRules[CD];
CreateRiemannDownRules[CD, Met];
CreateRicciRules[CD, Met];
CreateRicciScalarRules[CD, Met];
CreateEinsteinRules[CD, Met];

CreateEnMomRules[EnergyMomentum, Met, Density, Pressure, InternalEnergy, Velocity, BkgMetricS3];

MetricToStandard[expr_] := expr //. Join[StandardMetricRules[Met, BkgMetricS3], PPNMetricRules[Met, BkgMetricS3]];

End[]

EndPackage[]
