#!/usr/bin/env wolframscript
(* MPP manuscript v4.21: derivation checks (generated from MPP_checks_v4.21.ipynb). Run: wolframscript -file MPP_checks_v4.21.wls *)

(* Verification helpers for a manuscript notebook (Wolfram Language).
   check[name, expr, claimed, "Assumptions" -> ..., "TestPoint" -> {x -> 37/100, ...}]
     evaluates expr (held) and compares it with claimed. Give exact test points (37/100, not 0.37):
     they are raised to 50 digits before the 30-digit comparison. Verdicts:
       exact                              Simplify or FullSimplify of the difference gives 0
       numerically equal (20 digits)      |difference| < 10^-20 at TestPoint
       not evaluated (numerically equal)  Integrate/Sum/Limit/... left unevaluated, but equal at TestPoint
       not evaluated                      left unevaluated and no numerical value obtained
       undecided (no test point)          not simplified to 0 and no TestPoint given
       undecided (numerical evaluation failed)  the difference at TestPoint is not a number
       MISMATCH                           the difference at TestPoint is not small
       timeout                            expr took longer than $checkTimeLimit seconds
     The numerical tolerance is relative: 10^-20 times max(1, |claimed value at TestPoint|).
     True/False results (Resolve, Reduce) are compared directly: exact or MISMATCH.
   numcheck[name, value, "printed digits"]
     plain decimals only ("-0.0185", "12"), no exponent. Verdicts: matches printed (the printed
     string is the correctly rounded or the truncated value), MISMATCH, not numeric, bad printed string.
   Every result is stored in the association results and exported by exportResults[...].
   Name each check by the printed equation number, e.g. "(S48) J closed form, ab>0". *)

(* ---------- helpers ---------- *)
ClearAll["Global`*"];
$checkTimeLimit = 600;
results = <||>;
zeroQ[x_] := AllTrue[Flatten[{x}], # === 0 &];
stripCond[x_] := x /. ConditionalExpression[v_, _] :> v;

$unevaluatedHeads = Integrate | Sum | Product | Limit | Residue | MellinTransform | InverseMellinTransform |
  LaplaceTransform | InverseLaplaceTransform | FourierTransform | InverseFourierTransform |
  SeriesCoefficient | DSolve | DSolveValue | Expectation | NIntegrate | NSum;

SetAttributes[check, HoldAll];
Options[check] = {"Assumptions" -> True, "TestPoint" -> {}};
check[name_String, expr_, claimed_, OptionsPattern[]] := Module[
  {timC, resC, clC, assumC, tpC, verdC, numC, refC, scaleC, unevalQ},
  assumC = OptionValue["Assumptions"]; tpC = OptionValue["TestPoint"];
  {timC, resC} = AbsoluteTiming[TimeConstrained[expr, $checkTimeLimit, $Aborted]];
  resC = stripCond[resC]; clC = claimed;
  unevalQ = ! FreeQ[resC, $unevaluatedHeads];
  verdC = Which[
    resC === $Aborted, "timeout",
    BooleanQ[resC] && BooleanQ[clC], If[resC === clC, "exact", "MISMATCH"],
    unevalQ, "not evaluated",
    TrueQ[TimeConstrained[zeroQ[Simplify[resC - clC, assumC]], 120, False]], "exact",
    TrueQ[TimeConstrained[zeroQ[FullSimplify[resC - clC, assumC]], 300, False]], "exact",
    True, "check numerically"];
  numC = "n/a";
  If[MemberQ[{"check numerically", "not evaluated"}, verdC],
    numC = Quiet @ Check[TimeConstrained[N[(resC - clC) /. SetPrecision[tpC, 50], 30], $checkTimeLimit, "n/a"], "n/a"];
    refC = Quiet @ Check[TimeConstrained[N[clC /. SetPrecision[tpC, 50], 30], $checkTimeLimit, 1], 1];
    scaleC = Max[1, Abs[Select[Flatten[{refC}], NumericQ]]];   (* relative tolerance for large values *)
    verdC = Which[
      numC === "n/a" || ! AllTrue[Flatten[{numC}], NumericQ],
        Which[unevalQ, "not evaluated", tpC === {}, "undecided (no test point)",
          True, "undecided (numerical evaluation failed)"],
      AllTrue[Flatten[{numC}], Abs[#] < 10^-20 scaleC &],
        If[unevalQ, "not evaluated (numerically equal)", "numerically equal (20 digits)"],
      True, "MISMATCH"]];
  results[name] = <|"verdict" -> verdC, "result" -> ToString[resC, InputForm],
    "claimed" -> ToString[clC, InputForm], "difference_N30" -> ToString[numC, InputForm],
    "seconds" -> Round[timC, 0.01]|>;
  Print[name, "  ->  ", verdC, "   (", Round[timC, 0.01], " s)"];
  resC];

numcheck[name_String, val_, printed_String] := Module[{dec, p, v, verd},
  v = N[val, 40];
  verd = Which[
    ! StringMatchQ[printed, RegularExpression["-?[0-9]+(\\.[0-9]+)?"]], "bad printed string",
    ! NumericQ[v], "not numeric",
    True,
      dec = If[StringContainsQ[printed, "."], StringLength[Last[StringSplit[printed, "."]]], 0];
      p = ToExpression[printed <> "`40"];
      If[TrueQ[Abs[v - p] <= 10^-dec/2 || Abs[IntegerPart[v 10^dec] - p 10^dec] < 10^-20],
        "matches printed", "MISMATCH"]];
  results[name] = <|"verdict" -> verd,
    "result" -> ToString[N[val, 20], InputForm], "claimed" -> printed,
    "difference_N30" -> If[NumericQ[v] && NumericQ[p], ToString[N[v - p, 5], InputForm], "n/a"], "seconds" -> 0|>;
  Print[name, ":  ", N[val, 15], "  printed ", printed, "  ->  ", verd];
  verd === "matches printed"];

(* ---------- generic tools ---------- *)
(* Integral over the real line of a rational function f of x, by closing the contour in the
   upper half plane (the factor 2 Pi I is included). Valid when f decays at least like 1/x^2 and
   has no poles on the real axis. Give assumptions on the parameters, e.g. a > 0 && b > 0, so that
   the half plane of every pole can be decided. Returns $Failed with a message otherwise. *)
uhpIntegrate::undecided = "Cannot decide the half plane of the poles `1`. Give assumptions.";
uhpIntegrate::real = "Poles on the real axis `1`: the integral is not defined as written.";
uhpIntegrate::decay = "The integrand is not a rational function decaying at least like 1/x^2.";
uhpIntegrate[f_, x_, assum_ : True] := Catch[Module[{g, num, den, poles, up, real, undecided},
  g = Together[f]; num = Numerator[g]; den = Denominator[g];
  If[! PolynomialQ[num, x] || ! PolynomialQ[den, x] || Exponent[den, x] - Exponent[num, x] < 2,
    Message[uhpIntegrate::decay]; Throw[$Failed]];
  poles = DeleteDuplicates[x /. Solve[den == 0, x]];
  real = Select[poles, TrueQ[Simplify[Im[#] == 0, assum]] &];
  If[real =!= {}, Message[uhpIntegrate::real, real]; Throw[$Failed]];
  up = Select[poles, TrueQ[Simplify[Im[#] > 0, assum]] &];
  undecided = Select[poles, ! TrueQ[Simplify[Im[#] > 0, assum]] && ! TrueQ[Simplify[Im[#] < 0, assum]] &];
  If[undecided =!= {}, Message[uhpIntegrate::undecided, undecided]; Throw[$Failed]];
  2 Pi I Total[Residue[g, {x, #}] & /@ up]]];

(* Summary and export: put this in the last cell. *)
exportResults[manuscript_String, file_String : "wolfram_check_results.json"] := Module[{summary, out},
  summary = KeyValueMap[{#1, #2["verdict"], #2["seconds"]} &, results];
  If[summary === {}, Print["no checks were run"],
    Print[TableForm[summary, TableHeadings -> {None, {"check", "verdict", "s"}}]];
    Print["counts: ", Counts[summary[[All, 2]]]]];
  out = FileNameJoin[{Directory[], file}];
  Export[out, <|"version" -> $Version, "manuscript" -> manuscript, "results" -> results|>, "RawJSON"];
  Print["written: ", out]; out];

h[z_] := 1 - z ArcCot[z];                          (* Eqs. (10) and (S1): three-dimensional shape h_3(z) = 1 - z arctan(1/z), z > 0 *)
hp[z_] := -ArcCot[z] + z/(1 + z^2);                (* derivative h_3'(z) *)
g3[z_] := -1/(12 (1 + z^2)^2);                     (* text Sec. 2.9 and SI Sec. S14.6 (from SI Eq. (S79)): q^2 coefficient g_3(z) of L(q,qz) in three dimensions *)
R1[z_] := 4/(1 + z^2) ((1 + 3 z^2) - (2 + 3 z^2) z ArcCot[z]);   (* Eq. (36): R^(1)(z) *)
fL[r_] := 1/2 + (1 - r^2)/(4 r) Log[(1 + r)/(1 - r)];             (* Eq. (S1): f_L(r), 0 < r < 1 *)
alphaL = (4/(9 Pi))^(1/3);                         (* Eq. (3) at d = 3: alpha_3 = alpha_L *)
A0 = 3/(32 Pi^5);                                  (* SI Eq. (S18): kite prefactor A_0 *)
calA[d_] := 2 Pi^(d/2)/Gamma[d/2];                          (* text before (6): area A_d of the unit sphere S^(d-1) *)
calC[d_] := (4 Pi)^((d - 1)/2) Gamma[(d - 1)/2];             (* text before (3): Fourier transform C_d e^2/q^(d-1) of e^2/r *)
volV[d_] := Pi^(d/2)/Gamma[d/2 + 1];                         (* text before (3): volume V_d of the unit ball *)
alphaD[d_] := (2^(d - 1) Gamma[d/2 + 1]^2)^(-1/d);           (* (3): k_F a_B = 1/(alpha_d r_s) *)
calB[d_] := 2 calA[d] calC[d]/(2 Pi)^d;                      (* (6): q_TF^(d-1) = B_d/(k_F a_B) = B_d alpha_d r_s *)
E0d[d_, rs_] := 2/(alphaD[d] rs)^2;                          (* text Sec. 2.2: E_0 = hbar^2 k_F^2/m_e = 2/(alpha_d r_s)^2 Ry, with k_F a_B = 1/(alpha_d r_s) of (3) *)
qTFd[d_, rs_] := (calB[d] alphaD[d] rs)^(1/(d - 1));         (* (6): q_TF^(d-1) = B_d alpha_d r_s *)
nR[d_] := (d + 1)/(d - 1);                                   (* (28): ring resonance index N_r *)
hD[d_, z_] := 1 - Hypergeometric2F1[1, 1/2, d/2, -1/z^2];    (* (9): long-wavelength Lindhard shape h_d(z), z > 0 *)
h2D[z_] := 1 - z/Sqrt[1 + z^2];                              (* (10): h_2(z), z >= 0 *)
h4D[z_] := (Sqrt[1 + z^2] - z)^2;                            (* (10): h_4(z), z >= 0 *)
g2D[z_] := -z/(8 (1 + z^2)^(5/2));                           (* (S4): q^2 coefficient g_2(z) of the 2D L(q, qz), z >= 0 *)
lind2D[q_, z_] := 1 - (2/q) Re[Sqrt[(q/2 + I z)^2 - 1]];     (* text before (S4): 2D L(q, qz) at imaginary frequency, q > 0, z >= 0, principal Sqrt *)
lmxF[x_] := -x^2 Hypergeometric2F1[1, 2, 3, -x]/2;           (* (27)-(32): log(1 + x) - x without cancellation at small x *)
aGB = 2 (1 - Log[2])/Pi^2;                                    (* (30): Gell-Mann--Brueckner coefficient A in Ry *)
h2Dp[z_] := -1/(1 + z^2)^(3/2);                                 (* Eq. (10): derivative h_2'(z), z >= 0 *)
rVar[z_] := Sqrt[1 + z^2] - z;                                  (* Eq. (S2) and text before Eq. (38): r = s_z - z *)
aVar[z_] := z/Sqrt[1 + z^2];                                    (* Eq. (S2) and text before Eq. (38): a = z/s_z *)
R2D[z_] := (2/Pi) (h2D[z]^2 + 2 aVar[z]^2 (rVar[z] ArcTan[rVar[z]] - rVar[z]^2)/(1 + rVar[z]^2));   (* Eq. (38) first form = Eq. (S9): R_2(z) for z >= 0 (even in z) *)
aTwoL = 7/Pi + 4 Log[2]/Pi - 160/(9 Pi^2) - 16 Catalan/Pi^2;     (* Eqs. (40) and (59): a_2l in Ry *)
oddH[m_] := Sum[1/(2 j - 1), {j, 1, m}];                        (* SI Sec. S4: O_m = Sum_(j=1)^m 1/(2j-1) *)
bCoef[m_, z_] := (-I rVar[z])^Abs[m]/(I Sqrt[1 + z^2]);         (* SI Sec. S4: Fourier coefficients of b = (i z - cos theta)^(-1) *)

dvFreq17 = check["(17) on the real axis (|delta| = d)",
  Integrate[W^2/(W^2 + nu^2) (d^2 - nu^2)/(nu^2 + d^2)^2, {nu, 0, Infinity}, Assumptions -> W > 0 && d > 0]/Pi,
  W/(2 (W + d)^2), "Assumptions" -> W > 0 && d > 0, "TestPoint" -> {W -> 13/10, d -> 7/10}];
check["(17) by residues in the upper half plane, delta = -d and delta = +d",
  Module[{f}, f[del_] := W^2/(W^2 + nu^2)/(I nu - del)^2;
    {I (Residue[f[-d], {nu, I W}] + Residue[f[-d], {nu, I d}]), I Residue[f[d], {nu, I W}]}],
  {W/(2 (W + d)^2), W/(2 (W + d)^2)}, "Assumptions" -> W > 0 && d > 0, "TestPoint" -> {W -> 13/10, d -> 7/10}];

check["(S57) Re(i nu-delta)^-2 integrates to nu/(nu^2+delta^2) (antiderivative computed), absolute integral 1/|delta|",
  {Simplify[Integrate[ComplexExpand[Re[(I nu - del)^-2]], nu]],
   Integrate[Abs[ComplexExpand[Re[(I nu - d)^-2]]], {nu, 0, Infinity}, Assumptions -> d > 0]}, {nu/(nu^2 + del^2), 1/d},
  "Assumptions" -> d > 0, "TestPoint" -> {d -> 7/10, del -> 7/10, nu -> 1/3}];

dvJPlus = check["(S56) J closed form, ab>0",
  Integrate[(al be - nu^2)/((al^2 + nu^2) (be^2 + nu^2) (nu^2 + w^2)), {nu, 0, Infinity},
    Assumptions -> al > 0 && be > 0 && w > 0]/Pi,
  1/(2 w (al + w) (be + w)), "Assumptions" -> al > 0 && be > 0 && w > 0,
  "TestPoint" -> {al -> 1/3, be -> 2, w -> 7/10}];
dvJMinus = check["(S56) J closed form, ab<0",
  Integrate[(-al be - nu^2)/((al^2 + nu^2) (be^2 + nu^2) (nu^2 + w^2)), {nu, 0, Infinity},
    Assumptions -> al > 0 && be > 0 && w > 0]/Pi,
  -(al + be + 2 w)/(2 w (al + be) (al + w) (be + w)), "Assumptions" -> al > 0 && be > 0 && w > 0,
  "TestPoint" -> {al -> 1/3, be -> 2, w -> 7/10}];
check["(S54)-(S55) moments of 1/N(nu)",
  {Integrate[1/((al^2 + nu^2) (be^2 + nu^2) (w^2 + nu^2)), {nu, 0, Infinity}, Assumptions -> al > 0 && be > 0 && w > 0]/Pi,
   Integrate[nu^2/((al^2 + nu^2) (be^2 + nu^2) (w^2 + nu^2)), {nu, 0, Infinity}, Assumptions -> al > 0 && be > 0 && w > 0]/Pi},
  {(al + be + w)/(2 al be w (al + be) (be + w) (w + al)), 1/(2 (al + be) (be + w) (w + al))},
  "Assumptions" -> al > 0 && be > 0 && w > 0, "TestPoint" -> {al -> 1/3, be -> 2, w -> 7/10}];dvJInfMinus = Last[check["(S56),(18) limit w -> infinity of w^2 J, both signs, with J as computed in (S56) above",
  {Limit[w^2 dvJPlus, w -> Infinity, Assumptions -> al > 0 && be > 0], Limit[w^2 dvJMinus, w -> Infinity, Assumptions -> al > 0 && be > 0]},
  {0, -1/(al + be)}, "Assumptions" -> al > 0 && be > 0, "TestPoint" -> {al -> 1/3, be -> 2}]];
check["(S56),(17) W^2 J(W; delta, delta) with J as computed in (S56) above", W^2 dvJPlus /. {al -> d, be -> d, w -> W}, W/(2 (W + d)^2),
  "Assumptions" -> W > 0 && d > 0, "TestPoint" -> {W -> 13/10, d -> 7/10}];

check["(18) Re[1/((i nu-a)(i nu-b))] = (ab-nu^2)/((nu^2+a^2)(nu^2+b^2))",
  Simplify[ComplexExpand[Re[1/((I nu - a) (I nu - b))]]], (a b - nu^2)/((nu^2 + a^2) (nu^2 + b^2)),
  "TestPoint" -> {a -> 1, b -> -2, nu -> 3/10}];
check["(18) same signs: 0", Integrate[(a b - nu^2)/((nu^2 + a^2) (nu^2 + b^2)), {nu, 0, Infinity},
    Assumptions -> a > 0 && b > 0]/Pi, 0, "Assumptions" -> a > 0 && b > 0, "TestPoint" -> {a -> 1/3, b -> 2}];
check["(18) opposite signs: -1/(|a|+|b|)", Integrate[(-al be - nu^2)/((nu^2 + al^2) (nu^2 + be^2)), {nu, 0, Infinity},
    Assumptions -> al > 0 && be > 0]/Pi, -1/(al + be), "Assumptions" -> al > 0 && be > 0, "TestPoint" -> {al -> 1/3, be -> 2}];

check["(S81) Int Lambda^2 (c^2 - z^2)/((Lambda^2+z^2)(c^2+z^2)^2) dz = pi Lambda/(2(Lambda+c)^2)",
  Integrate[L^2 (c^2 - z^2)/((L^2 + z^2) (c^2 + z^2)^2), {z, 0, Infinity}, Assumptions -> L > 0 && c > 0],
  Pi L/(2 (L + c)^2), "Assumptions" -> L > 0 && c > 0, "TestPoint" -> {L -> 3, c -> 2/5}];
check["(S81) second equality: Pi Lambda^2 J(Lambda;c,c), with J of (S56) computed above (ab>0, alpha = beta = c, omega = Lambda)",
  Pi L^2 (dvJPlus /. {al -> c, be -> c, w -> L}), Pi L/(2 (L + c)^2), "Assumptions" -> L > 0 && c > 0, "TestPoint" -> {L -> 3, c -> 2/5}];

check["(S38) partial fractions of the sector denominators by Apart (text after (S38))",
  Apart[Simplify[ComplexExpand[Re[1/((I nu - a) (I nu - b))]]], nu], 1/(a - b) (-a/(nu^2 + a^2) + b/(nu^2 + b^2)),
  "TestPoint" -> {a -> 1/3, b -> 2, nu -> 7/10}];
check["(S39) the substitution nu = a Tan[u] absorbs the Lorentzian into the measure for any phi: integrand times Jacobian is phi(a Tan[u]), and u runs over (0, Pi/2) (nu = 0 at u = 0, 1/nu -> 0 as u -> Pi/2)",
  {Simplify[(a phi[nu]/(nu^2 + a^2) /. nu -> a Tan[u]) D[a Tan[u], u], a > 0], Limit[1/(a Tan[u]), u -> Pi/2, Direction -> "FromBelow", Assumptions -> a > 0], a Tan[0]},
  {phi[a Tan[u]], 0, 0}, "Assumptions" -> a > 0];
check["(S39) tangent substitution: H(a) for a Lorentzian line c^2/(c^2+nu^2) is c/(2(a+c))",
  Integrate[c^2/(c^2 + a^2 Tan[uu]^2), {uu, 0, Pi/2}, Assumptions -> a > 0 && c > 0]/Pi, c/(2 (a + c)),
  "Assumptions" -> a > 0 && c > 0, "TestPoint" -> {a -> 3/10, c -> 7/10}];

check["(3) alpha_d = 1/(k_F a_B r_s) from n_e V_d (r_s a_B)^d = 1 with n_e = 2 V_d k_F^d/(2 pi)^d, solved for xx = (k_F a_B r_s)^d",
  (xx /. First[Solve[(2 volV[d] xx/(2 Pi)^d) volV[d] == 1, xx]])^(-1/d), alphaD[d], "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(3) alpha_3 = alpha_L = (4/(9 pi))^(1/3)", alphaD[3], alphaL];
check["(3) alpha_2 = 1/Sqrt[2]", alphaD[2], 1/Sqrt[2]];

Module[{aB, ryd, asm, tp},
  aB = hb^2/(mE eC^2); ryd = mE eC^4/(2 hb^2);
  asm = hb > 0 && mE > 0 && eC > 0 && kF > 0 && rs > 0 && d > 1;
  tp = {hb -> 3/7, mE -> 5/4, eC -> 2/3, kF -> 7/5, rs -> 3/10, d -> 37/10};
  dvE0phys = check["(4) E_0 = hbar^2 k_F^2/m_e = 2 (k_F a_B)^2 Ry", (hb^2 kF^2/mE)/ryd, 2 (kF aB)^2,
    "Assumptions" -> asm, "TestPoint" -> tp];
  check["(text Sec. 2.2) coupling e^2/(hbar v_F) = 1/(k_F a_B)", eC^2/(hb (hb kF/mE)), 1/(kF aB),
    "Assumptions" -> asm, "TestPoint" -> tp];
  check["(text Sec. 2.2) E_0 = 2/(alpha_d r_s)^2 Ry with k_F a_B = 1/(alpha_d r_s) of (3)",
    dvE0phys /. kF -> 1/(alphaD[d] rs aB), 2/(alphaD[d] rs)^2, "Assumptions" -> asm, "TestPoint" -> tp];
  check["(text SI Sec. S1) loop factor C = v_0 k_F^d/E_0 = C_d alpha_d r_s with v_0 of (4)",
    ((calC[d] eC^2/kF^(d - 1)) kF^d/(hb^2 kF^2/mE)) /. kF -> 1/(alphaD[d] rs aB), calC[d] alphaD[d] rs,
    "Assumptions" -> asm, "TestPoint" -> tp]];

check["(text Sec. 2.2) d = 3, q >> 1: q^2 L -> 4/[3 (1 + 4 nu^2/q^4)] at fixed nu/q^2 = kap",
  Integrate[(1 - x^2)/2 Limit[q^2 (q x + q^2/2)/((q x + q^2/2)^2 + kap^2 q^4), q -> Infinity], {x, -1, 1}],
  4/(3 (1 + 4 kap^2)), "Assumptions" -> kap > 0, "TestPoint" -> {kap -> 37/100}];

check["(text Sec. 2.2) C_3 = 4 pi and C_2 = 2 pi", {calC[3], calC[2]}, {4 Pi, 2 Pi}];
check["(6) B_3 = 4/pi (q_TF^2 = 4/(pi a_B k_F)) and B_2 = 2", {calB[3], calB[2]}, {4/Pi, 2}];
check["(6) q_TF^(d-1) = C_d e^2 N(E_F)/k_F^(d-1) = B_d/(k_F a_B) (v_0 times the density of states)",
  calC[d] eC^2 (2 calA[d] kF^(d - 1)/((2 Pi)^d hb (hb kF/mE)))/kF^(d - 1), calB[d] mE eC^2/(hb^2 kF),
  "Assumptions" -> hb > 0 && mE > 0 && eC > 0 && kF > 0 && d > 1,
  "TestPoint" -> {hb -> 3/7, mE -> 5/4, eC -> 2/3, kF -> 7/5, d -> 37/10}];
check["(6) a_TF = B_d alpha_d at d = 3: a_TF = 4 alpha_L/pi", calB[3] alphaD[3], 4 alphaL/Pi];
check["(text Sec. 2.2) d = 2: q_TF = Sqrt[2] r_s", calB[2] alphaD[2] rs, Sqrt[2] rs,
  "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
check["(Tab. S2) A_2 = 2 pi and A_3 = 4 pi", {calA[2], calA[3]}, {2 Pi, 4 Pi}];
check["(text SI Sec. S1) Pi_0(q -> 0, 0) = -2 C A_d/(2 pi)^d with C = C_d alpha_d r_s gives q_TF^(d-1) = B_d alpha_d r_s",
  2 calC[d] alphaD[d] rs calA[d]/(2 Pi)^d, calB[d] alphaD[d] rs,
  "Assumptions" -> d > 1 && rs > 0, "TestPoint" -> {d -> 37/10, rs -> 3/10}];
check["(text SI Sec. S1) d = 3: C = C_3 alpha_L r_s = pi^2 q_TF^2, with r_s from (6) solved for q_TF", calC[3] alphaD[3] (rs /. First[Solve[qTF^2 == calB[3] alphaD[3] rs, rs]]), Pi^2 qTF^2,
  "Assumptions" -> qTF > 0, "TestPoint" -> {qTF -> 2/5}];

dvW6 = check["(7) Dyson equation W = lam v + lam v Pi0 W solved for W",
  wR /. First[Solve[wR == lamR vR + lamR vR pi0R wR, wR]], lamR vR/(1 - lamR vR pi0R),
  "TestPoint" -> {lamR -> 3/4, vR -> 5/3, pi0R -> -2/7}];
check["(7) equals lam v + (lam v) chi (lam v) with chi = Pi0/(1 - lam v Pi0) of (5)",
  lamR vR + lamR vR (pi0R/(1 - lamR vR pi0R)) lamR vR, lamR vR/(1 - lamR vR pi0R),
  "TestPoint" -> {lamR -> 3/4, vR -> 5/3, pi0R -> -2/7}];
check["(7) last form with v = q^(1-d) and Pi0 = -q_TF^(d-1) L of (6)",
  dvW6 /. {vR -> q^(1 - d), pi0R -> -qTF^(d - 1) lL},
  lamR/(q^(d - 1) + lamR qTF^(d - 1) lL),
  "Assumptions" -> q > 0 && qTF > 0 && lamR > 0 && lL > 0 && d > 1,
  "TestPoint" -> {q -> 3/7, qTF -> 2/5, lamR -> 3/4, lL -> 5/9, d -> 37/10}];

check["(8) c = -v Pi_0 = q_TF^(d-1) L/q^(d-1) with v = q^(1-d)", -(q^(1 - d)) (-(qTF^(d - 1) lL)), qTF^(d - 1) lL/q^(d - 1),
  "Assumptions" -> q > 0 && qTF > 0 && lL > 0 && d > 1, "TestPoint" -> {q -> 3/7, qTF -> 2/5, lL -> 5/9, d -> 37/10}];
check["(text SI Sec. S1) the weight (1 - mu^2)^((d-3)/2) on [-1, 1] has norm B(1/2, (d-1)/2)",
  Integrate[(1 - mu^2)^((d - 3)/2), {mu, -1, 1}, Assumptions -> d > 1], Beta[1/2, (d - 1)/2],
  "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(text SI Sec. S1) mu Re b = -mu^2/(mu^2 + z^2) with b = 1/(I z - mu)",
  ComplexExpand[mu Re[1/(I z - mu)]], -mu^2/(mu^2 + z^2),
  "Assumptions" -> Element[mu, Reals] && z > 0, "TestPoint" -> {mu -> 3/10, z -> 7/10}];
check["(9) Int_-1^1 P_d(mu) mu^2/(mu^2 + z^2) dmu = 1 - 2F1(1, 1/2; d/2; -1/z^2): Euler integral with the weight (1 - mu^2)^((d-3)/2)/B(1/2, (d-1)/2)",
  Integrate[(1 - mu^2)^((d - 3)/2) mu^2/(mu^2 + z^2), {mu, -1, 1}, Assumptions -> d > 1 && z > 0]/Beta[1/2, (d - 1)/2], hD[d, z],
  "Assumptions" -> d > 1 && z > 0, "TestPoint" -> {d -> 37/10, z -> 7/10}];
check["(9) the same at d = 5, where the integral is elementary",
  Integrate[(1 - mu^2) mu^2/(mu^2 + z^2), {mu, -1, 1}, Assumptions -> z > 0]/Beta[1/2, 2], hD[5, z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 13/10}];
check["(text SI Sec. S1) Euler's integral <(mu^2 + z^2)^-1> = z^-2 2F1(1, 1/2; d/2; -1/z^2)",
  Integrate[(1 - mu^2)^((d - 3)/2)/(mu^2 + z^2), {mu, -1, 1}, Assumptions -> d > 1 && z > 0]/Beta[1/2, (d - 1)/2], Hypergeometric2F1[1, 1/2, d/2, -1/z^2]/z^2,
  "Assumptions" -> d > 1 && z > 0, "TestPoint" -> {d -> 37/10, z -> 7/10}];
check["(text SI Sec. S1) d = 2: <(mu^2 + z^2)^-1> = 1/(z Sqrt[1 + z^2]), with mu = cos(theta) on the circle",
  Integrate[1/(Cos[th]^2 + z^2), {th, 0, Pi}, Assumptions -> z > 0]/Pi, 1/(z Sqrt[1 + z^2]),
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];
check["(text Sec. 2.2) h_d is even in z", hD[d, -z] - hD[d, z], 0];
check["(text Sec. 2.2) h_d(0) = 1 at d = 5 and d = 37/10 (limit z -> 0+)",
  {Limit[hD[5, z], z -> 0, Direction -> "FromAbove"], Limit[hD[37/10, z], z -> 0, Direction -> "FromAbove"]}, {1, 1}];
check["(text Sec. 2.2) h_d(0) = 1 for general d > 1: the z -> 0 limit of mu^2/(mu^2 + z^2) is 1 for mu != 0, so the average is the norm of the weight",
  Integrate[(1 - mu^2)^((d - 3)/2) Limit[mu^2/(mu^2 + z^2), z -> 0, Direction -> "FromAbove", Assumptions -> mu != 0], {mu, -1, 1}, Assumptions -> d > 1]/Beta[1/2, (d - 1)/2], 1,
  "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(text Sec. 2.2) h_d(z) = 1/(d z^2) + O(z^-4) as z -> Infinity: coefficients of w^0 and w^2 of h_d(1/w)",
  {SeriesCoefficient[hD[d, 1/w], {w, 0, 0}], SeriesCoefficient[hD[d, 1/w], {w, 0, 2}]}, {0, 1/d},
  "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];

Module[{lser, cfs},
  lser = Series[1 - (2/q) I Sqrt[1 + z^2] Sqrt[1 - (q^2/4 + I q z)/(1 + z^2)], {q, 0, 3}];
  cfs = Table[Simplify[ComplexExpand[Re[SeriesCoefficient[lser, n]]], z > 0], {n, -1, 3}];
  check["(S4) 2D: I Sqrt[1 + z^2] Sqrt[1 - (q^2/4 + I q z)/(1 + z^2)] squares to (q/2 + I z)^2 - 1",
    Expand[(I Sqrt[1 + z^2])^2 (1 - (q^2/4 + I q z)/(1 + z^2))], Expand[(q/2 + I z)^2 - 1],
    "Assumptions" -> z > 0 && q > 0, "TestPoint" -> {z -> 7/10, q -> 1/10}];
  check["(S4) 2D: this root is the principal Sqrt[(q/2 + I z)^2 - 1] at q = 1/10, z = 7/10",
    I Sqrt[1 + 49/100] Sqrt[1 - (1/400 + 7 I/100)/(1 + 49/100)], Sqrt[(1/20 + 7 I/10)^2 - 1]];
  check["(S4) 2D: Im[(q/2 + I z)^2 - 1] = q z > 0 and Im[1 - (q^2/4 + I q z)/(1 + z^2)] = -q z/(1 + z^2) < 0, so neither radicand meets the branch cut and the two continuous roots that agree at one point agree for all q, z > 0",
    {ComplexExpand[Im[(q/2 + I z)^2 - 1]], ComplexExpand[Im[1 - (q^2/4 + I q z)/(1 + z^2)]]}, {q z, -q z/(1 + z^2)},
    "Assumptions" -> z > 0 && q > 0, "TestPoint" -> {z -> 7/10, q -> 1/10}];
  check["(S4) 2D: coefficients of q^-1, q^0, ..., q^3 of L(q, qz) are {0, h_2, 0, g_2, 0}", cfs, {0, h2D[z], 0, g2D[z], 0},
    "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}]];
numcheck["(S4) g_2(0.3)", g2D[3/10], "-0.0302318708"];
check["(S4) g_2(1) = -1/(32 Sqrt[2])", g2D[1], -1/(32 Sqrt[2])];
numcheck["(S4) g_2(3)", g2D[3], "-0.0011858541"];
check["(text SI Sec. S1) z = 0: g_2(0) = 0 and L(q, 0) = 1 at q = 3/2 < 2", {g2D[0], lind2D[3/2, 0]}, {0, 1}];

check["(Tab. S2) d = 2: E_0 = 2/(alpha_2 r_s)^2 = 4/r_s^2 Ry", 2/(alphaD[2] rs)^2, 4/rs^2,
  "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
check["(Tab. S2) d = 3: E_0 = 2/(alpha_L r_s)^2 Ry = 32/(pi^2 q_TF^4) with r_s from (6) solved for q_TF",
  E0d[3, (rs /. First[Solve[qTF^2 == calB[3] alphaD[3] rs, rs]])], 32/(Pi^2 qTF^4), "Assumptions" -> qTF > 0, "TestPoint" -> {qTF -> 2/5}];
check["(Tab. S2) d = 2: (d/(8 pi)) E_0 q_TF^(d+1) = 2 Sqrt[2] r_s/pi Ry",
  (2/(8 Pi)) (2/(alphaD[2] rs)^2) (calB[2] alphaD[2] rs)^3, 2 Sqrt[2] rs/Pi, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
check["(Tab. S2) d = 3: (d/(8 pi)) E_0 q_TF^(d+1) = 12/pi^3 Ry",
  (3/(8 Pi)) (2/(alphaD[3] rs)^2) (calB[3] alphaD[3] rs)^2, 12/Pi^3, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];

check["(S1) h(z) = 1 - z arctan(1/z) equals h of the definitions cell", 1 - z ArcTan[1/z], h[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];
check["(text Sec. 2.2) h(z) = Int_0^1 a^2 da/(a^2 + z^2)", Integrate[a^2/(a^2 + z^2), {a, 0, 1}, Assumptions -> z > 0], h[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];

check["(S79)->(9) L(q,qz): coefficients of q^0..q^3 are {h, 0, g3, 0}",
  Module[{ser}, ser = Series[k/(2 q) Log[((k + q/2)^2 + z^2)/((k - q/2)^2 + z^2)], {q, 0, 3}];
    Table[Simplify[Integrate[SeriesCoefficient[ser, n], {k, 0, 1}, Assumptions -> z > 0], z > 0], {n, 0, 3}]],
  {h[z], 0, g3[z], 0}, "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];

check["(S19) Int_B |p-k|^-2 d^3k = 4 pi f_L(p) with the angular integral done first: Int_-1^1 dmu/(p^2 + k^2 - 2 p k mu), then the k integral",
  Integrate[2 Pi k^2 Integrate[1/(p^2 + k^2 - 2 p k mu), {mu, -1, 1}, Assumptions -> 0 < p < 1 && 0 < k < 1 && k != p], {k, 0, 1}, Assumptions -> 0 < p < 1],
  4 Pi fL[p], "Assumptions" -> 0 < p < 1, "TestPoint" -> {p -> 2/5}];
check["(S19) Int_B |p-k|^-2 d^3k = 4 pi f_L(p), 0<p<1, from the radial form",
  Integrate[2 Pi k/p Log[(k + p)/Abs[k - p]], {k, 0, 1}, Assumptions -> 0 < p < 1], 4 Pi fL[p],
  "Assumptions" -> 0 < p < 1, "TestPoint" -> {p -> 2/5}];
dvS91 = check["(S82) Int_0^1 r^2 f_L(r) dr = 1/4, so the same-ball integral is 16 pi^2/4 = 4 pi^2",
  {Integrate[r^2 fL[r], {r, 0, 1}], Integrate[r (1 - r^2) Log[(1 + r)/(1 - r)], {r, 0, 1}]}, {1/4, 1/3}];

check["(10) h_3 from (9) equals h of the definitions cell", hD[3, z], h[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];
check["(10) h_2 = 1 - z/Sqrt[1 + z^2] from (9) at d = 2", hD[2, z], h2D[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];
check["(10) h_4 = (Sqrt[1 + z^2] - z)^2 from (9) at d = 4", hD[4, z], h4D[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];
check["(10) h_2(0) = h_4(0) = 1", {h2D[0], h4D[0]}, {1, 1}];

dvRho = check["(19)/(S28) coupling average rho(c) = 2 Int lambda/(1+lambda c) = 2[c - log(1+c)]/c^2",
  2 Integrate[lam/(1 + lam c), {lam, 0, 1}, Assumptions -> c > 0], 2 (c - Log[1 + c])/c^2,
  "Assumptions" -> c > 0, "TestPoint" -> {c -> 7/10}];
check["(19) rho(c) = 1 - 2c/3 + c^2/2 + O(c^3), series of rho computed above", Normal[Series[dvRho, {c, 0, 2}]], 1 - 2 c/3 + c^2/2];

dvS93 = check["(S84) large-q Lindhard: coefficients of x^2, x^3, x^4 with x = 1/q (g4)",
  Module[{ser}, ser = Series[x/2 k Log[(vp + k x + k^2 x^2)/(vp - k x + k^2 x^2)], {x, 0, 4}];
    Simplify[Table[Integrate[SeriesCoefficient[ser, n], {k, 0, 1}, Assumptions -> vp > 0], {n, 2, 4}] /. vp -> u^2 + 1/4]],
  {4/(3 (1 + 4 u^2)), 0, 4/(3 (1 + 4 u^2)) (-4 (12 u^2 - 1)/(5 (1 + 4 u^2)^2))}, "TestPoint" -> {u -> 3/10}];
dvQs4 = check["(20) omega_p = q_*^2/2 and omega_p^2 = q_TF^2/3 solved for q_*^4: q_*^4 = 4 q_TF^2/3",
  qs4 /. First[Solve[qs4/4 == qTF2/3, qs4]], 4 qTF2/3, "Assumptions" -> qTF2 > 0, "TestPoint" -> {qTF2 -> 2}];
dvCHat = check["(20),(21) c = q_TF^2 L/q^2 with the coefficient of x^2 = q^-2 of L computed in (S84) and q_TF^2 from the q_*^4 above: c_hat = (q_*/q)^4/(1 + 4 u^2)",
  Simplify[(qTF2/q^2) (First[dvS93]/q^2) /. First[Solve[qs^4 == dvQs4, qTF2]], qs > 0 && q > 0], (qs/q)^4/(1 + 4 u^2),
  "Assumptions" -> qs > 0 && q > 0, "TestPoint" -> {qs -> 5/2, q -> 3, u -> 1/5}];

Module[{zr = stripCond[z /. First[Solve[r == Sqrt[1 + z^2] - z && r > 0, z]]], asm = 0 < r < 1, tp = {r -> 3/7}},
  check["(S2) r = Sqrt[1 + z^2] - z solved for z gives z = (1 - r^2)/(2 r); the solution is used below", zr, (1 - r^2)/(2 r), "Assumptions" -> asm, "TestPoint" -> tp];
  check["(S2) Sqrt[1 + z^2] = (1 + r^2)/(2 r)", Sqrt[1 + zr^2], (1 + r^2)/(2 r), "Assumptions" -> asm, "TestPoint" -> tp];
  check["(S2) a = z/Sqrt[1 + z^2] = (1 - r^2)/(1 + r^2)", zr/Sqrt[1 + zr^2], (1 - r^2)/(1 + r^2), "Assumptions" -> asm, "TestPoint" -> tp];
  check["(S2) h_2 = 2 r^2/(1 + r^2)", h2D[zr], 2 r^2/(1 + r^2), "Assumptions" -> asm, "TestPoint" -> tp];
  check["(S2) |dz| = (1 + r^2)/(2 r^2) dr", -D[zr, r], (1 + r^2)/(2 r^2), "Assumptions" -> asm, "TestPoint" -> tp];
  check["(text SI Sec. S1) h_4 = r^2", h4D[zr], r^2, "Assumptions" -> asm, "TestPoint" -> tp];
  dvS3h2 = check["(S3) 2 h_2^3 |dz/dr| = 8 r^4/(1 + r^2)^2", 2 h2D[zr]^3 (1 + r^2)/(2 r^2), 8 r^4/(1 + r^2)^2, "Assumptions" -> asm, "TestPoint" -> tp];
  dvS3h4 = check["(S3) 2 h_4^(5/3) |dz/dr| = r^(4/3) + r^(10/3)", 2 h4D[zr]^(5/3) (1 + r^2)/(2 r^2), r^(4/3) + r^(10/3), "Assumptions" -> asm, "TestPoint" -> tp]];
check["(S3) Int_0^1 of the r integrand computed above = 10 - 3 pi", Integrate[dvS3h2, {r, 0, 1}], 10 - 3 Pi];
dvIntH2cube = check["(S3) Int_R h_2^3 dz = 10 - 3 pi directly in z", 2 Integrate[h2D[z]^3, {z, 0, Infinity}], 10 - 3 Pi];
check["(S3) Int_0^1 of the r integrand of h_4^(5/3) computed above = 60/91", Integrate[dvS3h4, {r, 0, 1}], 60/91];
dvIntH4 = check["(S3) Int_R h_4^(5/3) dz = 60/91 directly in z", 2 Integrate[h4D[z]^(5/3), {z, 0, Infinity}], 60/91];
dvIntHsq = check["(text SI Sec. S1) d = 3: Int_R h^2 dz = (2 pi/3)(1 - log 2)", 2 Integrate[h[z]^2, {z, 0, Infinity}], 2 Pi (1 - Log[2])/3];

check["(S86) 1/(1+lambda c_hat) with c_hat as computed in (21), rewritten with nu = q^2 u: 1 - (lambda q_*^4/4)/(nu^2+omega^2), omega^2 = (q^4+lambda q_*^4)/4",
  Simplify[1/(1 + lam dvCHat) /. u -> nu/q^2, q > 0 && qs > 0 && lam > 0], 1 - (lam qs^4/4)/(nu^2 + (q^4 + lam qs^4)/4),
  "Assumptions" -> q > 0 && qs > 0 && lam > 0, "TestPoint" -> {nu -> 4/3, lam -> 1/2, qs -> 5, q -> 2}];
{dvLam21, dvTh21} = check["(21) Lambda and theta read off 1 - 1/(1 + lambda c_hat) (c_hat computed above) at nu = q z: theta Lambda^2/(Lambda^2 + z^2) with Lambda = (q/2) Sqrt[1 + lambda (q_*/q)^4], theta = lambda q_*^4/(q^4 + lambda q_*^4)",
  Module[{f = Together[1 - 1/(1 + lam dvCHat) /. u -> z/q], den, l2},
    den = Denominator[f]; l2 = Coefficient[den, z, 0]/Coefficient[den, z, 2];
    {Simplify[Sqrt[l2], q > 0 && qs > 0 && lam > 0], Simplify[Numerator[f]/(Coefficient[den, z, 2] l2), q > 0 && qs > 0 && lam > 0]}],
  {q/2 Sqrt[1 + lam (qs/q)^4], lam qs^4/(q^4 + lam qs^4)},
  "Assumptions" -> q > 0 && qs > 0 && lam > 0, "TestPoint" -> {q -> 2, qs -> 5, lam -> 1/2}];

check["(S69) Int_0^1 q^-s q^a log^m(1/q) dq = m!/(1+a-s)^(m+1) for general m >= 0",
  Integrate[q^(a - s) Log[1/q]^m, {q, 0, 1}, Assumptions -> a - s > -1 && m >= 0], m!/(1 + a - s)^(m + 1),
  "Assumptions" -> a - s > -1 && m >= 0, "TestPoint" -> {a -> 1, s -> 1/2, m -> 3}];
check["(S69) Int_0^1 q^-s q^a log^m(1/q) dq = m!/(1+a-s)^(m+1), m = 0..3",
  Table[Integrate[q^(a - s) Log[1/q]^m, {q, 0, 1}, Assumptions -> Re[a - s] > -1], {m, 0, 3}],
  Table[m!/(1 + a - s)^(m + 1), {m, 0, 3}], "TestPoint" -> {a -> 1, s -> 1/2}];
check["(S69) Int_1^inf t^(s-1-beta) dt = 1/(beta-s)",
  Integrate[t^(s - 1 - be), {t, 1, Infinity}, Assumptions -> be > s], 1/(be - s), "TestPoint" -> {be -> 2, s -> 1/2}];

check["(SI Prop. S2) toy case: k = exp(-q), g = 1/(1+t^2); coefficient of eps^2 log(eps) is (-1)^(0+1) kappa gamma = 1",
  Module[{ser}, ser = Normal[Series[Integrate[Exp[-q]/(1 + q^2/eps^2), {q, 0, Infinity}, Assumptions -> eps > 0], {eps, 0, 2}]];
    Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> LL, eps^2 LL]], 1];

check["(S70) -Res_{s=beta} (-1)^m m! kappa gamma eps^s/(s-beta)^(m+2) = (-1)^(m+1) kappa gamma eps^beta log^(m+1)(eps)/(m+1), m = 0..3",
  Table[-Residue[(-1)^m m! ka ga eps^s/(s - be)^(m + 2), {s, be}], {m, 0, 3}],
  Table[(-1)^(m + 1) ka ga eps^be Log[eps]^(m + 1)/(m + 1), {m, 0, 3}], "Assumptions" -> eps > 0,
  "TestPoint" -> {ka -> 2, ga -> 3, be -> 2, eps -> 1/10}];check["(S70) the same for general m >= 0: the residue of F(s)/(s - beta)^(m+2) is the coefficient of (s - beta)^(m+1) in F(s) = (-1)^m m! kappa gamma eps^s",
  -(-1)^m m! ka ga SeriesCoefficient[eps^s, {s, be, m + 1}, Assumptions -> m >= 0 && Element[m, Integers]],
  (-1)^(m + 1) ka ga eps^be Log[eps]^(m + 1)/(m + 1), "Assumptions" -> eps > 0 && m >= 0 && Element[m, Integers],
  "TestPoint" -> {ka -> 2, ga -> 3, be -> 2, eps -> 1/10, m -> 2}];

check["(S68) moment rule: 2 Int_0^1 lambda (lambda r)^p log^2(lambda r) dlambda",
  Simplify[PowerExpand[2 Integrate[lam (lam r)^p Log[lam r]^2, {lam, 0, 1}, Assumptions -> p > -2 && r > 0]], r > 0 && p > -2],
  2 r^p (Log[r]^2/(p + 2) - 2 Log[r]/(p + 2)^2 + 2/(p + 2)^3), "Assumptions" -> r > 0 && p > -2,
  "TestPoint" -> {r -> 3, p -> 1}];

check["(Tab. 1) low density: theta q/Lambda with Lambda = (q/2) Sqrt[1 + lambda (q_*/q)^4] and theta = lambda q_*^4/(q^4 + lambda q_*^4) of (21), q = q_* t: 2 lambda t^2/(t^4+lambda)^(3/2)",
  Simplify[(dvTh21 q/dvLam21) /. q -> t qs, t > 0 && lam > 0 && qs > 0], 2 lam t^2/(t^4 + lam)^(3/2),
  "Assumptions" -> t > 0 && lam > 0 && qs > 0, "TestPoint" -> {t -> 2/3, lam -> 1/2, qs -> 3}];
check["(S67) large-width term: 4 pi Int (Lambda^2/(Lambda^2+nu^2/q^2) - 1) A/nu^2 dnu/pi = -2 pi A/(q Lambda)",
  4 Pi Integrate[(L^2/(L^2 + nu^2/q^2) - 1) A/nu^2, {nu, 0, Infinity}, Assumptions -> L > 0 && q > 0]/Pi,
  -2 Pi A/(q L), "Assumptions" -> L > 0 && q > 0, "TestPoint" -> {L -> 3, q -> 4/5, A -> 1}];
dvMel25 = check["(25) Mellin transform of g(t) = 2 lambda t^2/(t^4+lambda)^(3/2)",
  MellinTransform[2 lam t^2/(t^4 + lam)^(3/2), t, s, Assumptions -> lam > 0],
  lam^(s/4) Gamma[(s + 2)/4] Gamma[(4 - s)/4]/Sqrt[Pi], "Assumptions" -> lam > 0, "TestPoint" -> {lam -> 7/10, s -> -3/2}];
check["(25) its residue at s = -2 is 2/sqrt(lambda)",
  Residue[dvMel25, {s, -2}], 2/Sqrt[lam],
  "Assumptions" -> lam > 0, "TestPoint" -> {lam -> 7/10}];

Module[{beta4 = (Zeta[4, 1/4] - Zeta[4, 3/4])/4^4},
  check["(2) Dirichlet beta(4) = Sum_(k>=0) (-1)^k/(2k+1)^4 = (zeta(4, 1/4) - zeta(4, 3/4))/4^4",
    Sum[(-1)^k/(2 k + 1)^4, {k, 0, Infinity}], beta4];
  numcheck["(2) a_0 = 2(log 2 - 1) + 2 G - 16 beta(4)/pi^2", 2 (Log[2] - 1) + 2 Catalan - 16 beta4/Pi^2, "-0.38499"]];
numcheck["(2) a_l = -2 Sqrt[2] [10/(3 pi) - 1]", -2 Sqrt[2] (10/(3 Pi) - 1), "-0.17263"];

check["(Fig. 1) open circles: s = 2 with N = 2 (d = 2) and N = 1 (d = 3), s = d + 1 with N = 2 (d = 3) and N = 3 (d = 2), s = 4 with N = 4 (d = 2) and N = 2 (d = 3)",
  {d /. First[Solve[2 == 2 (d - 1), d]], d /. First[Solve[2 == 1 (d - 1), d]], d /. First[Solve[d + 1 == 2 (d - 1), d]],
   d /. First[Solve[d + 1 == 3 (d - 1), d]], d /. First[Solve[4 == 4 (d - 1), d]], d /. First[Solve[4 == 2 (d - 1), d]]}, {2, 3, 3, 2, 2, 3}];
check["(Fig. 1) d = 2, s = 4 (insertion of order q^2, j = 0, m = 2): E_0 q_TF^(d+3) is proportional to r_s^3",
  Exponent[PowerExpand[(2/(alphaD[2] rs)^2) (calB[2] alphaD[2] rs)^5], rs], 3];

dvRoute = check["(S31b) route: Int dz/((z^2+a^2)(z^2+b^2)) = pi/(2ab(a+b))",
  Integrate[1/((z^2 + a^2) (z^2 + b^2)), {z, 0, Infinity}, Assumptions -> a > 0 && b > 0],
  Pi/(2 a b (a + b)), "Assumptions" -> a > 0 && b > 0, "TestPoint" -> {a -> 1/3, b -> 2}];
check["(S31b) route: Int Int ab/(a+b) = 2(1-log2)/3", Integrate[a b/(a + b), {a, 0, 1}, {b, 0, 1}], 2 (1 - Log[2])/3];
check["(S31b) route: with h = Int_0^1 a^2 da/(a^2 + z^2), Int_0^inf h^2 dz = Int Int a^2 b^2 (z integral computed above) da db = pi(1-log2)/3",
  Integrate[a^2 b^2 dvRoute, {a, 0, 1}, {b, 0, 1}], Pi (1 - Log[2])/3];
check["(S31b) Int_0^inf h^2 dz = pi(1-log2)/3 directly", Integrate[h[z]^2, {z, 0, Infinity}], Pi (1 - Log[2])/3];

check["(27) RPA coupling integral: Int_0^1 [c - c/(1 + lam c)] dlam = -[log(1 + c) - c]",
  Integrate[-c/(1 + lam c) + c, {lam, 0, 1}, Assumptions -> c > 0], -(Log[1 + c] - c),
  "Assumptions" -> c > 0, "TestPoint" -> {c -> 7/10}];
check["(27) prefactor: (1/2)(2 pi)^-d (2 pi)^-1 times the volume (2 pi)^d/(2 V_d) per electron times A_d = d/(8 pi)",
  (1/2) (2 Pi)^-d (2 Pi)^-1 ((2 Pi)^d/(2 volV[d])) calA[d], d/(8 Pi), "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(27) A_d = d V_d", calA[d], d volV[d], "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(27) nu = q z: phi(nu) dnu = q phi(q z) dz for any phi, and for an even integrand Int_0^Infinity = (1/2) Int_R",
  (phi[nu] /. nu -> q z) D[q z, z], q phi[q z], "Assumptions" -> q > 0];
check["(27) Int_0^Infinity dnu = (q/2) Int_R dz for an even integrand (test function 1/(1 + nu^2/q^2)^2)",
  Integrate[1/(1 + nu^2/q^2)^2, {nu, 0, Infinity}, Assumptions -> q > 0], (q/2) Integrate[1/(1 + z^2)^2, {z, -Infinity, Infinity}],
  "Assumptions" -> q > 0, "TestPoint" -> {q -> 3/7}];
check["(27) d = 3: (d/(8 pi)) E_0 = 12/(pi^3 q_TF^4) Ry, with E_0 of Sec. 2.2 and r_s from (6)", (3/(8 Pi)) E0d[3, (rs /. First[Solve[qTF^2 == calB[3] alphaD[3] rs, rs]])], 12/(Pi^3 qTF^4),
  "Assumptions" -> qTF > 0, "TestPoint" -> {qTF -> 2/5}];

dvD27 = check["(28) N_r = (d + 1)/(d - 1) solved for d: d = (N_r + 1)/(N_r - 1)", d /. First[Solve[nn == (d + 1)/(d - 1), d]], (nn + 1)/(nn - 1),
  "Assumptions" -> nn > 1, "TestPoint" -> {nn -> 5/2}];
check["(28) resonant dimensions 3, 2, 5/3, 3/2 for N_r = 2, 3, 4, 5", Table[dvD27 /. nn -> k, {k, 2, 5}], {3, 2, 5/3, 3/2}];
check["(text Sec. 2.8) kernel pole s = d + 1 meets the screening pole s = N (d - 1) at N = N_r", nn /. First[Solve[d + 1 == nn (d - 1), nn]], nR[d],
  "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(text Sec. 2.8) tail log(1 + x) - x = Sum_{N >= 2} (-1)^(N+1) x^N/N", Normal[Series[Log[1 + x] - x, {x, 0, 7}]], Sum[(-1)^(k + 1) x^k/k, {k, 2, 7}]];
check["(text Sec. 2.8) tail coefficients for general N >= 2: coefficient of x^N in log(1 + x) - x is (-1)^(N+1)/N",
  SeriesCoefficient[Log[1 + x] - x, {x, 0, n}, Assumptions -> Element[n, Integers] && n >= 2], (-1)^(n + 1)/n,
  "Assumptions" -> Element[n, Integers] && n >= 2, "TestPoint" -> {n -> 5}];
check["(text Sec. 2.8) helper lmxF for log(1 + x) - x of (27): -(x^2/2) 2F1(1, 2; 3; -x) (not printed)", lmxF[x], Log[1 + x] - x, "Assumptions" -> x > 0, "TestPoint" -> {x -> 7/10}];
check["(text Sec. 2.8) c of (8) with L -> h_d and t = q/q_TF is h_d/t^(d-1)", (qTF^(d - 1) hh/q^(d - 1)) /. q -> t qTF, hh/t^(d - 1),
  "Assumptions" -> qTF > 0 && t > 0 && hh > 0 && d > 1, "TestPoint" -> {qTF -> 2/5, t -> 3/7, hh -> 3/5, d -> 37/10}];

dvCoef28d3 = check["(29) d = 3: coefficient of eps^4 log(eps) in Int_0^1 q^3 [log(1 + h eps^2/q^2) - h eps^2/q^2] dq is (-1)^2 h^2/2",
  Module[{ser}, ser = Normal[Series[Integrate[q^3 (Log[1 + hh eps^2/q^2] - hh eps^2/q^2), {q, 0, 1},
        Assumptions -> 0 < eps < 1 && hh > 0, GenerateConditions -> False], {eps, 0, 4}]];
    Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> logEps, eps^4 logEps]],
  hh^2/2, "Assumptions" -> hh > 0, "TestPoint" -> {hh -> 3/5}];
dvCoef28d2 = check["(29) d = 2: coefficient of eps^3 log(eps) in Int_0^1 q^2 [log(1 + h eps/q) - h eps/q] dq is (-1)^3 h^3/3",
  Module[{ser}, ser = Normal[Series[Integrate[q^2 (Log[1 + hh eps/q] - hh eps/q), {q, 0, 1},
        Assumptions -> 0 < eps < 1 && hh > 0, GenerateConditions -> False], {eps, 0, 3}]];
    Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> logEps, eps^3 logEps]],
  -hh^3/3, "Assumptions" -> hh > 0, "TestPoint" -> {hh -> 3/5}];
check["(text Sec. 2.8) d = 2: the leading term is -h^2 eps^2/2 from the screening pole N = 2 at s = 2",
  Module[{ser}, ser = Normal[Series[Integrate[q^2 (Log[1 + hh eps/q] - hh eps/q), {q, 0, 1},
        Assumptions -> 0 < eps < 1 && hh > 0, GenerateConditions -> False], {eps, 0, 3}]];
    {Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> logEps, eps, 1], Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> logEps, eps, 2]}],
  {0, -hh^2/2}, "Assumptions" -> hh > 0, "TestPoint" -> {hh -> 3/5}];
check["(29) log q_TF = log r_s/(d - 1) + O(1) from q_TF^(d-1) = a_TF r_s of (6)",
  PowerExpand[Log[(aT rs)^(1/(d - 1))]] - Log[rs]/(d - 1), Log[aT]/(d - 1),
  "Assumptions" -> aT > 0 && rs > 0 && d > 1, "TestPoint" -> {aT -> 3/5, rs -> 3/10, d -> 37/10}];
check["(29) E_0 q_TF^(d+1) is proportional to r_s^((3-d)/(d-1)): exponent of r_s from the definitions", Simplify[Exponent[PowerExpand[E0d[d, rs] qTFd[d, rs]^(d + 1)], rs]], (3 - d)/(d - 1),
  "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];

dvJac5 = check["(S5) substitution u = h t^(1-d): t^d |dt/du| = h^(N_r) u^(-N_r - 1)/(d - 1)",
  Module[{tt = (hh/u)^(1/(d - 1))}, tt^d (-D[tt, u])], hh^nR[d] u^(-nR[d] - 1)/(d - 1),
  "Assumptions" -> hh > 0 && u > 0 && d > 3, "TestPoint" -> {hh -> 3/5, u -> 7/10, d -> 37/10}];
check["(S5) integration by parts with f = log(1 + u) - u and g' = u^(-N-1): g = -u^(-N)/N (computed), the boundary terms f g vanish at u -> 0 and u -> Infinity for 1 < N < 2, and -f' g = -u^(1-N)/(N (1 + u))",
  Module[{g = Integrate[u^(-nn - 1), u], f = Log[1 + u] - u},
    {g, Limit[f g, u -> 0, Direction -> "FromAbove", Assumptions -> 1 < nn < 2], Limit[f g, u -> Infinity, Assumptions -> 1 < nn < 2], Simplify[-D[f, u] g]}],
  {-u^(-nn)/nn, 0, 0, -u^(1 - nn)/(nn (1 + u))}, "Assumptions" -> u > 0 && 1 < nn < 2, "TestPoint" -> {u -> 7/10, nn -> 5/3}];
dvEulerB = check["(text SI Sec. S2) Int_0^Infinity u^(a-1)/(1 + u) du = pi/Sin[pi a] for 0 < a < 1",
  Integrate[u^(a - 1)/(1 + u), {u, 0, Infinity}, Assumptions -> 0 < a < 1], Pi/Sin[Pi a],
  "Assumptions" -> 0 < a < 1, "TestPoint" -> {a -> 37/100}];
check["(S5) the remaining integral -(1/N) Int u^(1-N)/(1 + u) du is the integral above with a = 2 - N: pi/(N Sin[pi N])",
  -(1/nn) (dvEulerB /. a -> 2 - nn), Pi/(nn Sin[Pi nn]), "Assumptions" -> 1 < nn < 2, "TestPoint" -> {nn -> 5/3}];
dvMellinU = check["(S5) Int_0^Infinity u^(-N-1) [log(1 + u) - u] du = pi/(N Sin[pi N]) for 1 < N < 2, directly",
  Integrate[u^(-nn - 1) (Log[1 + u] - u), {u, 0, Infinity}, Assumptions -> 1 < nn < 2], Pi/(nn Sin[Pi nn]),
  "Assumptions" -> 1 < nn < 2, "TestPoint" -> {nn -> 5/3}];
check["(S5) whole chain at d = 4: Int_0^Infinity t^4 [log(1 + h t^-3) - h t^-3] dt = h^(N_r) pi/((d - 1) N_r Sin[pi N_r]) with N_r = 5/3",
  Integrate[t^4 (Log[1 + hh t^-3] - hh t^-3), {t, 0, Infinity}, Assumptions -> hh > 0], hh^nR[4] Pi/((4 - 1) nR[4] Sin[Pi nR[4]]),
  "Assumptions" -> hh > 0, "TestPoint" -> {hh -> 3/5}];
check["(S5) whole chain at d = 7/2: Int_0^Infinity t^(7/2) [log(1 + h t^(-5/2)) - h t^(-5/2)] dt = h^(N_r) pi/((d - 1) N_r Sin[pi N_r]) with N_r = 9/5",
  Integrate[t^(7/2) (Log[1 + hh t^(-5/2)] - hh t^(-5/2)), {t, 0, Infinity}, Assumptions -> hh > 0], hh^nR[7/2] Pi/((7/2 - 1) nR[7/2] Sin[Pi nR[7/2]]),
  "Assumptions" -> hh > 0, "TestPoint" -> {hh -> 3/5}];

dvPref29 = check["(30) prefactor at d = 3: (d/(8 pi)) E_0 q_TF^(d+1) times the coefficient of eps^4 log(eps) computed for (29), divided by h^2, times 1/(d - 1) from log q_TF = log r_s/(d - 1): 3/pi^3",
  Simplify[(3/(8 Pi)) E0d[3, rs] qTFd[3, rs]^4 (dvCoef28d3/hh^2) (1/(3 - 1)), rs > 0], 3/Pi^3, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
check["(30) A = prefactor times Int_R h^2 dz (both computed above) = 2 (1 - log 2)/pi^2", dvPref29 dvIntHsq, 2 (1 - Log[2])/Pi^2];
numcheck["(30) Gell-Mann--Brueckner coefficient A = 0.06218", aGB, "0.06218"];

check["(31) d = 2: E_0 q_TF^3 = 8 Sqrt[2] r_s Ry", (2/(alphaD[2] rs)^2) (calB[2] alphaD[2] rs)^3, 8 Sqrt[2] rs,
  "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
dvPref30 = check["(31) prefactor at d = 2: (d/(8 pi)) E_0 q_TF^(d+1)/r_s times the coefficient of eps^3 log(eps) computed for (29), divided by h^3, times 1/(d - 1): -2 Sqrt[2]/(3 pi)",
  Simplify[(2/(8 Pi)) E0d[2, rs] qTFd[2, rs]^3/rs (dvCoef28d2/hh^3) (1/(2 - 1)), rs > 0], -2 Sqrt[2]/(3 Pi), "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
check["(31) a_l = prefactor times Int_R h_2^3 dz (both computed above) = -2 Sqrt[2] (10/(3 pi) - 1)", dvPref30 dvIntH2cube, -2 Sqrt[2] (10/(3 Pi) - 1)];
check["(text SI Sec. S2) d = 2: (2 Sqrt[2]/pi)(-1/3)(10 - 3 pi) = -2 Sqrt[2] [10/(3 pi) - 1]",
  (2 Sqrt[2]/Pi) (-1/3) (10 - 3 Pi), -2 Sqrt[2] (10/(3 Pi) - 1)];

check["(text SI Sec. S2) E_0 q_TF^(d+1) = 2 alpha_d^(N_r - 2) B_d^(N_r) r_s^(N_r - 2) Ry",
  (2/(alphaD[d] rs)^2) (calB[d] alphaD[d] rs)^nR[d], 2 alphaD[d]^(nR[d] - 2) calB[d]^nR[d] rs^(nR[d] - 2),
  "Assumptions" -> d > 3 && rs > 0, "TestPoint" -> {d -> 37/10, rs -> 3/10}];
dvA31 = check["(32) (d/(8 pi)) E_0 q_TF^(d+1) times 1/(d - 1) (Jacobian of (S5)) times the u integral computed in (S5) gives d alpha_d^(N_r - 2) B_d^(N_r) r_s^(N_r - 2)/(4 (d + 1) Sin[pi N_r]) (per unit Int h_d^(N_r) dz)",
  (d/(8 Pi)) E0d[d, rs] qTFd[d, rs]^(d + 1) Simplify[dvJac5/(hh^nR[d] u^(-nR[d] - 1)), hh > 0 && u > 0 && d > 3] (dvMellinU /. nn -> nR[d]),
  d alphaD[d]^(nR[d] - 2) calB[d]^nR[d] rs^(nR[d] - 2)/(4 (d + 1) Sin[Pi nR[d]]), "Assumptions" -> d > 3 && rs > 0, "TestPoint" -> {d -> 37/10, rs -> 3/10}];
check["(32) the power N_r - 2 of r_s is (3 - d)/(d - 1)", nR[d] - 2, (3 - d)/(d - 1), "Assumptions" -> d > 3, "TestPoint" -> {d -> 37/10}];
check["(text SI Sec. S2) d = 4: A_4 = 2 pi^2, C_4 = 4 pi^2, B_4 = 1, alpha_4 = 2^(-5/4), N_r = 5/3, Sin[5 pi/3] = -Sqrt[3]/2",
  {calA[4], calC[4], calB[4], alphaD[4], nR[4], Sin[5 Pi/3]}, {2 Pi^2, 4 Pi^2, 1, 2^(-5/4), 5/3, -Sqrt[3]/2}];
check["(32) A_4 = -24 2^(5/12)/(91 Sqrt[3]) Ry from (32) as computed above and Int h_4^(5/3) dz computed in (S3)",
  Simplify[PowerExpand[Simplify[dvA31 /. rs -> 1, d > 3]] /. d -> 4] dvIntH4, -24 2^(5/12)/(91 Sqrt[3])];
numcheck["(text Sec. 2.8) A_4 = -0.20325", -24 2^(5/12)/(91 Sqrt[3]), "-0.20325"];

check["(text SI Sec. S2) N_r = 2 - (d - 3)/(d - 1)", nR[d], 2 - (d - 3)/(d - 1), "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];
check["(text SI Sec. S2) Sin[pi N_r] = -(pi/2)(d - 3) + O((d - 3)^2)",
  {SeriesCoefficient[Sin[Pi nR[d]], {d, 3, 0}], SeriesCoefficient[Sin[Pi nR[d]], {d, 3, 1}]}, {0, -Pi/2}];
check["(text SI Sec. S2) remaining factors at d = 3: 3 B_3^2 Int_R h^2 dz/16 = (2/pi)(1 - log 2)",
  Limit[PowerExpand[Simplify[(dvA31 /. rs -> 1) Sin[Pi nR[d]], d > 3]], d -> 3] dvIntHsq, (2/Pi) (1 - Log[2])];
check["(text Sec. 2.8) (d - 3) A_d -> -2 A as d -> 3, with (32) as computed and Int_R h^2 dz computed in Sec. S1 (N_r -> 2)",
  Limit[(d - 3) PowerExpand[Simplify[dvA31 /. rs -> 1, d > 3]], d -> 3] dvIntHsq, -2 aGB];
check["(text Sec. 2.8) r_s^((3-d)/(d-1)) = 1 - (1/2)(d - 3) log r_s + O((d - 3)^2)",
  {SeriesCoefficient[rs^((3 - d)/(d - 1)), {d, 3, 0}], SeriesCoefficient[rs^((3 - d)/(d - 1)), {d, 3, 1}]}, {1, -Log[rs]/2},
  "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
numcheck["(text SI Sec. S2) -A = -0.0621814", -aGB, "-0.0621814"];
check["(text SI Sec. S2) explicit r_s^0 term of Eq. (21) of Schlesier et al.: residue of D Int R^2 du/(pi^3 (D - 3)) Ha at D = 3 is +2 A in Ry",
  2 Residue[dd dvIntHsq/(Pi^3 (dd - 3)), {dd, 3}], 2 aGB];

check["(text SI Sec. S2) B_d ~ 2 Sqrt[2]/Sqrt[pi d] as d -> Infinity", Limit[calB[d] Sqrt[Pi d]/(2 Sqrt[2]), d -> Infinity], 1];
check["(text SI Sec. S2) 1/alpha_d ~ d/e (Stirling)", Limit[1/(alphaD[d] d), d -> Infinity], 1/E];
check["(text SI Sec. S2) Sin[pi N_r] ~ -2 pi/d", Limit[d Sin[Pi nR[d]], d -> Infinity], -2 Pi];
check["(text SI Sec. S2) alpha_d^(N_r - 2) B_d^(N_r) ~ alpha_d^-1 B_d: (alpha_d B_d)^(2/(d-1)) -> 1",
  Limit[(alphaD[d] calB[d])^(2/(d - 1)), d -> Infinity], 1];
check["(text SI Sec. S2) Int_R x^2/(x^2 + y^2) dy = pi |x|",
  Integrate[xg^2/(xg^2 + y^2), {y, -Infinity, Infinity}, Assumptions -> Element[xg, Reals] && xg != 0], Pi Abs[xg],
  "Assumptions" -> Element[xg, Reals] && xg != 0, "TestPoint" -> {xg -> -7/10}];
check["(text SI Sec. S2) pi <|x|> = Sqrt[2 pi] for a standard normal x",
  Pi Expectation[Abs[xg], Distributed[xg, NormalDistribution[0, 1]]], Sqrt[2 Pi]];
dvAvgAbs = Last[check["(text SI Sec. S2) Int_R h_d dz = pi Int_-1^1 P_d(mu) |mu| dmu exactly (exchange the z and mu integrals and use Int_R mu^2/(mu^2 + z^2) dz = pi |mu|), and Int_-1^1 P_d(mu) |mu| dmu = 2/((d - 1) B(1/2, (d-1)/2))",
  {Integrate[mu^2/(mu^2 + z^2), {z, -Infinity, Infinity}, Assumptions -> mu > 0], 2 Integrate[mu (1 - mu^2)^((d - 3)/2), {mu, 0, 1}, Assumptions -> d > 1]/Beta[1/2, (d - 1)/2]},
  {Pi mu, 2/((d - 1) Beta[1/2, (d - 1)/2])}, "Assumptions" -> d > 1 && mu > 0, "TestPoint" -> {d -> 37/10, mu -> 3/10}]];
check["(text SI Sec. S2) Int_R h_d dz ~ Sqrt[2 pi/d]: Sqrt[d] pi Int P_d(mu) |mu| dmu -> Sqrt[2 pi] with that integral computed above",
  Limit[Sqrt[d] Pi dvAvgAbs, d -> Infinity], Sqrt[2 Pi]];
check["(text SI Sec. S2) A_d/d -> -1/(2 e pi): (32) as computed, divided by d r_s^(N_r - 2), times Int_R h_d dz = pi Int P_d(mu) |mu| dmu (h_d^(N_r) -> h_d since N_r -> 1)",
  Limit[(PowerExpand[Simplify[dvA31 /. rs -> 1, d > 3]]/d) Pi dvAvgAbs, d -> Infinity], -1/(2 E Pi)];
numcheck["(text Sec. 2.8) -1/(2 e pi) = -0.05855", -1/(2 E Pi), "-0.05855"];
numcheck["(text SI Sec. S2) -1/(2 e pi) = -0.058550", -1/(2 E Pi), "-0.058550"];
check["(text SI Sec. S2) large-d limit -D/(2 e pi^2) Ha of Schlesier et al. is pi/2 smaller in magnitude than -d/(4 e pi) Ha",
  (d/(4 E Pi))/(d/(2 E Pi^2)), Pi/2];

Module[{sigmaD, rSch},
  sigmaD[dd_] := Sum[(-1)^n/(n (n + nR[dd])) - (-1)^n/(n (nR[dd] - n)), {n, 1, Infinity}] - 1/nR[dd]^2;
  rSch[dd_, uu_] := Gamma[3/2] Gamma[(dd - 1)/2]/(2 uu^2 Gamma[dd/2 + 1]) Hypergeometric2F1[1, 3/2, dd/2 + 1, -1/uu^2];
  check["(text SI Sec. S2) series Sigma_D of Schlesier et al. = -pi/(N_r Sin[pi N_r]) at d = 4", sigmaD[4], -Pi/(nR[4] Sin[Pi nR[4]])];
  check["(text SI Sec. S2) series Sigma_D of Schlesier et al. = -pi/(N_r Sin[pi N_r]) at d = 5", sigmaD[5], -Pi/(nR[5] Sin[Pi nR[5]])];
  check["(text SI Sec. S2) series Sigma_D of Schlesier et al. = -pi/(N_r Sin[pi N_r]) at d = 6", sigmaD[6], -Pi/(nR[6] Sin[Pi nR[6]])];
  check["(text SI Sec. S2) series Sigma_D of Schlesier et al. = -pi/(N_r Sin[pi N_r]) at d = 9", sigmaD[9], -Pi/(nR[9] Sin[Pi nR[9]])];
  check["(text SI Sec. S2) series form of (9): the coefficient of y^(n+1) in 1 - 2F1(1, 1/2; c; y) equals that of -(y/(2c)) 2F1(1, 3/2; c + 1; y), so h_d(u) = 2F1(1, 3/2; d/2 + 1; -1/u^2)/(d u^2)",
    SeriesCoefficient[1 - Hypergeometric2F1[1, 1/2, cc, y], {y, 0, n + 1}, Assumptions -> Element[n, Integers] && n >= 0] -
      SeriesCoefficient[-(y/(2 cc)) Hypergeometric2F1[1, 3/2, cc + 1, y], {y, 0, n + 1}, Assumptions -> Element[n, Integers] && n >= 0], 0,
    "Assumptions" -> Element[n, Integers] && n >= 0 && cc > 0, "TestPoint" -> {n -> 3, cc -> 37/20}];
  check["(text SI Sec. S2) R_{(D-1)/2}(u, xi = 0) of Schlesier et al. = (pi B_d/4) h_d(u): with the series form of h_d the ratio is the constant Gamma(3/2) Gamma((d-1)/2) d/(2 Gamma(d/2 + 1)) = pi B_d/4",
    rSch[d, u]/(Hypergeometric2F1[1, 3/2, d/2 + 1, -1/u^2]/(d u^2)), Pi calB[d]/4,
    "Assumptions" -> d > 1 && u > 0, "TestPoint" -> {d -> 37/10, u -> 7/10}];
  check["(text SI Sec. S2) series Sigma_D of Schlesier et al. = -pi/(N_r Sin[pi N_r]) for general 1 < N < 2",
    Sum[(-1)^n/(n (n + nn)) - (-1)^n/(n (nn - n)), {n, 1, Infinity}] - 1/nn^2, -Pi/(nn Sin[Pi nn]),
    "Assumptions" -> 1 < nn < 2, "TestPoint" -> {nn -> 5/3}];
  check["(text SI Sec. S2) twice the leading coefficient of Eq. (21) of Schlesier et al. (Ha -> Ry) equals (32)",
    2 (-2 d (-Pi/(nR[d] Sin[Pi nR[d]]))/(Pi^3 (d - 1))) (Pi/(4 alphaD[d]))^((d - 3)/(d - 1)) (Pi calB[d]/4)^nR[d] ihN,
    d alphaD[d]^(nR[d] - 2) calB[d]^nR[d]/(4 (d + 1) Sin[Pi nR[d]]) ihN,
    "Assumptions" -> d > 3 && ihN > 0, "TestPoint" -> {d -> 37/10, ihN -> 3/5}]];

schlesierHa = <|4 -> -196/10000, 5 -> -285/10000, 6 -> -391/10000, 7 -> -509/10000, 8 -> -638/10000, 9 -> -773/10000|>;   (* Table I of Schlesier et al., Phys. Rev. B 102, 035123 (2020), xi = 0, in Hartree *)
numcheck["(Tab. 2) A_4 = -24 2^(5/12)/(91 Sqrt[3]) in closed form", -24 2^(5/12)/(91 Sqrt[3]), "-0.20325"];
numcheck["(Tab. 2) d = 4: Table I of Schlesier et al. converted to Ry", 2 schlesierHa[4], "-0.0392"];
numcheck["(Tab. 2) d = 4: ratio", (-24 2^(5/12)/(91 Sqrt[3]))/(2 schlesierHa[4]), "5.19"];
numcheck["(Tab. 2) d = 5: Table I of Schlesier et al. converted to Ry", 2 schlesierHa[5], "-0.0570"];
numcheck["(Tab. 2) d = 6: Table I of Schlesier et al. converted to Ry", 2 schlesierHa[6], "-0.0782"];
numcheck["(Tab. 2) d = 7: Table I of Schlesier et al. converted to Ry", 2 schlesierHa[7], "-0.1018"];
numcheck["(Tab. 2) d = 8: Table I of Schlesier et al. converted to Ry", 2 schlesierHa[8], "-0.1276"];
numcheck["(Tab. 2) d = 9: Table I of Schlesier et al. converted to Ry", 2 schlesierHa[9], "-0.1546"];

check["(33) first-order change of log(1+c) - c under c -> c + dc is dc [1/(1+c) - 1]",
  D[Log[1 + c + eps dc] - (c + eps dc), eps] /. eps -> 0, dc (1/(1 + c) - 1),
  "Assumptions" -> c > 0, "TestPoint" -> {c -> 37/100, dc -> 3/7}];
check["(33) closed-diagram weight: the coefficient of cf^(j-1) dc in (-1)^(j+1) c^j/j is j times that of cf^j, j = 2..6",
  Table[Coefficient[Expand[(-1)^(j + 1) (cf + dc)^j/j], cf^(j - 1) dc]/
    Coefficient[Expand[(-1)^(j + 1) (cf + dc)^j/j] /. dc -> 0, cf^j], {j, 2, 6}], Range[2, 6]];
check["(33) closed-diagram weight for general j >= 1: d/d(dc) of (-1)^(j+1) (c_f + dc)^j/j at dc = 0, divided by the coefficient (-1)^(j+1)/j of c_f^j, times c_f, is j",
  Simplify[cf (D[(-1)^(j + 1) (cf + dc)^j/j, dc] /. dc -> 0)/((-1)^(j + 1) cf^j/j), cf > 0], j, "Assumptions" -> cf > 0 && j >= 1, "TestPoint" -> {cf -> 3/7, j -> 5}];
check["(33) tail of the scaling function for general N >= 1: coefficient of x^N in 1/(1 + h x) - 1 is (-h)^N",
  SeriesCoefficient[1/(1 + hh x) - 1, {x, 0, n}, Assumptions -> Element[n, Integers] && n >= 1], (-hh)^n,
  "Assumptions" -> Element[n, Integers] && n >= 1, "TestPoint" -> {n -> 4, hh -> 3/5}];
check["(33) tail of the scaling function: 1/(1+c) - 1 = Sum_(N>=1) (-h)^N x^N with x = t^(1-d), through x^5",
  Normal[Series[1/(1 + hh x) - 1, {x, 0, 5}]], Sum[(-hh)^n x^n, {n, 1, 5}]];
check["(text SI Sec. S14.6) coupling average of the t^-2 tail: Int_0^1 lambda [(1 + lambda c)^-2 - 1] dlambda = -2c/3 + O(c^2)",
  Normal[Series[Integrate[lam ((1 + lam c)^-2 - 1), {lam, 0, 1}, Assumptions -> c > 0], {c, 0, 1}]], -2 c/3];

check["(text SI Sec. S14.6) window integral for the g3 term: coefficient of eps^6 log eps in eps^2 g Int_0^1 q^3 [-h eps^2/(q^2 + h eps^2)] dq is -h^2 g ((S70) with gamma = h^2, beta = 4)",
  Module[{ser}, ser = Normal[Series[eps^2 gs Integrate[q^3 (-hh eps^2/(q^2 + hh eps^2)), {q, 0, 1}, Assumptions -> eps > 0 && hh > 0], {eps, 0, 6}]];
    Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> LL, eps^6 LL]], -hh^2 gs, "TestPoint" -> {hh -> 3/5, gs -> -1/7}];
dvIntHg3 = check["(text SI Sec. S14.6) Int_R h^2 g3 dz = -Pi (Pi^2 - 6)/288", 2 Integrate[h[z]^2 g3[z], {z, 0, Infinity}], -Pi (Pi^2 - 6)/288];
dvC36b = check["(text SI Sec. S14.6) assembly of the g3 term: (d/(8 Pi)) E_0 q_TF^4 (definitions) times -Int h^2 g3 (computed above) times the coefficient of r_s log r_s in q_TF^2 log q_TF (definitions) = alphaL (Pi^2 - 6)/(12 Pi^3)",
  Simplify[(3/(8 Pi)) E0d[3, rs] qTFd[3, rs]^4, rs > 0] (-dvIntHg3) Coefficient[Expand[PowerExpand[qTFd[3, rs]^2 Log[qTFd[3, rs]]]] /. Log[rs] -> LL, rs LL],
  alphaL (Pi^2 - 6)/(12 Pi^3), "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
numcheck["(text Sec. 2.9) 0.00542", alphaL (Pi^2 - 6)/(12 Pi^3), "0.00542"];
numcheck["(text SI Sec. S3) C = 0.0184584246606864", alphaL (7 Pi^2/6 - 12 Log[2] - 1)/(2 Pi^3), "0.0184584246606864"];

check["(34),(S7) prefactor: C^2 A_d^2/(2 Pi)^(2d) = qTF^(d-1) alpha_d rs B_d/4 with C = C_d alpha_d rs, qTF^(d-1) = B_d alpha_d rs",
  (calC[d] ad rs)^2 calA[d]^2/(2 Pi)^(2 d), (calB[d] ad rs) ad rs calB[d]/4,
  "Assumptions" -> d > 1 && ad > 0 && rs > 0, "TestPoint" -> {d -> 37/10, ad -> 3/7, rs -> 2/5}];
check["(S7) qTF^(d-1) alpha_d rs = 2 C^2 A_d/((2 Pi)^d C_d)",
  (calB[d] ad rs) ad rs, 2 (calC[d] ad rs)^2 calA[d]/((2 Pi)^d calC[d]),
  "Assumptions" -> d > 1 && ad > 0 && rs > 0, "TestPoint" -> {d -> 37/10, ad -> 3/7, rs -> 2/5}];
check["(34),(S7) long-wavelength limit: g_p = 1/(i nu - q mu) with nu = q z is b/q, and (b - b')^2 = (mu - mu')^2 b^2 b'^2",
  {1/(I q z - q mu1) - (1/(I z - mu1))/q, (1/(I z - mu1) - 1/(I z - mu2))^2 - (mu1 - mu2)^2/((I z - mu1)^2 (I z - mu2)^2)},
  {0, 0}, "Assumptions" -> q > 0 && z > 0 && Element[mu1 | mu2, Reals],
  "TestPoint" -> {q -> 1/7, z -> 37/100, mu1 -> -1/3, mu2 -> 2/5}];
check["(S7) bound |mu| |b|^2 <= z^-2 for |mu| <= 1, z > 0 (squared form)",
  Resolve[ForAll[{mu1, z}, -1 <= mu1 <= 1 && z > 0, mu1^2 z^4 <= (mu1^2 + z^2)^2], Reals], True];

check["(33),(35) constants used: alpha_2 = 1/Sqrt[2], alpha_3 = alphaL, B_2 = 2, B_3 = 4/Pi",
  {alphaD[2], alphaD[3], calB[2], calB[3]}, {1/Sqrt[2], alphaL, 2, 4/Pi}];
check["(33),(35) prefactor (d/(8 Pi)) E0 qTF^(d+1) with E0 = 2/(alpha_d rs)^2 Ry: 12/Pi^3 (d = 3), 2 Sqrt[2] rs/Pi (d = 2)",
  Table[(d/(8 Pi)) (2/(alphaD[d] rs)^2) (calB[d] alphaD[d] rs)^((d + 1)/(d - 1)), {d, {3, 2}}],
  {12/Pi^3, 2 Sqrt[2] rs/Pi}, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 37/100}];
check["(text Sec. 2.9) E0 qTF^4 = 32/Pi^2 Ry (d = 3), (text Sec. 3.1) E0 alpha_2 rs qTF^3 = 8 rs^2 Ry (d = 2)",
  {(2/(alphaD[3] rs)^2) (calB[3] alphaD[3] rs)^2, (2/(alphaD[2] rs)^2) alphaD[2] rs (calB[2] alphaD[2] rs)^3},
  {32/Pi^2, 8 rs^2}, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 37/100}];

dvRes34 = check["(35) residue (S70) with m = 0: kernel kappa q and screening gamma t^-2 give -kappa gamma eps^2 log eps",
  -Residue[kk gg eps^s/(s - 2)^2, {s, 2}], -kk gg eps^2 Log[eps], "Assumptions" -> eps > 0];
check["(35) assembly for general N >= 1: the residue computed above with kappa = ar qTF^(d-1) R, gamma = (-h)^N, eps = qTF and log qTF = log rs/(d-1)",
  Simplify[(dvRes34 /. {kk -> ar qT^(d - 1) RR, gg -> (-hh)^n, eps -> qT}) /. Log[qT] -> Lrs/(d - 1), Element[n, Integers] && n >= 1 && hh > 0 && qT > 0],
  (-1)^(n + 1)/(d - 1) ar qT^(d + 1) hh^n RR Lrs,
  "Assumptions" -> d > 1 && qT > 0 && hh > 0 && Element[n, Integers] && n >= 1, "TestPoint" -> {d -> 37/10, qT -> 1/5, ar -> 2/7, RR -> 3/5, hh -> 4/9, Lrs -> -3, n -> 3}];
check["(35) log qTF = log rs/(d-1) + O(1): d log(qTF)/d rs = 1/((d-1) rs)",
  D[Log[(bb ad rs)^(1/(d - 1))], rs], 1/((d - 1) rs),
  "Assumptions" -> d > 1 && bb > 0 && ad > 0 && rs > 0, "TestPoint" -> {d -> 37/10, bb -> 2, ad -> 3/7, rs -> 2/5}];
check["(35) resonant dimensions N* (d-1) = 2: d = 3, 2, 5/3 for N* = 1, 2, 3",
  Table[d /. First[Solve[n (d - 1) == 2, d]], {n, 1, 3}], {3, 2, 5/3}];
check["(35) E0 rs qTF^(d+1) ~ rs^(-2 + 1 + (d+1)/(d-1)) = rs^(N*) at d = 1 + 2/N*",
  Simplify[Exponent[PowerExpand[E0d[d, rs] rs qTFd[d, rs]^(d + 1)], rs] /. d -> 1 + 2/n, n > 0], n, "Assumptions" -> n > 0, "TestPoint" -> {n -> 3}];

check["(36) R^(1) as printed, with arctan(1/z), agrees with the definition R1", R1[z],
  4/(1 + z^2) ((1 + 3 z^2) - (2 + 3 z^2) z ArcTan[1/z]), "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(S8) azimuthal average of 1/(2(1 - cos Theta)) at fixed mu, mu' is 1/(2|mu - mu'|)",
  Integrate[1/(2 (1 - mu1 mu2 - Sqrt[1 - mu1^2] Sqrt[1 - mu2^2] Cos[ph])), {ph, 0, 2 Pi}, Assumptions -> -1 < mu1 < mu2 < 1]/(2 Pi),
  1/(2 (mu2 - mu1)), "Assumptions" -> -1 < mu1 < mu2 < 1, "TestPoint" -> {mu1 -> -1/3, mu2 -> 2/5}];
check["(S8) Re[1/((I z - mu)^2 (I z - mu')^2)] = a(mu) a(mu') - b(mu) b(mu') with a = (mu^2 - z^2)/(mu^2 + z^2)^2, b = 2 mu z/(mu^2 + z^2)^2",
  ComplexExpand[Re[1/((I z - mu1)^2 (I z - mu2)^2)]],
  (mu1^2 - z^2)/(mu1^2 + z^2)^2 (mu2^2 - z^2)/(mu2^2 + z^2)^2 - (2 mu1 z/(mu1^2 + z^2)^2) (2 mu2 z/(mu2^2 + z^2)^2),
  "Assumptions" -> Element[mu1 | mu2, Reals] && z > 0, "TestPoint" -> {mu1 -> -1/3, mu2 -> 2/5, z -> 7/10}];
dvR3 = check["(36),(S8) R3 = R^(1)/(4 Pi): reduced double integral over -1 < mu2 < mu1 < 1, inner integral in mu2 first (rational), then mu1, symbolic in z",
  Module[{inner},
    inner = Integrate[mu1 mu2 (mu1 - mu2) ((mu1^2 - z^2)/(mu1^2 + z^2)^2 (mu2^2 - z^2)/(mu2^2 + z^2)^2 - 4 mu1 mu2 z^2/((mu1^2 + z^2)^2 (mu2^2 + z^2)^2)),
      {mu2, -1, mu1}, Assumptions -> z > 0 && -1 < mu1 < 1];
    FullSimplify[-(2/(8 Pi)) Integrate[inner, {mu1, -1, 1}, Assumptions -> z > 0], z > 0]],
  R1[z]/(4 Pi), "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(36),(S8) the closed form of R3 computed above, after FunctionExpand and ComplexExpand, minus R^(1)/(4 Pi) simplifies to zero",
  FullSimplify[ComplexExpand[FunctionExpand[dvR3]] - R1[z]/(4 Pi), z > 0], 0, "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(text SI Sec. S3) printed quadrature values R3(z), z = 0.1, 0.3, 1, 3, lie within 2.5 10^-17 of R^(1)/(4 Pi)",
  And @@ MapThread[Abs[N[R1[#1]/(4 Pi), 40] - #2] < 25/10^18 &, {{1/10, 3/10, 1, 3},
    {23049451763139286/10^17, 11645170899574007/10^17, 11619772367581344/10^18, 2451548733833656/10^19}}], True];
check["(text Sec. 2.9) Pi^(1) -> -(qTF^4/16) R^(1): -qTF^2 alpha_L r_s R3 of (34) with r_s from (6) and R3 = R^(1)/(4 Pi)",
  -qT^2 (alphaD[3] (rs /. First[Solve[qT^2 == calB[3] alphaD[3] rs, rs]])) RR/(4 Pi), -(qT^4/16) RR, "Assumptions" -> qT > 0, "TestPoint" -> {qT -> 2/5, RR -> 3/7}];
dvIntHR1 = check["(37) Int_R h R^(1) dz = Pi (Pi^2 - 12 log 2)/3", 2 Integrate[h[z] R1[z], {z, 0, Infinity}], Pi (Pi^2 - 12 Log[2])/3];
check["(text SI Sec. S3) Int_R h R3 dz = (Pi^2 - 12 log 2)/12", 2 Integrate[h[z] R1[z]/(4 Pi), {z, 0, Infinity}], (Pi^2 - 12 Log[2])/12];
dvPref36 = check["(37) prefactor from (35) with d = 3, N* = 1: (d/(8 Pi)) E_0 q_TF^4 alpha_3 (1/(d-1)) (1/(4 Pi)) from the definitions = 3 alphaL/(2 Pi^4)",
  Simplify[(3/(8 Pi)) E0d[3, rs] qTFd[3, rs]^4 alphaD[3] (1/(3 - 1)) (1/(4 Pi)), rs > 0], 3 alphaL/(2 Pi^4), "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3/10}];
dvC36 = check["(37) closed form: prefactor times Int_R h R^(1) dz (both computed above) = alphaL (Pi^2 - 12 log 2)/(2 Pi^3)",
  dvPref36 dvIntHR1, alphaL (Pi^2 - 12 Log[2])/(2 Pi^3)];
numcheck["(37) Carr-Maradudin coefficient 0.01304", alphaL (Pi^2 - 12 Log[2])/(2 Pi^3), "0.01304"];
check["(37),(1),(text Sec. 2.9) C = the two parts computed above, (37) and the g3 term of Sec. S14.6, = alphaL (7 Pi^2/6 - 12 log 2 - 1)/(2 Pi^3)",
  dvC36 + dvC36b, alphaL (7 Pi^2/6 - 12 Log[2] - 1)/(2 Pi^3)];

check["(37) window integral: coefficient of eps^2 log eps in Int_0^1 q [-h eps^2/(q^2 + h eps^2)] dq is h = -kappa gamma (kappa = 1, gamma = -h)",
  Module[{ser}, ser = Normal[Series[Integrate[q (-hh eps^2/(q^2 + hh eps^2)), {q, 0, 1}, Assumptions -> eps > 0 && hh > 0], {eps, 0, 2}]];
    Coefficient[Expand[PowerExpand[ser]] /. Log[eps] -> LL, eps^2 LL]], hh, "TestPoint" -> {hh -> 3/5}];
check["(text SI Sec. S14.6) (12/(Pi^3 qTF^4)) times (S70) for kappa = qTF^4 R^(1)/16, gamma = -h gives (3/(4 Pi^3)) h R^(1) qTF^2 log qTF",
  (3/(8 Pi)) E0d[3, rs /. First[Solve[qT^2 == calB[3] alphaD[3] rs, rs]]] (-(qT^4 RR/16) (-hh) qT^2 LqT), (3/(4 Pi^3)) hh RR qT^2 LqT,
  "Assumptions" -> qT > 0, "TestPoint" -> {qT -> 2/5, RR -> 3/7, hh -> 3/5, LqT -> -2}];
check["(text SI Sec. S14.6) qTF^2 log qTF = (2 alphaL/Pi) rs log rs + O(rs) with qTF^2 = 4 alphaL rs/Pi",
  Coefficient[Expand[PowerExpand[qTFd[3, rs]^2 Log[qTFd[3, rs]]]] /. Log[rs] -> LL, rs LL], 2 alphaL/Pi];
numcheck["(text SI Sec. S19) residue integral 0.0130393528", alphaL (Pi^2 - 12 Log[2])/(2 Pi^3), "0.0130393528"];
numcheck["(Tab. 4) Pi^(1) in rings, d = 3: 0.0130", alphaL (Pi^2 - 12 Log[2])/(2 Pi^3), "0.0130"];

check["(S2) r-substitution: z = (1-r^2)/(2r) gives r = s_z - z, s_z = (1+r^2)/(2r), a = (1-r^2)/(1+r^2), h2 = 2r^2/(1+r^2), |dz/dr| = (1+r^2)/(2r^2)",
  With[{zr = stripCond[z /. First[Solve[r == Sqrt[1 + z^2] - z && r > 0, z]]]}, {rVar[zr], Sqrt[1 + zr^2], aVar[zr], h2D[zr], -D[zr, r]}],
  {r, (1 + r^2)/(2 r), (1 - r^2)/(1 + r^2), 2 r^2/(1 + r^2), (1 + r^2)/(2 r^2)},
  "Assumptions" -> 0 < r < 1, "TestPoint" -> {r -> 37/100}];
check["(10) derivative used below: h2Dp = d h2/dz", D[h2D[z], z], h2Dp[z], "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(38),(S9) the two forms of R2 agree, and h2 + z h2' = h2^2 - 2 a^2 r^2/(1 + r^2)",
  {R2D[z] - (2/Pi) (h2D[z] + z h2Dp[z] + 2 aVar[z]^2 rVar[z] ArcTan[rVar[z]]/(1 + rVar[z]^2)), h2D[z] + z h2Dp[z]},
  {0, h2D[z]^2 - 2 aVar[z]^2 rVar[z]^2/(1 + rVar[z]^2)}, "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(S9) h2 + z h2' = 1 - z/s_z - z/s_z^3",
  h2D[z] + z h2Dp[z], 1 - z/Sqrt[1 + z^2] - z/(1 + z^2)^(3/2), "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];

check["(text Sec. 3.1) kernel Fourier coefficients relative to K0: <(cos m phi - 1)/(2|sin(phi/2)|)>_phi = -(2/Pi) O_m, m = 1..4",
  Table[Integrate[(Cos[m ph] - 1)/(2 Sin[ph/2]), {ph, 0, Pi}]/Pi, {m, 1, 4}], Table[-(2/Pi) oddH[m], {m, 1, 4}]];
check["(text Sec. 3.1) kernel Fourier coefficients for general m >= 1: cos(m phi) - cos((m-1) phi) = -2 sin((2m-1) phi/2) sin(phi/2), so K_m - K_(m-1) = -(1/Pi) Int_0^Pi sin((2m-1) phi/2) dphi = -2/(Pi (2m-1)) = -(2/Pi)(O_m - O_(m-1)); with K_0 = 0 this gives K_m = -(2/Pi) O_m",
  {Simplify[ExpToTrig[Expand[TrigToExp[Cos[m ph] - Cos[(m - 1) ph] + 2 Sin[(2 m - 1) ph/2] Sin[ph/2]]]]],
   Simplify[-Integrate[Sin[(2 m - 1) ph/2], {ph, 0, Pi}]/Pi, Element[m, Integers] && m >= 1],
   FullSimplify[-(2/Pi) (oddH[m] - oddH[m - 1]), Element[m, Integers] && m >= 1]},
  {0, -2/(Pi (2 m - 1)), -2/(Pi (2 m - 1))}, "Assumptions" -> Element[m, Integers] && m >= 1, "TestPoint" -> {m -> 5, ph -> 7/10}];
dvRootsB = w /. Solve[I z - (w + 1/w)/2 == 0, w];   (* SI Sec. S4: poles of b in w = e^(i theta), computed *)
check["(text SI Sec. S4) Fourier coefficients of b for general m by residues: with w = e^(i theta), the poles of b = 1/(i z - (w + 1/w)/2) computed by Solve satisfy w1 + w2 = 2 I z and w1 w2 = 1, and b = -2 w/((w - w1)(w - w2))",
  {Simplify[Total[dvRootsB]], Simplify[Times @@ dvRootsB], Simplify[1/(I z - (w + 1/w)/2) + 2 w/((w - dvRootsB[[1]]) (w - dvRootsB[[2]]))]},
  {2 I z, 1, 0}, "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100, w -> 1/3 + I/5}];
dvW1 = First[Select[dvRootsB, TrueQ[Resolve[ForAll[z, z > 0, ComplexExpand[Abs[#]^2] < 1], Reals]] &]];   (* the pole inside the unit circle *)
dvW2 = First[Select[dvRootsB, TrueQ[Resolve[ForAll[z, z > 0, ComplexExpand[Abs[#]^2] > 1], Reals]] &]];   (* the pole outside *)
check["(text SI Sec. S4) the pole of b inside the unit circle, selected above, is w1 = -I r", Simplify[dvW1 + I rVar[z], z > 0], 0,
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(text SI Sec. S4) <b e^(i m theta)>_theta = (1/(2 Pi I)) Contour w^(m-1) b dw is the residue at the inner pole for m >= 0 (b even in theta, so e^(-i m theta) gives the same): (-i r)^m/(i s_z)",
  Simplify[Residue[w^(m - 1) (-2 w/((w - dvW1) (w - dvW2))), {w, dvW1}], z > 0 && Element[m, Integers] && m >= 0],
  (-I rVar[z])^m/(I Sqrt[1 + z^2]), "Assumptions" -> z > 0 && Element[m, Integers] && m >= 0, "TestPoint" -> {z -> 37/100, m -> 3}];
check["(text SI Sec. S4) Fourier coefficients of b: <b e^(-i m theta)>_theta = (-i r)^|m|/(i s_z), m = 0..3, symbolic in z",
  Table[Integrate[Exp[-I m th]/(I z - Cos[th]), {th, 0, 2 Pi}, Assumptions -> z > 0]/(2 Pi), {m, 0, 3}], Table[bCoef[m, z], {m, 0, 3}],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(text SI Sec. S4) F1_m = (b_(m-1) + b_(m+1))/2 = (-i r)^(m-1) (1 - r^2)/(2 i s_z) and F1_m^2 = (-1)^m r^(2m) a^2, m = 1..4",
  Table[{(bCoef[m - 1, z] + bCoef[m + 1, z])/2 - (-I rVar[z])^(m - 1) (1 - rVar[z]^2)/(2 I Sqrt[1 + z^2]),
     ((bCoef[m - 1, z] + bCoef[m + 1, z])/2)^2 - (-1)^m rVar[z]^(2 m) aVar[z]^2}, {m, 1, 4}],
  Table[{0, 0}, {m, 1, 4}], "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
dvF1sq = check["(text SI Sec. S4) F1_m^2 for general m >= 1 from the coefficients of b: (-1)^m r^(2m) a^2",
  Simplify[((bCoef[m - 1, z] + bCoef[m + 1, z])/2)^2, Element[m, Integers] && m >= 1 && z > 0], (-1)^m rVar[z]^(2 m) aVar[z]^2,
  "Assumptions" -> Element[m, Integers] && m >= 1 && z > 0, "TestPoint" -> {z -> 37/100, m -> 5}];
dvF1ratio = check["(text SI Sec. S4) F1_(m+1)^2/F1_m^2 from F1_m^2 computed above is -r^2, independent of m, so F1_m^2 = F1_1^2 (-r^2)^(m-1)",
  Simplify[(dvF1sq /. m -> m + 1)/dvF1sq, Element[m, Integers] && m >= 1 && z > 0], -rVar[z]^2, "Assumptions" -> z > 0 && Element[m, Integers] && m >= 1, "TestPoint" -> {z -> 37/100, m -> 5}];
dvF11 = check["(text SI Sec. S4) F1_1^2 from F1_m^2 computed above: -r^2 a^2", Simplify[dvF1sq /. m -> 1, z > 0], -rVar[z]^2 aVar[z]^2,
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(text SI Sec. S4) F1_m = (-i r)^(m-1) (1 - r^2)/(2 i s_z) and F1_m^2 = (-1)^m r^(2m) a^2 for general m >= 1",
  {Simplify[(bCoef[m - 1, z] + bCoef[m + 1, z])/2 - (-I rVar[z])^(m - 1) (1 - rVar[z]^2)/(2 I Sqrt[1 + z^2]), Element[m, Integers] && m >= 1 && z > 0],
   Simplify[((bCoef[m - 1, z] + bCoef[m + 1, z])/2)^2 - (-1)^m rVar[z]^(2 m) aVar[z]^2, Element[m, Integers] && m >= 1 && z > 0]},
  {0, 0}, "Assumptions" -> Element[m, Integers] && m >= 1 && z > 0, "TestPoint" -> {z -> 37/100, m -> 5}];
dvReF21 = check["(text SI Sec. S4) Re F2_1 = Re <cos(theta) b^2 e^(-i theta)>_theta = h2 + z h2', symbolic in z",
  Simplify[ComplexExpand[Re[Integrate[Cos[th] Exp[-I th]/(I z - Cos[th])^2, {th, 0, 2 Pi}, Assumptions -> z > 0]/(2 Pi)]], z > 0],
  h2D[z] + z h2Dp[z], "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(text SI Sec. S4) generating function Sum_(m>=1) O_m x^m = Sqrt[x] ArcTanh[Sqrt[x]]/(1 - x): coefficients of x^1..x^8",
  Table[SeriesCoefficient[Sqrt[x] ArcTanh[Sqrt[x]]/(1 - x), {x, 0, k}], {k, 1, 8}], Table[oddH[k], {k, 1, 8}]];
dvGF = First[check["(text SI Sec. S4) generating function derived: Sum_(j>=1) x^j/(2j - 1) = Sqrt[x] ArcTanh[Sqrt[x]], and O_m - O_(m-1) = 1/(2m - 1), so Sum_(m>=1) O_m x^m = (1/(1-x)) Sum_(j>=1) x^j/(2j-1) (partial sums)",
  {Sum[x^j/(2 j - 1), {j, 1, Infinity}], FullSimplify[oddH[m] - oddH[m - 1], Element[m, Integers] && m >= 1]},
  {Sqrt[x] ArcTanh[Sqrt[x]], 1/(2 m - 1)}, "Assumptions" -> 0 < x < 1 && Element[m, Integers] && m >= 1, "TestPoint" -> {x -> 1/3, m -> 4}]]/(1 - x);
dvGFr = check["(text SI Sec. S4) at x = -r^2 the generating function computed above is -r ArcTan[r]/(1 + r^2)",
  Simplify[dvGF /. x -> -r^2, 0 < r < 1], -r ArcTan[r]/(1 + r^2), "Assumptions" -> 0 < r < 1, "TestPoint" -> {r -> 37/100}];
check["(38),(S9) Fourier assembly R2 = Pi^-1 Sum O_|m| T_m: self-energy terms T_(+-1) = Re F2_1 (computed above), vertex terms T_m = -2 F1_m^2 = -2 (F1_1^2/(-r^2)) (-r^2)^m (both computed above) summed over m >= 1 with the generating function at x = -r^2 (computed above), both signs of m",
  (1/Pi) (2 oddH[1] dvReF21 + 2 (-2) (dvF11/dvF1ratio) (dvGFr /. r -> rVar[z])), R2D[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];

check["(text SI Sec. S4) printed quadrature value R2(1) = 0.04977088289807023 lies within 7 10^-18 of (S9)",
  Abs[N[R2D[1], 40] - 4977088289807023/10^17] < 7/10^18, True];
numcheck["(Tab. S3) caption R2(0.3) = 0.3180676", R2D[3/10], "0.3180676"];
numcheck["(Tab. S3) caption R2(1) = 0.0497709", R2D[1], "0.0497709"];

check["(38)(i) static limit R2(0) = 2/Pi as z -> 0+", Limit[R2D[z], z -> 0, Direction -> "FromAbove"], 2/Pi];
check["(38)(i) Pi0 + Pi1 -> -qTF (h2(0) + (rs/Sqrt[2]) R2(0)) = -qTF (1 + Sqrt[2] rs/Pi) at nu = 0",
  -qT (h2D[0] + (rs/Sqrt[2]) R2D[0]), -qT (1 + Sqrt[2] rs/Pi)];
check["(38)(i) independent: Hartree-Fock compressibility of the 2D gas (n = 1/(Pi rs^2)) gives kappa0/kappa = 1 - Sqrt[2] rs/Pi, and 1/(1 - Sqrt[2] rs/Pi) = 1 + Sqrt[2] rs/Pi + O(rs^2)",
  {(D[nn (Pi nn - 8 Sqrt[2] Sqrt[Pi nn]/(3 Pi)), {nn, 2}]/D[nn Pi nn, {nn, 2}]) /. nn -> 1/(Pi rs^2),
   Normal[Series[1/(1 - Sqrt[2] rs/Pi), {rs, 0, 1}]]},
  {1 - Sqrt[2] rs/Pi, 1 + Sqrt[2] rs/Pi}, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 37/100}];
dvR2ser = check["(38)(ii),(S9) large-z expansion: coefficients of z^0..z^-8 are {0, 0, 0, 0, 5/(12 Pi), 0, -11/(20 Pi), 0, 1283/(2240 Pi)} (no z^-2 term)",
  Module[{ser}, ser = Normal[Series[Simplify[R2D[1/w], w > 0], {w, 0, 9}]]; Table[Coefficient[ser, w, k], {k, 0, 8}]],
  {0, 0, 0, 0, 5/(12 Pi), 0, -11/(20 Pi), 0, 1283/(2240 Pi)}];
check["(S9) remainder after the three terms is O(z^-10): the coefficient of z^-9 vanishes and that of z^-10 is a rational multiple of 1/Pi",
  Module[{ser}, ser = Normal[Series[Simplify[R2D[1/w], w > 0], {w, 0, 11}]]; {Coefficient[ser, w, 9], Element[Pi Coefficient[ser, w, 10], Rationals]}], {0, True}];
check["(38)(ii) -v Pi1 = rs^2 R2/q (Pi1 = -qTF (rs/Sqrt[2]) R2, qTF = Sqrt[2] rs, v = 1/q); z -> -I omega/q gives 5 rs^2 q^3/(12 Pi omega^4) + 11 rs^2 q^5/(20 Pi omega^6)",
  Expand[(1/q) qTFd[2, rs] alphaD[2] rs (dvR2ser[[5]]/z^4 + dvR2ser[[7]]/z^6) /. z -> -I om/q],
  5 rs^2 q^3/(12 Pi om^4) + 11 rs^2 q^5/(20 Pi om^6)];

check["(text SI Sec. S6) Eq. (25) of the reference divided by P(0) = -2/Pi^2 gives Delta P/P(0) = x^2 (log x + beta)/96",
  (-x^2 (Log[x] + bet)/(48 Pi^2))/(-2/Pi^2), x^2 (Log[x] + bet)/96];
check["(text SI Sec. S6) Delta P of Table I of the reference (-9.168 10^-5, -4.263 10^-4, -1.0418 10^-3, -1.9655 10^-3) divided by P(0) gives the quoted 4.524 10^-4, 2.104 10^-3, 5.141 10^-3, 9.699 10^-3",
  And @@ MapThread[Abs[#1/(-2/Pi^2) - #2] <= #3 &, {{-9168/10^8, -4263/10^7, -10418/10^7, -19655/10^7},
     {4524/10^7, 2104/10^6, 5141/10^6, 9699/10^6}, {5/10^8, 5/10^7, 5/10^7, 5/10^7}}], True];
check["(text SI Sec. S6) the static values 4.523 10^-4, 2.104 10^-3, 5.141 10^-3, 9.699 10^-3 agree with the reference (4.524 10^-4, 2.104 10^-3, 5.141 10^-3, 9.699 10^-3) to within one unit of the fourth digit",
  And @@ MapThread[Abs[#1 - #2] <= #3 &, {{4523/10^7, 2104/10^6, 5141/10^6, 9699/10^6}, {4524/10^7, 2104/10^6, 5141/10^6, 9699/10^6}, {1/10^7, 1/10^6, 1/10^6, 1/10^6}}], True];
check["(text SI Sec. S6) the printed z = 0.3 entries of Tab. S3 decrease by factors 3.62 and 3.66 when q is halved (about 3.6; an O(q) term would give 2)",
  3.5 < (12186/10^8)/(33626/10^9) < 3.7 && 3.5 < (33626/10^9)/(91990/10^10) < 3.7, True];

check["(text Sec. 3.1 (iii)) Ward coefficient qTF^(d-1) alpha_d rs Pi^-1 [h_d + z h_d'] for d = 2 is qTF (rs/Sqrt[2]) Pi^-1 [h2 + z h2'], static value rs^2/Pi = qTF rs/(Sqrt[2] Pi)",
  {calB[2] alphaD[2] rs alphaD[2] rs (h2D[z] + z h2Dp[z])/Pi, (calB[2] alphaD[2] rs alphaD[2] rs (h2D[z] + z h2Dp[z])/Pi) /. z -> 0},
  {(Sqrt[2] rs) (rs/Sqrt[2]) (h2D[z] + z h2Dp[z])/Pi, (Sqrt[2] rs) rs/(Sqrt[2] Pi)},
  "Assumptions" -> rs > 0 && z > 0, "TestPoint" -> {rs -> 2/5, z -> 37/100}];
numcheck["(text SI Sec. S4) Pi^-1 (h2 + z h2') at z = 10^-4: 0.318246", (h2D[1/10000] + h2Dp[1/10000]/10000)/Pi, "0.318246"];
numcheck["(text SI Sec. S4) Pi^-1 (h2 + z h2') at z = 0.3: 0.1429308", (h2D[3/10] + 3/10 h2Dp[3/10])/Pi, "0.1429308"];
numcheck["(text SI Sec. S4) Pi^-1 (h2 + z h2') at z = 1: -0.01930873", (h2D[1] + h2Dp[1])/Pi, "-0.01930873"];
check["(text SI Sec. S4) at phi = 0, (1/(2 Pi)) Int dtheta mu mu' Re(b b') = (1/(2 Pi)) Int dtheta cos^2 theta Re b^2 = h2 + z h2', symbolic in z",
  Simplify[ComplexExpand[Re[Integrate[Cos[th]^2/(I z - Cos[th])^2, {th, 0, 2 Pi}, Assumptions -> z > 0]/(2 Pi)]], z > 0], h2D[z] + z h2Dp[z],
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 3/10}];

check["(text Sec. 3.1) scaling function 1/(1+c) - 1 for d = 2 has the tail -h2/t + h2^2/t^2 - ...",
  Normal[Series[1/(1 + hh/t) - 1, {t, Infinity, 2}]], -hh/t + hh^2/t^2];
check["(39) from (35) with d = 2, N* = 2: (1/(4 Pi)) E0 alpha_2 rs qTF^3 (-1)^3/(d-1) = -(2/Pi) rs^2 (E0 alpha_2 rs qTF^3 = 8 rs^2 Ry)",
  Simplify[(2/(8 Pi)) E0d[2, rs] alphaD[2] rs qTFd[2, rs]^3 (-1)^3/(2 - 1), rs > 0], -(2/Pi) rs^2, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 2/5}];
check["(text SI Sec. S5) prefactor (1/(4 Pi)) E0 = 1/(Pi rs^2) Ry for d = 2",
  (2/(8 Pi)) (2/(alphaD[2] rs)^2), 1/(Pi rs^2), "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 2/5}];

dvRes10 = check["(S10) triple pole: residue rule (S70) with m = 1 gives +kappa gamma/2 eps^2 log^2 eps",
  -Residue[(-1)^1 1! kk gg eps^s/(s - 2)^3, {s, 2}], kk gg eps^2 Log[eps]^2/2, "Assumptions" -> eps > 0];
dvS10 = check["(S10) assembly for general N >= 1: the triple-pole residue computed above with kappa = ar qTF^(d-1) Pi^-1 (h + z h'), gamma = (-h)^N, eps = qTF and log qTF = log rs/(d-1) gives the coefficient (-1)^N/(2 Pi (d-1)^2)",
  Simplify[(dvRes10 /. {kk -> ar qT^(d - 1) HH/Pi, gg -> (-hh)^n, eps -> qT}) /. Log[qT] -> Lrs/(d - 1), Element[n, Integers] && n >= 1 && hh > 0 && qT > 0],
  (-1)^n/(2 Pi (d - 1)^2) ar qT^(d + 1) HH hh^n Lrs^2,
  "Assumptions" -> d > 1 && qT > 0 && hh > 0 && Element[n, Integers] && n >= 1, "TestPoint" -> {d -> 37/10, qT -> 1/5, ar -> 2/7, HH -> 3/5, hh -> 4/9, Lrs -> -3, n -> 2}];
check["(S10),(text Sec. 3.4) integration by parts: z h^N h' = z d(h^(N+1))/dz/(N+1)", z f[z]^n f'[z], z D[f[z]^(n + 1), z]/(n + 1)];
check["(S10) boundary terms z h^(N+1) vanish as z -> Infinity: d = 2 (N* = 2), d = 3 (N* = 1)",
  {Limit[z h2D[z]^3, z -> Infinity], Limit[z h[z]^2, z -> Infinity]}, {0, 0}];
dvVtx2 = check["(58),(S10),(S3) d = 2: Int_R h2^2 (h2 + z h2') dz = (2/3) Int_R h2^3 dz, with Int_R h2^3 dz = 10 - 3 Pi",
  {2 Integrate[h2D[z]^2 (h2D[z] + z h2Dp[z]), {z, 0, Infinity}], 2 Integrate[h2D[z]^3, {z, 0, Infinity}]}, {(2/3) (10 - 3 Pi), 10 - 3 Pi}];
dvVtx3 = check["(text Sec. 3.4),(S10) d = 3: Int_R h (h + z h') dz = (1/2) Int_R h^2 dz = (Pi/3)(1 - log 2)",
  2 Integrate[h[z] (h[z] + z hp[z]), {z, 0, Infinity}], (Pi/3) (1 - Log[2])];
check["(58),(S10) d = 2: (d/(8 Pi)) E0 alpha_2 rs qTF^3 = 2 rs^2/Pi, and (2 rs^2/Pi) ((-1)^2/(2 Pi)) (2/3)(10 - 3 Pi) = (2/(3 Pi^2))(10 - 3 Pi) rs^2",
  {Simplify[(2/(8 Pi)) E0d[2, rs] alphaD[2] rs qTFd[2, rs]^3, rs > 0],
   Simplify[(2/(8 Pi)) E0d[2, rs] alphaD[2] rs qTFd[2, rs]^3, rs > 0] Simplify[(dvS10/(ar qT^(d + 1) HH hh^n Lrs^2)) /. {n -> 2, d -> 2}] First[dvVtx2]},
  {2 rs^2/Pi, (2/(3 Pi^2)) (10 - 3 Pi) rs^2}, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 2/5}];
numcheck["(58) 0.03885", (2/(3 Pi^2)) (10 - 3 Pi), "0.03885"];
numcheck["(Tab. 4) Pi^(1)_vtx in rings, d = 2: 0.0389", (2/(3 Pi^2)) (10 - 3 Pi), "0.0389"];
check["(text Sec. 3.4),(S10) d = 3: (3/(8 Pi)) E0 alphaL rs qTF^4 = 12 alphaL rs/Pi^3, and (12 alphaL rs/Pi^3) ((-1)/(8 Pi)) (Pi/3)(1 - log 2) = 3 c_hd rs with c_hd of (50)",
  {Simplify[(3/(8 Pi)) E0d[3, rs] alphaD[3] rs qTFd[3, rs]^4, rs > 0],
   Simplify[(3/(8 Pi)) E0d[3, rs] alphaD[3] rs qTFd[3, rs]^4, rs > 0] Simplify[(dvS10/(ar qT^(d + 1) HH hh^n Lrs^2)) /. {n -> 1, d -> 3}] dvVtx3},
  {12 alphaL rs/Pi^3, 3 (-alphaL (1 - Log[2])/(6 Pi^3)) rs}, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 2/5}];
check["(text Sec. 3.4) N* + 1 = (d+1)/(d-1) = N_r for N* = 2/(d-1)", 2/(d - 1) + 1, (d + 1)/(d - 1), "Assumptions" -> d > 1, "TestPoint" -> {d -> 37/10}];

check["(S11) r-substitution of (39): -(1/Pi) 2 h2^2 R2 |dz/dr| = (-32 r^6/(1+r^2)^3 + 16 r^4 (1-r^2)^2/(1+r^2)^4 - 16 r^3 (1-r^2)^2 ArcTan[r]/(1+r^2)^4)/Pi^2",
  With[{zr = stripCond[z /. First[Solve[r == Sqrt[1 + z^2] - z && r > 0, z]]]}, -(1/Pi) 2 h2D[zr]^2 R2D[zr] (1 + r^2)/(2 r^2)],
  (-32 r^6/(1 + r^2)^3 + 16 r^4 (1 - r^2)^2/(1 + r^2)^4 - 16 r^3 (1 - r^2)^2 ArcTan[r]/(1 + r^2)^4)/Pi^2,
  "Assumptions" -> 0 < r < 1, "TestPoint" -> {r -> 37/100}];
dvA = check["(S11) A = 3/2 - 15 Pi/32", Integrate[r^6/(1 + r^2)^3, {r, 0, 1}], 3/2 - 15 Pi/32];
dvB = check["(S11) B = 13/6 - 11 Pi/16", Integrate[r^4 (1 - r^2)^2/(1 + r^2)^4, {r, 0, 1}], 13/6 - 11 Pi/16];
dvUint = check["(S11) u = 1 + r^2 (r = Sqrt[u - 1]) maps r^3 (1 - r^2)^2/(1 + r^2)^4 dr to (u - 1)(2 - u)^2/(2 u^4) du",
  Simplify[(r^3 (1 - r^2)^2/(1 + r^2)^4 /. r -> Sqrt[u - 1]) D[Sqrt[u - 1], u], u > 1], (u - 1) (2 - u)^2/(2 u^4),
  "Assumptions" -> u > 1, "TestPoint" -> {u -> 3/2}];
dvAntiC = check["(S11) antiderivative of the rational factor in u = 1 + r^2, computed: log(u)/2 + 5/(2u) - 2/u^2 + 2/(3u^3)",
  Integrate[Expand[dvUint], u], Log[u]/2 + 5/(2 u) - 2/u^2 + 2/(3 u^3), "Assumptions" -> u > 1, "TestPoint" -> {u -> 3/2}];
check["(S11) one integration by parts with the antiderivative F computed above: C = [F ArcTan[r]]_0^1 - Int_0^1 F/(1 + r^2) dr, which contains Int_0^1 log(1 + r^2)/(1 + r^2) dr",
  Module[{F = dvAntiC /. u -> 1 + r^2}, (F ArcTan[r] /. r -> 1) - (F ArcTan[r] /. r -> 0) - Integrate[F/(1 + r^2), {r, 0, 1}]],
  -5/18 - Pi Log[2]/8 + Pi/32 + Catalan/2];
check["(S11) Int_0^1 log(1 + r^2)/(1 + r^2) dr = (Pi/2) log 2 - G", Integrate[Log[1 + r^2]/(1 + r^2), {r, 0, 1}], Pi Log[2]/2 - Catalan];
dvC = check["(S11) C = -5/18 - Pi log 2/8 + Pi/32 + G/2, directly", Integrate[r^3 (1 - r^2)^2 ArcTan[r]/(1 + r^2)^4, {r, 0, 1}], -5/18 - Pi Log[2]/8 + Pi/32 + Catalan/2];
check["(S11)->(40) 2 (-32 A + 16 B - 16 C)/Pi^2 with A, B, C computed above is a_2l of (40)",
  2 (-32 dvA + 16 dvB - 16 dvC)/Pi^2, aTwoL];

check["(SI Prop. S1),(text Sec. 3.1),(text SI Sec. S7) one soft line: E0 rs^j qTF^(d+1+m) ~ rs^(j + (3-d+m)/(d-1)), = j + m/2 (d = 3), j + 1 + m (d = 2)",
  {Simplify[-2 + j + (d + 1 + m)/(d - 1) - (j + (3 - d + m)/(d - 1))], (j + (3 - d + m)/(d - 1)) /. d -> 3, (j + (3 - d + m)/(d - 1)) /. d -> 2},
  {0, j + m/2, j + 1 + m}];
check["(text Sec. 3.1),(text SI Sec. S7) orders: d = 3: Pi1 (j=1, m=0), g3 (0, 2) -> 1; d = 2: Pi1 (1, 0), linear Lindhard term (0, 1) -> 2; g2 (0, 2), second-order polarization and two Pi1 (2, 0) -> 3; q^2 log q of Pi1 (1, 2) -> 4; rings (0, 0): 0 (d = 3), 1 (d = 2)",
  Function[{dd, jj, mm}, jj + (3 - dd + mm)/(dd - 1)] @@@ {{3, 1, 0}, {3, 0, 2}, {2, 1, 0}, {2, 0, 1}, {2, 0, 2}, {2, 2, 0}, {2, 1, 2}, {3, 0, 0}, {2, 0, 0}},
  {1, 1, 2, 2, 3, 3, 4, 0, 1}];
check["(text SI Sec. S7) several soft lines: kF exponent of kF^-d [kF Q^(d+1)]^L [kF^(2-d)]^I Prod_l [kF^(d+2-2l)]^(V_l), Q = tau kF, Sum_l l V_l = 2 I, L = I - V + 1, is 2",
  Expand[-d + L + L (d + 1) + (2 - d) II + (d + 2) V - 2 (2 II) /. L -> II - V + 1], 2];
check["(text SI Sec. S7) Gamma_l = O(g_F kF^(-2(l-2))) with g_F ~ kF^(d-2) is kF^(d+2-2l) (exponents)", (d - 2) - 2 (l - 2), d + 2 - 2 l];
check["(SI Prop. S1),(text SI Sec. S7) L soft loops: kF^2 tau^(L(d+1)) ~ rs^(L(d+1)/(d-1) - 2) = rs^(3L-2) (d = 2), rs^(2L-2) (d = 3); L = 1, 2 in d = 2 give rs^1, rs^4",
  {(L (d + 1)/(d - 1) - 2) /. d -> 2, (L (d + 1)/(d - 1) - 2) /. d -> 3, Table[(LL (d + 1)/(d - 1) - 2) /. d -> 2, {LL, 1, 2}]},
  {3 L - 2, 2 L - 2, {1, 4}}];
check["(text Sec. 3.1),(text SI Sec. S7) unsymmetrized loops: two three-point vertices, each larger by tau^-1 ~ rs^(-1/(d-1)), move L = 2 to rs^2 (d = 2) and rs^1 (d = 3)",
  {(2 (d + 1)/(d - 1) - 2 - 2/(d - 1)) /. d -> 2, (2 (d + 1)/(d - 1) - 2 - 2/(d - 1)) /. d -> 3}, {2, 1}];
check["(text Sec. 2.9) g3 term: q^3 qTF^2 g3 h^2/t^4 = q^5 g3 h^2 (qTF/q)^6 with t = q/qTF (the case j = 0, m = 2 of the one-line counting)",
  q^3 qT^2 gs hh^2/(q/qT)^4, q^5 gs hh^2 (qT/q)^6, "Assumptions" -> q > 0 && qT > 0, "TestPoint" -> {q -> 1/3, qT -> 1/7, gs -> -1/5, hh -> 2/3}];

check["(S12),(Tab. S4) zeta = 1: a_l(1)/a_l(0) = 1/(4 Sqrt[2]): the t_- term h2(z/t_-) -> 0, the remaining Int_0^Infinity h2(z/Sqrt[2])^3 dz is computed, and the zeta = 0 denominator Int_0^Infinity (2 h2)^3 dz = 4 Int_R h2^3 dz (computed in (S3))",
  {Limit[h2D[z/tm], tm -> 0, Direction -> "FromAbove", Assumptions -> z > 0], Integrate[h2D[z/Sqrt[2]]^3, {z, 0, Infinity}]/(4 dvIntH2cube)}, {0, 1/(4 Sqrt[2])},
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 37/100}];
check["(S12),(Tab. S4),(text Sec. 3.1) zeta = 1: a_2l(1)/a_2l(0) = 1/8: the terms with t_- = Sqrt[1 - zeta] -> 0 vanish (R2(z/t)/t -> 0, h2(z/t) -> 0), z = Sqrt[2] u maps R2(z/Sqrt[2]) h2(z/Sqrt[2])^2 dz/Sqrt[2] to R2(u) h2(u)^2 du, and at zeta = 0 the integrand (2 R2)(2 h2)^2 is 8 R2 h2^2",
  {Limit[R2D[z/tm]/tm, tm -> 0, Direction -> "FromAbove", Assumptions -> z > 0], Limit[h2D[z/tm], tm -> 0, Direction -> "FromAbove", Assumptions -> z > 0],
   Simplify[((R2D[z/Sqrt[2]] h2D[z/Sqrt[2]]^2/Sqrt[2]) /. z -> Sqrt[2] u) D[Sqrt[2] u, u], u > 0], Simplify[(R2D[u] h2D[u]^2)/((2 R2D[u]) (2 h2D[u])^2), u > 0]},
  {0, 0, R2D[u] h2D[u]^2, 1/8}, "Assumptions" -> z > 0 && u > 0, "TestPoint" -> {z -> 37/100, u -> 37/100}];

numcheck["(40) a_2l = -0.1754613739208", aTwoL, "-0.1754613739208"];
numcheck["(40) a_2l/2 in Hartree = -0.0877306869604", aTwoL/2, "-0.0877306869604"];
check["(59) the closed form printed in (59) is that of (40)", 7/Pi + 4 Log[2]/Pi - 160/(9 Pi^2) - 16 Catalan/Pi^2, aTwoL];
numcheck["(59),(text Sec. 4) a_2l = -0.17546", aTwoL, "-0.17546"];
numcheck["(abstract) 7/Pi + 4 log 2/Pi - 160/(9 Pi^2) - 16 G/Pi^2 = -0.17546", 7/Pi + 4 Log[2]/Pi - 160/(9 Pi^2) - 16 Catalan/Pi^2, "-0.17546"];
check["(text Sec. 3.1) a_2l rs^2 log rs at rs = 0.1 is 4.0 10^-3 Ry to the printed digits", Abs[N[aTwoL (1/100) Log[1/10]] - 4/1000] < 5/10^5, True];
check["(text Sec. 3.1) at rs = 0.1 the term is about a tenth of a_l rs log rs (a_l of (2)): ratio between 0.09 and 0.11",
  0.09 < N[aTwoL (1/100) Log[1/10]/(-2 Sqrt[2] (10/(3 Pi) - 1) (1/10) Log[1/10])] < 0.11, True];

numcheck["(43) e_2x = log2/3 - 3 zeta(3)/(2 pi^2), value printed in SI Sec. S19", Log[2]/3 - 3 Zeta[3]/(2 Pi^2), "0.0483583"];

dvPiVtxPref = check["(S20) Pi_vtx = -2 C^2/(2 pi)^6 X: 2 C^2/(2 pi)^6 = qTF^4/(32 pi^2) with C = C_3 alpha_L r_s and r_s from (6)",
  2 (calC[3] alphaD[3] (rs /. First[Solve[qTF^2 == calB[3] alphaD[3] rs, rs]]))^2/(2 Pi)^6, qTF^4/(32 Pi^2), "Assumptions" -> qTF > 0, "TestPoint" -> {qTF -> 2}];
check["(41),(S18) A0: Y_D = -12 Pi_D/(S_D pi^3 qTF^4) with S_D = 4 and Pi_vtx = -(prefactor computed in (S20)) X",
  (-12/(4 Pi^3 qTF^4)) (-dvPiVtxPref), A0, "Assumptions" -> qTF > 0, "TestPoint" -> {qTF -> 2}];
check["(52) -12/(S_D pi^3 qTF^4) with S_D = 4 is -3/(pi^3 qTF^4)", -12/(4 Pi^3 qTF^4), -3/(Pi^3 qTF^4), "TestPoint" -> {qTF -> 2}];
check["(S43) BRM prefactor 2 A0 alphaL rs/pi^2 = 3 alphaL rs/(16 pi^7)", 2 A0 alphaL rs/Pi^2, 3 alphaL rs/(16 Pi^7), "TestPoint" -> {rs -> 2}];

check["(45) e_kite = (2/rs^2) Int_0^rs rho F(rho): F = rs^-1/4 gives (8/7) rs^-1/4",
  2/rs^2 Integrate[rho rho^(-1/4), {rho, 0, rs}, Assumptions -> rs > 0], 8/7 rs^(-1/4), "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3}];
check["(45) F = e + (rs/2) e' for e = rs^p gives (1 + p/2) rs^p, so F/e = 7/8 at p = -1/4",
  Simplify[(rs^p + rs/2 D[rs^p, rs])/rs^p], 1 + p/2, "Assumptions" -> rs > 0, "TestPoint" -> {rs -> 3, p -> -1/4}];
check["(45),(S68) F = c rs log^2 rs gives e_kite = (2/3) c rs log^2 rs + lower powers of log rs (factor 3/2)",
  Coefficient[Expand[PowerExpand[2/rs^2 Integrate[rho c rho Log[rho]^2, {rho, 0, rs}, Assumptions -> rs > 0]]], rs Log[rs]^2],
  2 c/3, "TestPoint" -> {c -> 1}];

check["(S34) partial fractions by Apart: c^2(c^2-z^2)/(c^2+z^2)^2 = 1 - 3z^2/(c^2+z^2) + 2z^4/(c^2+z^2)^2",
  Apart[c^2 (c^2 - z^2)/(c^2 + z^2)^2, c], 1 - 3 z^2/(c^2 + z^2) + 2 z^4/(c^2 + z^2)^2, "TestPoint" -> {c -> 1/3, z -> 7/10}];
dvAx = check["(S30),(47) a_x = 2 pi Int dOmega_p c_p^2 (c_p^2 - z^2)/(c_p^2 + z^2)^2 with dOmega_p = 2 pi dc_p, -1 < c_p < 1: 8 pi^2 (h + z h')",
  2 Pi 2 Pi Integrate[c^2 (c^2 - z^2)/(c^2 + z^2)^2, {c, -1, 1}, Assumptions -> z > 0], 8 Pi^2 (h[z] + z hp[z]),
  "Assumptions" -> z > 0, "TestPoint" -> {z -> 7/10}];
check["(47) Int_0^Z a_x dz = 8 pi^2 Z h(Z) -> 0: the bare kernel has no q log(1/q)",
  Module[{F = Integrate[dvAx, {z, 0, Z}, Assumptions -> Z > 0]}, {F, Limit[F, Z -> Infinity]}],
  {8 Pi^2 Z h[Z], 0}, "Assumptions" -> Z > 0, "TestPoint" -> {Z -> 5}];

dvKappa = check["(S26) kappa(Lambda) = 4 pi A0 2 pi Int dOmega_p c_p^2 Int_0^Infinity (dz/pi) Lambda^2 (c_p^2 - z^2)/((Lambda^2 + z^2)(c_p^2 + z^2)^2), with dOmega_p = 2 pi dc_p and the integrand even in c_p: z integral by (S81), then c integral",
  4 Pi A0 2 Pi 2 Pi 2 Integrate[c^2 Integrate[L^2 (c^2 - z^2)/((L^2 + z^2) (c^2 + z^2)^2), {z, 0, Infinity}, Assumptions -> L > 0 && c > 0]/Pi, {c, 0, 1}, Assumptions -> L > 0],
  3 L/(2 Pi^2) (1 - 2 L Log[1 + 1/L] + L/(1 + L)), "Assumptions" -> L > 0, "TestPoint" -> {L -> 3/10}];
check["(S26) limits of kappa computed above: kappa ~ 3 Lambda/(2pi^2) at 0 and ~ 1/(2 pi^2 Lambda) at infinity, so kappa(infinity) = 0 (text after (S27))",
  {Limit[dvKappa/L, L -> 0, Direction -> "FromAbove"], Limit[L dvKappa, L -> Infinity], Limit[dvKappa, L -> Infinity]}, {3/(2 Pi^2), 1/(2 Pi^2), 0}];

dvCalA = check["(S31) calA = Int_0^inf h a_x dz/pi with a_x computed in (S30) = 4 pi^2 (1-log2)/3", Integrate[h[z] dvAx, {z, 0, Infinity}]/Pi,
  4 Pi^2 (1 - Log[2])/3];
dvPrefWin = check["(S29) e_kite - e_2x = K[v (rho - 1)]: A0 times d^3q = 4 pi q^2 dq times v = 1/q^2 times the coefficient of c in rho(c) computed from (19), with c = q_TF^2 L/q^2, gives -(8 pi A0/3) q_TF^2 L/q^2 per dq",
  A0 4 Pi q^2 (1/q^2) SeriesCoefficient[dvRho, {c, 0, 1}] (qTF^2 lL/q^2), -(8 Pi A0/3) qTF^2 lL/q^2, "TestPoint" -> {q -> 1/3, qTF -> 1/5, lL -> 2/3}];
check["(S29) prefactor 8 pi A0/3 = 1/(4 pi^4)", 8 Pi A0/3, 1/(4 Pi^4)];
dvWin = check["(S32) window Int_{sqrt tau}^1 log(1/q)/q dq = log^2(tau)/8",
  Integrate[Log[1/q]/q, {q, Sqrt[tau], 1}, Assumptions -> 0 < tau < 1], Log[tau]^2/8,
  "Assumptions" -> 0 < tau < 1, "TestPoint" -> {tau -> 1/100}];
check["(50) c_hd = coefficient of r_s log^2 r_s in (prefactor of (S29)) calA tau (window integral of (S32)) with tau = q_TF^2 = B_3 alpha_L r_s, all computed above: -alphaL(1-log2)/(6 pi^3)",
  Coefficient[Expand[PowerExpand[(dvPrefWin q^2/(qTF^2 lL)) dvCalA tau dvWin /. tau -> calB[3] alphaD[3] rs]] /. Log[rs] -> LL, rs LL^2],
  -alphaL (1 - Log[2])/(6 Pi^3)];
numcheck["(50) c_hd", -alphaL (1 - Log[2])/(6 Pi^3), "-0.0008594457"];
check["(S78) second line: 2 Pi Int dOmega M_h(c,c) = calA, with M_Lambda(c,c) = c^2 Lambda^2 J(Lambda;c,c) [(S77) with J of (S56), computed], M_h = Int_0^1 M_Lambda dLambda, dOmega = 2 Pi dc on -1 < c < 1 with M even in c, and calA of (S31) (computed)",
  2 Pi 2 Pi 2 Integrate[Integrate[c^2 L^2 (dvJPlus /. {al -> c, be -> c, w -> L}), {L, 0, 1}, Assumptions -> c > 0], {c, 0, 1}], dvCalA];

check["(S83),(48) p_loc: Int w^2(1/4-u^2)/((u^2+w^2)(u^2+1/4)^2) du = 2 pi w/(1+2w)^2",
  Integrate[w^2 (1/4 - u^2)/((u^2 + w^2) (u^2 + 1/4)^2), {u, 0, Infinity}, Assumptions -> w > 0],
  2 Pi w/(1 + 2 w)^2, "Assumptions" -> w > 0, "TestPoint" -> {w -> 4/5}];
dvELam = check["(48) same-ball term of the radial kernel at large q: A0 (d^3q = 4 pi q^2) v(q) = 1/q^2, two same-ball regions with the Coulomb integral 16 pi^2 Int r^2 f_L dr = 4 pi^2 ((S82), computed), and the frequency integral (17) (computed) with W = q Lambda = q^2 w, |delta| = q^2/2: 3 p_loc(w)/(pi^3 q^2)",
  Simplify[A0 4 Pi q^2 (1/q^2) 2 (16 Pi^2 First[dvS91]) (dvFreq17 /. {W -> q^2 w, d -> q^2/2}), q > 0 && w > 0], 3 (2 Pi w/(1 + 2 w)^2)/(Pi^3 q^2),
  "Assumptions" -> q > 0 && w > 0, "TestPoint" -> {q -> 3, w -> 4/5}];
check["(48) prefactor 32 pi^2 A0 = 3/pi^3", 32 Pi^2 A0, 3/Pi^3];
check["(S92) (Lambda q)^2 J(Lambda q; a, a) at a = q^2/2, Lambda = w q, with J as computed in (S56), is p_loc/(pi q^2)",
  Simplify[(w q^2)^2 (dvJPlus /. {al -> q^2/2, be -> q^2/2, w -> w q^2}), w > 0 && q > 0], (2 Pi w/(1 + 2 w)^2)/(Pi q^2),
  "Assumptions" -> w > 0 && q > 0, "TestPoint" -> {w -> 4/5, q -> 3}];

dvA1int = check["(51) integrand of A1 = q_* e_kite: e_kite = -2 Int lambda dlambda Int theta_lambda e_Lambda dq with e_Lambda computed in (48), theta = lambda/(x^4 + lambda), w = Sqrt[1 + lambda/x^4]/2 ((21) with t = x), q = q_* x, dq = q_* dx",
  Simplify[((-2 lam dvTh21 (dvELam /. w -> dvLam21/q)) /. q -> qs x) qs qs, x > 0 && lam > 0 && qs > 0],
  -(12/Pi^2) lam (lam/(x^4 + lam)) (Sqrt[1 + lam/x^4]/2)/(x^2 (1 + Sqrt[1 + lam/x^4])^2), "Assumptions" -> x > 0 && lam > 0 && qs > 0,
  "TestPoint" -> {x -> 2/3, lam -> 1/2, qs -> 3}];
dvSep = check["(S90) x = lambda^(1/4) y separates the integrand of A1 computed above: lambda^(3/4) times -(12/pi^2)/(2 Sqrt[y^4 + 1] (y^2 + Sqrt[y^4 + 1])^2)",
  Simplify[(dvA1int /. x -> lam^(1/4) y) D[lam^(1/4) y, y], y > 0 && lam > 0],
  lam^(3/4) (-(12/Pi^2))/(2 Sqrt[y^4 + 1] (y^2 + Sqrt[y^4 + 1])^2), "Assumptions" -> y > 0 && lam > 0, "TestPoint" -> {y -> 2/3, lam -> 1/2}];
dvLam34 = check["(S90) Int_0^1 lambda^(3/4) dlambda = 4/7", Integrate[lam^(3/4), {lam, 0, 1}], 4/7];
dvCoefY = check["(S90) A1 = -(48/(7 pi^2)) Int_0^Infinity dy/(2 Sqrt[y^4 + 1] (y^2 + Sqrt[y^4 + 1])^2): the y-independent factor of the separated integrand above times Int lambda^(3/4)",
  Simplify[dvLam34 dvSep 2 Sqrt[y^4 + 1] (y^2 + Sqrt[y^4 + 1])^2/lam^(3/4), y > 0 && lam > 0], -48/(7 Pi^2), "Assumptions" -> y > 0 && lam > 0, "TestPoint" -> {y -> 2/3, lam -> 1/2}];
dvTauInt = check["(S91) y^2 = tau/Sqrt[1 - tau^2] turns the integrand of (S90) into (1/4) tau^(-1/2) (1 - tau)^(1/4) (1 + tau)^(-7/4) on 0 < tau < 1",
  Simplify[(1/(2 Sqrt[y^4 + 1] (y^2 + Sqrt[y^4 + 1])^2) /. y -> Sqrt[tau/Sqrt[1 - tau^2]]) D[Sqrt[tau/Sqrt[1 - tau^2]], tau], 0 < tau < 1],
  (1/4) tau^(-1/2) (1 - tau)^(1/4) (1 + tau)^(-7/4), "Assumptions" -> 0 < tau < 1, "TestPoint" -> {tau -> 3/10}];
dvSigInt = check["(S91) sigma = 2 tau/(1 + tau) turns the tau integrand computed above into 2^(-1/2) (1/4) sigma^(-1/2) (1 - sigma)^(1/4) on 0 < sigma < 1",
  Simplify[(dvTauInt /. tau -> sig/(2 - sig)) D[sig/(2 - sig), sig], 0 < sig < 1], 2^(-1/2) (1/4) sig^(-1/2) (1 - sig)^(1/4),
  "Assumptions" -> 0 < sig < 1, "TestPoint" -> {sig -> 3/10}];
dvIy = check["(S91) the sigma integral computed: (1/4) 2^(-1/2) B(1/2, 5/4)", Integrate[dvSigInt, {sig, 0, 1}], (1/4) 2^(-1/2) Beta[1/2, 5/4]];
check["(S91),(51) A1 = (factor computed in (S90)) times the y integral computed above = -(12/(7 sqrt2 pi^2)) B(1/2, 5/4) = -4 Gamma(1/4)/(7 sqrt2 pi^(3/2) Gamma(3/4))",
  {dvCoefY dvIy, dvCoefY dvIy}, {-12/(7 Sqrt[2] Pi^2) Beta[1/2, 5/4], -4 Gamma[1/4]/(7 Sqrt[2] Pi^(3/2) Gamma[3/4])}];
numcheck["(51) A1 (truncated in the manuscript)", -4 Gamma[1/4]/(7 Sqrt[2] Pi^(3/2) Gamma[3/4]), "-0.2146937609"];
numcheck["(51) q_*^-1 = (3 pi/(16 alphaL))^(1/4) rs^-1/4", (3 Pi/(16 alphaL))^(1/4), "1.0311348"];

numcheck["(Tab. 3) e_kite leading coefficient", (3 Pi/(16 alphaL))^(1/4) (-4 Gamma[1/4]/(7 Sqrt[2] Pi^(3/2) Gamma[3/4])), "-0.2213782293"];
numcheck["(Tab. 3) F = (7/8) e_kite", 7/8 (3 Pi/(16 alphaL))^(1/4) (-4 Gamma[1/4]/(7 Sqrt[2] Pi^(3/2) Gamma[3/4])), "-0.1937059506"];

check["(54) d/dPi [Pi/(1 - lambda v Pi)] = 1/eps^2", D[Pp/(1 - lam v Pp), Pp], 1/(1 - lam v Pp)^2,
  "TestPoint" -> {Pp -> -1/3, lam -> 1/2, v -> 2}];
dvFsum = check["(26),(text SI Sec. S17.4) free f-sum identity Int Delta_p delta_p d^3p = 4 pi q^2/3 (second ball shifted by p -> p - q)",
  Integrate[(r cc q + q^2/2) 2 Pi r^2, {r, 0, 1}, {cc, -1, 1}] - Integrate[(r cc q - q^2/2) 2 Pi r^2, {r, 0, 1}, {cc, -1, 1}],
  4 Pi q^2/3, "TestPoint" -> {q -> 2}];
check["(26) free response at large nu: nu^2 delta/(delta^2 + nu^2) -> delta, so nu^2 L(q, i nu) -> (1/(4 pi)) Int Delta_p delta_p d^3p (computed above), Pi_0 = -q_TF^2 L gives -q_TF^2 q^2/(3 nu^2), and q_TF^2 = B_3/(k_F a_B) = 4/(pi k_F a_B) by (6)",
  {Limit[nu^2 del/(del^2 + nu^2), nu -> Infinity], -qTF^2 dvFsum/(4 Pi), calB[3]/kFaB},
  {del, -qTF^2 q^2/3, 4/(Pi kFaB)}, "TestPoint" -> {del -> 3/7, qTF -> 2/5, q -> 1/2, kFaB -> 3}];
check["(26) assembled: Pi_0(q, i nu) ~ -q_TF^2 (1/(4 pi)) (Int Delta_p delta_p d^3p computed above)/nu^2 = -q_TF^2 q^2/(3 nu^2), the free value of chi_lambda",
  -qTF^2 (dvFsum/(4 Pi))/nu^2, -qTF^2 q^2/(3 nu^2), "TestPoint" -> {qTF -> 2/5, q -> 1/2, nu -> 3}];

check["(55),(S23) g_p^2 + g_k^2 - 2 g_p g_k = (g_p - g_k)^2 and delta_p - delta_k = q.(p-k)",
  {Expand[gp^2 + gk^2 - 2 gp gk - (gp - gk)^2],
   Expand[(px qx + py qy + pz qz + (qx^2 + qy^2 + qz^2)/2) - (kx qx + ky qy + kz qz + (qx^2 + qy^2 + qz^2)/2) -
     ((px - kx) qx + (py - ky) qy + (pz - kz) qz)]}, {0, 0}];

check["(S33) Fock cusp f_L(1+u) - f_L(1) = -(u/2) log(2/u) + O(u^2 log u), u > 0 (f_L(1) = 1/2)",
  Module[{ser}, ser = Normal[Series[(1 - (1 + uu)^2)/(4 (1 + uu)) Log[(2 + uu)/uu], {uu, 0, 1}, Assumptions -> uu > 0]];
    Simplify[PowerExpand[ser + uu/2 Log[2/uu]], uu > 0]], 0, "Assumptions" -> uu > 0, "TestPoint" -> {uu -> 1/1000}];
check["(56) Pi_vtx log coefficient: -(prefactor computed in (S20)) times a_x/(h + z h') with a_x computed in (S30) = -qTF^4/4",
  -dvPiVtxPref Simplify[dvAx/(h[z] + z hp[z]), z > 0], -qTF^4/4, "Assumptions" -> qTF > 0 && z > 0, "TestPoint" -> {qTF -> 2, z -> 7/10}];
check["(S35) Int_0^inf h Re(iz-a)^-2 dz/pi = [log(1+1/a) - 1/(1+a)]/2",
  Integrate[h[z] (a^2 - z^2)/(a^2 + z^2)^2, {z, 0, Infinity}, Assumptions -> a > 0]/Pi,
  (Log[1 + 1/a] - 1/(1 + a))/2, "Assumptions" -> a > 0, "TestPoint" -> {a -> 7/10}];

dvESig = check["(S94) self-energy radial kernel at large q: with f_L(s) = O(q^-2) dropped, the s integral of s/[Lambda q + (s^2 - r^2)/2]^2 over (q - r, q + r) is computed, its leading term at large q with Lambda = w q is taken, and the r integral with f_L (cf. (S82)) gives -3 p_loc(w)/(pi^3 q^2)",
  Module[{sInt, lead},
    sInt = Integrate[s/(L q + (s^2 - r^2)/2)^2, {s, q - r, q + r}, Assumptions -> q > 2 && 0 < r < 1 && L > 0];
    lead = Normal[Series[sInt /. L -> w q, {q, Infinity, 3}]];
    Simplify[-(3 w q/Pi^2) Integrate[r fL[r] lead, {r, 0, 1}], q > 0 && w > 0]],
  -3 (2 Pi w/(1 + 2 w)^2)/(Pi^3 q^2), "Assumptions" -> q > 2 && w > 0, "TestPoint" -> {w -> 4/5, q -> 3}];
check["(57) the self-energy pair cancels the kite at low density: the leading radial kernels computed in (48) and (S94) sum to zero, and both enter e_kite and e_SE with the same weight -2 Int lambda dlambda theta_lambda dq, so e_SE = -A1/q_*",
  Simplify[dvELam + dvESig], 0];

exportResults["MPP v4.21 (main text Eqs. (1)-(59), SI Eqs. (S1)-(S104))", "wolfram_check_results_v4.21.json"];

