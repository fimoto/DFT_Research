(* Symbolic and numerical checks for the reader-oriented JPA manuscript.
   Not executed in a Wolfram kernel by the distributor. No Python dependency.
   Run in a fresh kernel: Get["/full/path/kernels_checks.wl"]
   Each result retains its computed residual; unresolved CAS work is not PASS.
   Notation revision 2026-09-28: qTFsq denotes q_TF^2; the local variable m
   in J4 is the width ratio s_Lambda. Response coefficients use the common
   inverse-interaction unit; F4 retains the matching energy coefficient. *)
Begin["JPAReview`"];
ClearAll["JPAReview`*"];
results = {}; secondsPerCheck = 120;
SetAttributes[verify, HoldAll];
verify[id_, kind_, expression_, assumptions_:True] := Module[{value, status, elapsed},
 {elapsed,value}=AbsoluteTiming[TimeConstrained[
   FullSimplify[expression, Assumptions->assumptions],secondsPerCheck,$Aborted]];
 status=Which[TrueQ[value===0]||TrueQ[value===True],"PASS",
   value===$Aborted,"TIMEOUT",TrueQ[value===False],"FAIL",
   NumberQ[value],"FAIL",True,"UNRESOLVED"];
 AppendTo[results,<|"Check"->id,"Kind"->kind,"Status"->status,
   "Seconds"->elapsed,"Result"->ToString[value,InputForm]|>];
 Print[id," : ",status];
];
SetAttributes[numeric,HoldAll];
numeric[id_,expression_,target_,tolerance_:10^-15] := Module[{v,t,e,status},
 v=TimeConstrained[Check[expression,$Failed],secondsPerCheck,$Aborted];
 t=N[target,40];
 e=If[NumberQ[v],Abs[v-t]/Max[1,Abs[t]],Infinity];
 status=Which[v===$Aborted,"TIMEOUT",v===$Failed,"UNRESOLVED",
   !NumberQ[v],"UNRESOLVED",TrueQ[e<tolerance],"PASS",True,"FAIL"];
 AppendTo[results,<|"Check"->id,"Kind"->"numerical (not rigorous enclosure)",
  "Status"->status,"Result"->ToString[v,InputForm],"ScaledError"->ToString[e,InputForm],
  "Tolerance"->ToString[tolerance,InputForm]|>];Print[id," : ",status];
];

(* 1. Normalization and elementary transformations: current stable TeX source labels are used in result strings.
   Diagram symmetry factor and contraction are physical inputs. *)
pref=3/(32 Pi^5); d0=Pi^3 Log[2]/18-Pi Zeta[3]/4;
verify["N1 prefactor [app-eq:A0]","algebra",2*2*(4 Pi)^2*(2 Pi)^(-9)*3 Pi^2/4-pref];
verify["N2 Onsager conversion [app-eq:Onsager]","algebra",6 d0/Pi^3-(Log[2]/3-3 Zeta[3]/(2 Pi^2))];
verify["N3 ball Fourier transform","integral",
 4 Pi/rho Integrate[p Sin[p rho],{p,0,1},Assumptions->rho>0]-
 4 Pi (Sin[rho]-rho Cos[rho])/rho^3,rho>0];
verify["N4 Abel Coulomb transform","integral",
 Integrate[Exp[-eta r] Sin[k r]/k,{r,0,Infinity},Assumptions->eta>0&&k>0]-1/(k^2+eta^2),eta>0&&k>0];
verify["N5 Frullani regulated logarithm","integral",
 Integrate[Exp[-eta q](Cos[a q]-Cos[b q])/q,{q,0,Infinity},
 Assumptions->eta>0&&a>0&&b>0]-Log[(eta^2+b^2)/(eta^2+a^2)]/2,eta>0&&a>0&&b>0];
numeric["N6 D0 defining reduced integral",NIntegrate[
 Pi/3 Log[(u+1)/(u-1)]*(-Log[u](2/(u+1)-Log[1+1/u])+PolyLog[2,-1/u]+2 Log[1+1/u]),
 {u,1,2,10,Infinity},WorkingPrecision->40,AccuracyGoal->22,PrecisionGoal->22],d0];

(* 2. Frequency kernel: derive directly from rational integrands, not J. *)
Do[With[{sg=sign},
 verify[If[sg==1,"J1 same sign [app-app-eq:Jclosed]","J2 opposite signs [app-app-eq:Jclosed]"],"integral",
 Integrate[(sg a b-x^2)/((x^2+a^2)(x^2+b^2)(x^2+w^2)),
 {x,0,Infinity},Assumptions->a>0&&b>0&&w>0,GenerateConditions->True]/Pi-
 If[sg==1,1/(2 w (a+w)(b+w)),-(a+b+2 w)/(2 w (a+b)(a+w)(b+w))],a>0&&b>0&&w>0]
],{sign,{1,-1}}];
verify["J3 static same-sign cancellation [app-eq:Wdef]","integral",
 Integrate[(a b-x^2)/((x^2+a^2)(x^2+b^2)),{x,0,Infinity},Assumptions->a>0&&b>0],a>0&&b>0];
verify["J4 large-q profile p_loc [app-eq:largeq-block]","integral",
 Integrate[m^2(1/4-x^2)/((x^2+m^2)(x^2+1/4)^2),{x,0,Infinity},Assumptions->m>0]-2 Pi m/(1+2 m)^2,m>0];

(* 3. Mellin and large-pole coefficients. Profile endpoint identification
   and contour displacement remain analytic arguments in Appendices C,E,F. *)
verify["M1 Mellin kernel after Laplace integration [app-app-eq:Ks]","integral",
 Gamma[s+1] Integrate[t^(-s-1)/(1+t^2),{t,0,Infinity},Assumptions->-1<Re[s]<0,GenerateConditions->True]+
 Pi Gamma[s+1]/(2 Sin[Pi s/2]),-1<Re[s]<0];
km[s_]:=-Pi Gamma[s+1]/(2 Sin[Pi s/2]);
Do[With[{pt=row[[1]],ord=row[[2]],val=row[[3]]},
 verify["M2 Laurent at "<>ToString[pt],"limit",Limit[e^ord km[pt+e],e->0]-val]],
 {row,{{0,1,-1},{2,1,2},{4,1,-24},{-1,1,Pi/2},{-2,2,-1},{-3,1,-Pi/4},{-4,2,1/6}}}];
iw[s_]:=Pi Gamma[2-s] Gamma[1+s/2]/(2*2^(2-s)*Gamma[(3-s)/2]^2 Gamma[3-s/2]);
verify["M3 Weber-Schafheitlin I(0)","conditional algebra",iw[0]-1/4];
verify["M4 Weber-Schafheitlin derivative","conditional algebra",(D[iw[s],s]/.s->0)-(7/4-EulerGamma-Log[2])/4];
verify["M5 residue gives large-width correction [app-eq:largeMu]","conditional residue",
 -Residue[km[s]*(-1/(2 (s-1)^2)-(EulerGamma/2-3/8)/(s-1))*L^(-s),{s,1}]-
 (Pi Log[L]/4-Pi/16)/L,L>0];
verify["M6 real moment [app-app-sec:layers]","integral",-Pi/2 Integrate[w Exp[-w],{w,0,Infinity}]+Pi/2];
verify["M7 logarithmic real moment","integral",-Pi/2 Integrate[w Log[w] Exp[-w],{w,0,Infinity}]+Pi(1-EulerGamma)/2];

(* 4. Low density. These test the asymptotic dispersion integral, not
   the reduction of the complete RPA loops to that integral. *)
sLambda[x_,l_]:=Sqrt[1+l/x^4]/2;
kernel[x_,l_]:=l^2 sLambda[x,l]/(x^2(x^4+l)(1+2 sLambda[x,l])^2);
verify["L1 x=lambda^(1/4)t Jacobian [app-eq:A1-dispersion]","substitution",
 kernel[l^(1/4)t,l]*l^(1/4)-l^(3/4)*kernel[t,1],l>0&&t>0];
verify["L2 coupling integral","integral",Integrate[l^(3/4),{l,0,1}]-4/7];
verify["L3 recoil substitution [app-eq:A1-beta]","substitution",
 kernel[(s^2-1)^(-1/4),1]*s/(2(s^2-1)^(5/4))-(s-1)^(1/4)/(4(s+1)^(7/4)),s>1];
verify["L4 remaining improper integral","integral",
 Integrate[(s-1)^(1/4)/(s+1)^(7/4),{s,1,Infinity}]-Beta[5/4,1/2]/Sqrt[2]];
verify["L5 Gamma coefficient [app-eq:A1-closed]","algebra",
 FunctionExpand[-48/(7 Pi^2)*Beta[5/4,1/2]/(4 Sqrt[2])]+4 Gamma[1/4]/(7 Sqrt[2] Pi^(3/2) Gamma[3/4])];
numeric["L6 direct dispersion integral [app-eq:A1-dispersion]",
 -12/Pi^2 NIntegrate[kernel[x,l],{l,0,1},{x,0,1,Infinity},WorkingPrecision->40,
 AccuracyGoal->18,PrecisionGoal->18,MaxRecursion->30],
 -4 Gamma[1/4]/(7 Sqrt[2] Pi^(3/2) Gamma[3/4]),10^-12];

(* 5. High density and exchange insertion. Relative diagram sign and
   energy weight are input assumptions, not inferred by a CAS. *)
h[z_]:=1-z ArcTan[1/z]; alpha[z_]:=1-2 z ArcTan[1/z]+z^2/(1+z^2);
verify["H1 shell angular integral [app-eq:Az]","integral",
 Integrate[c^2(c^2-z^2)/(c^2+z^2)^2,{c,0,1},Assumptions->z>0]-alpha[z],z>0];
verify["H2 h+zhprime","algebra",h[z]+z D[h[z],z]-alpha[z],z>0];
(* h(z)=int_0^1 t^2/(t^2+z^2) dt. Integrate frequency first;
   the positive square integral becomes pi/2 int int a b/(a+b). *)
verify["H3 h integral representation","integral",Integrate[t^2/(t^2+z^2),{t,0,1},Assumptions->z>0]-h[z],z>0];
verify["H4 h squared frequency integral","integral",
 Integrate[a^2 b^2/((z^2+a^2)(z^2+b^2)),{z,0,Infinity},Assumptions->a>0&&b>0]-Pi a b/(2(a+b)),a>0&&b>0];
verify["H5 calA [app-eq:calA]","integral",
 2 Pi^2 Integrate[a b/(a+b),{a,0,1},{b,0,1}]-4 Pi^2(1-Log[2])/3];
verify["H6 window logarithm [app-eq:window]","integral",
 Integrate[Log[1/q]/q,{q,Sqrt[qTFsq],1},Assumptions->0<qTFsq<1]-Log[qTFsq]^2/8,0<qTFsq<1];
verify["H7 dimensionless coupling screening factor","integral",
 2 Integrate[l/(1+l c),{l,0,1},Assumptions->c>0]-2(c-Log[1+c])/c^2,c>0];
verify["H8 Clog2 normalization [app-eq:Clog2]","conditional algebra",
 -Pi pref/3*(4 aL/Pi)*(4 Pi^2(1-Log[2])/3)+aL(1-Log[2])/(6 Pi^3)];
Do[With[{sg=sign},verify["F1 Fock endpoint sign "<>ToString[sg],"series",
 Normal[Series[1/2+(1-(1+sg u)^2)/(4(1+sg u))*Log[(2+sg u)/u]-1/2,{u,0,1}]]+
 sg u Log[2/u]/2,u>0]],{sign,{1,-1}}];
verify["F2 derivative of a shifted resolvent","algebra",
 (D[1/(I z-d-e shift),e]/.e->0)-shift/(I z-d)^2];
verify["F3 whole-shell finite Fock term","integral",
 Integrate[(r+c)Log[2/(r+c)]-r Log[2/(-r)],{r,-c,0},Assumptions->c>0]-c^2(Log[2/c]+1/2),c>0];
verify["F4 selected energy normalization","conditional algebra",
 (-6/(Pi^3 qTFsq^2))*(-qTFsq^2/(32 Pi^2))-2 pref,qTFsq>0];

(* 6. Static coefficients from the full fixed-rectangle integrand (D.4).
   No Taylor coefficients are inserted into the left-hand calculation. *)
lt= q t/2+Sqrt[1-q^2(1-t^2)/4]; ls=q s/2+Sqrt[1-q^2(1-s^2)/4];
xx=lt u; yy=ls v;
hh=1+q xx yy(t+s)/(xx+yy)+q^2/4*((xx t-yy s)^2/(xx+yy)^2-1);
ff=lt^2 ls^2 u v/(xx+yy)^2/Sqrt[hh];
verify["S1 generic static Taylor coefficient [app-app-eq:staticavg]","series and integral",
 Integrate[FullSimplify[Normal[Series[ff,{q,0,2}]],Assumptions->u>0&&v>0&&Element[{t,s},Reals]],
 {t,-1,1},{s,-1,1},Assumptions->u>0&&v>0]/4-
 (u v/(u+v)^2+q^2 u^2 v^2(3 u v-3 u-3 v-5)/(12(u+v)^4)),u>0&&v>0];
verify["S2 c0 [app-app-eq:c0]","integral",Pi/4 Integrate[r/(1+r)^2,{r,0,1}]-Pi(2 Log[2]-1)/8];
verify["S3 c2 [app-app-eq:c2]","integral",-Pi/96 Integrate[r^2(r+14)/(1+r)^4,{r,0,1}]-(Pi/1152-Pi Log[2]/96)];
verify["S4 parity of exact fixed-domain integrand","algebra",
 (ff/.{q->-q,t->-t,s->-s})-ff,u>0&&v>0&&-1<t<1&&-1<s<1&&-1/4<q<1/4];
verify["S5 harmonic-mean bound","quantified inequality",
 Resolve[ForAll[{a,b},Implies[a>=0&&b>=0&&a+b>0,a b/(a+b)<=(a+b)/4]],Reals]];
verify["S6 exact resolvent modulus bound","quantified inequality",
 Resolve[ForAll[{a,b,z},Implies[-1<=a<=1&&-1<=b<=1&&z>=0,
 (a^2+z^2)(b^2+z^2)>=a^2 b^2(1+z^2)^2]],Reals]];
verify["S7 generic logarithmic recurrence","algebra",
 w^(n+1)/(w+qTFsq)-(w^n-qTFsq*w^n/(w+qTFsq)),w>0&&qTFsq>0&&Element[n,Integers]&&n>=0];
verify["S8 static low-density profile integral","integral",
 Integrate[1/(1+t^4),{t,0,Infinity}]-Pi/(2 Sqrt[2])];

(* 7. Frequency-first shell ingredients (B.8-B.10). The geometric
   shell remainder and removal of regulators are not certified here. *)
verify["B1 endpoint L1 mass [app-app-eq:lindhard-shell]","integral",
 Integrate[a-(1-(a-q/2)^2)/(2 q),{a,1-q/2,1},Assumptions->0<q<1]+
 Integrate[(1-(a-q/2)^2)/(2 q),{a,1,1+q/2},Assumptions->0<q<1]-q/4,0<q<1];
verify["B2 diagonal bound [app-app-eq:weighted-shell]","quantified inequality",
 Resolve[ForAll[{a,l},Implies[a>0&&l>0,a^2 l/(2(l+a)^2)<=a/8]],Reals]];
verify["B3 Fock integrated frequency kernel","integral",
 Integrate[l/(2(l+a)^2),{l,0,1},Assumptions->a>0]-(Log[1+1/a]-1/(1+a))/2,a>0];
verify["B4 overlap finite part","integral",
 -2 Pi^3 Integrate[r*Pi(4+r)(2-r)^2/12,{r,0,2}]+16 Pi^4/15];

summary=Counts[Lookup[results,"Status"]];
Print["Wolfram kernel: ",$Version];Print[summary];
Print["PASS certifies only the stated identity or numerical tolerance. TIMEOUT/UNRESOLVED require inspection; no synthetic PASS is used."];
Print["Full RPA asymptotic limits, contour shifts, diagram identification and all-order analytic arguments remain in the manuscript."];
(* Inspect JPAReview`results or export it yourself, e.g.
   Export["jpa_wolfram_results.json", JPAReview`results, "RawJSON"]. *)
End[];
