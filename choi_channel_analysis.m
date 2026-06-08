(* ::Package:: *)

(* ::Title:: *)
(*Setup*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*Get["C:\\Users\\Miguel\\Github\\Choi_echo\\QMB_offline.wl"];*)


(* ::Input:: *)
(*LaunchKernels[10];*)


(* ::Title:: *)
(*Choi Channel Analysis*)


(* ::Chapter::Closed:: *)
(*Definitions*)


(* ::Text:: *)
(*All reusable functions. No model parameters appear here.*)
(*Run this chapter once per session before evaluating anything downstream.*)


(*
(*No Mathematica built-in symbol is used as a local variable name anywhere in*)
(*this chapter. In particular, the Choi matrix argument is always choiMat,    *)
(*never D (reserved as the derivative operator).                              *)








(* ::Section::Closed:: *)
(*BLOCK 0 -- Kernel and packing utilities*)


(* ::Input:: *)
(*$CompileTarget = "WVM";*)
(*Needs["Developer`"];*)
(*ClearAll[PackedC];*)
(*PackedC[x_] := Developer`ToPackedArray[N[x]];*)


(* ::Section::Closed:: *)
(*BLOCK 1 -- Spectral unitary*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* Returns Function[t, U(t)] where                                  *)*)
(*(*   U(t) = P . Diag(e^{-i lambda t}) . P^dagger                   *)*)
(*(* Column-scaling implementation: avoids DiagonalMatrix allocation. *)*)
(*(*                                                                   *)*)
(*(* Inputs:                                                           *)*)
(*(*   allvecs -- D x D, rows are eigenvectors                         *)*)
(*(*   allvals -- length-D eigenvalue vector                           *)*)
(*(* ================================================================ *)*)
(*ClearAll[SpectralUnitary];*)
(*SpectralUnitary[allvecs_, allvals_] :=*)
(*  Module[{Pmat  = PackedC@Transpose[allvecs],*)
(*          Pinv  = PackedC@Conjugate[allvecs]},*)
(*    Function[{t},*)
(*      With[{ph = PackedC@Exp[-I allvals t]},*)
(*        PackedC@(Transpose[Transpose[Pmat]*ph] . Pinv)]]];*)


(* ::Section::Closed:: *)
(*BLOCK 2 -- Phi projector*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* |Phi+><Phi+| on R x S, where |Phi+> = dS^{-1/2} Sum_i |i>|i>.  *)*)
(*(* Precomputed once per session; independent of t and rhoE.         *)*)
(*(* ================================================================ *)*)
(*ClearAll[MakePhiProjector];*)
(*MakePhiProjector[dS_Integer?Positive] :=*)
(*  Module[{phi = PackedC@(Flatten[IdentityMatrix[dS]]/Sqrt[dS])},*)
(*    PackedC@Dyad[phi]];*)


(* ::Section::Closed:: *)
(*BLOCK 3 -- Choi channel builder*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* Constructs the time-dependent Choi matrix of the dynamical map   *)*)
(*(* E_t acting on subsystem S, via the system-environment            *)*)
(*(* representation (paper Eq. 2):                                    *)*)
(*(*                                                                   *)*)
(*(*   choiMat(t) = Tr_E[ (I_R x U(t)) (|Phi+><Phi+| x rhoE)         *)*)
(*(*                        (I_R x U(t))^dagger ]                     *)*)
(*(*                                                                   *)*)
(*(* where U(t) acts on S x E and I_R is the identity on the          *)*)
(*(* reference copy R of S.                                            *)*)
(*(*                                                                   *)*)
(*(* The returned object is a Function[t, choiMat(t)], i.e. a closure  *)*)
(*(* that holds rhoE, projPhi, and the spectral unitary internally.   *)*)
(*(*                                                                   *)*)
(*(* Options:                                                          *)*)
(*(*   "rhoE" -> dE x dE matrix  (environment density matrix)         *)*)
(*(*   "psi"  -> dS*dE vector    (global pure state; rhoE extracted)  *)*)
(*(*                                                                   *)*)
(*(* NOTE: use "rhoE" directly when the environment is already a pure  *)*)
(*(* product state (rhoE = |psiE><psiE|). "psi" is for cases where    *)*)
(*(* the environment state must be obtained by tracing a global state. *)*)
(*(* ================================================================ *)*)
(*ClearAll[ChoiChannelF];*)
(*Options[ChoiChannelF] = {"rhoE" -> None, "psi" -> None};*)
(**)
(*ChoiChannelF[dS_Integer?Positive, dE_Integer?Positive,*)
(*             unitaryF_Function, opts:OptionsPattern[]] :=*)
(*  Module[{rhoEval = OptionValue["rhoE"],*)
(*          psiVal  = OptionValue["psi"],*)
(*          projPhi, idRef},*)
(*    (* --- validate inputs --- *)*)
(*    If[rhoEval === None && psiVal === None,*)
(*      Message[ChoiChannelF::envspec]; Return[$Failed]];*)
(*    If[rhoEval =!= None && psiVal =!= None,*)
(*      Message[ChoiChannelF::envboth]; Return[$Failed]];*)
(*    (* --- build rhoE --- *)*)
(*    rhoEval = If[rhoEval === None,*)
(*      (* extract environment state from global pure state *)*)
(*      PackedC@Chop@MatrixPartialTrace[Dyad[psiVal], 1, {dS, dE}],*)
(*      PackedC@rhoEval];*)
(*    (* --- precompute time-independent objects --- *)*)
(*    projPhi = MakePhiProjector[dS];   (* |Phi+><Phi+|, dS^2 x dS^2 *)*)
(*    idRef   = IdentityMatrix[dS, SparseArray];   (* I_R, dS x dS *)*)
(*    (* --- return closure --- *)*)
(*    Function[{t},*)
(*      Module[{Ut = unitaryF[t], URt, OmegaMat, fullMat},*)
(*        URt      = KroneckerProduct[idRef, Ut];          (* I_R x U(t) *)*)
(*        OmegaMat = KroneckerProduct[projPhi, rhoEval];   (* |Phi+><Phi+| x rhoE *)*)
(*        fullMat  = URt . OmegaMat . ConjugateTranspose[URt];*)
(*        (* Trace out E: subsystem 3 of {dS, dS, dE} *)*)
(*        PackedC@Chop@MatrixPartialTrace[fullMat, 3, {dS, dS, dE}]]]];*)
(**)
(*ChoiChannelF::envspec =*)
(*  "Provide either \"rhoE\" -> (dE x dE matrix) or \"psi\" -> (dS*dE vector).";*)
(*ChoiChannelF::envboth =*)
(*  "Provide only one of \"rhoE\" or \"psi\", not both.";*)


(* ::Section::Closed:: *)
(*BLOCK 4 -- Choi state observables*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* All scalar quantities derived from the dS^2 x dS^2 Choi matrix. *)*)
(*(*                                                                   *)*)
(*(* Argument name: choiMat (never D, which is the derivative         *)*)
(*(* operator in Mathematica -- using D causes shadowing warnings      *)*)
(*(* that flood parallel links and kill worker kernels).               *)*)
(*(*                                                                   *)*)
(*(* ChoiPurity    : Tr[choiMat^2]. Range [1/dS^2, 1].               *)*)
(*(* ChoiVNEntropy : -Tr[choiMat log choiMat]. Range [0, log(dS^2)].  *)*)
(*(* ChoiEVals     : sorted descending eigenvalues of choiMat.         *)*)
(*(* ChoiAllObs    : single aggregator call returning an Association.  *)*)
(*(* ================================================================ *)*)
(*ClearAll[ChoiPurity, ChoiVNEntropy, ChoiEVals, ChoiAllObs];*)
(**)
(*ChoiPurity[choiMat_] := Chop@Re@Tr[choiMat . choiMat];*)
(**)
(*ChoiVNEntropy[choiMat_] :=*)
(*  Module[{ev = Select[Re[Eigenvalues[N[choiMat]]], # > 1.*^-14 &]},*)
(*    ev = ev/Total[ev];*)
(*    Chop@(-ev . Log[ev])];*)
(**)
(*ChoiEVals[choiMat_] :=*)
(*  Sort[Select[Re[Eigenvalues[N[choiMat]]], # > 1.*^-14 &], Greater];*)
(**)
(*ChoiAllObs[choiMat_] :=*)
(*  Module[{ev = Sort[Select[Re[Eigenvalues[N[choiMat]]], # > 1.*^-14 &], Greater],*)
(*          evNorm},*)
(*    evNorm = ev/Total[ev];*)
(*    <|"Purity"      -> Chop@Re@Tr[choiMat . choiMat],*)
(*      "VNEntropy"   -> Chop@(-evNorm . Log[evNorm]),*)
(*      "Eigenvalues" -> ev,*)
(*      "Rank"        -> Length[ev]|>];*)


(* ::Section::Closed:: *)
(*BLOCK 5 -- Time-series driver*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* ChoiSeries: evaluates choiF at each t in tlist and computes      *)*)
(*(* observables. Returns a List of Associations keyed by "Time".     *)*)
(*(*                                                                   *)*)
(*(* "Observables" option:                                             *)*)
(*(*   "Purity"  -- only Tr[choiMat^2]; no Eigenvalues call           *)*)
(*(*   "Entropy" -- purity + von Neumann entropy                       *)*)
(*(*   "All"     -- all of the above plus eigenvalue list and rank     *)*)
(*(*                                                                   *)*)
(*(* "Parallel" option: True (default) uses ParallelMap with           *)*)
(*(* DistributedContexts -> None (workers already loaded upstream).   *)*)
(*(* ================================================================ *)*)
(*ClearAll[ChoiSeries];*)
(*Options[ChoiSeries] = {"Observables" -> "All", "Parallel" -> True};*)
(**)
(*ChoiSeries[choiFn_Function, tlist_List, OptionsPattern[]] :=*)
(*  Module[{obsOpt = OptionValue["Observables"],*)
(*          parOpt = TrueQ[OptionValue["Parallel"]],*)
(*          stepFn},*)
(*    stepFn = Function[{t},*)
(*      Module[{choiMat = choiFn[t], obsResult},*)
(*        obsResult = Switch[obsOpt,*)
(*          "Purity",*)
(*            <|"Purity" -> ChoiPurity[choiMat]|>,*)
(*          "Entropy",*)
(*            Module[{ev = Sort[Select[Re[Eigenvalues[N[choiMat]]],*)
(*                                    # > 1.*^-14 &], Greater]},*)
(*              ev = ev/Total[ev];*)
(*              <|"Purity"    -> Chop@Re@Tr[choiMat . choiMat],*)
(*                "VNEntropy" -> Chop@(-ev . Log[ev])|>],*)
(*          _, (* "All" *)*)
(*            ChoiAllObs[choiMat]];*)
(*        Join[<|"Time" -> t|>, obsResult]]];*)
(*    If[parOpt,*)
(*      ParallelMap[stepFn, tlist,*)
(*        Method -> "CoarsestGrained",*)
(*        DistributedContexts -> None],*)
(*      Map[stepFn, tlist]]];*)


(* ::Section::Closed:: *)
(*BLOCK 6 -- Result extraction*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* ExtractTS  : extracts {time, value} pairs for a scalar key.      *)*)
(*(* ExtractEVals: pads eigenvalue lists to uniform length nEV and    *)*)
(*(*   returns a ragged-safe matrix for plotting.                      *)*)
(*(*   Missing eigenvalues (numerical zero) are filled with 0.        *)*)
(*(* ================================================================ *)*)
(*ClearAll[ExtractTS, ExtractEVals];*)
(**)
(*ExtractTS[res_List, key_String] :=*)
(*  {#["Time"], #[key]} & /@ res;*)
(**)
(*ExtractEVals[res_List, nKeep_Integer] :=*)
(*  (* Returns a nKeep-element list of {{t, lambda_k} ...} series.    *)*)
(*  Table[*)
(*    {#["Time"],*)
(*     If[Length[#["Eigenvalues"]] >= k, #["Eigenvalues"][[k]], 0.]} & /@ res,*)
(*    {k, 1, nKeep}];*)


(* ::Chapter::Closed:: *)
(*Model and Bipartition*)


(* ::Text:: *)
(*This is the only chapter that changes between runs.*)
(*All downstream chapters depend on the symbols defined here and run unchanged.*)


(* ::Section:: *)
(*Physical parameters*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* ISING MODEL: H = Sum_i (hx sx_i + hz sz_i) - J Sum_i sz_i sz_{i+1} *)*)
(*(* Open boundary conditions.                                        *)*)
(*(* To use a different model: replace Ha = ... in the Hamiltonian    *)*)
(*(* chapter with your preferred 2^L x 2^L matrix.                   *)*)
(*(* ================================================================ *)*)
(*J  = 1;      (* ZZ coupling *)*)
(*hx = 1;      (* transverse field *)*)
(*hz = 0.5;    (* longitudinal field; hz != 0 and J != 0 breaks integrability *)*)
(*L  = 6;      (* chain length *)*)


(* ::Section:: *)
(*Bipartition*)


(* ::Input:: *)
(*(* Subsystem S: first La spins (left end of the open chain).        *)*)
(*(* Environment E: remaining Lb spins.                               *)*)
(*(* Paper convention: La = 1 (single probe spin).                    *)*)
(*La = 1;*)
(*Lb = L - La;*)
(*dS = 2^La;   (* subsystem Hilbert space dimension *)*)
(*dE = 2^Lb;   (* environment Hilbert space dimension *)*)


(* ::Section::Closed:: *)
(*Time grid*)


(* ::Input:: *)
(*tMax  = 50.;*)
(*dt    = 0.1;*)
(*tlist = Range[0., tMax, dt];*)
(*Print["Time steps: ", Length[tlist]];*)


(* ::Chapter::Closed:: *)
(*Hamiltonian*)


(* ::Text:: *)
(*Build and diagonalize the full chain Hamiltonian.*)
(*Parity block-diagonalization is used for the Ising model to reduce*)
(*the matrix size by roughly half before the Eigenvalues call.*)
(*For a different model, replace only the "Build" section.*)


(* ::Section:: *)
(*Build*)


?LeaSpinChainHamiltonian


LeaSpinChainHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]


(* ::Input:: *)
(*Ha = IsingHamiltonian[hx, hz, J, L, BoundaryConditions -> "Open"];*)
(*Print["Ha dimensions: ", Dimensions[Ha]];*)


Ha=LeaSpinChainHamiltonian[5.,1.,0,2.,L,5];
Print["Ha dimensions: ", Dimensions[Ha]];


(* ::Section:: *)
(*Diagonalize with parity block structure*)


(* ::Input:: *)
(*system = BlockDiagonalize[Ha, Symmetry -> "Parity"];*)
(**)
(*{evEven, vecEven} = Transpose[Sort[Transpose[Quiet[Eigensystem[system[[1]]]]]]];*)
(*{evOdd,  vecOdd}  = Transpose[Sort[Transpose[Quiet[Eigensystem[system[[2]]]]]]];*)
(**)
(*{evenMap, oddMap} = buildMaps[L];*)
(*vecEvenFull = Table[ExpandBlockVector[v, evenMap, "even", L], {v, vecEven}];*)
(*vecOddFull  = Table[ExpandBlockVector[v, oddMap,  "odd",  L], {v, vecOdd}];*)
(**)
(*allVecs = Join[vecEvenFull, vecOddFull];*)
(*allVals = Join[evEven, evOdd];*)
(*sorted  = SortBy[Transpose[{allVals, allVecs}], First];*)
(*allVecs = sorted[[All, 2]];*)
(*allVals = sorted[[All, 1]];*)
(**)
(*Print["Levels: ", Length[allVals],*)
(*      "  E range: [", N@Min[allVals], ", ", N@Max[allVals], "]"];*)
(*Print["<r~> = ", N@MeanLevelSpacingRatio[allVals]];*)


{allVals,allVecs}=Transpose[Sort[Transpose[Eigensystem[Ha]]]];
Print["Levels: ", Length[allVals],
      "  E range: [", N@Min[allVals], ", ", N@Max[allVals], "]"];
Print["<r~> = ", N@MeanLevelSpacingRatio[allVals]];


(* ::Chapter::Closed:: *)
(*Environment Initial State*)


(* ::Text:: *)
(*Defines rhoE = |psiE><psiE|, the initial pure state of the environment.*)
(*The Choi matrix depends on this choice via paper Eq. 2.*)
(**)
(*"rhoE" is passed directly to ChoiChannelF because psiE lives entirely*)
(*in H_E. EnvStateFromPsi (which traces a global S x E state over S)*)
(*is only needed when the environment state is not already explicit.*)
(**)
(*The Haar-averaged Choi echo (paper Eq. 7) is recovered by repeating*)
(*the computation below over many independent seeds and averaging.*)
(*For a single realization, fix the seed for reproducibility.*)


(* ::Section:: *)
(*Random product state*)


(* ::Input:: *)
(*SeedRandom[42];   (* remove for true randomness *)*)
(*psiE = RandomChainProductState[Lb];*)
(*rhoE = PackedC@Dyad[psiE];*)
(*Print["rhoE: ", Dimensions[rhoE], "  Tr = ", N@Tr[rhoE], "  (should be {", dE, ",", dE, "} and 1)"];*)


(* ::Chapter::Closed:: *)
(*Dynamical Map: Choi Matrix*)


(* ::Text:: *)
(*Constructs the Function[t, choiMat(t)] closure and distributes all*)
(*required objects to parallel workers before the time loop.*)


(* ::Section:: *)
(*Build unitary and channel function*)


(* ::Input:: *)
(*unitaryFn = SpectralUnitary[allVecs, allVals];*)
(*choiFn    = ChoiChannelF[dS, dE, unitaryFn, "rhoE" -> rhoE];*)
(**)
(*(* Sanity check at t = 0.                                           *)*)
(*(* U(0) = identity => E_0 is the identity channel => choiMat(0)    *)*)
(*(* equals |Phi+><Phi+|, a rank-1 projector with Tr[choiMat^2] = 1. *)*)
(*choi0 = choiFn[0.];*)
(*Print["t=0  Tr[choiMat^2] = ", ChoiPurity[choi0], "  (should be 1)"];*)
(*Print["t=0  Tr[choiMat]   = ", N@Re@Tr[choi0],    "  (should be 1)"];*)


(* ::Section:: *)
(*Distribute to parallel kernels*)


(* ::Input:: *)
(*(* Distribute every symbol accessed inside the parallel time loop.  *)*)
(*(* DistributedContexts -> None in ChoiSeries prevents a redundant   *)*)
(*(* automatic redistribution on each ParallelMap call.               *)*)
(*DistributeDefinitions[*)
(*  choiFn,*)
(*  PackedC, MakePhiProjector, ChoiChannelF,*)
(*  ChoiPurity, ChoiVNEntropy, ChoiEVals, ChoiAllObs,*)
(*  MatrixPartialTrace, Dyad*)
(*];*)


(* ::Chapter:: *)
(*Time Evolution*)


(* ::Text:: *)
(*Evaluates choiMat(t) and its observables over tlist.*)
(*Use "Purity" for the fastest single-quantity sweep;*)
(*use "All" when eigenvalue spectra are needed.*)


(* ::Section:: *)
(*Compute*)


(* ::Input:: *)
(*AbsoluteTiming[*)
(*  results = ChoiSeries[choiFn, tlist,*)
(*    "Observables" -> "All",*)
(*    "Parallel"    -> True];*)
(*]*)


(* ::Section:: *)
(*Quick check*)


(* ::Input:: *)
(*Print["Results length: ", Length[results]];*)
(*Print["First entry: ", results[[1]]];*)
(*Print["Last entry:  ", results[[-1]]];*)


(* ::Chapter:: *)
(*Export*)


(* ::Text:: *)
(*Exports the full results list and parameter record.*)
(*The filename encodes the run parameters for unambiguous identification.*)


(* ::Input:: *)
(*runTag = "L" <> ToString[L] <>*)
(*         "_La" <> ToString[La] <>*)
(*         "_hx" <> ToString[hx] <>*)
(*         "_hz" <> ToString[hz] <>*)
(*         "_J"  <> ToString[J]  <>*)
(*         "_T"  <> ToString[tMax];*)
(**)
(*(* Full results list (Associations, one per time step) *)*)
(*Export["choi_results_" <> runTag <> ".mx", results];*)
(**)
(*(* Plain numeric table: {t, purity, entropy} for external use *)*)
(*numericTable = Table[*)
(*  {results[[k]]["Time"],*)
(*   results[[k]]["Purity"],*)
(*   results[[k]]["VNEntropy"]},*)
(*  {k, Length[results]}];*)
(*Export["choi_numeric_" <> runTag <> ".csv", numericTable];*)
(**)
(*(* Parameter record *)*)
(*params = <|"L" -> L, "La" -> La, "Lb" -> Lb, "dS" -> dS, "dE" -> dE,*)
(*            "J" -> J, "hx" -> hx, "hz" -> hz,*)
(*            "tMax" -> tMax, "dt" -> dt,*)
(*            "rMean" -> N@MeanLevelSpacingRatio[allVals]|>;*)
(*Export["choi_params_" <> runTag <> ".mx", params];*)
(**)
(*Print["Exported: ", runTag];*)


(* ::Chapter:: *)
(*Results*)


(* ::Text:: *)
(*All plots use linear axes for t to avoid the t=0 singularity on log scales,*)
(*except the eigenvalue plot which uses a log-linear axis because the early*)
(*dynamics span several orders of magnitude.*)
(*The first 4 eigenvalues are plotted regardless of dS to keep the figure*)
(*readable for larger bipartitions.*)


(* ::Section:: *)
(*Choi purity*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* Choi purity Tr[choiMat(t)^2]                                    *)*)
(*(*                                                                   *)*)
(*(* Bounds:                                                           *)*)
(*(*   1        at t = 0  (unitary channel, identity map)             *)*)
(*(*   1/dS^2   saturation (completely depolarizing channel)          *)*)
(*(* ================================================================ *)*)
(*purityTS = ExtractTS[results, "Purity"];*)
(**)
(*plotLabel = "L = " <> ToString[L] <>*)
(*            ",  La = " <> ToString[La] <>*)
(*            ",  hx = " <> ToString[hx] <>*)
(*            ",  hz = " <> ToString[hz] <>*)
(*            ",  J = "  <> ToString[J];*)
(**)
(*plotPurity = ListLinePlot[purityTS,*)
(*  Joined      -> True,*)
(*  PlotTheme   -> "Detailed",*)
(*  PlotLabel   -> Style["Choi Purity  " <> plotLabel, 20, Black],*)
(*  FrameStyle  -> Directive[Black, FontSize -> 18],*)
(*  FrameLabel  -> {Style["t", 16], Style["Tr[\[ScriptCapitalD](t)\[CenterDot]\[ScriptCapitalD](t)]", 16]},*)
(*  PlotStyle   -> Directive[Black, Thick],*)
(*  PlotRange   -> All,*)
(*  ImageSize   -> 900];*)
(**)
(*Show[plotPurity,*)
(*  Plot[1/dS^2, {t, 0, tMax},*)
(*    PlotStyle   -> Directive[Red, Dashed, Thick],*)
(*    PlotLegends -> {Style["1/dS\[CenterDot]dS = " <> ToString[N[1/dS^2]], 15, Red]}]]*)


(* ::Section:: *)
(*Export purity plot*)


(* ::Input:: *)
(*Export["plot_purity_L_8_La_4_" <> runTag <> ".png", %, ImageResolution -> 200];*)


(* ::Section::Closed:: *)
(*Von Neumann entropy of the Choi state*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* Channel entropy S(choiMat(t)) = -Tr[choiMat log choiMat]        *)*)
(*(*                                                                   *)*)
(*(* Bounds:                                                           *)*)
(*(*   0          at t = 0  (unitary channel)                         *)*)
(*(*   log(dS^2)  saturation (completely depolarizing channel)        *)*)
(*(* ================================================================ *)*)
(*entropyTS = ExtractTS[results, "VNEntropy"];*)
(**)
(*plotEntropy = ListLinePlot[entropyTS,*)
(*  Joined      -> True,*)
(*  PlotTheme   -> "Detailed",*)
(*  PlotLabel   -> Style["Choi Channel Entropy  " <> plotLabel, 20, Black],*)
(*  FrameStyle  -> Directive[Black, FontSize -> 18],*)
(*  FrameLabel  -> {Style["t", 16], Style["S(\[ScriptCapitalD](t))", 16]},*)
(*  PlotStyle   -> Directive[Black, Thick],*)
(*  PlotRange   -> All,*)
(*  ImageSize   -> 900];*)
(**)
(*Show[plotEntropy,*)
(*  Plot[Log[dS^2], {t, 0, tMax},*)
(*    PlotStyle   -> Directive[Red, Dashed, Thick],*)
(*    PlotLegends -> {Style["log(dS\[CenterDot]dS) = " <> ToString[N@Log[dS^2]], 15, Red]}]]*)


(* ::Section::Closed:: *)
(*Export entropy plot*)


(* ::Input:: *)
(*Export["plot_entropy_" <> runTag <> ".png", %, ImageResolution -> 200];*)


(* ::Section::Closed:: *)
(*Eigenvalue spectrum of the Choi state*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* Eigenvalues of choiMat(t) as a function of time.                *)*)
(*(*                                                                   *)*)
(*(* choiMat(t) is a valid density matrix (CPTP via Choi-Jamiolkowski):*)
(*(*   eigenvalues are real, non-negative, sum to 1.                  *)*)
(*(*   rank = number of Kraus operators in the minimal dilation.      *)*)
(*(*                                                                   *)*)
(*(* At t = 0: choiMat = |Phi+><Phi+| => rank 1, lambda_1 = 1.       *)*)
(*(*                                                                   *)*)
(*(* Only the first 4 eigenvalues are plotted to keep the figure      *)*)
(*(* readable for La > 1 (where dS^2 > 4).                            *)*)
(*(* ================================================================ *)*)


nPlot    = 4;   (* fixed: always show the 4 largest eigenvalues *)
evSeries = ExtractEVals[results, nPlot];

evColors = {Black, Red, Blue, Darker[Green]};
evStyles = Directive[#, Thick] & /@ evColors;
evLabels = Table[Style["\[Lambda]" <> ToString[k], 15], {k, 1, nPlot}];

ListLogLinearPlot[evSeries,
  Joined       -> True,
  PlotTheme    -> "Detailed",
  PlotLabel    -> Style["Choi Eigenvalues (top 4)  " <> plotLabel, 20, Black],
  FrameStyle   -> Directive[Black, FontSize -> 18],
  FrameLabel   -> {Style["t", 16], Style["\[Lambda]k(t)", 16]},
  PlotStyle    -> evStyles,
  PlotRange    -> All,
  ImageSize    -> 900,
  PlotLegends  -> Placed[LineLegend[evColors, evLabels,
                    LegendLayout -> "Column"], Right]]


(* ::Section::Closed:: *)
(*Export eigenvalue plot*)


(* ::Input:: *)
(*Export["plot_eigenvalues_" <> runTag <> ".png", %, ImageResolution -> 200];*)


(* ::Section:: *)
(*Saturation values*)


(* ================================================================ *)
(* Time-averaged saturation (paper Eq. 8 analogue for T = tMax).   *)
(* Uses the trapezoidal rule for both purity and entropy.           *)
(*                                                                   *)
(* decoherence[] from QMB_offline.wl computes the trapezoidal       *)
(* integral of a purity time series and divides by tMax.            *)
(* ================================================================ *)
purityVals  = #["Purity"]   & /@ results;
entropyVals = #["VNEntropy"] & /@ results;

timeAvgPurity = decoherence[purityVals, tlist];

timeAvgEntropy =
  (1/tMax)*Total[
    Table[(entropyVals[[k]] + entropyVals[[k+1]])/2*(tlist[[k+1]] - tlist[[k]]),
          {k, Length[tlist]-50, Length[tlist] - 1}]];

Print["\n--- Saturation summary ---"];
Print["<Tr[D^2]>_T          = ", N@timeAvgPurity];
Print["Depolarizing  1/dS^2 = ", N[1/dS^2]];
Print["<S(D)>_T             = ", N@timeAvgEntropy];
Print["Depolarizing  log(dS^2) = ", N@Log[dS^2]];
Print["<r~>                 = ", N@MeanLevelSpacingRatio[allVals]];


(* ================================================================ *)
(* Time-averaged saturation (paper Eq. 8 analogue for T = tMax).   *)
(* Uses the trapezoidal rule for both purity and entropy.           *)
(*                                                                   *)
(* decoherence[] from QMB_offline.wl computes the trapezoidal       *)
(* integral of a purity time series and divides by tMax.            *)
(* ================================================================ *)
purityVals  = #["Purity"]   & /@ results;
entropyVals = #["VNEntropy"] & /@ results;

timeAvgPurity = decoherence[purityVals, tlist];

timeAvgEntropy =
  (1/tMax)*Total[
    Table[(entropyVals[[k]] + entropyVals[[k+1]])/2*(tlist[[k+1]] - tlist[[k]]),
          {k, Length[tlist]-50, Length[tlist] - 1}]];

Print["\n--- Saturation summary ---"];
Print["<Tr[D^2]>_T          = ", N@timeAvgPurity];
Print["Depolarizing  1/dS^2 = ", N[1/dS^2]];
Print["<S(D)>_T             = ", N@timeAvgEntropy];
Print["Depolarizing  log(dS^2) = ", N@Log[dS^2]];
Print["<r~>                 = ", N@MeanLevelSpacingRatio[allVals]];

