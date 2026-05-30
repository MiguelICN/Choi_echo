(* ::Package:: *)

(* ::Title::Closed:: *)
(*Setup*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*Get["C:\\Users\\Miguel\\Github\\Choi_echo\\QMB_offline.wl"];*)


(* ::Input:: *)
(*LaunchKernels[10];*)


(* ::Title:: *)
(*Ising chaometer*)


(* ::Input:: *)
(*J=1;*)
(*h=1.05;*)
(*L=4;*)
(*g=0.5;*)


(* ::Input:: *)
(*visited=ConstantArray[False,2^L,SparseArray];*)
(*evenMap={};*)
(*oddMap={};*)
(*Do[*)
(*If[!visited[[i+1]],*)
(*j=revIndex[i,L];*)
(*visited[[i+1]]=True;*)
(*visited[[j+1]]=True;*)
(*If[i==j,AppendTo[evenMap,{i}],AppendTo[evenMap,{i,j}];*)
(*AppendTo[oddMap,{i,j}];*)
(*];]*)
(*,{i,0,(2^L)-1}];*)


(* ::Input:: *)
(*Ha=IsingHamiltonian[h,g,J,L,BoundaryConditions->"Open"];*)
(*system=BlockDiagonalize[Ha, Symmetry->"Parity"];*)
(*{eveneigenval,eveneigenvec}=Transpose[Sort[Transpose[Quiet[Eigensystem[system[[1]]]]]]];*)
(*{oddeigenval,oddeigenvec}=Transpose[Sort[Transpose[Quiet[Eigensystem[system[[2]]]]]]];*)
(*neweven=Table[ExpandBlockVector[i,evenMap,"even",L],{i,eveneigenvec}];*)
(*newodd=Table[ExpandBlockVector[i,oddMap,"odd",L],{i,oddeigenvec}];*)
(*allvecs=Join[neweven,newodd];*)
(*allvals=Join[eveneigenval,oddeigenval];*)
(*newall=SortBy[Transpose[{allvals,allvecs}],First];*)
(*allvecs=newall[[All,2]];*)
(*allvals=newall[[All,1]];*)


(* ::Input:: *)
(*La=2;*)
(*Lb=L-La;*)
(*sdiageigenvec=ParallelTable[sdiag[L,La,Lb,i],{i,allvecs},DistributedContexts->Full];*)


(* ::Input:: *)
(*EntanglementEntropySVD[psi_, dA_, dB_] :=*)
(*  Module[{mat, sv2},*)
(*    mat = ArrayReshape[psi, {dA, dB}];*)
(*    sv2 = SingularValueList[mat]^2;*)
(*    (* Hard threshold at machine epsilon level *)*)
(*    sv2 = Select[sv2, # > 1.*^-14 &];*)
(*    (* Renormalize to absorb discarded numerical noise *)*)
(*    sv2 = sv2 / Total[sv2];*)
(*    -sv2 . Log[sv2]*)
(*  ]*)


(* ::Input:: *)
(*La=2;*)
(*Lb=L-La;*)
(*sdiageigenvec=ParallelTable[EntanglementEntropySVD[i,2^La,2^Lb],{i,allvecs},DistributedContexts->Full];*)


(* ::Input:: *)
(*test=Show[{ListPlot[Transpose[{allvals,sdiageigenvec}],PlotTheme->"Detailed",Axes->True,GridLines->None,PlotLabel->Style["Half chain Entanglement, Open Ising, L = "<>ToString[L]<>", hz = "<>ToString[g]<>"",35,Black],PlotStyle->{Directive[Gray,Opacity[0.5]]},PlotRange->All,ImageSize->1200,FrameStyle->Directive[Black,FontSize->30],FrameLabel->{Style["E"],Style["S(E)"]}],Plot[pageEntropy[La,Lb],{x,Min[allvals],Max[allvals]},PlotStyle->{Directive[Black,Dashed]},PlotLegends->{Style["Page value = "<>ToString[N[pageEntropy[La,Lb]]],20,Black]}]},PlotRange->All];*)


(* ::Input:: *)
(*Export["pagecurveising_L_4.png",test,ImageResolution->300]*)


(* ::Input:: *)
(*productstate2=Table[RandomChainProductState[L],20];*)
(*(*productstate2=Table[Flatten[KroneckerProduct[haar[1],haar[L-1]]],50];*)*)
(*iprs=ParallelTable[Quiet[Total[Abs[Conjugate[allvecs] . i]^4]],{i,productstate2},DistributedContexts->Full];*)
(*energies=ParallelTable[Chop[Conjugate[i] . Ha . i],{i,productstate2},DistributedContexts->Full];*)


(* ::Input:: *)
(*data=SortBy[Transpose[{iprs,energies,productstate2}],First];*)
(*productstate=data[[1;;5]][[All,3]];*)


(* ::Input:: *)
(*tlist=Table[i,{i,0,200,0.05}];*)
(*tlist//Length*)


(* ::Input:: *)
(*fakesvn=Table[ParallelTable[{t,pur[L,La,Lb,i,allvals,allvecs,t]},{t,tlist},DistributedContexts->Full],{i,productstate}];*)


(* ::Input:: *)
(*(*Step 1:Determine the overall energy range*)*)
(*sorted=data[[1;;5]][[All,1]];*)
(*{Emin,Emax}={Min[sorted],Max[sorted]};*)
(*(*Step 2:Define a function to map an energy value to a temperature-like color*)*)
(*colorMapping[e_]:=ColorData["FruitPunchColors"][Rescale[e,{Emin,Emax}]];*)
(*(*plot*)*)
(*plots=Table[ListLogLogPlot[fakesvn[[i]],PlotTheme->"Detailed",PlotLabel->Style["L = "<>ToString[L]<>", La = "<>ToString[La]<>", J = "<>ToString[J]<>", hz = "<>ToString[g]<>", \!\(\*SubscriptBox[\(\[Psi]\), \(E\)]\) = RPS",30,Black],Joined->True,Axes->True,GridLines->None,FrameStyle->Directive[Black,FontSize->30],FrameLabel->{Style["t"],Style["Svn(t)"]},PlotRange->All,ImageSize->700,PlotStyle->Directive[colorMapping[sorted[[i]]]]],{i,Length[productstate]}];*)
(*Legended[Show[{plots,LogLogPlot[pageEntropy[La,Lb],{x,10^-1,Max[tlist]},PlotStyle->{Directive[Black,Dashed]}]}],BarLegend[{"FruitPunchColors",{Emin,Emax}},LabelStyle->{FontSize->12}]]*)


(* ::Title:: *)
(*Choi matrix ising*)


(* ::Chapter:: *)
(*Test*)


(* ::Input:: *)
(*i=1;*)
(*\[Psi]=productstate[[i]];*)
(*rho=Chop[MatrixPartialTrace[Dyad[\[Psi]],1,{2^La,2^Lb}]];*)


(* ::Input:: *)
(*{L,La,Lb}*)
(*\[Psi]//Dimensions*)
(*rho//Dimensions*)


(* ::Input:: *)
(*A=KroneckerProduct[IdentityMatrix[2^(Lb)]/2^(La),rho];*)
(*P=Transpose[allvecs];*)
(*Pinv=Conjugate[allvecs];*)


(* ::Input:: *)
(*A//Dimensions*)
(*allvecs//Dimensions*)


(* ::Input:: *)
(*ClearAll[U];*)
(*U[t_]:=P . DiagonalMatrix[Exp[-I*allvals*t]] . Pinv;*)


(* ::Input:: *)
(*AbsoluteTiming[chois=ParallelTable[ChoiMatrix[U[t],A,La,Lb],{t,tlist}];*)
(*truepurity=Chop[Map[Tr,Map[# . #&,chois]]];]*)


(* ::Input:: *)
(*chois[[1]]//Dimensions*)


(* ::Input:: *)
(*plot=ListLogLogPlot[Transpose[{tlist,truepurity}],PlotTheme->"Detailed",PlotLabel->Style["L = "<>ToString[L]<>", La = "<>ToString[La]<>", J = "<>ToString[J]<>", hz = "<>ToString[g]<>", RPS \!\(\*SubscriptBox[\(\[Psi]\), \(E\)]\)",30,Black],Joined->True,Axes->True,GridLines->None,FrameStyle->Directive[Black,FontSize->30],FrameLabel->{Style["t"],Style["Choi's Purity"]},PlotRange->All,ImageSize->1000,PlotStyle->Directive[colorMapping[sorted[[i]]],Thick]];*)
(*final=Show[{plot,LogLogPlot[1/((2^La)^2),{x,10^-1,tlist[[-1]]},PlotStyle->Directive[Red,Dashed]]}]*)


(* ::Input:: *)
(*Transpose[{tlist,truepurity}][[1;;50]]*)


(* ::Chapter:: *)
(*Improved code c (BEST CODE SO FAR)*)


(* ::Input:: *)
(*Needs["Developer`"];*)
(*ClearAll[PackC];*)
(*PackC[x_]:=Developer`ToPackedArray[x];*)
(*(*1) Environment state from a global pure state|\[Psi]\:27e9 on S\[CircleTimes]E*)*)
(*ClearAll[EnvStateFromPsi];*)
(*EnvStateFromPsi[psi_,m_Integer?Positive,n_Integer?Positive]:=PackC@Chop@MatrixPartialTrace[Dyad[psi],1,{m,n}];*)
(*(*2) Fast spectral unitary builder:U(t)=V.diag(e^{-i \[Lambda] t}).V^\[Dagger],implemented via column scaling (no DiagonalMatrix allocation)*)*)
(*ClearAll[SpectralUnitary];*)
(*SpectralUnitary[allvecs_,allvals_]:=Module[{P=PackC@Transpose[allvecs],Pinv=PackC@Conjugate[allvecs]},Function[{t},With[{ph=PackC@Exp[-I allvals t]},PackC@(Transpose[Transpose[P]*ph] . Pinv)]]];*)
(*ClearAll[PhiProjector];*)
(*PhiProjector[m_Integer?Positive]:=Module[{phi=PackC@(Flatten[IdentityMatrix[m]]/Sqrt[m])},PackC@Dyad[phi]];*)
(*ClearAll[ChoiSysBuilder];*)
(*Options[ChoiSysBuilder]={"rhoE"->None,"psi"->None};*)
(*ChoiSysBuilder[m_Integer?Positive,n_Integer?Positive,unitaryF_Function,opts:OptionsPattern[]]:=Module[{rhoE=OptionValue["rhoE"],psi=OptionValue["psi"],projPhi,idm},If[rhoE===None&&psi===None,Message[ChoiSysBuilder::envspec];Return[$Failed]];*)
(*If[rhoE=!=None&&psi=!=None,Message[ChoiSysBuilder::envboth];Return[$Failed]];*)
(*rhoE=If[rhoE===None,EnvStateFromPsi[psi,m,n],rhoE]//PackC;*)
(*projPhi=PhiProjector[m];*)
(*idm=IdentityMatrix[m,SparseArray];(*sparse identity*)Function[{t},Module[{U=unitaryF[t],UR,Omega,X},UR=KroneckerProduct[idm,U];(*acts on R\[CircleTimes]S\[CircleTimes]E*)Omega=KroneckerProduct[projPhi,rhoE];(* |\[CapitalPhi]\:27e9\:27e8\[CapitalPhi]|\[CircleTimes]\[Rho]E*)X=UR . Omega . ConjugateTranspose[UR];*)
(*PackC@Chop@MatrixPartialTrace[X,3,{m,m,n}]]]];*)
(**)
(*ChoiSysBuilder::envspec="Provide either \"rhoE\" -> (n\[Times]n matrix) or \"psi\" -> (mn-vector).";*)
(*ChoiSysBuilder::envboth="Provide only one of \"rhoE\" or \"psi\", not both.";*)
(**)
(*ClearAll[ProcessPurity];*)
(*Options[ProcessPurity]={"Parallel"->True};*)
(*ProcessPurity[choiF_Function,tlist_List,OptionsPattern[]]:=Module[{par=TrueQ[OptionValue["Parallel"]]},If[par,ParallelMap[(Chop@Re@Tr[# . #])&@*choiF,tlist,Method->"CoarsestGrained"],Map[(Chop@Re@Tr[# . #])&@*choiF,tlist]]];*)


(* ::Input:: *)
(*m=2^La;*)
(*n=2^Lb;*)
(*i=1;*)
(*\[Psi]=productstate[[i]];*)
(*(*Option A:environment given as a pure global state|\[Psi]\:27e9*)*)
(*rhoE=EnvStateFromPsi[\[Psi],m,n];*)
(*unitaryF=SpectralUnitary[allvecs,allvals];*)
(*choiF=ChoiSysBuilder[m,n,unitaryF,"rhoE"->rhoE];*)


(* ::Input:: *)
(*(*Purity over your tlist*)*)
(*AbsoluteTiming[purity=ProcessPurity[choiF,tlist];]*)


(* ::Input:: *)
(*plot=ListLogLogPlot[Transpose[{tlist,purity}],PlotTheme->"Detailed",PlotLabel->Style["L = "<>ToString[L]<>", La = "<>ToString[La]<>", J = "<>ToString[J]<>", hz = "<>ToString[g]<>", RPS \!\(\*SubscriptBox[\(\[Psi]\), \(E\)]\)",30,Black],Joined->True,Axes->True,GridLines->None,FrameStyle->Directive[Black,FontSize->30],FrameLabel->{Style["t"],Style["Choi's Purity"]},PlotRange->All,ImageSize->1000,PlotStyle->Directive[Black,Dashed]];*)
(*finalb=Show[{plot,LogLogPlot[1/((2^La)^2),{x,10^-1,tlist[[-1]]},PlotStyle->Directive[Red,Dashed]]}];*)
(*Show[{final,finalb}]*)
(*decoherence[purity,tlist]*)


(* ::Input:: *)
(*Transpose[{tlist,purity}][[1;;50]]*)
