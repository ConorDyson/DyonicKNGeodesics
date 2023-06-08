(* ::Package:: *)

(* ::Title:: *)
(*DionicKNGeoOrbit subpackage of DionicKNGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["DionicKNGeodesics`DionicKNGeoOrbit`"];


DionicKNGeoOrbit::usage = "For generic Orbits and returns a DionicKNGeoOrbitFunction[..] which stores the orbital trajectory and parameters.";
DionicKNGeoOrbitFunction::usage = "DionicKNGeoOrbitFunction[a,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";
Begin["`Private`"];


(* ::Section:: *)
(*DionicKN*)


(* ::Subsection:: *)
(*Real Roots Generic Orbit (Mino)*)


DionicKNGeoOrbit::r0outofboundsGen = "Intial radius `1` is not between the outer radial root `2`, and inner radial root `3`.";
DionicKNGeoOrbit::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";
DionicKNGeoOrbit::TooMuchENergy = "the total value of \!\(\*SuperscriptBox[\(\[ScriptCapitalO]\), \(2\)]\)+\!\(\*SuperscriptBox[\(a\), \(2\)]\) = `1` is greater than one your are breakig cosmic censorship, stop this!"


Options[DionicKNGeoRealOrbitMino] = {"Roots"-> Automatic};
DionicKNGeoRealOrbitMino[a_, \[ScriptCapitalO]_ , \[CapitalOmega]_, \[Delta]_ , \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalK]_, initCoords_, OptionsPattern[Options[DionicKNGeoRealOrbitMino]]] := Module[
	{consts,assoc,M,Mino,t0,r0,\[Phi]0,kr,hR,hP,hM,\[CapitalUpsilon]t,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]tr,\[CapitalUpsilon]tz,\[CapitalUpsilon]z,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,t,r,z,z0,\[Theta]0, \[Theta], J,\[Phi], r1,r2,r3,r4,MCMP,MCPP,MCPM,MCMM, Minoz, MinozFunc, MinoRFunc,m ,z1,z2,z3,z4,\[ScriptCapitalQ], RM, RP,ROOTS,RealRoots,ComplexRoots, kz, Z1, Z2, MinoR},
	
	M=1;
	m=1;
	J = 1-\[ScriptCapitalE]^2;
	\[ScriptCapitalQ] = \[ScriptCapitalK];
		
	RM = M-Sqrt[M^2-a^2-\[ScriptCapitalO]^2];
	RP = M+Sqrt[M^2-a^2-\[ScriptCapitalO]^2];
	consts = <|"a"->a,"\[ScriptCapitalO]"->\[ScriptCapitalO],"\[Delta]"->\[Delta], "\[CapitalOmega]"->\[CapitalOmega], "\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalK]"->\[ScriptCapitalK]|>;
	
	If[\[ScriptCapitalO]^2+a^2 >1,
		Message[DionicKNGeoPlunge::TooMuchENergy,\[ScriptCapitalO]^2+a^2];
		Return[$Failed]
		];
		
	 
	rROOTS =r/.Solve[(\[Delta]  r-m (a^2+r^2) \[ScriptCapitalE]+a m \[ScriptCapitalL])^2/m^2-(r^2+(-a \[ScriptCapitalE]+\[ScriptCapitalL])^2+\[ScriptCapitalK])(r-RP)(r-RM)==0,r];
	{r1,r2,r3,r4}= rROOTS;
	kr = Sqrt[((r1-r2) (r3-r4))/((r1-r3) (r2-r4))];
	Ar = -(1-\[ScriptCapitalE]^2);
	AMr[\[Lambda]_]:=JacobiAmplitude[Sqrt[-Ar(r1-r3) (r2-r4) ]/2 \[Lambda],kr^2];

	MinoR[r_]:=(-2  EllipticF[ArcSin[ Sqrt[((r1-r3)(r-r4))/((r3-r4)(r1-r))]],kr^2])/Sqrt[-Ar(r1-r3) (r2-r4) ];
	r[\[Lambda]_]:=((r1-r3) r4+r1 (r3-r4) Sin[AMr[\[Lambda]]]^2)/(r1-r3+(r3-r4) Sin[AMr[\[Lambda]]]^2);
	RINT\[Lambda][\[Lambda]_]:= 2  /Sqrt[-Ar(r1-r3) (r2-r4) ] (r1 EllipticF[AMr[\[Lambda]],kr^2]+(-r1+r4) EllipticPi[(-r3+r4)/(r1-r3),AMr[\[Lambda]],kr^2]);
	R2INT\[Lambda][\[Lambda]_]:=-(-1/(Sqrt[-Ar (r1-r3) (r2-r4)] (kr^2 (r1-r3)+r3-r4)))((r1-r3) (r1-r4) (r3-r4) EllipticE[AMr[\[Lambda]],kr^2]+(kr^2 (r1-r3)+r3-r4) (r1^2-r3 r4+r1 (r3+r4)) EllipticF[AMr[\[Lambda]],kr^2]-(r1-r4) ((r3-r4) (2 r1+r3+r4)+kr^2 (r1-r3) (r1+r3+2 r4)) EllipticPi[(-r3+r4)/(r1-r3),AMr[\[Lambda]],kr^2]+((r1-r3) (r1-r4) (r3-r4)^2 Cos[AMr[\[Lambda]]] Sin[AMr[\[Lambda]]] Sqrt[1-kr^2 Sin[AMr[\[Lambda]]]^2])/(r1-r3+(r3-r4) Sin[AMr[\[Lambda]]]^2));
	RMINT\[Lambda][\[Lambda]_]:= (-2 ((RM-r4) EllipticF[AMr[\[Lambda]],kr^2]+(-r1+r4) EllipticPi[-(((RM-r1) (r3-r4))/((r1-r3) (RM-r4))),AMr[\[Lambda]],kr^2]))/((RM-r1) (RM-r4) Sqrt[Ar (r1-r3) (-r2+r4)]);

	RPINT\[Lambda][\[Lambda]_]:= (-2 ((RP-r4) EllipticF[AMr[\[Lambda]],kr^2]+(-r1+r4) EllipticPi[-(((RP-r1) (r3-r4))/((r1-r3) (RP-r4))),AMr[\[Lambda]],kr^2]))/((RP-r1) (RP-r4) Sqrt[Ar (r1-r3) (-r2+r4)]);
	
	
	ZROOTS =z/.Solve[\[ScriptCapitalK]-z^2 \[ScriptCapitalK]-z 1/m^2 (\[CapitalOmega]^2 z+2 a m \[CapitalOmega] (-1+z^2) \[ScriptCapitalE]+a^2 m^2 z (-1+z^2) (-1+\[ScriptCapitalE]^2)+2 m \[CapitalOmega] \[ScriptCapitalL]+m^2 z \[ScriptCapitalL]^2)==0,z];
	
	{z4,z3,z2,z1}= ZROOTS;
	kz =Sqrt[((z1-z4) (z2-z3))/((z1-z3) (z2-z4))];
	Az = a^2 (1-\[ScriptCapitalE]^2);
	AMZ[\[Lambda]_]:=JacobiAmplitude[ 1/2 Sqrt[Az (z1-z3) (z2-z4)] \[Lambda],kz^2];
	
	Minoz[z_]:=(2 EllipticF[ArcSin[ Sqrt[((z2-z4)(z-z3))/((z2-z3)(z-z4))]],kz^2])/Sqrt[Az (z1-z3) (z2-z4)];
	z[\[Lambda]_]:=(z3 (-z2+z4)+(z2-z3) z4 Sin[AMZ[\[Lambda]]]^2)/(-z2+z4+(z2-z3) Sin[AMZ[\[Lambda]]]^2);
	ZINT\[Lambda][\[Lambda]_]:= 2 /Sqrt[Az (z1-z3) (z2-z4)] (z4 EllipticF[AMZ[\[Lambda]],kz^2]+(z3-z4) EllipticPi[(z2-z3)/(z2-z4),AMZ[\[Lambda]],kz^2]);
	Z2INT\[Lambda][\[Lambda]_]:=((-((z2-z3) (z2-z4) (z3-z4) EllipticE[AMZ[\[Lambda]],kz^2])-((-1+kz^2) z2+z3-kz^2 z4) (z2 (z3-z4)-z4 (z3+z4)) EllipticF[AMZ[\[Lambda]],kz^2]+(z3-z4) ((-1+kz^2) z2^2+z3^2+2 z2 (kz^2 z3-z4)-2 (-1+kz^2) z3 z4-kz^2 z4^2) EllipticPi[(z2-z3)/(z2-z4),AMZ[\[Lambda]],kz^2]+-(z2-z3) Sin[AMZ[\[Lambda]]] ((z2-z3) (z2-z4) (z3-z4) Cos[AMZ[\[Lambda]]] )/(-z2+z4+(z2-z3) Sin[AMZ[\[Lambda]]]^2) Sqrt[(1-kz^2 Sin[AMZ[\[Lambda]]]^2)])/(Sqrt[Az (z1-z3) (z2-z4)] ((-1+kz^2) z2+z3-kz^2 z4)));
	ZM1INT\[Lambda][\[Lambda]_]:= 2 /Sqrt[Az (z1-z3) (z2-z4)] 1/((1-z3) (1-z4)) ((-1+z3) EllipticF[AMZ[\[Lambda]],kz^2]+(-z3+z4) EllipticPi[((z2-z3) (1-z4))/((1-z3) (z2-z4)),AMZ[\[Lambda]],kz^2]);
	ZP1INT\[Lambda][\[Lambda]_]:= 2 /Sqrt[Az (z1-z3) (z2-z4)] 1/((1+z3) (1+z4)) ((1+z3) EllipticF[AMZ[\[Lambda]],kz^2]+(-z3+z4) EllipticPi[((z2-z3) (-1-z4))/((-1-z3) (z2-z4)),AMZ[\[Lambda]],kz^2]);


		
	z0 = Cos[\[Theta]0];
	
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,r3+0.01,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0<=r3||r0>=r4,
		Message[DionicKNGeoOrbit::r0outofboundsGen,r0,r4,r3];
		Return[$Failed]
		];
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[DionicKNGeoOrbit::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
	
	
	
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[DionicKNGeoOrbit::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
		
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r]- MinoR[r0]],Listable];
	
	If[Abs[Arg]==\[Pi]/2, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	
	If[\[ScriptCapitalK]==0, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	
	


	tr[\[Lambda]_]:=  (2 a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE]+RM RP \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]-a \[ScriptCapitalL]-(RM \[Delta])/m-(RP \[Delta])/m)\[Lambda]+RINT\[Lambda][\[Lambda]] ((RM+RP) \[ScriptCapitalE]-\[Delta]/m)+R2INT\[Lambda][\[Lambda]] \[ScriptCapitalE]-((a^2+RP^2) (a^2 m \[ScriptCapitalE]-a m \[ScriptCapitalL]+RP (m RP \[ScriptCapitalE]-\[Delta])))/(m (RM-RP)) RPINT\[Lambda][\[Lambda]]+((a^2+RM^2) (a^2 m \[ScriptCapitalE]-a m \[ScriptCapitalL]+RM (m RM \[ScriptCapitalE]-\[Delta])))/(m (RM-RP)) RMINT\[Lambda][\[Lambda]] ;
	\[Phi]r[\[Lambda]_]:= (((a (a^2 m \[ScriptCapitalE]-a m \[ScriptCapitalL]+RM (m RM \[ScriptCapitalE]-\[Delta])))/(m (RM-RP)))RMINT\[Lambda][\[Lambda]]+ (((a (-a^2 m \[ScriptCapitalE]+a m \[ScriptCapitalL]+RP (-m RP \[ScriptCapitalE]+\[Delta])))/(m (RM-RP)) ) )RPINT\[Lambda][\[Lambda]])+a \[ScriptCapitalE] \[Lambda]  ;
	tz[\[Lambda]_]:= + a^2 Z2INT\[Lambda][\[Lambda]] \[ScriptCapitalE]+(a \[CapitalOmega])/m ZINT\[Lambda][\[Lambda]]-a (a  \[ScriptCapitalE]- \[ScriptCapitalL])\[Lambda];
	\[Phi]z[\[Lambda]_]:=ZP1INT\[Lambda][\[Lambda]]((m \[ScriptCapitalL]-\[CapitalOmega])/(2 m)) -ZM1INT\[Lambda][\[Lambda]]( (m \[ScriptCapitalL]+\[CapitalOmega])/(2 m)) -a \[ScriptCapitalE] \[Lambda];



(*Note: If including Precission Error May occur in the case of evaluating R2INT\[Lambda][0] due to a bug with mathematica in
 evaluating terms of the form EllipticE[0, num`n], num is any number and n is less thna machine precision*)
	t=Function[{Global`\[Lambda]}, Evaluate[  tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ Minoz[z0]]-tr[MinoR[r0]]-tz[Minoz[z0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda]+ MinoR[r0]]], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + Minoz[z0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[ \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ Minoz[z0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[Minoz[z0]] + \[Phi]0], Listable];
	
	assoc = Association[
		"a" -> a,
		"Parametrization"->"Mino", 
		"ChargeParameters"->{\[ScriptCapitalO] , \[CapitalOmega], \[Delta]},
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalK], 
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {r4,r3,r2,r1},
		"PolarRoots"-> {z1,z2,z3,z3},
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(r\)]\)"-> kr,"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialMinoTime"-> MinoRFunc,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Inner\)], \(\)]\)"-> Abs[MinoR[R1]] - MinoR[r0](* Time outgoing to inner root*),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Outer\)], \(\)]\)"-> Abs[MinoR[R2]] - MinoR[r0](* Time outgoing to outer root*)|>,

		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	DionicKNGeoOrbitFunction[a ,\[ScriptCapitalO], \[CapitalOmega], \[Delta], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalK], assoc]
]


(* ::Section:: *)
(*DionicKNGeoOrbit and DionicKNGeoOrbitFuction*)


DionicKNGeoOrbit::fourcomplexroots = "This function can only currently solve plunging geodesics in the paramater space where we have four finite real roots or two finite real roots and two finite complex roots.";
DionicKNGeoOrbit::highenergy = "This code does not currently support plunging solutions with \[ScriptCapitalE]>1, \[ScriptCapitalE]=`1` has been selected ";
DionicKNGeoOrbit::negativecarter = "This code does not currently support plunging solutions with Q<0, Q=`1` has been selected ";
DionicKNGeoOrbit::rIoutbounds = "Given ISSO radius `1` is not in the allowed range between `2` and `3` for spin `4`.";
DionicKNGeoOrbit::Inclinationoutofbounds = "Given inclination angle `1` is not between -\!\(\*FractionBox[\(\[Pi]\), \(2\)]\) and \!\(\*FractionBox[\(\[Pi]\), \(2\)]\)."
DionicKNGeoOrbit::aoutofbounds1 = "The Orbits packages does not currently support the extremal a=1 case"
DionicKNGeoOrbit::aoutofbounds0 = "The Orbits packages does not currently support the zero spin a=0 case"
DionicKNGeoOrbit::NoBoundMotion = "The Paramaters you provided do not correspond to bound motion"


DionicKNGeoOrbit[a_:0.7 ,\[ScriptCapitalO]_:0.3, \[CapitalOmega]_:3, \[Delta]_:3,\[ScriptCapitalE]_:0.8,\[ScriptCapitalL]_:0.3,\[ScriptCapitalK]_:3, initCoords:{_,_,_,_}:Automatic,OptionsPattern[]]:=Module[
	{param, method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
		
	If[a==1,
		Message[DionicKNGeoOrbit::aoutofbounds1];
		Return[$Failed]
		];
	
	If[\[ScriptCapitalE]>=1,
		Message[DionicKNGeoOrbit::highenergy,\[ScriptCapitalE]];
		Return[$Failed]
		];
		
	If[\[ScriptCapitalK]<0,
		Message[DionicKNGeoOrbit::negativecarter,\[ScriptCapitalK]];
		Return[$Failed]
		];

	ROOTS = r/.NSolve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL] - \[Delta] r)^2-(r^2-2r+a^2+\[ScriptCapitalO]^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalK])==0,r];

	RealRoots = Select[ROOTS,PossibleZeroQ@Im[#]&];
	ComplexRoots = Select[ROOTS,Not@PossibleZeroQ@Im[#]&];

	If[Length[ComplexRoots]== 2, Message[DionicKNGeoOrbit::NoBoundMotion]; Return[$Failed]];
	If[Length[ComplexRoots]== 0, Return[DionicKNGeoRealOrbitMino[a,\[ScriptCapitalO], \[CapitalOmega], \[Delta], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalK], initCoords,"Roots"->ROOTS]]];
	];


DionicKNGeoOrbitFunction /:
 MakeBoxes[kgof:DionicKNGeoOrbitFunction[a_,\[ScriptCapitalO]_,\[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalK]_,assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalO]: ", \[ScriptCapitalO]}], "  ",
                  BoxForm`SummaryItem[{"\[CapitalOmega]: ", \[CapitalOmega]}], "  ",
                  BoxForm`SummaryItem[{"\[Delta]: ", \[Delta]}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalE]: ", \[ScriptCapitalE]}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalL]: ", \[ScriptCapitalL]}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalK]: ", \[ScriptCapitalK]}]}],
             BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}]};
  extended = {BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}],
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["AngularMomentum"]}],
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["CarterConstant"]}]};
  BoxForm`ArrangeSummaryBox[
    DionicKNGeoOrbitFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];


DionicKNGeoOrbitFunction[a_,\[ScriptCapitalO]_,\[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalK]_,assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
DionicKNGeoOrbitFunction[a_,\[ScriptCapitalO]_,\[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalK]_,assoc_][y_?StringQ] := assoc[y]
Keys[g_DionicKNGeoOrbitFunction]^:=Keys[g[[5]]]


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
