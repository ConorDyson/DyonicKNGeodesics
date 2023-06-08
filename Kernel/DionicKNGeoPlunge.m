(* ::Package:: *)

(* ::Title:: *)
(*DionicKNGeoOrbit subpackage of DionicKNGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["DionicKNGeodesics`DionicKNGeoPlunge`"];


DionicKNGeoPlunge::usage = "Takes either DionicKNGeoPlunge[a, {\[ScriptCapitalO],\[CapitalOmega],\[Delta], \[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]}]";
DionicKNGeoPlungeFunction::usage = "DionicKNGeoPlungeFunction[...] an object for storing the trajectory and orbital parameters in the assoc Association.";
Begin["`Private`"];


(* ::Section:: *)
(*DionicKN*)


(* ::Subsection:: *)
(*Real Roots Generic Plunge (Mino)*)


DionicKNGeoPlunge::r0outofboundsGen = "Intial radius `1` is not between the outer radial root `2`, and inner radial root `3`.";
DionicKNGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";
DionicKNGeoPlunge::TooMuchENergy = "the total value of \!\(\*SuperscriptBox[\(\[ScriptCapitalO]\), \(2\)]\)+\!\(\*SuperscriptBox[\(a\), \(2\)]\) = `1` is greater than one your are breakig cosmic censorship, stop this!"


Options[DionicKNGeoRealPlungeMino] = {"Roots"-> Automatic};
DionicKNGeoRealPlungeMino[a_, \[ScriptCapitalO]_ , \[CapitalOmega]_,\[Delta]_, \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalK]_, initCoords_, OptionsPattern[Options[DionicKNGeoRealPlungeMino]]] := Module[
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
		
	
	 
	rROOTS =r/.Solve[(\[Delta] r-m (a^2+r^2) \[ScriptCapitalE]+a m \[ScriptCapitalL])^2/m^2-(r^2+(-a \[ScriptCapitalE]+\[ScriptCapitalL])^2+\[ScriptCapitalK])(r-RP)(r-RM)==0,r];
	{r4,r3,r2,r1}= rROOTS;
	kr =Sqrt[((r1-r2) (r3-r4))/((r1-r3) (r2-r4))];
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
					{t0,r0,\[Theta]0,\[Phi]0}={0,r4+0.01,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0<=r4||r0>=r3,
		Message[DionicKNGeoPlunge::r0outofboundsGen,r0,r3,r4];
		Return[$Failed]
		];
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[DionicKNGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
	
	
	
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[DionicKNGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
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
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalK], 
		"ChargeParameters"->{\[ScriptCapitalO] , \[CapitalOmega], \[Delta]},
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {r1,r2,r3,r4},
		"PolarRoots"-> {z1,z2,z3,z3},
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(r\)]\)"-> kr,"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Inner\)], \(\)]\)"-> Abs[MinoR[R1]] - MinoR[r0](* Time outgoing to inner root*),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Outer\)], \(\)]\)"-> Abs[MinoR[R2]] - MinoR[r0](* Time outgoing to outer root*)|>,

		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	DionicKNGeoPlungeFunction[a,\[ScriptCapitalO], \[CapitalOmega],\[Delta], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalK], assoc]
]


(* ::Subsection:: *)
(*Complex Roots Generic Plunge (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


DionicKNGeoPlunge::r0outofboundsGen = "Intial radius `1` is not between the outer radial root `2`, and inner radial root `3`.";
DionicKNGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";


Options[DionicKNGeoComplexPlungeMino]= {"Roots"->Automatic};
DionicKNGeoComplexPlungeMino[a_ , \[ScriptCapitalO]_ , \[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalK]_, initCoords_, OptionsPattern[Options[DionicKNGeoComplexPlungeMino]]] := Module[{consts,MCPP, MCPM, MCMP , MCMM , assoc,MinoR,\[ScriptCapitalQ],m,z1,z2,z3,z4,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],t0,r0,\[Phi]0,M,D1M,D1P,D2M,D2P,Minoz,MinozFunc,MinoRFunc,\[Theta]0,z0,e,b,c,d,A,B,chi,kr,p2,f, t, r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots,kz,Z1, Z2, AMR,CNR,SNR,DNR,AMZ},
	
	M=1;
	m=1;
	J = 1-\[ScriptCapitalE]^2;
	\[ScriptCapitalQ] = \[ScriptCapitalK];
		
	RM = M-Sqrt[M^2-a^2-\[ScriptCapitalO]^2];
	RP = M+Sqrt[M^2-a^2-\[ScriptCapitalO]^2];
	consts = <|"a"->a,"\[ScriptCapitalO]"->\[ScriptCapitalO],"\[Delta]"->\[Delta], "\[CapitalOmega]"->\[CapitalOmega], "\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalK]"->\[ScriptCapitalK]|>;
	
		
	If[\[ScriptCapitalO]^2+a^2 >1,
		Message[DionicKNGeoPlunge::TooMuchENergy,a^2+\[ScriptCapitalO]^2];
		Return[$Failed]
		];
		
	
	ROOTS = r/.Solve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL] - \[Delta] r)^2-(r^2-2M*r+a^2+\[ScriptCapitalO]^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r];
	RealRoots = Select[ROOTS,Im[#]==0&];
	ComplexRoots = Select[ROOTS,Im[#]!=0&];
	R1= RealRoots[[1]];
	R2= RealRoots[[2]];
	R3=ComplexRoots[[1]];
	R4= ComplexRoots[[2]];
	
	If[a==0, {Z1,Z2}={0,\[ScriptCapitalL]}, {Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]}];
	
	e = R2;
	b = R1;
	c = Re[R3];
	d = Abs[Im[R4]];
	A = Sqrt[(e-c)^2+d^2];
	B  = Sqrt[(b-c)^2+d^2];
	chi = (A*b+e*B)/(A+B);
	kr = Sqrt[((e-b)^2-(A-B)^2)/(4*A*B)];
	p2 = b*A^2+e*B^2-(e+b)*A*B;
	f = (4 A B)/(A-B)^2;
	kz = a*Sqrt[J](Z1/Z2);
	
	\[CapitalUpsilon]r = 2 \[Pi] Sqrt[A B (1-\[ScriptCapitalE]^2) ]/(4 EllipticK[kr^2]);
	\[CapitalUpsilon]\[Theta] = \[Pi]/2 EllipticK[kz^2]/Z2;
	
	
	D1M=Sqrt[ -f];
	D2M=Sqrt[(4 A B (RM-b) (e-RM))]/(A( RM-b)+B(e- RM));

	D1P=Sqrt[ -f];
	D2P=Sqrt[(4 A B (RP-b) (e-RP))]/(A(RP- b)+B( e- RP));
	
	
	AMR[\[Lambda]_]:=JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2];
	SNR[\[Lambda]_]:=JacobiSN[Sqrt[A B J] \[Lambda],kr^2];
	CNR[\[Lambda]_]:=JacobiCN[Sqrt[A B J] \[Lambda],kr^2];
	DNR[\[Lambda]_]:=JacobiDN[Sqrt[A B J] \[Lambda],kr^2];

	AMZ[\[Lambda]_]:=JacobiAmplitude[Z2 \[Lambda],kz^2];
	
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,R1,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0<R1||r0>R2,
		Message[DionicKNGeoPlunge::r0outofboundsGen,r0,R2,R1];
		Return[$Failed]
		];
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[DionicKNGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
		
			
	z0 = Cos[\[Theta]0];
	
	Minoz[z_] :=InverseJacobiSN[z/Z1,kz^2]/Z2;
	If[Abs[Arg]==\[Pi]/2, Minoz[z_] :=0];
	If[a==0, Minoz[z_] :=0];
	If[\[ScriptCapitalQ]==0, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];


	MinoR[x_]:= -1/Sqrt[J A B] EllipticF[(\[Pi]/2 -ArcSin[(B(e-x)-A(x-b))/(   B(e-x)+A(x-b)  )  ]),kr^2];
	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r] -MinoR[r0] ],Listable];
	
	
	r[\[Lambda]_] := ((A-B) (A b-B e) SNR[\[Lambda]]^2+2 (A B (b+e)+A B (b-e) CNR[\[Lambda]]))/(4 A B+(A-B)^2 SNR[\[Lambda]]^2);

	z[\[Lambda]_]:= Z1*JacobiSN[Z2 \[Lambda],kz^2];
	If[kz==0, z[\[Lambda]_]:= Z1*Sin[Z2*\[Lambda]]];

(*Integrals*)
	
	RINT\[Lambda][\[Lambda]_] := ((A b-B e)/(A-B) \[Lambda]-1/ Sqrt[ J] ArcTan[(e-b)/(2 Sqrt[A B])  SNR[\[Lambda]]/Sqrt[1-kr^2 (SNR[\[Lambda]])^2]]+((A+B) (e-b))/(2 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]);
	R2INT\[Lambda][\[Lambda]_]:=\[Lambda] /(A-B) (A b^2-B e^2)+ Sqrt[A B ]/Sqrt[ J] (EllipticE[AMR[\[Lambda]],kr^2])-((A+B) (A^2+2 b^2-B^2-2 e^2))/(4 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]-(Sqrt[A B ] (A+B-(A-B)CNR[\[Lambda]]))/((A-B) Sqrt[ J]) (SNR[\[Lambda]] DNR[\[Lambda]])/(f+(SNR[\[Lambda]])^2)+ (A^2+2 b^2-B^2-2 e^2)/(4 (e-b) Sqrt[ J])ArcTan[(f-(1+2 f kr^2) SNR[\[Lambda]]^2),2 SNR[\[Lambda]] DNR[\[Lambda]]Sqrt[f (1+f kr^2)]](*ArcTan[(2 SNR[\[Lambda]] DNR[\[Lambda]]Sqrt[f (1+f kr^2)])/(f-(1+2 f kr^2) SNR[\[Lambda]]^2)]*);
	If[a!=0,RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))+((e-b) (A (b-RM)+B (e-RM)))/(2 Sqrt[A B J] (b-RM) (-e+RM) (A (b-RM)-B (e-RM))) EllipticPi[1/D2M^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]-  Sqrt[(e-b)]/(Sqrt[ J] Sqrt[ (RM-b) (e-RM)] Sqrt[ (A^2 (RM-b)-(e-RM) (b^2-B^2+e RM-b (e+RM)))]) 1/4 (Log[((D2M  Sqrt[1-D2M^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)/((D2M  Sqrt[1-D2M^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)])];
	If[a==0,RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))];
	RPINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RP)+B (-e+RP))+((e-b) (A (b-RP)+B (e-RP)))/(2 Sqrt[A B J] (b-RP) (-e+RP) (A (b-RP)+B (-e+RP))) EllipticPi[1/D2P^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]- Sqrt[(e-b)]/(Sqrt[ J] Sqrt[(RP-b) (e-RP)] Sqrt[ (A^2 (RP-b)-(e-RP) (b^2-B^2+e RP-b (e+RP)))]) 1/4 (Log[((D2P  Sqrt[1-D2P^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)/((D2P  Sqrt[1-D2P^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)]);


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

	MCPP = -Abs[MinoR[RP]] -MinoR[r0]+(2\[Pi])/\[CapitalUpsilon]r;
	MCMP = -Abs[MinoR[RM]] -MinoR[r0]+(2\[Pi])/\[CapitalUpsilon]r;
	MCMM = Abs[MinoR[RM]] - MinoR[r0];
	MCPM = Abs[MinoR[RP]] - MinoR[r0];
	assoc = Association[
		"a" -> a,
		"Parametrization"->"Mino", 
		"ChargeParameters"->{\[ScriptCapitalO] , \[CapitalOmega], \[Delta]},
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalK], 
		"ConstantsOfMotion" -> consts,
		"Periods"-> <|"\!\(\*SubscriptBox[\(\[Lambda]\), \(r\)]\)"-> 2\[Pi]/\[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[Lambda]\), \(\[Theta]\)]\)" -> 2\[Pi]/\[CapitalUpsilon]\[Theta]|>,
		"RadialRoots"-> {R1,R2,R3,R4},
		"PolarRoots"-> {z1,z2,z3,z3},
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(r\)]\)"-> kr,"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Inner\)], \(\)]\)"-> Abs[MinoR[R1]] - MinoR[r0](* Time outgoing to inner root*),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Outer\)], \(\)]\)"-> Abs[MinoR[R2]] - MinoR[r0](* Time outgoing to outer root*)|>,

		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	DionicKNGeoPlungeFunction[a,\[ScriptCapitalO], \[CapitalOmega],\[Delta], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalK], assoc]
]


(* ::Subsection::Closed:: *)
(*Real Roots Generic Bound (Mino)*)


DionicKNGeoPlunge::r0outofboundsGen = "Intial radius `1` is not between the outer radial root `2`, and inner radial root `3`.";
DionicKNGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";
DionicKNGeoPlunge::TooMuchENergy = "the total value of \!\(\*SuperscriptBox[\(\[ScriptCapitalO]\), \(2\)]\)+\!\(\*SuperscriptBox[\(a\), \(2\)]\) = `1` is greater than 1 your are breaking cosmic censorship, stop this!"


(*Options[DionicKNGeoRealBoundMino] = {"Roots"-> Automatic};
DionicKNGeoRealBoundMino[a_, \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalK]_, \[ScriptCapitalO]_ , \[CapitalOmega]_, initCoords_, OptionsPattern[Options[DionicKNGeoRealBoundMino]]] := Module[
	{consts,assoc,M,Mino,t0,r0,\[Phi]0,kr,hR,hP,hM,\[CapitalUpsilon]t,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]tr,\[CapitalUpsilon]tz,\[CapitalUpsilon]z,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,t,r,z,z0,\[Theta]0, \[Theta], J,\[Phi], R1,R2,R3,R4,MCMP,MCPP,MCPM,MCMM, Minoz, MinozFunc, MinoRFunc,m , RM, RP,ROOTS,RealRoots,ComplexRoots, kz, Z1, Z2, MinoR},
	
	M=1;
	m=1;
	J = 1-\[ScriptCapitalE]^2;
	\[Delta] = 0;
	\[ScriptCapitalQ] = \[ScriptCapitalK];
		
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	consts = <|"\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalK]"->\[ScriptCapitalK]|>;
	
	If[\[ScriptCapitalO]^2+a^2 >1,
		Message[DionicKNGeoPlunge::TooMuchENergy,\[ScriptCapitalO]^2+a^2];
		Return[$Failed]
		];
		
	
	 
	rROOTS =r/.Solve[(\[Delta]  r-m (a^2+r^2) \[ScriptCapitalE]+a m \[ScriptCapitalL])^2/m^2-(r^2+(-a \[ScriptCapitalE]+\[ScriptCapitalL])^2+\[ScriptCapitalK])(r-RP)(r-RM)==0,r];
	{r4,r3,r2,r1}= rROOTS;
	kr =Sqrt[((r1-r2) (r3-r4))/((r1-r3) (r2-r4))];
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



	
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,R4,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0<R4||r0>R3,
		Message[DionicKNGeoPlunge::r0outofboundsGen,r0,R3,R4];
		Return[$Failed]
		];
		
		
	z0 = Cos[\[Theta]0];
	
	
	
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[DionicKNGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
		
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r]- MinoR[r0]],Listable];
	
	If[Abs[Arg]==\[Pi]/2, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	
	If[\[ScriptCapitalQ]==0, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	
	
	

	
	

	\[Phi]r[\[Lambda]_] := RPINT\[Lambda][\[Lambda]]((a (-a^2 m \[ScriptCapitalE]+a m \[ScriptCapitalL]+RP (-m RP \[ScriptCapitalE]+\[Delta])))/(m (RM-RP)) ) + RMINT\[Lambda][\[Lambda]]((a (a^2 m \[ScriptCapitalE]-a m \[ScriptCapitalL]+RM (m RM \[ScriptCapitalE]-\[Delta])))/(m (RM-RP))) +a \[ScriptCapitalE] \[Lambda];
	\[Phi]z[\[Lambda]_]:= ZP1INT\[Lambda][\[Lambda]]((m \[ScriptCapitalL]-\[CapitalOmega])/(2 m)) -ZM1INT\[Lambda][\[Lambda]]( (m \[ScriptCapitalL]+\[CapitalOmega])/(2 m)) -a \[ScriptCapitalE] \[Lambda];
	


	
	r[\[Lambda]_]:=(R1(R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-R4(R3-R1))/((R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-(R3-R1));
	z[\[Lambda]_]:=(z3 (-z2+z4)+(z2-z3) z4 Sin[AMZ[\[Lambda]]]^2)/(-z2+z4+(z2-z3) Sin[AMZ[\[Lambda]]]^2);
	
	
	tr[\[Lambda]_] := (2 a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE]+RM RP \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]-a \[ScriptCapitalL]-(RM \[Delta])/m-(RP \[Delta])/m)\[Lambda]+RINT\[Lambda][\[Lambda]] ((RM+RP) \[ScriptCapitalE]-\[Delta]/m)+R2INT\[Lambda][\[Lambda]] \[ScriptCapitalE]-((a^2+RP^2) (a^2 m \[ScriptCapitalE]-a m \[ScriptCapitalL]+RP (m RP \[ScriptCapitalE]-\[Delta])))/(m (RM-RP)) RPINT\[Lambda][\[Lambda]]+((a^2+RM^2) (a^2 m \[ScriptCapitalE]-a m \[ScriptCapitalL]+RM (m RM \[ScriptCapitalE]-\[Delta])))/(m (RM-RP)) RMINT\[Lambda][\[Lambda]] ;
	tz[\[Lambda]_]:= + a^2 Z2INT\[Lambda][\[Lambda]] \[ScriptCapitalE]+(a \[CapitalOmega])/m ZINT\[Lambda][\[Lambda]]-a (a  \[ScriptCapitalE]- \[ScriptCapitalL])\[Lambda];

	t=Function[{Global`\[Lambda]}, Evaluate[ tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ Minoz[z0]]-tr[MinoR[r0]]-tz[Minoz[z0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda] + MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + Minoz[z0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[\[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ Minoz[z0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[Minoz[z0]] + \[Phi]0], Listable];
	

	
	assoc = Association[
		"a" -> a,
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency"-> \[CapitalUpsilon]r,
		"PolarFrequency"-> \[CapitalUpsilon]z,
		"Frequencies"-> <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"->\[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]z |>,
		"RadialPeriod"-> 2\[Pi]/\[CapitalUpsilon]r,
		"PolarPeriod"-> 2\[Pi]/\[CapitalUpsilon]z,
		"Periods"-> <|"\!\(\*SubscriptBox[\(\[Lambda]\), \(r\)]\)"-> 2\[Pi]/\[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[Lambda]\), \(\[Theta]\)]\)" -> 2\[Pi]/\[CapitalUpsilon]z|>,
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(r\)]\)"-> kr,"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialRoots"-> {R4,R3,R2,R1},
		"RadialMinoTime"-> MinoRFunc,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Inner\)], \(\)]\)"-> Abs[MinoR[R4]] - MinoR[r0](* Time outgoing to inner root*),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Outer\)], \(\)]\)"-> Abs[MinoR[R3]] - MinoR[r0](* Time outgoing to outer root*)|>,
		"PolarRoots"-> {z1,z2,z3,z4},
		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	DionicKNGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]
*)


(* ::Section:: *)
(*DionicKNGeoOrbit and DionicKNGeoOrbitFuction*)


DionicKNGeoPlunge::fourcomplexroots = "This function can only currently solve plunging geodesics in the paramater space where we have four finite real roots or two finite real roots and two finite complex roots.";
DionicKNGeoPlunge::highenergy = "This code does not currently support plunging solutions with \[ScriptCapitalE]>1, \[ScriptCapitalE]=`1` has been selected ";
DionicKNGeoPlunge::negativecarter = "This code does not currently support plunging solutions with Q<0, Q=`1` has been selected ";
DionicKNGeoPlunge::rIoutbounds = "Given ISSO radius `1` is not in the allowed range between `2` and `3` for spin `4`.";
DionicKNGeoPlunge::Inclinationoutofbounds = "Given inclination angle `1` is not between -\!\(\*FractionBox[\(\[Pi]\), \(2\)]\) and \!\(\*FractionBox[\(\[Pi]\), \(2\)]\)."
DionicKNGeoPlunge::aoutofbounds1 = "The plunges packages does not currently support the extremal a=1 case"
DionicKNGeoPlunge::aoutofbounds0 = "The plunges packages does not currently support the zero spin a=0 case"


DionicKNGeoPlunge[a_:0.7 ,\[ScriptCapitalO]_:0.3, \[CapitalOmega]_:3, \[Delta]_:3,\[ScriptCapitalE]_:0.8,\[ScriptCapitalL]_:0.3,\[ScriptCapitalK]_:3, initCoords:{_,_,_,_}:Automatic,OptionsPattern[]]:=Module[
	{param, method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
		
	If[a==1,
		Message[DionicKNGeoPlunge::aoutofbounds1];
		Return[$Failed]
		];
	
	If[\[ScriptCapitalE]>=1,
		Message[DionicKNGeoPlunge::highenergy,\[ScriptCapitalE]];
		Return[$Failed]
		];
		
	If[\[ScriptCapitalK]<0,
		Message[DionicKNGeoPlunge::negativecarter,\[ScriptCapitalK]];
		Return[$Failed]
		];

	ROOTS = r/.NSolve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL]-\[Delta] r)^2-(r^2-2r+a^2+\[ScriptCapitalO]^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalK])==0,r];

	RealRoots = Select[ROOTS,PossibleZeroQ@Im[#]&];
	ComplexRoots = Select[ROOTS,Not@PossibleZeroQ@Im[#]&];

	If[Length[ComplexRoots]== 2, Return[DionicKNGeoComplexPlungeMino[a,\[ScriptCapitalO], \[CapitalOmega] ,\[Delta], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalK], initCoords,"Roots"->ROOTS]]];
	If[Length[ComplexRoots]== 0, Return[DionicKNGeoRealPlungeMino[a,\[ScriptCapitalO], \[CapitalOmega] ,\[Delta], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalK], initCoords,"Roots"->ROOTS]]];
	];


DionicKNGeoPlungeFunction /:
 MakeBoxes[kgof:DionicKNGeoPlungeFunction[a_,\[ScriptCapitalO]_,\[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalK]_,assoc_], form:(StandardForm|TraditionalForm)] :=
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
    DionicKNGeoPlungeFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];


DionicKNGeoPlungeFunction[a_,\[ScriptCapitalO]_,\[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalK]_,assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
DionicKNGeoPlungeFunction[a_,\[ScriptCapitalO]_,\[CapitalOmega]_,\[Delta]_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalK]_,assoc_][y_?StringQ] := assoc[y]
Keys[g_DionicKNGeoPlungeFunction]^:=Keys[g[[5]]]


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
