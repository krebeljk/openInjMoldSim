 //update all tensor calculations 
	L = fvc::grad(U);
	W = 1.0/2.0*( L.T() - L);
	D = 1.0/2.0*( L + L.T());

	WA2 = -(W & A2);
	A2W = ( A2& W); 
	DA2 = lambda*(D &  A2);
	A2D = lambda*(A2 & D) ;
	magnitudeD = pow(0.5*(D && D.T()),0.5);

	dyadA2D = A2 && D;
	dyadUnitTensorD = UnitTensor && D;

	B = A2 & A2;
	trD = tr(D);
	dyadUnitD = UnitTensor && D;
	dyadBD = B && D;
	DB = D & B;
	BD = B & D;

	
//Fiber Interaction coefficient
//Bay
scalar CI;
if ( FiberVolfrac*r > 1.3){
	CI = 0.0184*Foam::exp(-0.7148* FiberVolfrac*r);
}
else{
//Phan-Thien(2002)
	CI = 0.03*(1.0-Foam::exp(-0.224* FiberVolfrac*r));
}
//-------------------------------------------------------
//hybrid closure for A4 tensor
//Fiber Closure Equations
if(FiberClosure=="hybrid"){

	ClosureQuadrat = dyadA2D  * A2;
	ClosureLinear = 	 -1.0/35.0*( dyadUnitTensorD * UnitTensor + 2.0*D  ) 
					+ 1.0/7.0*( dyadUnitTensorD * A2 + 2.0*(A2 & D) 
					+ dyadA2D * UnitTensor + 2.0*(D & A2));
					
	fHybrid  = 1.0-27.0*det(A2);				

	ClosureA4 = fHybrid*ClosureQuadrat + (1.0-fHybrid)*ClosureLinear;

	SourceInteract = 2.0*CI*magnitudeD*(UnitTensor - 3.0*A2);

	SourceA2 = - WA2 - A2W  + DA2 + A2D - 2.0*lambda * ClosureA4 + SourceInteract;
}

//-------------------------------------------------------
//IBOF-5 closure of Chung and Kwon
if(FiberClosure=="IBOF"){
//Calculate A2 Invariants 
A2II = 0.5 - 0.5*(A2.component(tensor::XX)*A2.component(tensor::XX) + 2.0*A2.component(tensor::XY)*A2.component(tensor::XY) + 2.0*A2.component(tensor::XZ)*A2.component(tensor::XZ)+
2.0*A2.component(tensor::YZ)*A2.component(tensor::YZ) + A2.component(tensor::YY)*A2.component(tensor::YY) 
+ A2.component(tensor::ZZ)*A2.component(tensor::ZZ) );

A2III = A2.component(tensor::XX)*
(A2.component(tensor::YY)*A2.component(tensor::ZZ) - 
A2.component(tensor::YZ)*A2.component(tensor::YZ)) + A2.component(tensor::XY)*(A2.component(tensor::YZ)*A2.component(tensor::XZ)-
A2.component(tensor::XY)*A2.component(tensor::ZZ))+
A2.component(tensor::XZ)*(A2.component(tensor::XY)*A2.component(tensor::YZ)-
A2.component(tensor::YY)*A2.component(tensor::XZ));

beta3 = c3[0] + c3[1]*A2II + c3[2]*A2II*A2II + c3[3]*A2III
+ c3[4]*A2III*A2III + c3[5]*A2II*A2III + c3[6]* A2II*A2II*A2III+
c3[7]*A2II*A2III*A2III + c3[8]*A2II*A2II*A2II + c3[9]*A2III*A2III*A2III+
c3[10]*A2II*A2II*A2II*A2III + c3[11]*A2II*A2II*A2III*A2III+
c3[12]*A2II*A2III*A2III*A2III +c3[13]*A2II*A2II*A2II*A2II+
c3[14]*A2III*A2III*A2III*A2III + c3[15]*A2II*A2II*A2II*A2II*A2III+
c3[16]*A2II*A2II*A2II*A2III*A2III + c3[17]*A2II*A2II*A2III*A2III*A2III+
c3[18]*A2II*A2III*A2III*A2III*A2III+c3[19]*A2II*A2II*A2II*A2II*A2II+
c3[20]*A2III*A2III*A2III*A2III*A2III;

beta4 = c4[0] + c4[1]*A2II + c4[2]*A2II*A2II + c4[3]*A2III
+ c4[4]*A2III*A2III + c4[5]*A2II*A2III + c4[6]* A2II*A2II*A2III+
c4[7]*A2II*A2III*A2III + c4[8]*A2II*A2II*A2II + c4[9]*A2III*A2III*A2III+
c4[10]*A2II*A2II*A2II*A2III + c4[11]*A2II*A2II*A2III*A2III+
c4[12]*A2II*A2III*A2III*A2III +c4[13]*A2II*A2II*A2II*A2II+
c4[14]*A2III*A2III*A2III*A2III + c4[15]*A2II*A2II*A2II*A2II*A2III+
c4[16]*A2II*A2II*A2II*A2III*A2III + c4[17]*A2II*A2II*A2III*A2III*A2III+
c4[18]*A2II*A2III*A2III*A2III*A2III+c4[19]*A2II*A2II*A2II*A2II*A2II+
c4[20]*A2III*A2III*A2III*A2III*A2III;

beta6 = c6[0] + c6[1]*A2II + c6[2]*A2II*A2II + c6[3]*A2III
+ c6[4]*A2III*A2III + c6[5]*A2II*A2III + c6[6]*A2II*A2II*A2III+
c6[7]*A2II*A2III*A2III + c6[8]*A2II*A2II*A2II + c6[9]*A2III*A2III*A2III+
c6[10]*A2II*A2II*A2II*A2III + c6[11]*A2II*A2II*A2III*A2III+
c6[12]*A2II*A2III*A2III*A2III + c6[13]*A2II*A2II*A2II*A2II+
c6[14]*A2III*A2III*A2III*A2III + c6[15]*A2II*A2II*A2II*A2II*A2III+
c6[16]*A2II*A2II*A2II*A2III*A2III + c6[17]*A2II*A2II*A2III*A2III*A2III+
c6[18]*A2II*A2III*A2III*A2III*A2III + c6[19]*A2II*A2II*A2II*A2II*A2II+
c6[20]*A2III*A2III*A2III*A2III*A2III;


beta1 = 3.0/5.0*(-1.0/7.0+1.0/5.0*beta3*(1.0/7.0+4.0/7.0*A2II+8.0/3.0*A2III) - beta4*(1.0/5.0-8.0/15.0*A2II - 14.0/15.0*A2III) - beta6*(1.0/35.0-24.0/105.0*A2III-4.0/35.0*A2II + 16.0/15.0*A2II*A2III + 8.0/35.0*A2II*A2II));
beta2 = 6.0/7.0*(1.0-1.0/5.0*beta3*(1.0+4.0*A2II)+7.0/5.0*beta4*(1.0/6.0-A2II)-beta6*(-1.0/5.0 + 2.0/3.0*A2III + 4.0/5.0*A2II - 8.0/5.0*A2II*A2II));
beta5 = -4.0/5.0*beta3 - 7.0/5.0*beta4 - 6.0/5.0*beta6*(1.0-4.0/3.0*A2II);

  ClosureA4 = beta1/3.0*( trD*UnitTensor + 2.0*D)+
beta2/6.0*( UnitTensor * (dyadA2D) + A2 * (dyadUnitD) + 4.0*A2D )
+ beta3/3.0*( A2 *dyadA2D + 2.0*(A2& DA2) )+
beta4/6.0*(UnitTensor*dyadBD + B*(dyadUnitD) + 2.0*(B & D.T()) + 2.0* BD)
+ beta5/6.0*(A2*dyadBD + B*dyadA2D + 4.0*(A2 & DB)) 
+ beta6/3.0*(B*dyadBD + 2.0*(B & DB) );

ClosureA4_imp = 2.0*lambda*(beta2/6.0*dyadUnitD + beta3/3.0*dyadA2D + beta5/6.0*dyadBD);

ClosureA4_exp = beta1/3.0*( trD*UnitTensor + 2.0*D)+
		beta2/6.0*( UnitTensor * (dyadA2D) +  4.0*A2D )
		+ beta3/3.0*( 2.0*(A2& DA2) )+
		beta4/6.0*(UnitTensor*dyadBD + B*(dyadUnitD) + 2.0*(B & D.T()) + 2.0* BD)
		+ beta5/6.0*(B*dyadA2D + 4.0*(A2 & DB)) 
		+ beta6/3.0*(B*dyadBD + 2.0*(B & DB) );

SourceInteract_exp = 2.0*CI*magnitudeD*UnitTensor;

SourceA2 = - WA2 - A2W  + lambda*DA2 + lambda*A2D - 2.0*lambda * ClosureA4_exp + SourceInteract_exp; //Vorzeichen von W andersrum als Advani definiert!!
}
