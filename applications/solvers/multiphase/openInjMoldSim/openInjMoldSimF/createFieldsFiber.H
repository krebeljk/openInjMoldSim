    Info<< "Reading fiberProperties\n" << endl;
    IOdictionary fiberProperties
    (
        IOobject
        (
            "fiberProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Switch FiberSuspension(
	  fiberProperties.lookup("FiberSuspension")
    );
    
    if(FiberSuspension)
    {
      Info << "Fiber suspension model on " << endl;
    }
      
      
	scalar mu1(
		readScalar(fiberProperties.lookup("mu1"))
	); 
	 
	scalar mu2(
		readScalar(fiberProperties.lookup("mu2"))
	); 

	scalar r(
		readScalar(fiberProperties.lookup("r"))
	); 

	scalar lambda = (r*r-1.0)/(r*r+1.0);
	
	scalar FiberVolfrac(
		readScalar(fiberProperties.lookup("FiberVolfrac"))
	); 
	
	word FiberClosure(fiberProperties.lookup("FiberClosure"));	
	if (FiberClosure=="hybrid"){
		Info<< "Fiber Closure Model: hybrid " << endl;
	}
	if (FiberClosure=="IBOF"){
		Info<< "Fiber Closure Model: IBOF " << endl;
	}
	
	tensor ut (1,0,0,0,1,0,0,0,1);
	tensor zerotensor (0,0,0,0,0,0,0,0,0);
  
	volTensorField 	UnitTensor 
    	(
        IOobject
        (
            "UnitTensor",
            runTime.timeName(),
            mesh,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        	mesh,
		ut
    	);


//VelocityGradient Tensor
	volTensorField L
    (
        IOobject
        (
            "L",
            runTime.timeName(),
            mesh,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	fvc::grad(U) //normal convention L=grad(U)^T
    );

//Rate of Deformation Tensor
	volTensorField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0/2.0*( L + L.T())
    );

//vorticity tensor
	volTensorField W
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	1.0/2.0*( L.T() - L) //ACHTUNG: DEFINITION OF L andersrum
    );

volSymmTensorField StressFiber
	(
		IOobject
		(
		"StressFiber",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
		),
        mesh,
        dimensionedSymmTensor("zero", dimensionSet(1,-1,-2,0,0,0,0), symmTensor::zero)
	);

	volSymmTensorField A2 
    	(
        IOobject
        (
            "A2",
            runTime.timeName(),
            mesh,
	    IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
		mesh
    	);

	volScalarField fHybrid
	( 
        IOobject
        (
            "fHybrid",
            runTime.timeName(),
            mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
		1.0-27.0*det(A2)
    	);

//------------------------------------------------------------------------------------
	volTensorField WA2 = -(W & A2);
	volTensorField A2W = ( A2& W); 
	volTensorField DA2 = lambda*(D &  A2);
	volTensorField A2D = lambda*(A2 & D) ;
	volScalarField magnitudeD = pow(0.5*(D && D.T()),0.5);

	volScalarField dyadA2D = A2 && D;
	volScalarField dyadUnitTensorD = UnitTensor && D;

	volTensorField B = A2 & A2;
	volScalarField trD = tr(D);
	volScalarField dyadUnitD = UnitTensor && D;
	volScalarField dyadBD = B && D;
	volTensorField DB = D & B;
	volTensorField BD = B & D;
//------------------------------------------------------------------------------------
//define additional field variables for hybrid closure model
 	volSymmTensorField ClosureQuadrat = dyadA2D  * A2;
	
	volTensorField ClosureLinear = 	 -1.0/35.0*( dyadUnitTensorD * UnitTensor + 2.0*D  ) 
					+ 1.0/7.0*( dyadUnitTensorD * A2 + 2.0*(A2 & D) 
					+ dyadA2D * UnitTensor + 2.0*(D & A2));
	volTensorField ClosureA4 = fHybrid*ClosureQuadrat + (1.0-fHybrid)*ClosureLinear;
	
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

	volTensorField SourceInteract = 2.0*CI*magnitudeD*(UnitTensor - 3.0*A2);
	volTensorField SourceA2 = - WA2 - A2W  + DA2 + A2D - 2.0*lambda * ClosureA4 + SourceInteract;
//------------------------------------------------------------------------------------
	/*
	volTensorField SourceA2 
	(
        IOobject
        (
            "SourceA2",
            runTime.timeName(),
            mesh,
	    IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh,
	dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0), 0.0) 
    );
*/
//Fiber Closure Equations IBOF
//Chung and Kwon 2002, Journal of Rheology

//Calculate A2 Invariants 
	
volScalarField A2II = 0.5 - 0.5*(A2.component(tensor::XX)*A2.component(tensor::XX) + 2.0*A2.component(tensor::XY)*A2.component(tensor::XY) + 2.0*A2.component(tensor::XZ)*A2.component(tensor::XZ)+
2.0*A2.component(tensor::YZ)*A2.component(tensor::YZ) + A2.component(tensor::YY)*A2.component(tensor::YY) 
+ A2.component(tensor::ZZ)*A2.component(tensor::ZZ) );

volScalarField A2III = A2.component(tensor::XX)*
(A2.component(tensor::YY)*A2.component(tensor::ZZ) - 
A2.component(tensor::YZ)*A2.component(tensor::YZ)) + A2.component(tensor::XY)*(A2.component(tensor::YZ)*A2.component(tensor::XZ)-
A2.component(tensor::XY)*A2.component(tensor::ZZ))+
A2.component(tensor::XZ)*(A2.component(tensor::XY)*A2.component(tensor::YZ)-
A2.component(tensor::YY)*A2.component(tensor::XZ));

List<double> c3(21);
List<double> c4(21);
List<double> c6(21);

c3[0] = 0.24940908165786e2;
c3[1] = -0.435101153160329e3;
c3[2] = 0.372389335663877e4;
c3[3] = 0.703443657916476e4;
c3[4] = 0.823995187366106e6;
c3[5] = -0.133931929894245e6;
c3[6] = 0.880683515327916e6;
c3[7] = -0.991630690741981e7;
c3[8] = -0.159392396237307e5;
c3[9] = 0.800970026849796e7;
c3[10] = -0.237010458689252e7;
c3[11] = 0.379010599355267e8;
c3[12] = -0.337010820273821e8;
c3[13] = 0.322219416256417e5;
c3[14] = -0.257258805870567e9;
c3[15] = 0.214419090344474e7;
c3[16] = -0.449275591851490e8;
c3[17] = -0.213133920223355e8;
c3[18] = 0.157076702372204e10;
c3[19] = -0.232153488525298e5;
c3[20] = -0.395769398304473e10;

c4[0] = -0.497217790110754e0;
c4[1] = 0.234980797511405e2;
c4[2] = -0.391044251397838e3;
c4[3] = 0.153965820593506e3;
c4[4] = 0.152772950743819e6;
c4[5] = -0.213755248785646e4;
c4[6] = -0.400138947092812e4;
c4[7] = -0.185949305922308e7;
c4[8] = 0.296004865275814e4;
c4[9] = 0.247717810054366e7;
c4[10] =0.101013983339062e6;
c4[11] =0.732341494213578e7;
c4[12] =-0.147919027644202e8;
c4[13] =-0.104092072189767e5;
c4[14] = -0.635149929624336e8;
c4[15] = -0.247435106210237e6;
c4[16] = -0.902980378929272e7;
c4[17] = 0.724969796807399e7;
c4[18] = 0.487093452892595e9;
c4[19] = 0.138088690964946e5;
c4[20] = -0.160162178614234e10;

c6[0] = 0.234146291570999e2;
c6[1] = -0.412048043372534e3;
c6[2] = 0.319553200392089e4;
c6[3] = 0.573259594331015e4;
c6[4] = -0.485212803064813e5;
c6[5] = -0.605006113515592e5;
c6[6] = -0.477173740017567e5;
c6[7] = 0.599066486689836e7;
c6[8] = -0.110656935176569e5;
c6[9] = -0.460543580680696e8;
c6[10] =0.203042960322874e7;
c6[11] =-0.556606156734835e8;
c6[12] =0.567424911007837e9;
c6[13] =0.128967058686204e5;
c6[14] =-0.152752854956514e10;
c6[15] =-0.499321746092534e7;
c6[16] =0.132124828143333e9;
c6[17] =-0.162359994620983e10;
c6[18] =0.792526849882218e10;
c6[19] =0.466767581292985e4;
c6[20] = -0.128050778279459e11;

volScalarField beta3 = c3[0] + c3[1]*A2II + c3[2]*A2II*A2II + c3[3]*A2III
+ c3[4]*A2III*A2III + c3[5]*A2II*A2III + c3[6]* A2II*A2II*A2III+
c3[7]*A2II*A2III*A2III + c3[8]*A2II*A2II*A2II + c3[9]*A2III*A2III*A2III+
c3[10]*A2II*A2II*A2II*A2III + c3[11]*A2II*A2II*A2III*A2III+
c3[12]*A2II*A2III*A2III*A2III +c3[13]*A2II*A2II*A2II*A2II+
c3[14]*A2III*A2III*A2III*A2III + c3[15]*A2II*A2II*A2II*A2II*A2III+
c3[16]*A2II*A2II*A2II*A2III*A2III + c3[17]*A2II*A2II*A2III*A2III*A2III+
c3[18]*A2II*A2III*A2III*A2III*A2III+c3[19]*A2II*A2II*A2II*A2II*A2II+
c3[20]*A2III*A2III*A2III*A2III*A2III;

volScalarField beta4 = c4[0] + c4[1]*A2II + c4[2]*A2II*A2II + c4[3]*A2III
+ c4[4]*A2III*A2III + c4[5]*A2II*A2III + c4[6]* A2II*A2II*A2III+
c4[7]*A2II*A2III*A2III + c4[8]*A2II*A2II*A2II + c4[9]*A2III*A2III*A2III+
c4[10]*A2II*A2II*A2II*A2III + c4[11]*A2II*A2II*A2III*A2III+
c4[12]*A2II*A2III*A2III*A2III +c4[13]*A2II*A2II*A2II*A2II+
c4[14]*A2III*A2III*A2III*A2III + c4[15]*A2II*A2II*A2II*A2II*A2III+
c4[16]*A2II*A2II*A2II*A2III*A2III + c4[17]*A2II*A2II*A2III*A2III*A2III+
c4[18]*A2II*A2III*A2III*A2III*A2III+c4[19]*A2II*A2II*A2II*A2II*A2II+
c4[20]*A2III*A2III*A2III*A2III*A2III;

volScalarField beta6 = c6[0] + c6[1]*A2II + c6[2]*A2II*A2II + c6[3]*A2III
+ c6[4]*A2III*A2III + c6[5]*A2II*A2III + c6[6]*A2II*A2II*A2III+
c6[7]*A2II*A2III*A2III + c6[8]*A2II*A2II*A2II + c6[9]*A2III*A2III*A2III+
c6[10]*A2II*A2II*A2II*A2III + c6[11]*A2II*A2II*A2III*A2III+
c6[12]*A2II*A2III*A2III*A2III + c6[13]*A2II*A2II*A2II*A2II+
c6[14]*A2III*A2III*A2III*A2III + c6[15]*A2II*A2II*A2II*A2II*A2III+
c6[16]*A2II*A2II*A2II*A2III*A2III + c6[17]*A2II*A2II*A2III*A2III*A2III+
c6[18]*A2II*A2III*A2III*A2III*A2III + c6[19]*A2II*A2II*A2II*A2II*A2II+
c6[20]*A2III*A2III*A2III*A2III*A2III;


volScalarField beta1 = 3.0/5.0*(-1.0/7.0+1.0/5.0*beta3*(1.0/7.0+4.0/7.0*A2II+8.0/3.0*A2III) - beta4*(1.0/5.0-8.0/15.0*A2II - 14.0/15.0*A2III) - beta6*(1.0/35.0-24.0/105.0*A2III-4.0/35.0*A2II + 16.0/15.0*A2II*A2III + 8.0/35.0*A2II*A2II));

volScalarField beta2 = 6.0/7.0*(1.0-1.0/5.0*beta3*(1.0+4.0*A2II)+7.0/5.0*beta4*(1.0/6.0-A2II)-beta6*(-1.0/5.0 + 2.0/3.0*A2III + 4.0/5.0*A2II - 8.0/5.0*A2II*A2II));

volScalarField beta5 = -4.0/5.0*beta3 - 7.0/5.0*beta4 - 6.0/5.0*beta6*(1.0-4.0/3.0*A2II);

  ClosureA4 = beta1/3.0*( trD*UnitTensor + 2.0*D)+
beta2/6.0*( UnitTensor * (dyadA2D) + A2 * (dyadUnitD) + 4.0*A2D )
+ beta3/3.0*( A2 *dyadA2D + 2.0*(A2& DA2) )+
beta4/6.0*(UnitTensor*dyadBD + B*(dyadUnitD) + 2.0*(B & D.T()) + 2.0* BD)
+ beta5/6.0*(A2*dyadBD + B*dyadA2D + 4.0*(A2 & DB)) 
+ beta6/3.0*(B*dyadBD + 2.0*(B & DB) );

volScalarField ClosureA4_imp = 2.0*lambda*(beta2/6.0*dyadUnitD + beta3/3.0*dyadA2D + beta5/6.0*dyadBD);

volTensorField ClosureA4_exp = beta1/3.0*( trD*UnitTensor + 2.0*D)+
		beta2/6.0*( UnitTensor * (dyadA2D) +  4.0*A2D )
		+ beta3/3.0*( 2.0*(A2& DA2) )+
		beta4/6.0*(UnitTensor*dyadBD + B*(dyadUnitD) + 2.0*(B & D.T()) + 2.0* BD)
		+ beta5/6.0*(B*dyadA2D + 4.0*(A2 & DB)) 
		+ beta6/3.0*(B*dyadBD + 2.0*(B & DB) );

volTensorField SourceInteract_exp = 2.0*CI*magnitudeD*UnitTensor;

SourceA2 = - WA2 - A2W  + lambda*DA2 + lambda*A2D - 2.0*lambda * ClosureA4_exp + SourceInteract_exp; //Vorzeichen von W andersrum als Advani definiert!!
	  
volScalarField Sp_imp1 = 6.0*CI*(pow(0.5*(D && D.T()),0.5));
