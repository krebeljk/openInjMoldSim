//A2 equation
Sp_imp1 = 6.0*CI*(pow(0.5*(D && D.T()),0.5));

if(FiberClosure=="IBOF"){			
        fvSymmTensorMatrix A2Eqn
	(
			fvm::ddt(A2) + fvm::div(phi, A2) + fvm::SuSp(-fvc::div(phi),A2) +
			fvm::Sp(Sp_imp1,A2) + fvm::Sp(ClosureA4_imp,A2)
	);
	solve(A2Eqn == symm(SourceA2) * pos(alpha1-0.9));
}
else{
	fvSymmTensorMatrix A2Eqn
			(
			fvm::ddt(A2) + fvm::div(phi, A2) + fvm::SuSp(-fvc::div(phi),A2));
	solve(A2Eqn == symm(SourceA2) * pos(alpha1-0.9));

}
