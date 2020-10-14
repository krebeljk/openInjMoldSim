//A2 equation
Sp_imp1 = 6.0*CI*(pow(0.5*(D && D.T()),0.5));

if(FiberClosure=="IBOF"){			
        fvTensorMatrix A2Eqn
	(
			fvm::ddt(A2) + fvm::div(phi, A2) + fvm::SuSp(-fvc::div(phi),A2) +
			fvm::Sp(Sp_imp1,A2) + fvm::Sp(ClosureA4_imp,A2)
	);
	solve(A2Eqn == SourceA2 * pos(alpha1-0.9));
}
else{
	fvTensorMatrix A2Eqn
			(
			fvm::ddt(A2) + fvm::div(phi, A2) + fvm::SuSp(-fvc::div(phi),A2));
	solve(A2Eqn == SourceA2 * pos(alpha1-0.9));

}

	//make symmetric...prevent numerical roundoff errors
	A2Tmp = symm(A2);

if(CorrectA2) //switch set in fiberProperties
{
      for(label celli=0; celli<U.size(); celli++)
      {
	  A2[celli] = A2Tmp[celli];
      }
}
