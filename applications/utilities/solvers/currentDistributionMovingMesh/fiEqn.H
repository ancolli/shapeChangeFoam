{
//potential field
    solve
    (
       - fvm::laplacian(keff, fi)
    );


Info<< "Min/max fi:" << min(fi).value() << ' '
    << max(fi).value() << endl;

//check if wall_anode exist
string patchNameNeeded = "wall_anode";
scalar patchOK = 0;

    forAll(mesh.boundaryMesh(), p_ith)
    {     
      string patchName = mesh.boundaryMesh()[p_ith].name();	
      if (patchNameNeeded == patchName)  // checks if "wall_anode" is in the current patchName
      {
        patchOK = 1;
      }
    }  

if (patchOK==1)  // checks if "wall_anode" is in the current patchName
      {
// if wall_anode exists then: Calculates flux on the anode

        label patchi = mesh.boundaryMesh().findPatchID("wall_anode");   
        Flux.boundaryFieldRef()[patchi] = keff.boundaryFieldRef()[patchi]*fi.boundaryFieldRef()[patchi].snGrad();
			    
// Calculates average Flux
	//scalar area = gSum(mesh.magSf().boundaryField()[patchi]);
        scalar sumField = 0.0;

        //if (area > 0)
        {
        	sumField = gSum
                (
                	mesh.magSf().boundaryField()[patchi]
                      * Flux.boundaryField()[patchi]
                );// /area;
        }

        Info<< "Total current = " << sumField << " A "<< endl;

// reading and writing probeLocation (reference electrode)
        label cellNo = mesh.findCell(vector(controlProperties.lookup("refLocation")));
	scalar Vref = fi[cellNo];
        Results.set("Vref", Vref);
// calculate cell potential
        
        scalar Vcell_old = 0.0;
        scalar var = 0.0;
        
        scalar mult = readScalar(Results.lookup("mult")); //scalar mult = 0.0;
        scalar Vcell = readScalar(Results.lookup("Vcell")); 

        scalar control = readScalar(controlProperties.lookup("control"));  
        scalar I = readScalar(controlProperties.lookup("I"));
        scalar maxPot = readScalar(controlProperties.lookup("maxPot"));
        scalar alpha = readScalar(controlProperties.lookup("alpha"));
        scalar tol_I = readScalar(controlProperties.lookup("tol_I"));

        scalar ba = readScalar(transportProperties.lookup("ba"));   
        
        Vcell=(maxPot < 0.0) ? -mag(Vcell) : mag(Vcell); // 
// potentiostatic mode
        scalar maxvar = readScalar(controlProperties.lookup("maxvar"));
        scalar tol_ref = readScalar(controlProperties.lookup("tol_ref"));
	scalar ref = readScalar(controlProperties.lookup("ref"));
	
        var=(ref < 0.0) ? -ref-Vref : ref-(Vcell-Vref); // <0 cathodic overpotential, >0 anodic overpotential      
        var=(var > 0.0) ? min(var,maxvar) : max(var,-maxvar);
      
        var=(mag(var/ref) > tol_ref) ? var : 0.0; // limito ref a x%
        var = (control==1) ? var : ((control==2) ? 0.0 : 0.0); //
        
// galvanostatic   
        mult = (mag(sumField-I)/mag(I) > tol_I) ? (sumField-I)/(I/ba) : 0.0; // limito corriente a x%	
// overrelaxation 
        Vcell_old = Vcell;
        
        mult=(maxPot < 0.0) ? -mult : mult; // 
        
	Vcell = (control==1) ? Vcell+var : ((control==2) ? maxPot : Vcell-mult); //Vcell = Vcell-mult;//
        Vcell = Vcell_old+alpha*(Vcell-Vcell_old);
	//Vcell = max(min(Vcell, maxPot),ref); // limited potentials
        
        // writing Vcell, I and Vref  
        Results.set("Vcell", Vcell);
        Results.set("I", sumField); 
        Results.set("Vref", (ref < 0.0) ? -Vref : (Vcell-Vref));//
        Results.set("mult", mult);
        Results.set("var", var); 
               
   }
   else
   {
	FatalErrorIn
        (
            " Foam::patchProbes::findElements(const fvMesh&)"
        )   << " No patch name: wall_anode found "
            << " the patch name: wall_anode should exist "
            << endl
            << exit(FatalError);     
   }
}
    

