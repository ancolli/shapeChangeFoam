codeInclude
#{
    #include "addToRunTimeSelectionTable.H"
    #include "volFields.H"
    #include "pointPatchFields.H"
    #include "PrimitivePatchInterpolation.H"

#};

code
#{
    vectorField p(this->patch().localPoints());
    vectorField displacement(p.size(),vector(0,0,0));
    
    const polyMesh& mesh = this->internalField().mesh()();
    const label& patchID = this->patch().index();
    const vectorField faceNorm = -mesh.boundaryMesh()[patchID].faceNormals(); 
    //vectorField faceNorm(p.size(),vector(0,1,0));
    
    const volScalarField& field = this->db().objectRegistry::lookupObject<volScalarField>("fi");
    const volScalarField& keff = this->db().objectRegistry::lookupObject<volScalarField>("keff");
    
    scalarField Flux = -(63.5/8960/1000/96485/2.0)*keff.boundaryField()[patchID]*field.boundaryField()[patchID].snGrad(); 
    


    PrimitivePatchInterpolation<primitivePatch> patchInterpolator (field.mesh().boundaryMesh()[patchID]); 
    
    vectorField FluxPointValues = patchInterpolator.faceToPointInterpolate(Flux*faceNorm); 
    
    scalar MaxBump = max(p.component(0));
    scalar MinBump = min(p.component(0));
        
    scalar MaxParallel = returnReduce(MaxBump, maxOp<scalar>());
    scalar MinParallel = returnReduce(MinBump, minOp<scalar>());
   
    displacement = FluxPointValues; 
    
    Info << "displacement = " <<  max(mag(displacement)) << endl;
    
    forAll(p, idx) 
    {
    	displacement[idx].x() = (p[idx].x() == MinParallel || p[idx].x() == MaxParallel) ? vector(0,0,0).x() : FluxPointValues[idx].x();
    }      
    operator==(displacement);   
#};
   
codeOptions
#{
    -I$(WM_PROJECT_USER_DIR)/applications/lnInclude
#};
