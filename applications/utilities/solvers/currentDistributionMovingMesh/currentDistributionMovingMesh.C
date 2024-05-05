/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    currentDistribution

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids and the Laplace Equiation."
    );
    
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
//    #include "createUfIfPresent.H"
//    #include "CourantNo.H"
//    #include "setInitialDeltaT.H"
	
    //turbulence->validate();
    
   /* if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }*/


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
     // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    scalar mult = 0.0;
    scalar var = 0.0;

    while (runTime.run())
    {      
        //#include "readDyMControls.H"

        /*if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }*/

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
        	do {
                	#include "fiEqn.H"
                	mult = readScalar(Results.lookup("mult"));
                	var = readScalar(Results.lookup("var"));	                	 
                }
                while (mag(mult) > 1e-5 || mag(var) > 1e-5);
                
                
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    //MRF.update();

                    /*if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        //#include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        //#include "meshCourantNo.H"
                    }*/
                }
            }
            
            
            
            /*solve
            (
                - fvm::laplacian(kf, fi)
            );*/

            //#include "UEqn.H"           
            //#include "CEqn.H"

            // --- Pressure corrector loop
            /*while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }*/
        }
        
        runTime.write();

        runTime.printExecutionTime(Info);
        
    }
    
    runTime.writeNow();

    Info<< "End\n" << endl;

    return 0;

}

// ************************************************************************* //
