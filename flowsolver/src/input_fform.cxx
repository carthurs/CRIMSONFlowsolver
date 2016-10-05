#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string.h>

#include "CInput.h"
#include "common_c.h"

using namespace std;
void print_error_code(int ierr);

int input_fform() {

	int ierr = 0;
	int i, j;
	char* path_to_config = 0, *inpfilename_env = 0;
	char complete_filename[256], inpfname[256];

	try {
		// get the input file stream
		path_to_config = getenv("PHASTA_CONFIG");
		if (path_to_config)
			strcpy(complete_filename, path_to_config);
		else
			strcpy(complete_filename, ".");
		strcat(complete_filename, "/input.config");
		inpfilename_env = getenv("PHASTA_INPFILE");
		if (inpfilename_env)
			strcpy(inpfname, inpfilename_env);
		else
			strcpy(inpfname, "solver.inp");
		printf("\n Complete Filename: %s \n", complete_filename);
		printf("\n Local Config: %s \n\n", inpfname);
		string def(complete_filename);
		CInput inp(inpfname, def);

		// Disabled Features

		conpar.iALE = inp.GetValue("iALE");
		conpar.icoord = inp.GetValue("icoord");
		conpar.irs = inp.GetValue("irs");
		conpar.iexec = inp.GetValue("iexec");
		timpar.ntseq = inp.GetValue("ntseq");
		solpar.imap = inp.GetValue("imap");

		// ALE
		if( (string) inp.GetValue("Arbitrary Lagrangian Eulerian description") == "True")
		{
			aleFlags.aleOn = 1; 
		}			
		else
		{
			aleFlags.aleOn = 0; 
		} 

		if( (string) inp.GetValue("Rigid body motion") == "True")
		{
			aleFlags.rigidOn = 1; 
		}			
		else
		{
			aleFlags.rigidOn = 0; 
		} 
		
		
		// Solution Control Keywords

		if ((string) inp.GetValue("Equation of State") == "Incompressible")
			matdat.matflg[0][0] = -1;
		if ((string) inp.GetValue("Equation of State") == "Compressible")
			matdat.matflg[0][0] = 0;
		inpdat.Delt[0] = inp.GetValue("Time Step Size");
		inpdat.nstep[0] = inp.GetValue("Number of Timesteps");
		if ((string) inp.GetValue("Viscous Control") == "Viscous")
			conpar.navier = 1;
		else
			conpar.navier = 0;

//		if ((string) inp.GetValue("Turbulence Model") == "No-Model") {
//			turbvari.irans = 0;
//			turbvari.iles = 0;
//		} else if ((string) inp.GetValue("Turbulence Model") == "LES") {
//			turbvari.iles = 1;
//			turbvari.irans = 0;
//		} else if ((string) inp.GetValue("Turbulence Model") == "RANS-SA") {
//			turbvari.iles = 0;
//			turbvari.irans = -1;
//		} else if ((string) inp.GetValue("Turbulence Model") == "RANS") {
//			turbvari.iles = 0;
//			turbvari.irans = -1; // assume S-A for backward compatibility
//		} else if ((string) inp.GetValue("Turbulence Model") == "RANS-KE") {
//			turbvari.iles = 0;
//			turbvari.irans = -2;
//		} else if ((string) inp.GetValue("Turbulence Model") == "DES") {
//			turbvari.iles = 1;
//			turbvari.irans = -1;
//		} else {
//			cout
//			<< " Turbulence Model: Only Legal Values ( No-Model, LES, RANS-SA, RANS-KE, DES )";
//			cout << endl;
//			exit(1);
//		}

		//if (turbvari.iles * turbvari.irans != 0)
			//turbvar.eles = inp.GetValue("DES Edge Length");

		int solflow, solheat, solscalr, ilset;
		((string) inp.GetValue("Solve Flow") == "True") ?
				solflow = 1 : solflow = 0;
		((string) inp.GetValue("Solve Heat") == "True") ?
				solheat = 1 : solheat = 0;
		//for compressible solheat= False so
		if ((string) inp.GetValue("Equation of State") == "Compressible")
			solheat = 0;
		ilset = (int) inp.GetValue("Solve Level Set");
		solscalr = (int) inp.GetValue("Solve Scalars");
		solscalr += ilset;
//		if (turbvari.irans == -1)
//			solscalr++;
//		if (turbvari.irans == -2)
//			solscalr = solscalr + 2;
		if (solscalr > 4) {
			cout << " Only Four Scalars are supported \n";
			cout << " Please reduce number of scalars \n";
			exit(1);
		}
		inpdat.impl[0] = 10 * solflow + solscalr * 100 + solheat;

		levlset.iLSet = ilset;
		if (ilset > 0) {
			levlset.epsilon_ls = inp.GetValue(
					"Number of Elements Across Interface");
			levlset.epsilon_lsd = inp.GetValue(
					"Number of Elements Across Interface for Redistancing");
			levlset.dtlset = inp.GetValue("Pseudo Time step for Redistancing");
			levlset.iExpLSSclr2 = inp.GetValue(
					"Explicit Solve for Redistance Field");
			levlset.iExpLSSclr1 = inp.GetValue(
					"Explicit Solve for Scalar 1 Field");
			if ((string) inp.GetValue("Apply Volume Constraint") == "True") {
				levlset.ivconstraint = 1;
			} else if ((string) inp.GetValue("Apply Volume Constraint")
					== "False") {
				levlset.ivconstraint = 0;
			} else {
				cout
				<< "Apply Volume Constraint: Only Legal Values (True, False) ";
				cout << endl;
				exit(1);
			}
		}

		vector<double> vec;

		// OUTPUT CONTROL KEY WORDS.

		conpar.necho = inp.GetValue("Verbosity Level");
		outpar.ntout = inp.GetValue("Number of Timesteps between Restarts");
		if ((string) inp.GetValue("Print Statistics") == "True")
			outpar.ioform = 2;
		else
			outpar.ioform = 1;

		if ((string) inp.GetValue("Print Wall Fluxes") == "True")
			outpar.iowflux = 1;
		else
			outpar.iowflux = 0;

		if ((string) inp.GetValue("Print FieldView") == "True")
			outpar.iofieldv = 1;
		else
			outpar.iofieldv = 0;

		if ((string) inp.GetValue("Print ybar") == "True")
			outpar.ioybar = 1;
		else
			outpar.ioybar = 0;

		strcpy(outpar.iotype,
				((string) inp.GetValue("Data Block Format")).c_str());
		//strcpy( phasta_iotype , ((string)inp.GetValue("Data Block Format")).c_str());
		//turbvar.sonfathvar = inp.GetValue("Number of Father Nodes");

		if ((string) inp.GetValue("Print Residual at End of Step") == "True")
			genpar.lstres = 1;
		else
			genpar.lstres = 0;

		if ((string) inp.GetValue("Print Error Indicators") == "True")
			turbvar.ierrcalc = 1;
		else
			turbvar.ierrcalc = 0;

		if ((string) inp.GetValue("Print Velocity Hessian") == "True")
			turbvar.ihessian = 1;
		else
			turbvar.ihessian = 0;

		if (turbvar.ierrcalc == 1)
			turbvari.ierrsmooth = inp.GetValue(
					"Number of Error Smoothing Iterations");

		int nsrfCM = inp.GetValue("Number of Force Surfaces");
		if (nsrfCM > 0) {
			vector<int> ivec = inp.GetValue(
					"Surface ID's for Force Calculation");
			for (i = 0; i < MAXSURF + 1; i++)
				aerfrc.nsrflist[i] = 0;
			for (i = 0; i < nsrfCM; i++) {
				aerfrc.nsrflist[ivec[i]] = 1;
				//        cout <<"surface in force list "<< ivec[i] << endl;
			}
			ivec.erase(ivec.begin(), ivec.end());
		}

		aerfrc.isrfIM = inp.GetValue("Surface ID for Integrated Mass");
		//Limiting
//		vec = inp.GetValue("Limit u1");
//		for (i = 0; i < 3; i++) {
//			turbvar.ylimit[0][i] = vec[i];
//		}
//		vec.erase(vec.begin(), vec.end());
//
//		vec = inp.GetValue("Limit u2");
//		for (i = 0; i < 3; i++) {
//			turbvar.ylimit[1][i] = vec[i];
//		}
//		vec.erase(vec.begin(), vec.end());
//
//		vec = inp.GetValue("Limit u3");
//		for (i = 0; i < 3; i++) {
//			turbvar.ylimit[2][i] = vec[i];
//		}
//		vec.erase(vec.begin(), vec.end());
//
//		vec = inp.GetValue("Limit Pressure");
//		for (i = 0; i < 3; i++) {
//			turbvar.ylimit[3][i] = vec[i];
//		}
//		vec.erase(vec.begin(), vec.end());
//
//		vec = inp.GetValue("Limit Temperature");
//		for (i = 0; i < 3; i++) {
//			turbvar.ylimit[4][i] = vec[i];
//		}
//		vec.erase(vec.begin(), vec.end());

		//Material Properties Keywords
		matdat.nummat = levlset.iLSet + 1;
		if ((string) inp.GetValue("Shear Law") == "Constant Viscosity")
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][1] = 0;

		if ((string) inp.GetValue("Bulk Viscosity Law")
				== "Constant Bulk Viscosity")
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][2] = 0;

		mmatpar.pr = inp.GetValue("Prandtl Number");

		if ((string) inp.GetValue("Conductivity Law")
				== "Constant Conductivity")
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][3] = 0;

		vec = inp.GetValue("Density");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][0][0] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Viscosity");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][1][0] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		//      vec = inp.GetValue("Specific Heat");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][2][0] = 0;
		}
		//      vec.erase(vec.begin(),vec.end());

		vec = inp.GetValue("Thermal Conductivity");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][3][0] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Scalar Diffusivity");
		for (i = 0; i < solscalr; i++) {
			sclrs.scdiff[i] = vec[i];
		}
		vec.erase(vec.begin(), vec.end());

		//if ((string) inp.GetValue("Zero Mean Pressure") == "True")
		//	turbvar.pzero = 1;

		//turbvar.rmutarget = inp.GetValue("Target Viscosity For Step NSTEP");

		if ((string) inp.GetValue("Body Force Option") == "None") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 0;
		} else if ((string) inp.GetValue("Body Force Option") == "Vector") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 1;
		} else if ((string) inp.GetValue("Body Force Option")
				== "User e3source.f") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 3;
		} else if ((string) inp.GetValue("Body Force Option") == "Boussinesq") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 2;
		} else if ((string) inp.GetValue("Body Force Option")
				== "Cooling Analytic") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 4;
		} else if ((string) inp.GetValue("Body Force Option")
				== "Cooling Initial Condition") {
			for (i = 0; i < levlset.iLSet + 1; i++)
				matdat.matflg[i][4] = 5;
		}

		// the following block of stuff is common to all cooling type sponges.
		// Specific stuff belongs in the conditionals above

//		if (matdat.matflg[0][4] >= 4) {
//			spongevar.betamax = inp.GetValue(
//					"Maximum Value of Sponge Parameter");
//			spongevar.zinsponge = inp.GetValue(
//					"Inflow Cooling Sponge Ends at z");
//			spongevar.zoutsponge = inp.GetValue(
//					"Outflow Cooling Sponge Begins at z");
//			spongevar.radsponge = inp.GetValue(
//					"Radial Cooling Sponge Begins at r");
//			spongevar.grthosponge = inp.GetValue(
//					"Sponge Growth Coefficient Outflow");
//			spongevar.grthisponge = inp.GetValue(
//					"Sponge Growth Coefficient Inflow");
//
//			spongevar.spongecontinuity = 0;
//			spongevar.spongemomentum1 = 0;
//			spongevar.spongemomentum2 = 0;
//			spongevar.spongemomentum3 = 0;
//			spongevar.spongeenergy = 0;
//
//			if ((string) inp.GetValue("Sponge for Continuity Equation")
//					== "True")
//				spongevar.spongecontinuity = 1;
//			if ((string) inp.GetValue("Sponge for x Momentum Equation")
//					== "True")
//				spongevar.spongemomentum1 = 1;
//			if ((string) inp.GetValue("Sponge for y Momentum Equation")
//					== "True")
//				spongevar.spongemomentum2 = 1;
//			if ((string) inp.GetValue("Sponge for z Momentum Equation")
//					== "True")
//				spongevar.spongemomentum3 = 1;
//			if ((string) inp.GetValue("Sponge for Energy Equation") == "True")
//				spongevar.spongeenergy = 1;
//
//		}

		vec = inp.GetValue("Body Force");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][4][0] = vec[0 + i * 3];
			matdat.datmat[i][4][1] = vec[1 + i * 3];
			matdat.datmat[i][4][2] = vec[2 + i * 3];
		}
		vec.erase(vec.begin(), vec.end());

		vec = inp.GetValue("Body Force Pressure Gradient");
		for (i = 0; i < levlset.iLSet + 1; i++) {
			matdat.datmat[i][6][0] = vec[0 + i * 3];
			matdat.datmat[i][6][1] = vec[1 + i * 3];
			matdat.datmat[i][6][2] = vec[2 + i * 3];
		}
		vec.erase(vec.begin(), vec.end());

		if ((string) inp.GetValue("Surface Tension Option") == "No") {
			genpar.isurf = 0;
		} else if ((string) inp.GetValue("Surface Tension Option") == "Yes") {
			genpar.isurf = 1;
		} else {
			cout << " Surface Tension: Only Legal Values (Yes, No) ";
			cout << endl;
			exit(1);
		}
		if (genpar.isurf > 0) {
			genpar.Bo = inp.GetValue("Bond Number");
		}

		genpar.EntropyPressure = inp.GetValue(
				"Entropy Form of Pressure Constraint on Weight Space");

		if ((string) inp.GetValue("Rotating Frame of Reference") == "True") {
			matdat.matflg[0][5] = 1;
			vec = inp.GetValue("Rotating Frame of Reference Rotation Rate");
			matdat.datmat[0][5][0] = vec[0];
			matdat.datmat[0][5][1] = vec[1];
			matdat.datmat[0][5][2] = vec[2];
			vec.erase(vec.begin(), vec.end());
		} else {
			matdat.matflg[0][5] = 0;
			matdat.datmat[0][5][0] = 0.;
			matdat.datmat[0][5][1] = 0.;
			matdat.datmat[0][5][2] = 0.;
		}

		//Linear Solver parameters
		inpdat.memLSFlag=0;
		if ((string) inp.GetValue("Solver Type")
				== "ACUSIM with P Projection") {
			incomp.iprjFlag = 0;
			incomp.ipresPrjFlag = 1;
		} else if ((string) inp.GetValue("Solver Type") == "ACUSIM") {
			incomp.iprjFlag = 0;
			incomp.ipresPrjFlag = 0;
		} else if ((string) inp.GetValue("Solver Type")
				== "ACUSIM with Velocity Projection") {
			incomp.iprjFlag = 1;
			incomp.ipresPrjFlag = 0;
		} else if ((string) inp.GetValue("Solver Type")
				== "ACUSIM with Full Projection") {
			incomp.iprjFlag = 1;
			incomp.ipresPrjFlag = 1;
		} else if ((string) inp.GetValue("Solver Type")
				== "GMRES Matrix Free") {
			inpdat.impl[0] += 10 * solflow;
		} else if ((string) inp.GetValue("Solver Type") == "GMRES EBE") {
			inpdat.impl[0] += 20 * solflow;
		} else if( (string)inp.GetValue("Solver Type") =="memLS"){
		    inpdat.memLSFlag=1;
		    printf("**********************************************************************\n"); 
		    printf("*** IF USING MEMLS ENSURE THE SETTINGS IN INPUT.CONFIG ARE CHANGED ***\n");		    
		    printf("**********************************************************************\n"); 
		}
		//GMRES sparse is assumed default and has the value of 10, MFG 20,
		// EBE 30

		//    inpdat.niter[0] = inp.GetValue("Number of Solves per Time Step");
		solpar.nGMRES = inp.GetValue("Number of GMRES Sweeps per Solve");
		solpar.Kspace = inp.GetValue(
				"Number of Krylov Vectors per GMRES Sweep");
		inpdat.LHSupd[0] = inp.GetValue(
				"Number of Solves per Left-hand-side Formation");
		inpdat.epstol[0] = inp.GetValue("Tolerance on Momentum Equations");
		inpdat.epstol[6] = inp.GetValue("Tolerance on Continuity Equations");
		inpdat.epstol[7] = inp.GetValue("Tolerance on memLS NS Solver");
		incomp.prestol = inp.GetValue(
				"Tolerance on ACUSIM Pressure Projection");
		incomp.minIters = inp.GetValue(
				"Minimum Number of Iterations per Nonlinear Iteration");
		incomp.maxIters = inp.GetValue(
				"Maximum Number of Iterations per Nonlinear Iteration");
		inpdat.deltol[0][0] = inp.GetValue("Velocity Delta Ratio");
		inpdat.deltol[1][0] = inp.GetValue("Pressure Delta Ratio");
		incomp.nPrjs = inp.GetValue("Number of Velocity Projection Vectors");
		incomp.nPresPrjs = inp.GetValue(
				"Number of Pressure Projection Vectors");
		incomp.iverbose = inp.GetValue("ACUSIM Verbosity Level");

		if (solheat == 1) {
			inpdat.epstol[1] = inp.GetValue("Temperature Solver Tolerance");
			inpdat.LHSupd[1] =
					inp.GetValue(
							"Number of Solves of Temperature per Left-hand-side Formation");
		}

		// The following is where you should put any inputs that are able to
		// input differently for each scalar.  It is a little tedious in the code
		// but it should make the solver.inp easier to understand. Note this will
		// require some care with regression tests.

		if (solscalr > 0) {
			inpdat.epstol[2] = inp.GetValue("Scalar 1 Solver Tolerance");
			inpdat.LHSupd[2] =
					inp.GetValue(
							"Number of Solves of Scalar 1 per Left-hand-side Formation");

//			vec = inp.GetValue("Limit Scalar 1");
//			for (i = 0; i < 3; i++) {
//				turbvar.ylimit[5][i] = vec[i];
//			}
//			vec.erase(vec.begin(), vec.end());
		}

		if (solscalr > 1) {
			inpdat.epstol[3] = inp.GetValue("Scalar 2 Solver Tolerance");
			inpdat.LHSupd[3] =
					inp.GetValue(
							"Number of Solves of Scalar 2 per Left-hand-side Formation");

//			vec = inp.GetValue("Limit Scalar 2");
//			for (i = 0; i < 3; i++) {
//				turbvar.ylimit[6][i] = vec[i];
//			}
//			vec.erase(vec.begin(), vec.end());
		}

		if (solscalr > 2) {
			inpdat.epstol[4] = inp.GetValue("Scalar 3 Solver Tolerance");
			inpdat.LHSupd[4] =
					inp.GetValue(
							"Number of Solves of Scalar 3 per Left-hand-side Formation");

//			vec = inp.GetValue("Limit Scalar 3");
//			for (i = 0; i < 3; i++) {
//				turbvar.ylimit[7][i] = vec[i];
//			}
//			vec.erase(vec.begin(), vec.end());
		}

		if (solscalr > 3) {
			inpdat.epstol[5] = inp.GetValue("Scalar 4 Solver Tolerance");
			inpdat.LHSupd[5] =
					inp.GetValue(
							"Number of Solves of Scalar 4 per Left-hand-side Formation");

//			vec = inp.GetValue("Limit Scalar 4");
//			for (i = 0; i < 3; i++) {
//				turbvar.ylimit[8][i] = vec[i];
//			}
//			vec.erase(vec.begin(), vec.end());
		}

		// DISCRETIZATION CONTROL

		genpar.ipord = inp.GetValue("Basis Function Order");
		if ((string) inp.GetValue("Time Integration Rule") == "First Order")
			inpdat.rhoinf[0] = -1;
		else
			inpdat.rhoinf[0] = (double) inp.GetValue(
					"Time Integration Rho Infinity");
		if ((string) inp.GetValue("Predictor at Start of Step")
				== "Same Velocity")
			genpar.ipred = 1;
		if ((string) inp.GetValue("Predictor at Start of Step")
				== "Zero Acceleration")
			genpar.ipred = 2;
		if ((string) inp.GetValue("Predictor at Start of Step")
				== "Same Acceleration")
			genpar.ipred = 3;
		if ((string) inp.GetValue("Predictor at Start of Step") == "Same Delta")
			genpar.ipred = 4;

		if ((string) inp.GetValue("Weak Form") == "Galerkin")
			solpar.ivart = 1;
		if ((string) inp.GetValue("Weak Form") == "SUPG")
			solpar.ivart = 2;

		if ((string) inp.GetValue("Flow Advection Form") == "Convective")
			solpar.iconvflow = 2;
		else if ((string) inp.GetValue("Flow Advection Form") == "Conservative")
			solpar.iconvflow = 1;
		if ((string) inp.GetValue("Scalar Advection Form") == "Convective")
			solpar.iconvsclr = 2;
		else if ((string) inp.GetValue("Scalar Advection Form")
				== "Conservative")
			solpar.iconvsclr = 1;
		if ((string) inp.GetValue("Use Conservative Scalar Convection Velocity")
				== "True")
			sclrs.consrv_sclr_conv_vel = 1;
		else if ((string) inp.GetValue(
				"Use Conservative Scalar Convection Velocity") == "False")
			sclrs.consrv_sclr_conv_vel = 0;
		// TAU INPUT
		if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Shakib")
			genpar.itau = 0;
		else if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Franca")
			genpar.itau = 1;
		else if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Jansen(dev)")
			genpar.itau = 2;
		else if ((string) inp.GetValue("Tau Matrix") == "Diagonal-Compressible")
			genpar.itau = 3;
		else if ((string) inp.GetValue("Tau Matrix") == "Matrix-Mallet")
			genpar.itau = 10;
		else if ((string) inp.GetValue("Tau Matrix") == "Matrix-Modal")
			genpar.itau = 11;

		genpar.dtsfct = inp.GetValue("Tau Time Constant");
		genpar.taucfct = inp.GetValue("Tau C Scale Factor");

		// FLOW DISCONTINUITY CAPTURING

		if ((string) inp.GetValue("Discontinuity Capturing") == "Off")
			solpar.iDC = 0;
		else if ((string) inp.GetValue("Discontinuity Capturing")
				== "DC-mallet")
			solpar.iDC = 1;
		else if ((string) inp.GetValue("Discontinuity Capturing")
				== "DC-quadratic")
			solpar.iDC = 2;
		else if ((string) inp.GetValue("Discontinuity Capturing")
				== "DC-minimum")
			solpar.iDC = 3;
		else {
			cout << "Condition not defined for Discontinuity Capturing \n ";
			exit(1);
		}

		// SCALAR DISCONTINUITY CAPTURING

		vector<int> ivec = inp.GetValue("Scalar Discontinuity Capturing");
		for (i = 0; i < 2; i++)
			solpar.idcsclr[i] = ivec[i];
		ivec.erase(ivec.begin(), ivec.end());

		//        if((string)inp.GetValue("Scalar Discontinuity Capturing") == "No") solpar.idcsclr = 0;
		//      else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "1") solpar.idcsclr = 1;
		//   else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "2") solpar.idcsclr = 2;
		//   else {
		//        cout<< "Condition not defined for Scalar Discontinuity Capturing \n ";
		//        exit(1);
		//      }
		if ((string) inp.GetValue("Include Viscous Correction in Stabilization")
				== "True") {
			if (genpar.ipord == 1)
				genpar.idiff = 1;
			else
				genpar.idiff = 2;
		} else {
			genpar.idiff = 0;
		}

		timdat.flmpl = inp.GetValue("Lumped Mass Fraction on Left-hand-side");
		timdat.flmpr = inp.GetValue("Lumped Mass Fraction on Right-hand-side");

		if ((string) inp.GetValue("Dump CFL") == "True")
			timdat.iCFLworst = 1;

		intdat.intg[0][0] = inp.GetValue("Quadrature Rule on Interior");
		intdat.intg[0][1] = inp.GetValue("Quadrature Rule on Boundary");
		genpar.ibksiz = inp.GetValue("Number of Elements Per Block");

		((string) inp.GetValue("Turn Off Source Terms for Scalars") == "True") ?
				sclrs.nosource = 1 : sclrs.nosource = 0;
		sclrs.tdecay = inp.GetValue("Decay Multiplier for Scalars");

		// TURBULENCE MODELING PARAMETER
//		int tpturb = turbvari.iles - turbvari.irans;
//		int ifrule;
//		if (tpturb != 0) {
//
//			turbvari.nohomog = inp.GetValue("Number of Homogenous Directions");
//
//			if ((string) inp.GetValue("Turbulence Wall Model Type")
//					== "Slip Velocity")
//				turbvar.itwmod = 1;
//			else if ((string) inp.GetValue("Turbulence Wall Model Type")
//					== "Effective Viscosity")
//				turbvar.itwmod = 2;
//			else
//				turbvar.itwmod = 0;
//			if (turbvari.irans < 0)
//				turbvar.itwmod = turbvar.itwmod * (-1);
//			ifrule = inp.GetValue("Velocity Averaging Steps");
//			turbvar.wtavei = (ifrule > 0) ? 1.0 / ifrule : -1.0 / ifrule;
//
//			if (turbvari.iles == 1) {
//
//				if ((string) inp.GetValue("Dynamic Model Type") == "Bardina")
//					turbvari.iles += 10;
//				else if ((string) inp.GetValue("Dynamic Model Type")
//						== "Projection")
//					turbvari.iles += 20;
//
//				ifrule = inp.GetValue("Filter Integration Rule");
//				turbvari.iles += ifrule - 1;
//				ifrule = inp.GetValue("Dynamic Model Averaging Steps");
//				turbvar.dtavei = (ifrule > 0) ? 1.0 / ifrule : -1.0 / ifrule;
//				turbvar.fwr1 = inp.GetValue("Filter Width Ratio");
//				turbvar.flump = inp.GetValue("Lumping Factor for Filter");
//
//				if ((string) inp.GetValue("Model Statistics") == "True") {
//					turbvari.modlstats = 1;
//				} else {
//					turbvari.modlstats = 0;
//				}
//
//				if ((string) inp.GetValue("Double Filter") == "True") {
//					turbvari.i2filt = 1;
//				} else {
//					turbvari.i2filt = 0;
//				}
//
//				if ((string) inp.GetValue("Model/SUPG Dissipation") == "True") {
//					turbvari.idis = 1;
//				} else {
//					turbvari.idis = 0;
//				}
//
//				if ((string) inp.GetValue("Dynamic Model Type") == "Standard") {
//
//					if ((string) inp.GetValue("Dynamic Sub-Model Type")
//							== "None")
//						turbvari.isubmod = 0;
//					else if ((string) inp.GetValue("Dynamic Sub-Model Type")
//							== "DFWR")
//						turbvari.isubmod = 1;
//					else if ((string) inp.GetValue("Dynamic Sub-Model Type")
//							== "SUPG")
//						turbvari.isubmod = 2;
//				} else if ((string) inp.GetValue("Dynamic Model Type")
//						== "Projection") {
//
//					if ((string) inp.GetValue("Projection Filter Type")
//							== "Linear")
//						turbvari.ifproj = 0;
//					else if ((string) inp.GetValue("Projection Filter Type")
//							== "Quadratic")
//						turbvari.ifproj = 1;
//
//					if ((string) inp.GetValue("Dynamic Sub-Model Type")
//							== "None")
//						turbvari.isubmod = 0;
//					else if ((string) inp.GetValue("Dynamic Sub-Model Type")
//							== "ConsistentProj")
//						turbvari.isubmod = 1;
//				}
//
//			}
//		}

		// SPEBC MODELING PARAMETERS

//		if ((spebcvr.irscale = inp.GetValue("SPEBC Model Active")) >= 0) {
//
//			//ifrule = inp.GetValue("Velocity Averaging Steps");
//			//turbvar.wtavei =
//			//		(ifrule > 0) ? 1.0 / ifrule : 1.0 / inpdat.nstep[0];
//			spebcvr.intpres = inp.GetValue("Interpolate Pressure");
//			spebcvr.plandist = inp.GetValue("Distance between Planes");
//			spebcvr.thetag = inp.GetValue("Theta Angle of Arc");
//			spebcvr.ds = inp.GetValue("Distance for Velocity Averaging");
//			spebcvr.tolerence = inp.GetValue("SPEBC Cylindrical Tolerance");
//			spebcvr.radcyl = inp.GetValue("Radius of recycle plane");
//			spebcvr.rbltin = inp.GetValue("Inlet Boundary Layer Thickness");
//			spebcvr.rvscal = inp.GetValue("Vertical Velocity Scale Factor");
//		}

		// CARDIOVASCULAR MODELING PARAMETERS
		if ((string) inp.GetValue("Time Varying Boundary Conditions From File")
				== "True")
			nomodule.itvn = 1;
		else
			nomodule.itvn = 0;
		if (nomodule.itvn == 1)
			nomodule.bcttimescale = inp.GetValue("BCT Time Scale Factor");

		nomodule.ipvsq = 0;
		if (nomodule.icardio = inp.GetValue("Number of Coupled Surfaces")) {
			if (nomodule.icardio > MAXSURF) {
				cout << "Number of Coupled Surfaces > MAXSURF \n";
				exit(1);
			}
			if ((string) inp.GetValue("Pressure Coupling") == "None")
				nomodule.ipvsq = 0;
			if ((string) inp.GetValue("Pressure Coupling") == "Explicit")
				nomodule.ipvsq = 1;
			if ((string) inp.GetValue("Pressure Coupling") == "Implicit")
				nomodule.ipvsq = 2;
			if ((string) inp.GetValue("Pressure Coupling") == "P-Implicit")
				nomodule.ipvsq = 3;
			
           /**********************************************************
            ***          Influx Stabilisation Coefficient          ***
            **********************************************************/
      
            nomodule.stabflux_coeff = inp.GetValue("Influx Coefficient"); 

     		/**********************************************************
      		***                Global Node Numbering               ***
      		**********************************************************/

      		nomodule.indsurf = int(0);
      		if ((string) inp.GetValue("Global Node Numbering") == "True")
      		{
          		nomodule.indsurf = int(1);
      		}

      		nomodule.geombcHasObservationFields = int(1);
      		if ((string) inp.GetValue("Geombc Has Observation Fields") == "False")
      		{
      			nomodule.geombcHasObservationFields = int(0);
      		}

      		nomodule.geombcHasNodeTags = int(1);
	      	if ((string) inp.GetValue("Geombc Has Node Tags") == "False")
	      	{
	      		nomodule.geombcHasNodeTags = int(0);
	      	}

		    /**********************************************************
		     ***               Heart Model Parameters               ***
		     **********************************************************/

		    if ((string) inp.GetValue("Heart Model") == "True")
		    {
	            nomodule.iheart = 1;        
                nomodule.heartparam[0] = inp.GetValue("Aortic Surface");

		        nomodule.heartparam[1] = inp.GetValue("Preload");
                nomodule.heartparam[2] = inp.GetValue("Mitral Valve");

                nomodule.heartparam[3] = inp.GetValue("End Diastolic Volume");
		        nomodule.heartparam[4] = inp.GetValue("Unstressed Volume");
                		       
		        nomodule.heartparam[5] = inp.GetValue("Maximum Elastance");
                nomodule.heartparam[6] = inp.GetValue("Time to Max Elastance");     
		        nomodule.heartparam[7] = inp.GetValue("Time to Relax");     
                nomodule.heartparam[8] = inp.GetValue("Period");		        

                nomodule.heartparam[9] = inp.GetValue("Ventricular Resistance");     
		        nomodule.heartparam[10] = inp.GetValue("Aortic Valve");
		        
		        		       		      
		        // backflow parameters        
		        if ((string) inp.GetValue("Backflow") == "True")
		        {
		           nomodule.heartparam[11] = 1;    
		           nomodule.heartparam[12] = inp.GetValue("Backflow Magnitude");
		           nomodule.heartparam[13] = inp.GetValue("Backflow Steepness");
		           nomodule.heartparam[14] = inp.GetValue("Backflow Closure");
		           nomodule.heartparam[15] = inp.GetValue("Backflow Time");           
		        }
		        else 
		        {
		           nomodule.heartparam[11] = 0;
		           nomodule.heartparam[12] = 0.0;
		           nomodule.heartparam[13] = 0.0;
		           nomodule.heartparam[14] = 0.0;
		           nomodule.heartparam[15] = 0.0;
		        }        
		      } 
		      else 
		      {
		        nomodule.iheart  = 0;
		      };  


     		/**********************************************************
      		***                                                     ***
      		**********************************************************/

			if ((string) inp.GetValue("Inflow Coupling") == "True")
				nomodule.incp = 1;
			else
				nomodule.incp = 0;
			if (nomodule.incp == 1) {
				nomodule.numINCPSrfs = inp.GetValue(
						"Number of Coupled Inflow Surfaces");
				ivec = inp.GetValue("List of Coupled Inflow Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistINCP[i] = 0;
				for (i = 0; i < nomodule.numINCPSrfs; i++) {
					nomodule.nsrflistINCP[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("Inflow Parameters From File")
						== "True")
					nomodule.incpfile = 1;
				else
					nomodule.incpfile = 0;
			}
			if (nomodule.numResistSrfs = inp.GetValue(
					"Number of Resistance Surfaces")) {
				ivec = inp.GetValue("List of Resistance Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistResist[i] = 0;
				for (i = 0; i < nomodule.numResistSrfs; i++) {
					nomodule.nsrflistResist[i + 1] = ivec[i];
				}
				vec = inp.GetValue("Resistance Values");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.ValueListResist[i] = 0;
				for (i = 0; i < nomodule.numResistSrfs; i++)
					nomodule.ValueListResist[i + 1] = vec[i];
				vec.erase(vec.begin(), vec.end());
			}
			if (nomodule.numImpSrfs = inp.GetValue(
					"Number of Impedance Surfaces")) {
				ivec = inp.GetValue("List of Impedance Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistImp[i] = 0;
				for (i = 0; i < nomodule.numImpSrfs; i++) {
					nomodule.nsrflistImp[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("Impedance From File") == "True")
					nomodule.impfile = 1;
				else
					nomodule.impfile = 0;
			}
			nomodule.ircrfile = 0; // value remains if RCR Values From File == False; changed below if True
			if (nomodule.numRCRSrfs = inp.GetValue("Number of RCR Surfaces")) {
				ivec = inp.GetValue("List of RCR Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistRCR[i] = 0;
				for (i = 0; i < nomodule.numRCRSrfs; i++) {
					nomodule.nsrflistRCR[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("RCR Values From File") == "True")
				{
					nomodule.ircrfile = 1;
				}
			}
//			if(nomodule.numTRCRSrfs=inp.GetValue("Number of Time-varying RCR Surfaces")){
//				ivec = inp.GetValue("List of Time-varying RCR Surfaces");
//				for(i=0;i<MAXSURF+1; i++)
//					nomodule.nsrflistTRCR[i] = 0;
//				for(i=0; i< nomodule.numTRCRSrfs; i++){
//					nomodule.nsrflistTRCR[i+1] = ivec[i];
//				}
//				if ( (string)inp.GetValue("Time-varying RCR Values From File") == "True")
//					nomodule.itrcrfile = 1;
//				else
//					nomodule.itrcrfile = 0;
////				if ( (string)inp.GetValue("Regulation of Flow") == "True")
////					nomodule.regflow = 1; else nomodule.regflow = 0;
////				if(nomodule.regflow){
////					nomodule.numRegSrfs=inp.GetValue("Number of Regulated Surfaces");
////					ivec = inp.GetValue("List of Regulated Surfaces");
////					for(i=0;i<MAXSURF+1; i++) nomodule.nsrflistReg[i] = 0;
////					for(i=0; i< nomodule.numRegSrfs; i++){
////						nomodule.nsrflistReg[i+1] = ivec[i];
////					}
////				}
//			}
			// Nan rcr ---------------------------
			if(grcrbccom.numGRCRSrfs=inp.GetValue("Number of experimental RCR Surfaces")){
				ivec = inp.GetValue("List of experimental RCR Surfaces");
				for(i=0;i<MAXSURF+1; i++)
					grcrbccom.nsrflistGRCR[i] = 0;
				for(i=0; i< grcrbccom.numGRCRSrfs; i++){
					grcrbccom.nsrflistGRCR[i+1] = ivec[i];
				}
				if ( (string)inp.GetValue("experimental RCR Values From File") == "True")
					grcrbccom.igrcrfile = 1;
				else
					grcrbccom.igrcrfile = 0;

			}
			// -----------------------------------
			if (nomodule.numCORSrfs = inp.GetValue(
					"Number of Coronary Surfaces")) {
				ivec = inp.GetValue("List of Coronary Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistCOR[i] = 0;
				for (i = 0; i < nomodule.numCORSrfs; i++) {
					nomodule.nsrflistCOR[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue("Coronary Values From File")
						== "True")
					nomodule.icorfile = 1;
				else
					nomodule.icorfile = 0;
			}

			if(nomodule.numControlledCoronarySrfs=inp.GetValue("Number of Controlled Coronary Surfaces"))
			{
				ivec = inp.GetValue("List of Controlled Coronary Surfaces");

				for(i=0; i<MAXSURF+1; i++) 
				{
					nomodule.indicesOfCoronarySurfaces[i] = 0;
				}

				for(i=0; i< nomodule.numControlledCoronarySrfs; i++)
				{
					nomodule.indicesOfCoronarySurfaces[i+1]=ivec[i];
				}

			}

			if (nomodule.numNetlistLPNSrfs = inp.GetValue("Number of Netlist LPN Surfaces"))
			{
				ivec = inp.GetValue("List of Netlist LPN Surfaces");

				for (i=0; i<MAXSURF+1; i++)
				{
					nomodule.indicesOfNetlistSurfaces[i] = 0;
				}

				for (i=0; i<nomodule.numNetlistLPNSrfs; i++)
				{
					nomodule.indicesOfNetlistSurfaces[i+1] = ivec[i];
				}
			}

			nomodule.numLoopClosingCircuits = inp.GetValue("Number of Loop Closing Netlist Circuits");

			if((string)inp.GetValue("Input prescribed HR and peak systolic pressure from file")=="True"){
		      nomodule.inputHRandSP = int(1);
		    }
		    else
		    {
		      nomodule.inputHRandSP = int(0);
		    }

		    if((string)inp.GetValue("Simulate in Purely Zero Dimensions")=="True"){
		      nomodule.pureZeroDSimulation = int(1);
		    }
		    else
		    {
		      nomodule.pureZeroDSimulation = int(0);
		    }

		    if (nomodule.num3DConnectedComponents = inp.GetValue("Number of Connected Components of 3D Domain"))
			{
				for (i=0; i<MAXSURF+1; i++)
				{
					nomodule.surfacesOfEachConnectedComponent[i] = 0;
				}

				bool domainHasMultipleConnectedComponents = (nomodule.num3DConnectedComponents > 1);
				if (domainHasMultipleConnectedComponents)
				{
					ivec = inp.GetValue("List of Surfaces In Each Connected Component Separated by -1s");
					int lengthOfConnectedComponentLineInSolverInp = nomodule.numNetlistLPNSrfs + nomodule.num3DConnectedComponents - 1;
					for (i=0; i<lengthOfConnectedComponentLineInSolverInp; i++)
					{
						nomodule.surfacesOfEachConnectedComponent[i+1] = ivec[i];
					}
				}
				else // If there's only one connected component, then this vector should just be the netlist surface indices (as the zeroD domain only supports Netlist boundary conditions currently.)
				{
					for (i=0; i<nomodule.numNetlistLPNSrfs; i++)
					{
						nomodule.surfacesOfEachConnectedComponent[i+1] = nomodule.indicesOfNetlistSurfaces[i+1];
					}	
				}
			}

			if (nomodule.numVisFluxSrfs = inp.GetValue(
					"Number of Surfaces which zero out in-plane tractions")) {
				ivec = inp.GetValue(
						"List of Surfaces which zero out in-plane tractions");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistVisFlux[i] = 0;
				for (i = 0; i < nomodule.numVisFluxSrfs; i++) {
					nomodule.nsrflistVisFlux[i + 1] = ivec[i];
				}
			}
			if (nomodule.numCalcSrfs = inp.GetValue(
					"Number of Surfaces which Output Pressure and Flow")) {
				ivec = inp.GetValue("List of Output Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistCalc[i] = 0;
				for (i = 0; i < nomodule.numCalcSrfs; i++) {
					nomodule.nsrflistCalc[i + 1] = ivec[i];
				}
			}
			if (nomodule.numDirCalcSrfs =
					inp.GetValue(
							"Number of Dirichlet Surfaces Which Output Pressure and Flow")) {
				ivec = inp.GetValue("List of Dirichlet Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistDirCalc[i] = 0;
				for (i = 0; i < nomodule.numDirCalcSrfs; i++) {
					nomodule.nsrflistDirCalc[i + 1] = ivec[i];
				}
			}
			if ((string) inp.GetValue("Lagrange Multipliers") == "True")
				nomodule.Lagrange = 1;
			else
				nomodule.Lagrange = 0;
			if (nomodule.Lagrange == 1) {
				nomodule.numLagrangeSrfs = inp.GetValue(
						"Number of Constrained Surfaces");
				ivec = inp.GetValue("List of Constrained Surfaces");
				for (i = 0; i < MAXSURF + 1; i++)
					nomodule.nsrflistLagrange[i] = 0;
				for (i = 0; i < nomodule.numLagrangeSrfs; i++) {
					nomodule.nsrflistLagrange[i + 1] = ivec[i];
				}
				if ((string) inp.GetValue(
						"Constrained Surface Information From File") == "True")
					nomodule.iLagfile = 1;
				else
					nomodule.iLagfile = 0;
			}
		}
		nomodule.rescontrol = 0;
		if ((string) inp.GetValue("Residual Control") == "True")
			nomodule.rescontrol = 1;
		if (nomodule.rescontrol == 1) {
			nomodule.ResCriteria = inp.GetValue("Residual Criteria");
			nomodule.MinNumIter = inp.GetValue("Minimum Required Iterations");
		}

		if ((string) inp.GetValue("Has masterController.py Control Script") == "True")
		{
			nomodule.hasMasterPythonControlScript = 1;
		}
		else
		{
			nomodule.hasMasterPythonControlScript = 0;	
		}

		nomodule.ideformwall = 0;
		if ((string) inp.GetValue("Deformable Wall") == "True") {
			nomodule.ideformwall = 1;
			nomodule.rhovw = inp.GetValue("Density of Vessel Wall");
			nomodule.rnuvw = inp.GetValue("Poisson Ratio of Vessel Wall");
			nomodule.rshearconstantvw = inp.GetValue(
					"Shear Constant of Vessel Wall");
			nomodule.nProps = inp.GetValue(
					"Number of Wall Properties per Node");

			if ((string) inp.GetValue("Wall Mass Matrix for LHS") == "True")
				nomodule.iwallmassfactor = 1;
			else
				nomodule.iwallmassfactor = 0;
			if ((string) inp.GetValue("Wall Stiffness Matrix for LHS")
					== "True")
				nomodule.iwallstiffactor = 1;
			else
				nomodule.iwallstiffactor = 0;

			// new "boundary element tags"
			if ((string) inp.GetValue("Use Boundary Element Tags") == "True") {
				nomodule.iuseBET = 1;
			}
			else
                nomodule.iuseBET = 0;

			// SWB vessel wall properties
			if ((string) inp.GetValue("Use SWB File") == "True") {
				nomodule.iUseSWB = 1;

				// only use the thickness value in SWB
				// this will allow the use of the wall regions only for stiffness
				nomodule.iUseSWBthickonly = 0;
				if ((string) inp.GetValue("Use SWB Wall Thickness Only") == "True")
					nomodule.iUseSWBthickonly = 1;
			}
			else
				nomodule.iUseSWB = 0;

			// wall regions
			if (nomodule.iUseSWB == 0 || (nomodule.iUseSWB == 1 && nomodule.iUseSWBthickonly == 1 )) {
				if (nomodule.numWallRegions = inp.GetValue("Number of Wall Regions")) {

					cout << "Number of Wall Regions " << nomodule.numWallRegions << endl;

					ivec = inp.GetValue("List of Wall Region Surfaces");
					for (i = 0; i < MAXREGIONS + 1; i++)
						nomodule.nsrflistWallRegions[i] = 0;
					for (i = 0; i < nomodule.numWallRegions; i++) {
						nomodule.nsrflistWallRegions[i + 1] = ivec[i];
					}

					nomodule.nWallETagID = inp.GetValue("Wall Elastic Modulus Region Tag ID");

					vec = inp.GetValue("Wall Elastic Modulus Values");
					for (i = 0; i < MAXREGIONS + 1; i++)
						nomodule.ValueListWallE[i] = 0;
					for (i = 0; i < nomodule.numWallRegions; i++)
						//nomodule.ValueListWallE[ nomodule.nsrflistWallRegions[i + 1] ] = vec[i];
						nomodule.ValueListWallE[ i + 1 ] = vec[i];
					vec.erase(vec.begin(), vec.end());

					nomodule.nWallhTagID = inp.GetValue("Wall Thickness Region Tag ID");

					vec = inp.GetValue("Wall Thickness Values");
					for (i = 0; i < MAXREGIONS + 1; i++)
						nomodule.ValueListWallh[i] = 0;
					for (i = 0; i < nomodule.numWallRegions; i++)
						//nomodule.ValueListWallh[ nomodule.nsrflistWallRegions[i + 1] ] = vec[i];
						nomodule.ValueListWallh[ i + 1 ] = vec[i];
					vec.erase(vec.begin(), vec.end());
				}

				// default values
				nomodule.evw = inp.GetValue("Young Mod of Vessel Wall");
				nomodule.thicknessvw = inp.GetValue("Thickness of Vessel Wall");

			}

			if ((string) inp.GetValue("Wall Damping Term") == "True") {
				nomodule.iwalldamp = 1;
				nomodule.tissSuppDampCoeff = inp.GetValue("Damping Coefficient for Tissue Support");
			} else
				nomodule.iwalldamp = 0;

			if ((string) inp.GetValue("Wall External Support Term") == "True") {
				nomodule.iwallsupp = 1;
				nomodule.tissSuppStiffCoeff = inp.GetValue("Stiffness Coefficient for Tissue Support");
			} else
				nomodule.iwallsupp = 0;

			if ((string) inp.GetValue("Axial Tethering Damping Term") == "True") {
				nomodule.iringdamp = 1;
				nomodule.tissSuppRingDampCoeff = inp.GetValue("Axial Tethering Damping Coefficient");
			} else
				nomodule.iringdamp = 0;

			if ((string) inp.GetValue("Axial Tethering Stiffness Term") == "True") {
				nomodule.iringsupp = 1;
				nomodule.tissSuppRingStiffCoeff = inp.GetValue("Axial Tethering Stiffness Coefficient");
			} else
				nomodule.iringsupp = 0;

			if ((string) inp.GetValue("Measure Distance to Wall Data") == "True") {
				nomodule.imeasdist = 1;
			} else
				nomodule.imeasdist = 0;

			if ((string) inp.GetValue("Wall State Filter Term") == "True") {
				nomodule.imeasdist = 1;
				nomodule.idistancenudge = 1;
				nomodule.stateFilterCoeff = inp.GetValue("Wall State Filter Coefficient");
			} else {
				nomodule.idistancenudge = 0;
			}

			if ((string) inp.GetValue("Use Reference Displacements") == "True") {
				nomodule.iinitialprestress = 1;

				if ((string) inp.GetValue("Update Reference Displacements with Average") == "True")
					nomodule.iupdateprestress = 1;
				else
					nomodule.iupdateprestress = 0;

			}
			else {
				nomodule.iinitialprestress = 0;
				nomodule.iupdateprestress = 0;
			}
		}

		// Scaling Parameters Keywords

		outpar.ro = inp.GetValue("Density");
		outpar.vel = inp.GetValue("Velocity");
		outpar.press = inp.GetValue("Pressure");
		outpar.temper = inp.GetValue("Temperature");
		outpar.entrop = inp.GetValue("Entropy");

		// Step Sequencing

		ivec = inp.GetValue("Step Construction");
		sequence.seqsize = ivec.size();
		if (sequence.seqsize > 100 || sequence.seqsize < 2)
			cerr << "Sequence size must be between 2 and 100 " << endl;

		for (i = 0; i < sequence.seqsize; i++)
			sequence.stepseq[i] = ivec[i];

	} catch (exception &e) {
		cerr << endl << "Input exception: " << e.what() << endl << endl;
		ierr = 001;
		print_error_code(ierr);
		return ierr;
	}

	return ierr;

}

void print_error_code(int ierr) {
  /*
    Return Error codes:
    0xx         Input error
    1xx         Solution Control
    105         Turbulence Model not supported

    2xx         Material Properties

    3xx         Output Control

    4xx         Discretization Control

    5xx         Scaling Parameters

    6xx         Linear Solver Control
    601         linear solver type not supported
  */
  cout << endl << endl << "Input error detected: " << endl << endl;
  if ( ierr == 001 ) {
    cout << endl << "Input Directive not understood" << endl << endl;
  }
  if ( ierr == 105 ) {
    cout << endl << "Turbulence Model Not Supported" << endl << endl;
  }
  if ( ierr == 601 ) {
    cout << endl << "Linear Solver Type Not Supported" << endl << endl;
  }

}
