using System;

namespace VPSO.Problem
{
    public class Problem
    {
        public Problem()
        {
            Solution = new Position(Constants.DMax);
            SwarmSize = new SwarmSize(Constants.DMax);
        }
        public int Constraint;			// Number of constraints
        public double Epsilon; 		// Admissible error
        public double EpsilonConstraint;  // Admissible error for each constraint
        public int EvaluationMaximum; 			// Maximum number of fitness evaluations
        public int function; 		// Function code
        public double ObjectiveValue; 	// Objective value					
        public Position Solution; // Solution position (if known, just for tests)	
        public SwarmSize SwarmSize;			// Search space

        public static Problem problemDef(int functionCode, string fLandscape)
        {
            int d;
            int initDiff = 0; // For each function, you can define a specific 
            // initialisation space. In that case, you also have
            // to set initDiff to a non zero value
            Problem pb = new Problem(); // Initialised just for my stupid compiler
            int scanNb;
            float z = 0.0f;

            int nAtoms; // For Lennard-Jones problem
            double[] lennard_jones = new[]
                                         {-1, -3, -6, -9.103852, -12.71, -16.505384,-19.821489,-24.113360,-28.422532,
                                          -32.77,-37.97,-44.33,-47.84,-52.32};

            pb.function = functionCode;
            pb.Epsilon = 0.00000;		// Acceptable error. Defalut value
            pb.ObjectiveValue = 0;       // Objective value. Default value
            pb.Constraint = 0; 				// Number of constraints. Default value
            pb.SwarmSize.valueNb = 0;


            // Define the solution point, for test
            // NEEDED when Parameters.stop = 2 
            // i.e. when stop criterion is distance_to_solution < epsilon
            for (d = 0; d < 30; d++)
            {
                pb.Solution.x[d] = 0;
            }


            // ------------------ Search space
            switch (functionCode)
            {
                case 100: // CEC 2005 F1
                    pb.SwarmSize.D = 30;//30; 
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100;
                        pb.SwarmSize.max[d] = 100;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.EvaluationMaximum = pb.SwarmSize.D * 10000;
                    pb.EvaluationMaximum = 1000;
                    pb.Epsilon = 0.000001;	//Acceptable error
                    pb.ObjectiveValue = -450;       // Objective value
                    break;

                case 102:		// Rosenbrock. CEC 2005 F6
                    pb.SwarmSize.D = 10;	// 10

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100; pb.SwarmSize.max[d] = 100;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = pb.SwarmSize.D * 10000;
                    pb.Epsilon = 0.01;	// Acceptable error
                    pb.ObjectiveValue = 390;
                    //pb.EvaluationMaximum=100000; //pb.epsilon=0;
                    break;

                case 103:// CEC 2005 F9, Rastrigin
                    pb.SwarmSize.D = 30;	// 30 
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -5;
                        pb.SwarmSize.max[d] = 5;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.Epsilon = 0; //0.01; // 0.01;	// Acceptable error
                    pb.ObjectiveValue = -330;       // Objective value
                    pb.EvaluationMaximum = pb.SwarmSize.D * 10000;
                    pb.EvaluationMaximum = 100000;
                    break;

                case 104:// CEC 2005 F2  Schwefel
                    pb.SwarmSize.D = 30;
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100;
                        pb.SwarmSize.max[d] = 100;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.Epsilon = 0.00001;	// Acceptable error
                    pb.ObjectiveValue = -450;       // Objective value
                    pb.EvaluationMaximum = pb.SwarmSize.D * 10000;
                    pb.EvaluationMaximum = 1000;
                    break;

                case 105:// CEC 2005 F7  Griewank (NON rotated)
                    pb.SwarmSize.D = 10;	 // 10 
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -600;
                        pb.SwarmSize.max[d] = 600;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.Epsilon = 0.01;	//Acceptable error
                    pb.ObjectiveValue = -180;       // Objective value
                    pb.EvaluationMaximum = pb.SwarmSize.D * 10000;
                    break;

                case 106:// CEC 2005 F8 Ackley (NON rotated)
                    pb.SwarmSize.D = 10;
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -32;
                        pb.SwarmSize.max[d] = 32;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.Epsilon = 0.0001;	// Acceptable error
                    pb.ObjectiveValue = -140;       // Objective value
                    pb.EvaluationMaximum = pb.SwarmSize.D * 10000;

                    break;

                case 0:			// Parabola
                    pb.SwarmSize.D = 30;// 30 	// Dimension

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100; // -100
                        pb.SwarmSize.max[d] = 100;	// 100
                        pb.SwarmSize.q.Q[d] = 0;	// granularity/quantum/step   1 => integer  
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 6000;// 6000	// Max number of evaluations for each run
                    pb.Epsilon = 0.0001;	//0.0001 Acceptable error
                    pb.ObjectiveValue = 0;       // Objective value
                    break;

                case 1:		// Griewank
                    pb.SwarmSize.D = 10;

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100;
                        pb.SwarmSize.max[d] = 100;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 400000;
                    pb.Epsilon = 0.05;	// Acceptable error
                    pb.ObjectiveValue = 0.000;       // Objective value
                    break;

                case 2:		// Rosenbrock
                    pb.SwarmSize.D = 30;

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -30; pb.SwarmSize.max[d] = 30;

                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 300000; // 40000	
                    pb.Epsilon = 0.0001;	//0.0001 Acceptable error
                    pb.ObjectiveValue = 0;       // Objective value
                    break;

                case 3:		// Rastrigin
                    pb.SwarmSize.D = 30;
                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -5.12; //-10; 
                        pb.SwarmSize.max[d] = 5.12; //10; 	  
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.EvaluationMaximum = 200000; // 40000
                    pb.Epsilon = 0.0001;	//0.0001 Acceptable error
                    pb.ObjectiveValue = 0;       // Objective value
                    break;


                case 4:		// Tripod
                    pb.SwarmSize.D = 2;	// 2

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100; // -100
                        pb.SwarmSize.max[d] = 100; // 100
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.Epsilon = 0.0001;
                    pb.EvaluationMaximum = 10000; //10000; 
                    pb.ObjectiveValue = 0; // Objective value
                    break;

                case 5: // Ackley
                    pb.SwarmSize.D = 10;
                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -32; // -32
                        pb.SwarmSize.max[d] = 32; // 32
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 3200;
                    pb.Epsilon = 0.000;
                    pb.ObjectiveValue = 0;
                    break;

                case 6: // Center-bias test function
                    pb.SwarmSize.D = 1;	// Dimension <=30

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -100; // 100
                        pb.SwarmSize.max[d] = 100; // 100
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.minS[d] = -105;
                        pb.SwarmSize.maxS[d] = 105;
                    }

                    pb.EvaluationMaximum = 100; //
                    pb.Epsilon = 0.00;
                    break;

                    // Pressure vessel (confinement method)
                case 7: //  penalty method
                    // Solutions with quantisation 0.0625
                    // (1.125, 0.625, 58.2901, 43.6927) => 7197.729
                    // (1.125, 0.625, 55.8592, 57.7315) => 7197.729
                    // If no granularity => min = 6059.7143
                    pb.Constraint = 3;
                    pb.SwarmSize.D = 4;

                    pb.SwarmSize.min[0] = 1.125; pb.SwarmSize.max[0] = 12.5;
                    pb.SwarmSize.q.Q[0] = 0.0625;
                    pb.SwarmSize.min[1] = 0.625; pb.SwarmSize.max[1] = 12.5;
                    pb.SwarmSize.q.Q[1] = 0.0625;
                    pb.SwarmSize.min[2] = 0.00000001; pb.SwarmSize.max[2] = 240;
                    pb.SwarmSize.q.Q[2] = 0;
                    pb.SwarmSize.min[3] = 0.00000001; pb.SwarmSize.max[3] = 240;
                    pb.SwarmSize.q.Q[3] = 0;

                    /*
                            pb.SwarmSize.min[0] = 0.0625; pb.SwarmSize.max[0] = 99;
                            pb.SwarmSize.q.q[0] = 0.; //0.0625;
                            pb.SwarmSize.min[1] = 0.0625; pb.SwarmSize.max[1] = 99; 
                            pb.SwarmSize.q.q[1] = 0.; //0.0625; 

                            pb.SwarmSize.min[2] = 10; pb.SwarmSize.max[2] = 200; 
                            pb.SwarmSize.q.q[2] = 0;
                            pb.SwarmSize.min[3] = 10; pb.SwarmSize.max[3] = 200; 
                            pb.SwarmSize.q.q[3] = 0;
                    */
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.EvaluationMaximum = 50000;
                    pb.Epsilon = 0.00001; //0.0000000001;	
                    pb.ObjectiveValue = 7197.72893; //7197.7277771; 
                    break;

                case 8: // Compression spring
                    pb.Constraint = 4; // See confin.c
                    pb.SwarmSize.D = 3;

                    pb.SwarmSize.min[0] = 1; pb.SwarmSize.max[0] = 70; pb.SwarmSize.q.Q[0] = 1;
                    pb.SwarmSize.min[1] = 0.6; pb.SwarmSize.max[1] = 3; pb.SwarmSize.q.Q[1] = 0;
                    pb.SwarmSize.min[2] = 0.207; pb.SwarmSize.max[2] = 0.5; pb.SwarmSize.q.Q[2] = 0.001;

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    pb.EvaluationMaximum = 20000;
                    pb.Epsilon = 1.0e-10;
                    pb.ObjectiveValue = 2.6254214578;
                    break;

                case 9: // Gear train

                    pb.SwarmSize.D = 4;

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = 12;
                        pb.SwarmSize.max[d] = 60;
                        pb.SwarmSize.q.Q[d] = 1;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 20000;
                    pb.Epsilon = 1.0e-13;
                    pb.ObjectiveValue = 2.7e-12;
                    break;

                case 10: // Cellular phone
                    pb.SwarmSize.D = 2 * 10; //2*10  2*nb_of_stations

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = 0;
                        pb.SwarmSize.max[d] = 100;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 20000; // 20000
                    pb.Epsilon = 1e-9;
                    pb.ObjectiveValue = 0.005530517; // Best known result (2010-01-03)
                    // pb.epsilon=0; pb.ObjectiveValue=0;
                    break;

                case 11: // PAPR/OFDM
                    pb.SwarmSize.D = 16;
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -1;
                        pb.SwarmSize.max[d] = 1;
                        pb.SwarmSize.q.Q[d] = 1;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 10000;

                    pb.Epsilon = 0.000001;	//0.0001 Acceptable error
                    // pb.ObjectiveValue = 0.01328;  // Objective value for T1
                    pb.ObjectiveValue = 0;
                    break;
                case 12: // Schwefel
                    pb.SwarmSize.D = 30;

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -500; // -32
                        pb.SwarmSize.max[d] = 500; // 32
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 900000;
                    pb.Epsilon = 0.000;
                    pb.ObjectiveValue = 0;
                    break;

                case 13: // Cutting stock  (see also perf()) 
                    pb.SwarmSize.valueNb = 5;
                    pb.SwarmSize.D = 0; for (d = 0; d < pb.SwarmSize.valueNb; d++) pb.SwarmSize.D = pb.SwarmSize.D + Constants.pieceNb[d];
                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = Constants.valueList[0];
                        pb.SwarmSize.max[d] = Constants.valueList[4];
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 100000;
                    pb.Epsilon = 0.00001;
                    pb.ObjectiveValue = 940;

                    break;

                case 14: // Polygon, smallest perimeter
                    pb.SwarmSize.D = 13; // 2*(number of sides)-3 
                    // Polar coordinates
                    // The first point is fixed theta0=0, rho0=0
                    // For the second point, the angle theta1 is 0
                    // So the variables are
                    // rho1
                    // theta2 rho2
                    // theta3 rho3
                    // etc.

                    for (d = 0; d < pb.SwarmSize.D; d = d + 2) // rho
                    {
                        pb.SwarmSize.min[d] = Math.PI / (pb.SwarmSize.D + 3); //0;
                        pb.SwarmSize.max[d] = 1;
                    }

                    for (d = 1; d < pb.SwarmSize.D; d = d + 2) // theta
                    {
                        pb.SwarmSize.min[d] = 0; // 0
                        pb.SwarmSize.max[d] = Math.PI / 2; // 0
                    }


                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 100000;
                    pb.Epsilon = 0.000000;
                    pb.ObjectiveValue = 0;

                    break;


                    // Constrained g13, with Confinement method
                    // Optimum 0.0539498

                case 15: // with Penalty method
                    pb.Constraint = 3;
                    pb.SwarmSize.D = 5;
                    pb.EpsilonConstraint = 0.0001;

                    for (d = 0; d < 2; d++)
                    {
                        pb.SwarmSize.min[d] = -2.3;
                        pb.SwarmSize.max[d] = 2.3;
                    }


                    for (d = 2; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -3.2;
                        pb.SwarmSize.max[d] = 3.2;
                    }

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 100000; //340000 ; 
                    pb.Epsilon = 0.0000001;
                    pb.ObjectiveValue = 0; //0.0539498;
                    break;

                case 16: // G3 (constrained)
                    pb.SwarmSize.D = 10;
                    pb.ObjectiveValue = 0;
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = 0;
                        pb.SwarmSize.max[d] = 1;
                        pb.SwarmSize.q.Q[d] = 0;
                    }
                    pb.EvaluationMaximum = 340000; //340000;
                    pb.ObjectiveValue = 0;
                    pb.Epsilon = 1.0e-6;

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }


                    break;

                case 17: // Lennard-Jones
                    nAtoms = 5; // in {2, ..., 15}
                    pb.SwarmSize.D = 3 * nAtoms; pb.ObjectiveValue = lennard_jones[nAtoms - 2];
                    pb.EvaluationMaximum = 5000 + 3000 * nAtoms * (nAtoms - 1); // Empirical rule
                    pb.Epsilon = 1.0e-6;
                    // Note: with this acceptable error, nAtoms=10 seems to be the maximum
                    //       possible value for a non-null success rate  (5%)
                    // 			 with SPSO 2007

                    //pb.SwarmSize.D=3*21; pb.ObjectiveValue=-81.684;	
                    //pb.SwarmSize.D=3*27; pb.ObjectiveValue=-112.87358;
                    //pb.SwarmSize.D=3*38; pb.ObjectiveValue=-173.928427;

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -2;
                        pb.SwarmSize.max[d] = 2;
                        pb.SwarmSize.q.Q[d] = 0;
                    }

                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }
                    break;
                case 99: // Test

                    pb.SwarmSize.D = 2;	// Dimension

                    // Boundaries
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        pb.SwarmSize.min[d] = -1;
                        pb.SwarmSize.max[d] = 1;
                        pb.SwarmSize.q.Q[d] = 0;
                        pb.SwarmSize.maxS[d] = pb.SwarmSize.max[d];
                        pb.SwarmSize.minS[d] = pb.SwarmSize.min[d];
                    }

                    pb.EvaluationMaximum = 100000;
                    pb.ObjectiveValue = -1.4;
                    pb.Epsilon = 1.0e-6;
                    break;

                case 1000: // Landscape on the file fLandscape.txt
                    // WARNING. Just valid for 1D function
                    pb.SwarmSize.D = 1;	// Dimension
                    //TODO: fix Landscape values			
                    //scanNb=fscanf(fLandscape,"%i",&funct.N);

                    pb.SwarmSize.min[0] = Constants.Infinity;
                    pb.SwarmSize.max[0] = -Constants.Infinity;
                    pb.SwarmSize.q.Q[0] = 1;
                    pb.SwarmSize.maxS[0] = pb.SwarmSize.max[0];
                    pb.SwarmSize.minS[0] = pb.SwarmSize.min[0];

                    for (d = 0; d < Program.funct.N; d++)
                    {
                        //scanNb=fscanf(fLandscape,"%f",&z);funct.x[d]=z;
                        if (z < pb.SwarmSize.min[0]) pb.SwarmSize.min[0] = z;
                        if (z > pb.SwarmSize.max[0]) pb.SwarmSize.max[0] = z;

                        //scanNb=fscanf(fLandscape,"%f",&z);funct.fx[d]=z;
                    }
                    pb.EvaluationMaximum = 100;
                    break;

            } // End of 	switch (pb.function)

            pb.SwarmSize.q.Size = pb.SwarmSize.D;


            if (initDiff == 0) // If no specific initialisation space
            {
                for (d = 0; d < pb.SwarmSize.D; d++)
                {
                    pb.SwarmSize.maxInit[d] = pb.SwarmSize.max[d];
                    pb.SwarmSize.minInit[d] = pb.SwarmSize.min[d];
                }
            }

            return pb;
        }

        public static Fitness perf(Position x, Problem pb)
        {
            return FitnessEvaluator.Evaluate(x, pb);
        }
    };
}