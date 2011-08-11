using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO.Problem
{
    public abstract class ProblemBase
    {
        

        public int Constraint { get; set; }			// Number of constraints
        public double Epsilon { get; set; } 		// Admissible error
        public double EpsilonConstraint { get; set; }  // Admissible error for each constraint
        public int EvaluationMaximum { get; set; } 			// Maximum number of fitness evaluations
        public double ObjectiveValue { get; set; } 	// Objective value					
        public Position Solution { get; set; } // Solution position (if known, just for tests)	
        public SwarmSize SwarmSize { get; set; }			// Search space
        protected ProblemBase()
        {
            Solution = new Position(Constants.DMax);
            SwarmSize = new SwarmSize(Constants.DMax);
            int d;
            int initDiff = 0; // For each function, you can define a specific 
            // initialisation space. In that case, you also have
            // to set initDiff to a non zero value
            int scanNb;
            float z = 0.0f;

            int nAtoms; // For Lennard-Jones problem
            double[] lennard_jones = new[]
                                         {-1, -3, -6, -9.103852, -12.71, -16.505384,-19.821489,-24.113360,-28.422532,
                                          -32.77,-37.97,-44.33,-47.84,-52.32};
            Epsilon = 0.00000;		// Acceptable error. Defalut value
            ObjectiveValue = 0;       // Objective value. Default value
            Constraint = 0; 				// Number of constraints. Default value
            SwarmSize.valueNb = 0;


            // Define the solution point, for test
            // NEEDED when Parameters.stop = 2 
            // i.e. when stop criterion is distance_to_solution < epsilon
            for (d = 0; d < 30; d++)
            {
                Solution.x[d] = 0;
            }
        }

        public abstract Fitness Evaluate(Position x);
    }
}
