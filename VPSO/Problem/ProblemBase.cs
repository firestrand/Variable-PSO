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

            Epsilon = 0.00000;		// Acceptable error. Defalut value
            ObjectiveValue = 0;       // Objective value. Default value
            Constraint = 0; 				// Number of constraints. Default value
            SwarmSize.valueNb = 0;

            // Define the solution point, for test
            // NEEDED when Parameters.stop = 2 
            // i.e. when stop criterion is distance_to_solution < epsilon
            for (int d = 0; d < 30; d++)
            {
                Solution.x[d] = 0;
            }

            
        }

        public abstract Fitness Evaluate(Position x);
    }
}
