using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO.Problem
{
    public class CEC2005F1Circle : ProblemBase
    {
        public CEC2005F1Circle()
        {
            SwarmSize.D = 30;//30; 
            for (int d = 0; d < SwarmSize.D; d++)
            {
                SwarmSize.min[d] = -100;
                SwarmSize.max[d] = 100;
                SwarmSize.q.Q[d] = 0;
                SwarmSize.maxS[d] = SwarmSize.max[d];
                SwarmSize.minS[d] = SwarmSize.min[d];
                SwarmSize.maxInit[d] = SwarmSize.max[d];
                SwarmSize.minInit[d] = SwarmSize.min[d];
            }
            EvaluationMaximum = SwarmSize.D * 10000;
            EvaluationMaximum = 1000;
            Epsilon = 0.000001;	//Acceptable error
            ObjectiveValue = -450;       // Objective value

            SwarmSize.q.Size = SwarmSize.D;
        }
        public override Fitness Evaluate(Position x)
        {
            double value = FitnessEvaluator.EvaluateSphere(x);
            var fitness = new Fitness(Constants.fMax);
            fitness.f[0] = Math.Abs(value - ObjectiveValue);
            Program.nEval = Program.nEval + 1;
            return fitness;
        }
    }
}
