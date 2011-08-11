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
            var offset0 = new[]{ 
                                    -3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
                                    -8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000, 
                                    -1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
                                    6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001, 
                                    3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001, 
                                    -6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
                                };
            double f = -450;
            for (int d = 0; d < x.size; d++)
            {
                double xd = x.x[d];
                x.x[d] = xd - offset0[d];
                f = f + xd * xd;
            }
            var ff = new Fitness(Constants.fMax);
            ff.f[0] = Math.Abs(f - ObjectiveValue);
            return ff;
        }
    }
}
