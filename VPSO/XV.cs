using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public class XV // For move and confinement subprograms
    {
        public XV()
        {
            x=new Position(Constants.DMax);
            v=new Velocity(Constants.DMax);
        }

        public Position x;
        public Velocity v;

        public void Confinement(
                              Problem.IProblem pb)
        {
            // Confinement and evaluation
            // Note: the two are together, for depending on the clamping option
            // 				there may be no evaluation at all
            // Set to the bounds, and v to zero

            for (int d = 0; d < pb.SwarmSize.D; d++)
            {
                if (x.x[d] < pb.SwarmSize.minS[d])
                {
                    x.x[d] = pb.SwarmSize.minS[d];
                    v.V[d] = 0;
                }

                if (x.x[d] > pb.SwarmSize.maxS[d])
                {
                    x.x[d] = pb.SwarmSize.maxS[d];
                    v.V[d] = 0;
                }
            }

            if (pb.Constraint != 0)
            {
                Fitness ff = Position.Constraint(x, pb);

                for (int i = 1; i < ff.size; i++)
                {
                    if (ff.f[i] <= 0) continue;

                    // If a constraint is not respected the velocity is decreased
                    // and the position moved accordingly in th HOPE that the new position will be OK
                    // (but we don't check it)
                    for (int d = 0; d < pb.SwarmSize.D; d++)
                    {
                        v.V[d] = v.V[d]/2;
                        x.x[d] = x.x[d] - v.V[d];
                    }
                    
                }
            }
            x = Position.Quantis(x, pb.SwarmSize);

            // Evaluation
            x.f = pb.Evaluate(x);
        }
    };
}
