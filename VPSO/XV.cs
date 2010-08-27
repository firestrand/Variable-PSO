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

        public static XV confinement(XV xv0,
                              Problem pb)
        {
            // Confinement and evaluation
            // Note: the two are together, for depending on the clamping option
            // 				there may be no evaluation at all
            int d;
            Fitness ff; // Just for constraints
            int i;
            XV xv;

            xv = xv0;

            // Set to the bounds, and v to zero

            for (d = 0; d < pb.SwarmSize.D; d++)
            {
                if (xv.x.x[d] < pb.SwarmSize.minS[d])
                {
                    xv.x.x[d] = pb.SwarmSize.minS[d];
                    xv.v.V[d] = 0;
                }

                if (xv.x.x[d] > pb.SwarmSize.maxS[d])
                {
                    xv.x.x[d] = pb.SwarmSize.maxS[d];
                    xv.v.V[d] = 0;
                }
            }

            xv.x = Position.quantis(xv.x, pb.SwarmSize);

            if (pb.constraint == 0) goto evaluation;

            ff = Position.constraint(xv.x, pb.function, pb.epsConstr);

            for (i = 1; i < ff.size; i++)
            {
                if (ff.f[i] <= 0) continue;

                // If a constraint is not respected
                // the velocity is decreased
                // and the position moved accordingly
                // in th HOPE that the new position will be OK
                // (but we don't check it)
                for (d = 0; d < pb.SwarmSize.D; d++)
                {
                    xv.v.V[d] = xv.v.V[d] / 2;
                    xv.x.x[d] = xv.x.x[d] - xv.v.V[d];
                }
                xv.x = Position.quantis(xv.x, pb.SwarmSize);
            }

            evaluation:
            // Evaluation
            xv.x.f = Problem.perf(xv.x, pb);
            return xv;
        }
    };
}
