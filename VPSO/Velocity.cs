using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public class Velocity
    {
        public Velocity(int dMax)
        {
            V = new double[dMax];
            Size = 0;
        }
        public int Size;				// Number of dimensions D  
        public double[] V; 		// Components

        public static Velocity Initialize(Position pos, SwarmSize swarmSize)
        {
            int d;

            var vel = new Velocity(Constants.DMax) {Size = pos.size};

            // Half-diff  0.5*(alea-x)

            for (d = 0; d < vel.Size; d++)
            {
                vel.V[d] = (Alea.NextDouble(swarmSize.min[d], swarmSize.max[d]) - pos.x[d]) / 2;
            }

            return vel;
        }
    };
}
