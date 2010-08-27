using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public static class Tools
    {
        public static double distanceL(Position x1, Position x2, double L)
        {  // Distance between two positions
            // L = 2 => Euclidean	
            int d;
            double n;

            n = 0;

            for (d = 0; d < x1.size; d++)
                n = n + Math.Pow(Math.Abs(x1.x[d] - x2.x[d]), L);

            n = Math.Pow(n, 1 / L);
            return n;
        }

        public static int sign(double x)
        {    // Be careful: not the classical one 
            if (x <= 0) return -1;
            return 1;
        }
    }
}
