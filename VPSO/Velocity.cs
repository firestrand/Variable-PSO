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
    };
}
