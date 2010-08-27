using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public class SwarmSize
    {
        public SwarmSize(int dMax)
        {
            D = 0;
            max = new double[dMax];
            maxInit = new double[dMax];
            min = new double[dMax];
            minInit = new double[dMax];
            q = new Quantum(dMax);
            maxS = new double[dMax];
            minS = new double[dMax];
            valueNb = 0;
        }
        public int D;
        public double[] max;  	// Initialisation, and feasible space
        public double[] maxInit;
        public double[] min;
        public double[] minInit;
        public Quantum q;		// Quantisation step size. 0 => continuous problem
        public double[] maxS;	// Search space
        public double[] minS;
        public int valueNb; 		 	// if >0, the search space is a list of "valueNb"values.
    };
}
