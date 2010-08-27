using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public class Parameters
    {
        public Parameters()
        {
            w = 0;
            c1 = 0;
            c2 = 0;
            formula = 0;
            S = 0;
            spreadProba = 0;
        }
        public double w; 			// Inertial weight
        public double c1;			// Confidence coefficient
        public double c2;			// Idem		
        public int formula; 		// Formula for stagnation
        public int S;				// Initial swarm size 
        public double spreadProba; // probability threshold
    };
}
