using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public class Landscape
    {
        public Landscape(int nMax)
        {
            x = new double[nMax];
            fx = new double[nMax];
            N = 0;
        }
        public int N;
        public double[] x;//[NMax];
        public double[] fx;//[NMax];
    };
}
