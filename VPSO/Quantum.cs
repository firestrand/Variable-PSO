using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public class Quantum
    {
        public Quantum(int dMax)
        {
            Q = new double[dMax];
            Size = 0;
        }
        public double[] Q;
        public int Size;
    };
}
