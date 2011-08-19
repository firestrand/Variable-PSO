using System;
using VPSO.Problem;

namespace VPSO
{
    public class Position
    {
        public Position(int dMax)
        {
            x = new double[dMax];
            f = new Fitness(dMax);
            size = 0;
        }
        public Fitness f;				// Fitness value  + constraints <=0
        public int size;				// Number of dimensions  D
        public double[] x;		// Coordinates
        public Position Clone()
        {
            var retVal = new Position(x.Length);
            retVal.size = size;
            retVal.f = f.Clone();
            x.CopyTo(retVal.x,0);
            return retVal;
        }

        public static Position Quantis(Position x, SwarmSize swarmSize)
        {
            /*
             Quantisatition of a position
             Only values like x+k*q (k integer) are admissible 
             */
            Position quantx = x.Clone();

            for (int d = 0; d < x.size; d++)
            {
                if (swarmSize.q.Q[d] > Constants.Zero)	// Note that qd can't be < 0
                {
                    quantx.x[d] = swarmSize.q.Q[d] * Math.Floor(0.5 + x.x[d] / swarmSize.q.Q[d]);
                }
            }
            return quantx;
        }

        public static Position ValueAccept(Position x, int valueNb)
        {
            // valueList[] is a global variable (see main.h)
            // Move the position to the nearest acceptable value, 
            // for each dimension independently
            int d;
            int i;

            Position xv;
            xv = x;

            for (d = 0; d < x.size; d++)
            {
                if (x.x[d] <= Constants.valueList[0])
                {
                    xv.x[d] = Constants.valueList[0];
                    continue;
                }

                if (x.x[d] >= Constants.valueList[valueNb - 1])
                {
                    xv.x[d] = Constants.valueList[valueNb - 1];
                    continue;
                }


                for (i = 1; i < valueNb; i++)
                {
                    if (x.x[d] >= Constants.valueList[i - 1] && x.x[d] <= Constants.valueList[i])
                    {
                        if (x.x[d] - Constants.valueList[i - 1] < Constants.valueList[i] - x.x[d])
                            xv.x[d] = Constants.valueList[i - 1];
                        else
                            xv.x[d] = Constants.valueList[i];
                        break;
                    }
                }
            }
            return xv;
        }

        public static Position Discrete(Position pos0, IProblem pb)
        {
            if (pb.SwarmSize.valueNb > 0) // The search space is a list of values
            {
                return ValueAccept(pos0, pb.SwarmSize.valueNb);
            }

            // Quantisation
            Position pos = Quantis(pos0, pb.SwarmSize);
            return pos;

        }
        public static Fitness Constraint(Position x, IProblem pb)
        {
            throw new NotImplementedException();
        }
        public static Fitness Constraint(Position x, int functCode, double epsConstr)
        {
            // ff[0] is defined in perf()
            // Variables specific to Coil compressing spring
            const double fmax = 1000.0;
            const double fp = 300;
            double Cf;
            double K;
            double sp;
            double lf;

            const double S = 189000.0;
            const double lmax = 14.0;
            const double spm = 6.0;
            const double sw = 1.25;
            const double G = 11500000;
            Fitness ff = new Fitness(Constants.DMax) {size = 1};

            switch (functCode)
            {

                case 7:
                    ff.size = 4;

                    ff.f[1] = 0.0193 * x.x[2] - x.x[0];
                    ff.f[2] = 0.00954 * x.x[2] - x.x[1];
                    ff.f[3] = 750 * 1728 - Math.PI * x.x[2] * x.x[2] * (x.x[3] + (4.0 / 3) * x.x[2]);
                    break;

                case 8:
                    ff.size = 5;

                    Cf = 1 + 0.75 * x.x[2] / (x.x[1] - x.x[2]) + 0.615 * x.x[2] / x.x[1];
                    K = 0.125 * G * Math.Pow(x.x[2], 4) / (x.x[0] * x.x[1] * x.x[1] * x.x[1]);
                    sp = fp / K;
                    lf = fmax / K + 1.05 * (x.x[0] + 2) * x.x[2];

                    ff.f[1] = 8 * Cf * fmax * x.x[1] / (Math.PI * x.x[2] * x.x[2] * x.x[2]) - S;
                    ff.f[2] = lf - lmax;
                    ff.f[3] = sp - spm;
                    ff.f[4] = sw - (fmax - fp) / K;
                    break;

                case 15:
                    ff.size = 4;

                    ff.f[1] = Math.Abs(x.x[0] * x.x[0] + x.x[1] * x.x[1] + x.x[2] * x.x[2]
                                       + x.x[3] * x.x[3] + x.x[4] * x.x[4] - 10) - epsConstr; // Constraint h1<=eps
                    ff.f[2] = Math.Abs(x.x[1] * x.x[2] - 5 * x.x[3] * x.x[4]) - epsConstr; // Constraint h2<=eps;
                    ff.f[3] = Math.Abs(Math.Pow(x.x[0], 3) + Math.Pow(x.x[1], 3) + 1) - epsConstr; // Constraint h3<=eps
                    break;

            }

            return ff;
        }
        /// <summary>
        ///  Initialise a position
        /// Note: possible constraints are not checked here. 
        /// The position may be unfeasible	
        /// </summary>
        public static Position Initialize(SwarmSize swarmSize)
        {
            int d;
            var pos = new Position(Constants.DMax) {size = swarmSize.D};

            //  Random uniform
            for (d = 0; d < pos.size; d++)
            {
                pos.x[d] = Alea.NextDouble(swarmSize.min[d], swarmSize.max[d]);
            }
            if (swarmSize.valueNb > 0) // If only some values are acceptable
                pos = ValueAccept(pos, swarmSize.valueNb);

            return pos;
        }
    };
}