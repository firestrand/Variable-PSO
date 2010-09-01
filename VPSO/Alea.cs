using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using VPSO;

namespace VPSO
{
    public class Alea
    {
        static ulong kiss_x = 1;
        static ulong kiss_y = 2;
        static ulong kiss_z = 4;
        static ulong kiss_w = 8;
        static ulong kiss_carry = 0;
        static ulong kiss_k;
        static ulong kiss_m;
        public static void seed_rand_kiss(ulong seed)
        {
            kiss_x = seed | 1;
            kiss_y = seed | 2;
            kiss_z = seed | 4;
            kiss_w = seed | 8;
            kiss_carry = 0;
        }

        public static ulong rand_kiss()
        {
            kiss_x = kiss_x * 69069 + 1;
            kiss_y ^= kiss_y << 13;
            kiss_y ^= kiss_y >> 17;
            kiss_y ^= kiss_y << 5;
            kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
            kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
            kiss_z = kiss_w;
            kiss_w = kiss_m;
            kiss_carry = kiss_k >> 30;
            return kiss_x + kiss_y + kiss_w;
        }
        public static double NextDouble(double a, double b)
        {				// random number (uniform distribution) in  [a b]
            // randOption is a global parameter (see also Parameters.rand)
            double r;

            r = a + (double)rand_kiss() * (b - a) / Constants.RAND_MAX_KISS;

            return r;
        }
        public static int alea_integer(int a, int b)
        {				// Integer random number in [a b]
            int ir;
            double r;

            r = NextDouble(0, 1);
            ir = (int)(a + r * (b + 1 - a));

            if (ir > b) ir = b;
            return ir;
        }

        // ===========================================================
        public static double alea_normal(double mean, double std_dev)
        {
            /*
             Use the polar form of the Box-Muller transformation to obtain a pseudo
             random number from a Gaussian distribution 
             */
            double x1, x2, w, y1;
            // double y2;

            do
            {
                x1 = 2.0 * NextDouble(0, 1) - 1.0;
                x2 = 2.0 * NextDouble(0, 1) - 1.0;
                w = x1 * x1 + x2 * x2;
            }
            while (w >= 1.0);

            w = Math.Sqrt(-2.0 * Math.Log(w) / w);
            y1 = x1 * w;
            // y2 = x2 * w;
            if (NextDouble(0, 1) < 0.5) y1 = -y1;
            y1 = y1 * std_dev + mean;
            return y1;
        }

        // ===========================================================
        public static Velocity alea_sphere(int D, double radius1, double radius2,
                                    int randType, double randGamma)
        {
            /*  ******* Random point in a hypersphere ********
             Maurice Clerc 2003-07-11
            Last update 2010-02-15

            Put  a random point inside the hypersphere S (center 0 
              between radius1 and radius2
  
            */

            int j;
            double length;
            double r;
            Velocity v = new Velocity(D);

            v.Size = D;

            // ----------------------------------- Step 1.  Direction
            length = 0;
            for (j = 0; j < D; j++)
            {
                v.V[j] = alea_normal(0, 1);
                length = length + v.V[j] * v.V[j];
            }

            length = Math.Sqrt(length);

            if (length > 0 && radius2 > 0)
            {
                //----------------------------------- Step 2.   Random radius
                switch (randType)
                {
                    default:
                        r = NextDouble(0, 1);
                        break;

                    case 1:
                        r = Math.Abs(alea_normal(0, randGamma));
                        break;

                }

                for (j = 0; j < D; j++)
                {
                    v.V[j] = r * (radius2 - radius1) * v.V[j] / length;
                }
            }
            else
            {
                for (j = 0; j < D; j++) v.V[j] = 0;
            }
            return v;
        }

        //==================================================
        static int aleaInformant(double[] tab, int S)
        {
            /*
                Define at random an informant, but according to a non uniform distribution
                which is given by the table tab[]
                Note 1: tab doesn't need to be normalised 
            */

            int t;
            double[] tabCumul = new double[Constants.SMax];
            double r;

            tabCumul[0] = tab[0];
            for (t = 1; t < S; t++) tabCumul[t] = tabCumul[t - 1] + tab[t];

            r = NextDouble(0, tabCumul[t - 1]);

            for (t = 0; t < S; t++)
            {
                if (tabCumul[t] >= r) return t;
            }

            return S - 1;
        }

        //==================================================
        public static void aleaIndex(int[] index, int S)
        {
            int[] indexTemp = new int[Constants.SMax];
            int length;
            int rank;
            int s;
            int t;

            length = S;
            for (s = 0; s < S; s++) indexTemp[s] = s; //=index[s];

            for (s = 0; s < S; s++)
            {
                rank = alea_integer(0, length - 1);
                index[s] = indexTemp[rank];
                //printf("\nalea152 s %i rank %i",s,rank);//printf("  ");
                if (rank < length - 1)	// Compact
                {
                    for (t = rank; t < length - 1; t++)
                        indexTemp[t] = indexTemp[t + 1];
                }
                length = length - 1;
            }
        }
    }
}
