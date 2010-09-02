using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using VPSO;

namespace VPSO
{
    public class Alea
    {
        static readonly Random rand = new Random(1);
        //static ulong kiss_x = 1;
        //static ulong kiss_y = 2;
        //static ulong kiss_z = 4;
        //static ulong kiss_w = 8;
        //static ulong kiss_carry = 0;
        //static ulong kiss_k;
        //static ulong kiss_m;
        //public static void seed_rand_kiss(ulong seed)
        //{
        //    kiss_x = seed | 1;
        //    kiss_y = seed | 2;
        //    kiss_z = seed | 4;
        //    kiss_w = seed | 8;
        //    kiss_carry = 0;
        //}

        //public static ulong rand_kiss()
        //{
        //    kiss_x = kiss_x * 69069 + 1;
        //    kiss_y ^= kiss_y << 13;
        //    kiss_y ^= kiss_y >> 17;
        //    kiss_y ^= kiss_y << 5;
        //    kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
        //    kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
        //    kiss_z = kiss_w;
        //    kiss_w = kiss_m;
        //    kiss_carry = kiss_k >> 30;
        //    return kiss_x + kiss_y + kiss_w;
        //}
        public static double NextDouble(double a, double b)
        {				// random number (uniform distribution) in  [a b]
            // randOption is a global parameter (see also Parameters.rand)
            //double r;

            //r = a + (double)rand_kiss() * (b - a) / Constants.RAND_MAX_KISS;

            return a + rand.NextDouble() * (b - a);//r;
        }
        public static int NextInteger(int a, int b)
        {				// Integer random number in [a b]
            //int ir;
            //double r;

            //r = NextDouble(0, 1);
            //ir = (int)(a + r * (b + 1 - a));

            //if (ir > b) ir = b;
            return rand.Next(a,b);
        }

        // ===========================================================
        public static double NextNormal(double mean, double stdDev)
        {
            /*
             Use the polar form of the Box-Muller transformation to obtain a pseudo
             random number from a Gaussian distribution 
             */
            double x1, w, y1;
            // double y2;

            do
            {
                x1 = 2.0 * NextDouble(0, 1) - 1.0;
                double x2 = 2.0 * NextDouble(0, 1) - 1.0;
                w = x1 * x1 + x2 * x2;
            }
            while (w >= 1.0);

            w = Math.Sqrt(-2.0 * Math.Log(w) / w);
            // y2 = x2 * w;
            if (NextDouble(0, 1) < 0.5)
                y1 = -x1 * w;
            else
                y1 = x1 * w;
            y1 = y1 * stdDev + mean;
            return y1;
        }
        public static void Shuffle(int[] index, int s)
        {
            var indexTemp = new int[Constants.SMax];

            int length = s;
            for (int i = 0; i < s; i++) 
                indexTemp[i] = i; //=index[s];

            for (int i = 0; i < s; i++)
            {
                int rank = NextInteger(0, length - 1);
                index[i] = indexTemp[rank];
                //printf("\nalea152 s %i rank %i",s,rank);//printf("  ");
                if (rank < length - 1)	// Compact
                {
                    for (int t = rank; t < length - 1; t++)
                        indexTemp[t] = indexTemp[t + 1];
                }
                length = length - 1;
            }
        }
    }
}
