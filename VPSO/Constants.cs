using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace VPSO
{
    public static class Constants
    {
        public const UInt64 RAND_MAX_KISS = UInt64.MaxValue;
        public const int DMax = 32;			// Max number of dimensions of the search space
        public const int fMax = 6;			// Max number of constraints +1
        public const int RMax = 101;		// Max number of runs
        public const int SMax = 200;		// Max swarm size
        public const double Zero = 1.0e-30;	// To avoid numerical instabilities with some compilers
        public const double Infinity = 1.0e30;
        public const int MMax = 500; // Max size of a "memory" list (cf Lennard-Jones problem)
        public const int NMax = 1000;		// Max number of points when the landscape is read on a file
        public static float[] valueList = new float[] { 1380, 1820, 1930, 2050, 2150 };
        public static int[] pieceNb = new[] { 15, 6, 3, 12, 15 };
        public static float toCut = 5600;
    }
}
