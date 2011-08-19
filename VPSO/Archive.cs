using System;
using VPSO.Problem;

namespace VPSO
{
    public class Archive
    {
        public Archive(int mMax)
        {
            M = new Position[mMax];
            Rank = 0;
            Size = 0;
        }
        public Position[] M; // Memorise positions
        public int Rank;
        public int Size;

        public void memSave(Position P)
        {
            // Save a position
            // Is useful to generate a new particle in a promising area
            // The memPos list is a global variable
            M[Rank] = P;

            if (Size < M.Length - 1)
            {
                Size = Size + 1;
                Rank = Rank + 1;
            }
            else Rank = 0; // We re-use the memory cyclically 
        }

        public Position InitializeFar(IProblem pb)
        {
            // Try to find a new position that is "far" from all the memorised ones
            //Note: memPos is a global variable
            double[] coord = new double[Constants.MMax];
            double[] interv = new double[2];

            var xFar = new Position(Constants.DMax) { size = pb.SwarmSize.D };

            for (int d = 0; d < pb.SwarmSize.D; d++) // For each dimension
            {
                for (int n = 0; n < this.Size; n++) coord[n] = this.M[n].x[d]; // All the coordinates on
                // this dimension

                Array.Sort(coord); // Sort them
                // by increasing order

                // Find the biggest intervall
                interv[0] = coord[0];
                interv[1] = coord[1];
                double delta = interv[1] - interv[0];

                for (int n = 1; n < this.Size - 1; n++)
                {
                    if (coord[n + 1] - coord[n] < delta) continue;

                    interv[0] = coord[n];
                    interv[1] = coord[n + 1];
                    delta = interv[1] - interv[0];
                }



                // Particular case, xMax
                if (pb.SwarmSize.max[d] - coord[this.Size - 1] > delta)
                {
                    xFar.x[d] = pb.SwarmSize.max[d];
                    delta = pb.SwarmSize.max[d] - coord[this.Size - 1];
                }

                // Particular case, xMin
                if (coord[0] - pb.SwarmSize.min[d] > delta)
                    xFar.x[d] = pb.SwarmSize.min[d];
                else
                    xFar.x[d] = 0.5 * (interv[1] + interv[0]);// Take the middle
            }

            xFar = Position.Discrete(xFar, pb);
            xFar.f = pb.Evaluate(xFar);
            return xFar;

        }
    };
}