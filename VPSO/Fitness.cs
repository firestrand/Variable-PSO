namespace VPSO
{
    public class Fitness
    {
        public Fitness(int fMax)
        {
            f = new double[fMax];
            size = 0;
        }
        public int size;
        public double[] f;
        public Fitness Clone()
        {
            var retVal = new Fitness(f.Length);
            f.CopyTo(retVal.f,0);
            retVal.size = size;
            return retVal;
        }
        public double errorFC()
        {
            // Error, including the constraints
            double error = f[0];

            if (size == 1) return error;

            for (int n = 1; n < size; n++)
                if (f[n] > 0) error = error + f[n];

            return error;

        }
    };
}