namespace VPSO
{
    public class Result
    {
        public Result()
        {
            SW=new Swarm(Constants.SMax);
            error=new Fitness(Constants.fMax);
        }
        public double nEval; 			// Number of evaluations  
        public Swarm SW;			// Final swarm
        public Fitness error;				// Numerical result of the run
        public double convRate; 			// For information: how fast was the convergence
    };
}