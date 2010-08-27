namespace VPSO
{
    public class Swarm
    {
        public Swarm(int sMax)
        {
            P = new Position[sMax];
            V = new Velocity[sMax];
            X = new Position[sMax];
            S = 0;
            best = 0;
        }
        public int best; 				// Rank of the best particle
        public Position[] P;	// Previous best positions found by each particle
        public int S; 					// Swarm size 
        public Velocity[] V;	// Velocities
        public Position[] X;	// Positions
        //	int worst;				// Rank of the worst particle
    };
}