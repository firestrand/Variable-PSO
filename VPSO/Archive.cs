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
    };
}