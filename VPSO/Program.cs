using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using VPSO.Problem;

namespace VPSO
{
    public class Program
    {
        // Global variables
        public static int dLBest; 			// For information 
        public static Landscape funct; // When (x,f(x)) is read from a file
        public static int iter;
        public static Archive memPos = new Archive(Constants.MMax);
        public static double nEval;			// Number of fitness evaluations
        public static int run;


        // For Cutting stock
        // http://en.wikipedia.org/wiki/Cutting_stock_problem

        static void Main(string[] args)
        {
            Position bestBest = new Position(Constants.DMax);
            double convRateMean;
            int d;					// Current dimension
            double error;				// Current error
            double errorMean;			// Average error
            double errorMin = Constants.Infinity;	// Best result over all runs
            double[] errorMeanBest = new double[Constants.RMax];
            double evalMean;			// Mean number of evaluations
            int functionCode;

            int n;
            int nFailure;				// Number of unsuccessful runs
            double logProgressMean = 0;
            Parameters parameters = new Parameters();
            ProblemBase pb;
            int runMax;
            Result result;
            double successRate;
            double variance;


            // ----------------------------------------------- PROBLEM
            functionCode = 102; //TODO: Something seems wrong with Rosenbrockf6
            /* (see problemDef, cellular_phone, cec2005pb for precise definitions)
             0 Parabola (Sphere)
             1 Griewank
             2 Rosenbrock (Banana)
             3 Rastrigin
             4 Tripod (dimension 2)
             5 Ackley
             6 Center-bias test function
             7 Pressure vessel (penalty method)
             1007        "        (confinement method)
             8 Compression spring
             1008				"				 (confinement method)
             9 Gear train
             10 Cellular phone

             12 Schwefel
             13 Cutting stock
             14 convex small polygon, 4 edges, smallest perimeter
             15 constrained problem g13 (penalty method)
             16 constrained problem g3 (penalty method)
             17 Lennard-Jones problem
             1015 				"               (confinement method)

             CEC 2005 benchmark  (no more than 30D. See cec2005data.c)
             100 F1 (shifted Parabola/Sphere) 
             102 F6 (shifted Rosenbrock) 
             103 F9 (shifted Rastrigin) 
             104 F2 Schwefel 
             105 F7 Griewank  (NOT rotated)
             106 F8 Ackley  (NOT rotated) 

             99 Test

             */

            runMax = 100; // Numbers of runs


            // =========================================================== 

            pb = new CEC2005F1Circle();//Problem.Problem.problemDef(functionCode, null);//fLandscape); fLandscape is for when the landscape is read from a file.

            /* -----------------------------------------------------
             PARAMETERS and OPTIONS 
            */
            parameters.formula = 1; // Which formula to use to compute the number of
            // iterations between two checks. See SpreadIter()
            parameters.spreadProba = 1; //1; 
            // We want that the probability for a particle to be
            // informed my any other one after T time steps
            // is at least spreadProba. Then T is automatically
            // computed.
            // Particular case: set to zero => T=1
            // Suggested values:
            // formula   spreadProba
            //	1						1.0
            //  2, 3				0.9
            // with PSO_B, formula=1 and spredProba=0.55	

            // -------------------------------------------------------------AUTOMATIC VALUES
            // You may modify them, if you really want to ...
            parameters.S = (int)(40 + 2 * Math.Sqrt(pb.SwarmSize.D));	// Initial swarm size 
            //Parameters.S=pb.SwarmSize.D+1;

            // According to Clerc's Stagnation Analysis
            parameters.w = 1 / (2 * Math.Log((double)2)); // 0.721
            parameters.c1 = 0.5 + Math.Log((double)2); // 1.193
            parameters.c2 = parameters.c1;


            //---------------------------------------------------------- WARNINGS and ERRORS

            if (runMax > Constants.RMax)  // Maximum number of runs
            {
                runMax = Constants.RMax;
                Console.WriteLine("WARNING. No more than {0} runs", Constants.RMax);
            }

            if (parameters.S > Constants.SMax)
            {
                parameters.S = Constants.SMax;
                Console.WriteLine("WARNING. The swarm size can not be bigger than {0}", Constants.SMax);
            }

            if (parameters.spreadProba > 1)
            {
                throw new ArgumentException("Parameters.spreadProba must be smaller than 1");
            }


            // Display information about the problem and the parameters
            // Some information


            //---------------
            convRateMean = 0;
            errorMean = 0;
            evalMean = 0;
            nFailure = 0;

            //------------------------------------- RUNS	
            for (run = 0; run < runMax; run++)
            {
                //printf("\n run %i/%i",run+1,runMax); 
                // For some options, positions must be memorised as much a possible
                //fprintf(f_swarm,"run %i\n",run);

                result = PSO(parameters, pb, 0);
                error = result.error.errorFC();

                convRateMean = convRateMean + result.convRate;

                if (error > pb.Epsilon) // Failure
                {
                    nFailure = nFailure + 1;
                }

                // Result display
                //printf("\nmean %1.2e", errorMean/(run+1)); // Temporary error mean
                //printf(", %0.0f %% ", (double)100*(run+1-nFailure)/(run+1)); // Temporary success rate
                //printf (" %i/%i. Eval %f. Error %f ", run+1,runMax, result.nEval, error);
                //printf("  x(0)= %f",(double)result.SW.P[result.SW.best].x[0]);

                // Save result
                //fprintf( f_run, "\n%i %f %e ", run, result.nEval,  error );
                //for ( d = 0; d < pb.SwarmSize.D; d++ ) 
                //    fprintf( f_run, " %0.20f",  (double)result.SW.P[result.SW.best].x[d] );


                // Compute/store some statistical information
                if (run == 0)
                {
                    errorMin = error;
                    bestBest = result.SW.P[result.SW.best];
                }
                else if (error < errorMin)
                {
                    errorMin = error;
                    bestBest = result.SW.P[result.SW.best];
                }

                evalMean = evalMean + result.nEval;
                errorMean = errorMean + error;
                errorMeanBest[run] = error;
                logProgressMean = logProgressMean - Math.Log(error);

            }		// End loop on "run"


            // Display statistical information
            // Means
            convRateMean = convRateMean / (double)runMax;
            evalMean = evalMean / (double)runMax;
            errorMean = errorMean / (double)runMax;
            logProgressMean = logProgressMean / (double)runMax;

            Console.WriteLine("Conv. rate (mean)= {0}", convRateMean);
            Console.WriteLine("Eval. (mean)= {0}", evalMean);
            Console.WriteLine("Error (mean) = {0}, {1}", errorMean, errorMean);

            // Variance
            variance = 0;

            for (run = 0; run < runMax; run++)
                variance = variance + Math.Pow(errorMeanBest[run] - errorMean, 2);

            variance = Math.Sqrt(variance / runMax);
            Console.WriteLine("Std. dev. {0}", variance);
            Console.WriteLine("Log_progress (mean) = {0}", logProgressMean);

            // Success rate and minimum value
            Console.WriteLine("Failure(s) {0}", nFailure);
            successRate = 1 - nFailure / (double)runMax;
            Console.WriteLine("Success rate = {0}%", 100 * successRate);

            //if (run > 1)
            Console.WriteLine("Best min value = {0}", bestBest.f.f[0]);
            if (pb.Constraint > 0)
            {
                Console.WriteLine("{0} constraints (should be negative): ", pb.Constraint);
                for (n = 0; n < pb.Constraint; n++) Console.WriteLine("{0} ", bestBest.f.f[n + 1]);
            }

            // Display the best position
            Console.WriteLine("Position: ");
            for (d = 0; d < pb.SwarmSize.D; d++) Console.WriteLine(" {0}", (double)bestBest.x[d]);

            // Save
            /*
             fprintf(f_synth,"\n"); for (d=0;d<SwarmSize.D;d++) fprintf(f_synth,"%f ",
                 pb.offset[d]);
             fprintf(f_synth,"    %f %f %f %.0f%% %f",errorMean,variance,errorMin,
                 successRate,evalMean); 

             fprintf( f_synth, "\n%f %f %f %f %.0f%% %f ", shift,
                 errorMean, variance, errorMin, successRate, evalMean );
             */
            //fprintf (f_synth, "%f %f %.0f%% %f   ",
            //   errorMean, variance, 100*successRate, evalMean);		   

            // Save the best result
            //fprintf(f_synth,"\n ");
            //for(n=0;n<pb.constraint+1;n++) fprintf( f_synth, "%f  ", bestBest.f.f[n]);
            //fprintf(f_synth," X");
            //for ( d = 0; d < pb.SwarmSize.D; d++ ) fprintf( f_synth, " %f", (double) bestBest.x[d] );


            //-----------------------------------------------------
            // Repeat some information
            Console.WriteLine("Function {0}", functionCode); Console.WriteLine("Dimension {0}", pb.SwarmSize.D);
            Console.ReadLine();
            return; // End of main program
        }


        // ===========================================================

        //============================================================= COEFF_SC
        static double coeff_SC(int D)
        {

            //  The D-cube whose edge is a
            // and  the D-sphere whose radius is r=a*coeff_S_C
            // have the same volume
            double coeff;
            double c1 = 0, c2;
            int d;
            double d3;
            double x;

            int option = 0; // 0 => exact value
            // 1 => mean value

            if (D == 1) return 0.5; // The radius is just the half intervall

            x = (double)D;
            if ((2 * (int)(D / 2) == D) || option == 1) // D even
            {
                d3 = 1; for (d = 2; d <= D / 2; d++) d3 = d3 * d; // (D/2)!
                c1 = Math.Pow(d3, 1 / x) / Math.Sqrt(Math.PI);
                if (option == 0) return c1;
            }

            // D odd
            {
                d3 = 1; for (d = 2; d <= D; d++) d3 = d3 * d; // D!
                c2 = 0.5 * Math.Pow(d3, 1 / x) / Math.Pow(Math.PI, 0.5 - 0.5 / x);
                if (option == 0) return c2;
            }

            coeff = (c1 + c2) / 2;

            return coeff;
        }
        //================================================== KISS
        /*
         A good pseudo-random numbers generator

         The idea is to use simple, fast, individually promising
         generators to get a composite that will be fast, easy to code
         have a very long period and pass all the tests put to it.
         The three components of KISS are
         x(n)=a*x(n-1)+1 mod 2^32
         y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
         z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
         The y's are a shift register sequence on 32bit binary vectors
         period 2^32-1;
         The z's are a simple multiply-with-carry sequence with period
         2^63+2^32-1.  The period of KISS is thus
         2^32*(2^32-1)*(2^63+2^32-1) > 2^127
         */






        //================================================
        static int compareDoubles(double a, double b)
        {
            return a < b ? -1 : a > b ? 1 : 0;//TODO: Factor this out
        }


        // ===============================================
        static int betterThan(Fitness f1, Fitness f2)
        {
            int n, n1, n2;

            if (f1.size == 1) goto fitness; // No constraints (except the search space)

            // Criterion "Number of respected constraints
            n1 = 0;
            n2 = 0;

            for (n = 1; n < f1.size; n++)
            {
                if (f1.f[n] < 0) n1 = n1 + 1;
                if (f2.f[n] < 0) n2 = n2 + 1;
            }

            if (n1 > n2) return 1;
            if (n1 < n2) return 0;

            // Here, n1=n2

            fitness:
            // Criterion "Total Fitness"
            if (f1.errorFC() < f2.errorFC() - Constants.Zero) return 1;
            return 0;

        }

        //============================================================
        static double potentiel(Position X, Archive list)
        {
            double dist;
            int n;
            double pot = 0;

            for (n = 0; n < list.Size; n++)
            {
                dist = Tools.distanceL(X, list.M[n], 2);
                if (dist > Constants.Zero) pot = pot + 1.0 / dist;
                else return Constants.Infinity;
            }

            return pot;

        }
        static XV move(Result R, int s, int g, IProblem pb,
                    Parameters parameters)
        {
            double c1, c2;
            int d;
            Velocity GX = new Velocity(Constants.DMax);
            Velocity PX = new Velocity(Constants.DMax);

            XV xv = new XV();
            double w;

            xv.x.size = pb.SwarmSize.D;
            xv.v.Size = pb.SwarmSize.D;
            PX.Size = pb.SwarmSize.D;
            GX.Size = pb.SwarmSize.D;


            w = parameters.w;
            c1 = parameters.c1;
            c2 = parameters.c2;

            // Update the velocity (Exploration tendency)
            for (d = 0; d < pb.SwarmSize.D; d++)
            {
                xv.v.V[d] = w * R.SW.V[s].V[d];
            }

            // --------------------------------------------------------NORMAL PSO STRATEGY

            // Prepare Exploitation tendency
            for (d = 0; d < pb.SwarmSize.D; d++) // p-x
            {
                PX.V[d] = R.SW.P[s].x[d] - R.SW.X[s].x[d];
            }

            if (g != s)
            {
                for (d = 0; d < pb.SwarmSize.D; d++) // g-x
                {
                    GX.V[d] = R.SW.P[g].x[d] - R.SW.X[s].x[d];
                }
            }

            // Complete the velocity
            for (d = 0; d < pb.SwarmSize.D; d++)
            {
                xv.v.V[d] = xv.v.V[d] + Alea.NextDouble(0, c1) * PX.V[d];

            }

            if (g != s) // If the best neighbour is not the particle itself
            {
                for (d = 0; d < pb.SwarmSize.D; d++)
                {
                    xv.v.V[d] = xv.v.V[d] + Alea.NextDouble(0, c2) * GX.V[d];
                }
            }

            // Update the position
            for (d = 0; d < pb.SwarmSize.D; d++)
            {
                xv.x.x[d] = R.SW.X[s].x[d] + xv.v.V[d];
            }

            // NOTE that the position is not yet evaluated
            return xv;
        }

        //================================================


        //=================================================
        public static double lennard_jones(Position x)
        {
            /*
                This is for black-box optimisation. Therefore, we are not supposed to know
                that there are some symmetries. That is why the dimension of the problem is
                3*nb_of_atoms, as it could be 3*nb_of_atoms-6
            */
            int d;
            int dim = 3;
            double dist;
            double f;
            int i, j;
            int nPoints = x.size / dim;
            Position x1 = new Position(Constants.DMax);
            Position x2 = new Position(Constants.DMax);
            double zz;

            x1.size = dim; x2.size = dim;

            f = 0;
            for (i = 0; i < nPoints - 1; i++)
            {
                for (d = 0; d < dim; d++) x1.x[d] = x.x[3 * i + d];
                for (j = i + 1; j < nPoints; j++)
                {
                    for (d = 0; d < dim; d++) x2.x[d] = x.x[3 * j + d];

                    dist = Tools.distanceL(x1, x2, 2);

                    zz = Math.Pow(dist, -6);
                    f = f + zz * (zz - 1);
                }
            }
            f = 4 * f;
            //printf("\n %f",f);
            return f;
        }


        // ===========================================================
        static double normL(Velocity v, double L)
        {   // L-norm of a vector
            int d;
            double n;

            n = 0;

            for (d = 0; d < v.Size; d++)
                n = n + Math.Pow(Math.Abs(v.V[d]), L);

            n = Math.Pow(n, 1 / L);
            return n;
        }


        // ===========================================================

        static Result PSO(Parameters parameters, IProblem pb, int level)
        {
            int added; // For information 
            Fitness error = new Fitness(Constants.fMax);
            Fitness errorInit = new Fitness(Constants.fMax);  // Just for information
            Fitness errorPrev = new Fitness(Constants.fMax);
            double errorTot;

            int g;
            int sBest; // Rank in g of the best of the bests
            int improvTot; // Number of particles that have improved their previous best
            int[] index = new int[Constants.SMax];
            int initLinks;	// Flag to (re)init or not the information links
            int initLinkNb;
            int iterBegin;
            int[,] LINKS = new int[Constants.SMax, Constants.SMax];	// Information links
            int m;
            int moved;
            int n;
            int noStop;
            Result R = new Result();
            int removed; // For information

            int s0 = 0;
            int s, s1, s2;
            int spread;
            int stagnation = 0;
            int swarmMod;
            int sWorst;

            XV xvNorm = new XV();

            // -----------------------------------------------------
            // INITIALISATION

            R.SW.S = parameters.S; // Initial size of the swarm
            memPos.Rank = 0; 					// Rank (in M) where to memorise a new position
            memPos.Size = 0; 					// Number of memorised positions

            // Positions
            for (s = 0; s < R.SW.S; s++)
            {
                R.SW.X[s] = Position.Initialize(pb.SwarmSize);
                memPos.memSave(R.SW.X[s]); 		// Save the position
            }

            // Velocities
            for (s = 0; s < R.SW.S; s++)
            {
                R.SW.V[s] = Velocity.Initialize(R.SW.X[s], pb.SwarmSize);
            }

            // Discrete values
            // Note: may be removed if you are sure that the initialisation
            //        takes discretisation into account (or if there is none)
            for (s = 0; s < R.SW.S; s++)
            {
                R.SW.X[s] = Position.Discrete(R.SW.X[s], pb);
            }


            // Note: at this point no confinement is needed:
            // initialisation is supposed to be OK from this point of view
            // (but some constraints may be not respected)


            // First evaluations
            for (s = 0; s < R.SW.S; s++)
            {
                R.SW.X[s].f = pb.Evaluate(R.SW.X[s]);
                R.SW.P[s] = R.SW.X[s];	// Best position = current one
            }

            // Save the positions
            for (s = 0; s < R.SW.S; s++) memPos.memSave(R.SW.X[s]);

            // Find the best 
            R.SW.best = best(R.SW);
            error = R.SW.P[R.SW.best].f;

            // Display the best
            Console.WriteLine("Best value after init. {0} ", R.SW.P[R.SW.best].f.f[0]);
            if (pb.Constraint > 0)
            {
                Console.WriteLine("Constraints (should be < 0) ");
                for (n = 0; n < pb.Constraint; n++) Console.WriteLine("{0} ", error.f[n + 1]);
            }

            //fprintf(f_run,"\n Best value after init. %f ", errorPrev );
            //printf( "\n Position :\n" );
            //for ( d = 0; d < pb.SwarmSize.D; d++ ) printf( " %f", R.SW.P[R.SW.best].x[d] );

            initLinks = 1;		// So that information links will beinitialized
            initLinkNb = 0; // Count the number of iterations between two reinit of the links
            iter = 0;
            nEval = 0;
            noStop = 0;
            added = 0; removed = 0; // For information
            spread = spreadIter(parameters.spreadProba, R.SW.S, parameters.formula); // Number of iterations
            // needed to "spread" the information
            errorInit = error; // For information

            // ---------------------------------------------- ITERATIONS
            while (noStop == 0)
            {

                //printf("\niter %i",iter);
                //fprintf(f_run,"\niter %i",iter); 
                iter = iter + 1;
                errorPrev = error;
                Alea.Shuffle(index, R.SW.S); // Random numbering of the particles

                if (initLinks == 1)	// Bidirectional ring topology. Randomly built
                {
                    initLinks = 0;
                    initLinkNb = 0; // Count the number of iterations since the last re-init
                    // of the links

                    // Init to zero (no link)
                    for (s = 0; s < R.SW.S; s++)
                    {
                        for (m = 0; m < R.SW.S; m++) LINKS[m, s] = 0;
                    }

                    // Information links (bidirectional ring)
                    for (s = 0; s < R.SW.S - 1; s++)
                    {
                        LINKS[index[s], index[s + 1]] = 1;
                        LINKS[index[s + 1], index[s]] = 1;
                    }

                    LINKS[index[0], index[R.SW.S - 1]] = 1;
                    LINKS[index[R.SW.S - 1], index[0]] = 1;

                    // Each particle informs itself
                    for (m = 0; m < R.SW.S; m++) LINKS[m, m] = 1;
                }

                // Loop on particles, for move
                improvTot = 0;

                for (s0 = 0; s0 < R.SW.S; s0++)
                {
                    s = index[s0];

                    // Find the best informant
                    g = s;
                    for (m = 0; m < R.SW.S; m++)
                    {
                        if (m == s) continue;
                        if (LINKS[m, s] == 1 && betterThan(R.SW.P[m].f, R.SW.P[g].f) == 1)
                            g = m;
                    }

                    // Move	
                    xvNorm = move(R, s, g, pb, parameters);
                    xvNorm.x = Position.Discrete(xvNorm.x, pb);

                    // Confinement and evaluation
                    xvNorm.Confinement(pb);

                    // New position and new velocity
                    R.SW.X[s] = xvNorm.x;
                    R.SW.V[s] = xvNorm.v;

                    // Update the best previous position
                    if (betterThan(R.SW.X[s].f, R.SW.P[s].f) == 1) // Improvement of the previous best
                    {
                        R.SW.P[s] = R.SW.X[s].Clone();
                        improvTot++; // Increase the number of improvements during this iteration

                        // Memorise the improved position
                        memPos.memSave(R.SW.P[s]);

                        // Update the best of the bests
                        if (betterThan(R.SW.P[s].f, R.SW.P[R.SW.best].f) == 1) R.SW.best = s;
                    }

                    // Decide to stop or not
                    errorTot = R.SW.P[R.SW.best].f.errorFC();

                    if (errorTot > pb.Epsilon && nEval < pb.EvaluationMaximum)
                        noStop = 0;	// Won't stop
                    else // Failure
                    {
                        noStop = 1;	// Will stop
                        goto end;
                    }

                }			// End of "for (s0=0 ...  "	= end of the iteration (move)

                /*-------------------------------------------------- Adaptations
                 Rule 1:
                 Check every "spread" iterations after each re-init of the links
                 If no improvement of the global best
                 => reinit links before the next iteration

                 Rule 2:
                 if no improvement of the global best during "spread" iterations
                 => Try to add a particle (and initialise it in a non-searched area)
                 => re-init links before the next iteration
                 Note that the condition is slightly different from the one of Rule 1

                 Rule 3:
                 if "enough" local improvements during the iteration
                 => try to remove a particle (keep at least D+1 ones)

                 */

                // Rule 1 - Re-initializing the information links
                // Check if improvement since last re-init of the links
                initLinkNb = initLinkNb + 1; // Number of iterations since the last check 

                if (initLinkNb >= spread) // It's time to check
                {
                    initLinkNb = 0; // Reset to zero the number of iterations since the last check
                    // The swarm size may have been modified, so must be "spread"
                    spread = spreadIter(parameters.spreadProba, R.SW.S, parameters.formula);

                    if (betterThan(error, errorPrev) == 1)	// Improvement
                        initLinks = 0;	 // No need of structural adaptation								
                    else			// No improvement			
                        initLinks = 1;	// Information links will be	reinitialized	
                }
                else initLinks = 0;  // To early, no need to check 

                sWorst = worst(R.SW); // Rank of the worst particle, before any adaptation

                // Rule 2 - Adding a particle
                // Check global stagnation (improvement of the global best) 
                if (betterThan(R.SW.P[R.SW.best].f, errorPrev) == 1) stagnation = 0;	// Improvement
                else stagnation++; // No improvement during this iteration

                swarmMod = 0; // Information flag

                if (stagnation >= spread)  // Too many iterations without global improvement
                // =>  add a particle
                {
                    if (R.SW.S < Constants.SMax) // if not too many particles
                    {
                        s = R.SW.S;
                        R.SW.X[s] = memPos.InitializeFar(pb); // Init in a non-searched area
                        R.SW.X[s] = Position.Discrete(xvNorm.x, pb); // If discrete search space
                        R.SW.X[s].f = pb.Evaluate(R.SW.X[s]);	 // Evaluation			
                        R.SW.V[s] = Velocity.Initialize(R.SW.X[s], pb.SwarmSize); // Init velocity						
                        R.SW.P[s] = R.SW.X[s].Clone(); // Previous best = current position								
                        R.SW.S = R.SW.S + 1; // Increase the swarm size

                        //fprintf(f_swarm,"%i %i  %f\n",iter, R.SW.S,error.f[0]);
                        // Count the number of added particles (for information)
                        added++;
                        initLinks = 1; // Links will be reinitialised
                        stagnation = 0; // Reset the count for stagnation
                        swarmMod = 1; // A particle has been added
                        //printf("\n iter %i, added %i, S %i, spread %i",iter, added, R.SW.S, spread);
                    }
                }

                // Rule 3 - Removing a particle
                // If enough improvements of some particles, remove the worst 
                // (but keep at least D+1 particles)
                // NOTE: this is "the worst" without taking into account the particle
                // that has (possibly) been added
                // NOTE: it is perfectly possible to have a particle added
                // (because of no improvement of the global best)  AND
                // a particle removed (because enough _local_ improvements)

                if (R.SW.S > pb.SwarmSize.D + 1 && improvTot > 0.5 * R.SW.S)
                {
                    if ((swarmMod == 0 && sWorst < R.SW.S - 1) || swarmMod == 1)
                    // if the worst is not the last 
                    {
                        R.SW.P[sWorst] = R.SW.P[R.SW.S - 1]; // ... replace it by the last
                        R.SW.V[sWorst] = R.SW.V[R.SW.S - 1];
                        R.SW.X[sWorst] = R.SW.X[R.SW.S - 1];

                        // Compact the matrix of the links
                        for (s1 = 0; s1 < R.SW.S; s1++)  // For each line, compact the columns
                            for (s2 = sWorst; s2 < R.SW.S - 1; s2++) LINKS[s1, s2] = LINKS[s1, s2 + 1];

                        for (s2 = 0; s2 < R.SW.S - 1; s2++)	// For each column, compact the lines
                            for (s1 = sWorst; s1 < R.SW.S - 1; s1++) LINKS[s1, s2] = LINKS[s1 + 1, s2];
                    }
                    R.SW.S = R.SW.S - 1; // Decrease the swarm size
                    if (s < R.SW.best) R.SW.best = R.SW.best - 1; // The rank of the best may 
                    // have been modified
                    // Count the number of removed particles (for information)
                    removed++;
                    swarmMod = -1; // A particle has been remowed
                    //printf("\n iter %i, removed %i, S %i, spread %i",iter, removed, R.SW.S, spread);
                }

                // Save on a the result of the iteration on a file
                if (swarmMod != 0)
                {
                    //fprintf(f_swarm,"%i %i  %f\n",iter, R.SW.S,R.SW.P[R.SW.best].f.f[0]);	
                }

                // End of the iteration 
            end: ;

            } // End of "while (noStop==0)"


            // Convergence rate, just for information
            R.convRate = (errorInit.f[0] - R.SW.P[R.SW.best].f.f[0]) / errorInit.f[0];

            // Information about the evolution of the swarm size
            Console.WriteLine("{0} iterations, +{1} -{2} particles", iter, added, removed);

            // Final number of evaluations
            R.nEval = nEval;
            // Final fitness
            R.error = R.SW.P[R.SW.best].f;
            return R;
        }
        //==============================================================================
        static int best(Swarm SW)
        {
            // Find the rank of the best position
            // 	Remember that f is Math.Abs(fitness-ObjectiveValue)
            // 	We want to minimise it
            int s = 0;
            int best = 0;

            for (s = 1; s < SW.S; s++)
            {
                if (betterThan(SW.P[s].f, SW.P[best].f) == 1)
                    best = s;
            }
            return best;
        }
        //==============================================================================
        static int worst(Swarm SW)
        {
            // Find the rank of the worst position
            int s = 0;
            int worst = 0;

            for (s = 1; s < SW.S; s++)
            {
                if (betterThan(SW.P[worst].f, SW.P[s].f) == 1)
                    worst = s;
            }
            return worst;
        }
        //==============================================================================
        //==============================================================================
        static int spreadIter(double spreadProba, int S, int formula)
        {
            // Number of iterations to spread information
            switch (formula)
            {
                default: // 1
                    return (int)Math.Ceiling(0.5 + spreadProba * 0.5 * S);
                case 2:
                    return (int)Math.Ceiling(0.5 + Math.Log(1 - spreadProba) / Math.Log(1 - 2.0 / S));
                case 3:
                    return (int)Math.Ceiling(0.5 + Math.Log(1 - spreadProba) / Math.Log(1 - 3.0 / S));
            }
        }

    }
}