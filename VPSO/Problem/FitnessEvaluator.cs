using System;

namespace VPSO.Problem
{
    /// <summary>Computes benchmark fitness independently of problem configuration.</summary>
    internal static class FitnessEvaluator
    {
        private static readonly double[] offset_0 = new[]
                                    {
                                        -3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
                                        -8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000,
                                        -1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
                                        6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001,
                                        3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001,
                                        -6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
                                    };
        private static readonly double[] offset_2 = new[]
                                    {
                                        8.1023200e+001, -4.8395000e+001,  1.9231600e+001, -2.5231000e+000,  7.0433800e+001,
                                        4.7177400e+001, -7.8358000e+000, -8.6669300e+001,  5.7853200e+001, -9.9533000e+000,
                                        2.0777800e+001,  5.2548600e+001,  7.5926300e+001,  4.2877300e+001, -5.8272000e+001,
                                        -1.6972800e+001,  7.8384500e+001,  7.5042700e+001, -1.6151300e+001,  7.0856900e+001,
                                        -7.9579500e+001, -2.6483700e+001,  5.6369900e+001, -8.8224900e+001, -6.4999600e+001,
                                        -5.3502200e+001, -5.4230000e+001,  1.8682600e+001, -4.1006100e+001, -5.4213400e+001
                                    };
        private static readonly double[] offset_3 = new[]
                                    {
                                        1.9005000e+000, -1.5644000e+000, -9.7880000e-001, -2.2536000e+000,  2.4990000e+000,
                                        -3.2853000e+000,  9.7590000e-001, -3.6661000e+000,  9.8500000e-002, -3.2465000e+000,
                                        3.8060000e+000, -2.6834000e+000, -1.3701000e+000,  4.1821000e+000,  2.4856000e+000,
                                        -4.2237000e+000,  3.3653000e+000,  2.1532000e+000, -3.0929000e+000,  4.3105000e+000,
                                        -2.9861000e+000,  3.4936000e+000, -2.7289000e+000, -4.1266000e+000, -2.5900000e+000,
                                        1.3124000e+000, -1.7990000e+000, -1.1890000e+000, -1.0530000e-001, -3.1074000e+000
                                    };
        private static readonly double[] offset_4 = new[]
                                    {
                                        3.5626700e+001, -8.2912300e+001, -1.0642300e+001, -8.3581500e+001,  8.3155200e+001,
                                        4.7048000e+001, -8.9435900e+001, -2.7421900e+001,  7.6144800e+001, -3.9059500e+001,
                                        4.8885700e+001, -3.9828000e+000, -7.1924300e+001,  6.4194700e+001, -4.7733800e+001,
                                        -5.9896000e+000 ,-2.6282800e+001, -5.9181100e+001,  1.4602800e+001, -8.5478000e+001,
                                        -5.0490100e+001,  9.2400000e-001,  3.2397800e+001,  3.0238800e+001, -8.5094900e+001,
                                        6.0119700e+001, -3.6218300e+001, -8.5883000e+000, -5.1971000e+000,  8.1553100e+001
                                    };
        private static readonly double[] offset_5 = new[]
                                    {
                                        -2.7626840e+002, -1.1911000e+001, -5.7878840e+002, -2.8764860e+002, -8.4385800e+001,
                                        -2.2867530e+002, -4.5815160e+002, -2.0221450e+002, -1.0586420e+002, -9.6489800e+001,
                                        -3.9574680e+002, -5.7294980e+002, -2.7036410e+002, -5.6685430e+002, -1.5242040e+002,
                                        -5.8838190e+002, -2.8288920e+002, -4.8888650e+002, -3.4698170e+002, -4.5304470e+002,
                                        -5.0658570e+002, -4.7599870e+002, -3.6204920e+002, -2.3323670e+002, -4.9198640e+002,
                                        -5.4408980e+002, -7.3445600e+001, -5.2690110e+002, -5.0225610e+002, -5.3723530e+002
                                    };
        private static readonly double[] offset_6 = new[]
                                    {
                                        -1.6823000e+001,  1.4976900e+001,  6.1690000e+000,  9.5566000e+000,  1.9541700e+001,
                                        -1.7190000e+001, -1.8824800e+001,  8.5110000e-001, -1.5116200e+001,  1.0793400e+001,
                                        7.4091000e+000,  8.6171000e+000, -1.6564100e+001, -6.6800000e+000,  1.4543300e+001,
                                        7.0454000e+000, -1.8621500e+001,  1.4556100e+001, -1.1594200e+001, -1.9153100e+001,
                                        -4.7372000e+000,  9.2590000e-001,  1.3241200e+001, -5.2947000e+000,  1.8416000e+000,
                                        4.5618000e+000, -1.8890500e+001,  9.8008000e+000, -1.5426500e+001,  1.2722000e+000
                                    };

        // Preserve the legacy sphere formula and its in-place position shift.
        internal static double EvaluateSphere(Position position)
        {
            double value = -450;
            for (int d = 0; d < position.size; d++)
            {
                double coordinate = position.x[d];
                position.x[d] = coordinate - offset_0[d];
                value = value + coordinate * coordinate;
            }
            return value;
        }

        internal static Fitness Evaluate(Position x, Problem pb)
        {				// Evaluate the fitness value for the particle of rank s

            double c;
            int d;
            double DD;
            double dTheta;
            double dx1, dx2;
            int i, j, k;
            int n;
            double f = 0;
            Fitness ff = new Fitness(Constants.fMax);
            int grid;
            double p;
            double s11, s12, s21, s22;
            double sum1, sum2;
            double t0, tt, t1;
            double x1, x2, x3, x4;
            double xd;
            Position xs;
            double[,] xy = new double[8, 2];
            double y, y1, y2;
            double z, z2;
            double[,] zz = new double[8, 8];

            // Shifted Parabola/Sphere (CEC 2005 benchmark)

            // Shifted Rosenbrock (CEC 2005 benchmark)

            // Shifted Rastrigin (CEC 2005)


            // Shifted Schwefel (F2 CEC 2005)


            // Shifted Griewank (CEC 2005)


            // Shifted Ackley (CEC 2005)


            // Data for center-bias test function (code 6)
            double radius = 10; //
            int shift = 1;  // Set to zero in order to have the center on {0, ...,0}, to 1 else
            double[] center = new[]
                                  {
                                      80, -20,  1.9231600e+001, -2.5231000e+000,  7.0433800e+001,
                                      4.7177400e+001, -7.8358000e+000, -8.6669300e+001,  5.7853200e+001, -9.9533000e+000,
                                      2.0777800e+001,  5.2548600e+001,  7.5926300e+001,  4.2877300e+001, -5.8272000e+001,
                                      -1.6972800e+001,  7.8384500e+001,  7.5042700e+001, -1.6151300e+001,  7.0856900e+001,
                                      -7.9579500e+001, -2.6483700e+001,  5.6369900e+001, -8.8224900e+001, -6.4999600e+001,
                                      -5.3502200e+001, -5.4230000e+001,  1.8682600e+001, -4.1006100e+001, -5.4213400e+001
                                  };
            Position summit = new Position(Constants.DMax);

            Program.nEval = Program.nEval + 1; // Number of fitness evaluations (global variable)

            //---------------------------------------------------

            xs = x;
            ff.size = 1; // Default value (no constraint)

            switch (pb.function)
            {
                default:
                    throw new ArgumentException("Wrong code function");
                    break;
                    // CEC 2005 benchmark
                case 100: // Parabola (Sphere) CEC 2005
                    f = EvaluateSphere(xs);
                    break;

                case 102:  // Rosenbrock

                    for (d = 0; d < xs.size; d++)
                    {
                        xs.x[d] = xs.x[d] - offset_2[d];
                    }

                    f = 390;
                    for (d = 1; d < xs.size; d++)
                    {
                        tt = xs.x[d - 1] - 1;
                        f = f + tt * tt;

                        tt = xs.x[d - 1] * xs.x[d - 1] - xs.x[d];
                        f = f + 100 * tt * tt;
                    }
                    break;

                case 103: // Rastrigin
                    for (d = 0; d < xs.size; d++)
                    {
                        xs.x[d] = xs.x[d] - offset_3[d];
                    }
                    f = -330;
                    k = 10;

                    for (d = 0; d < xs.size; d++)
                    {
                        xd = xs.x[d];
                        f = f + xd * xd - k * Math.Cos(2 * Math.PI * xd);
                    }
                    f = f + xs.size * k;
                    break;

                case 104: // Schwefel (F2)
                    for (d = 0; d < xs.size; d++)
                    {
                        xs.x[d] = xs.x[d] - offset_4[d];
                    }

                    f = -450;
                    for (d = 0; d < xs.size; d++)
                    {
                        sum2 = 0.0;
                        for (k = 0; k <= d; k++)
                        {
                            sum2 += xs.x[k];
                        }
                        f += sum2 * sum2;
                    }
                    break;


                case 105: // Griewank. WARNING: in the CEC 2005 benchmark it is rotated
                    sum1 = 0.0;
                    sum2 = 1.0;
                    f = -180;
                    for (d = 0; d < xs.size; d++)
                    {
                        xd = xs.x[d] - offset_5[d];
                        sum1 += xd * xd;
                        sum2 *= Math.Cos(xd / Math.Sqrt(1.0 + d));
                    }
                    f = f + 1.0 + sum1 / 4000.0 - sum2;
                    break;

                case 106: // Ackley
                    f = -140;

                    sum1 = 0.0;
                    sum2 = 0.0;
                    for (d = 0; d < xs.size; d++)
                    {
                        xd = xs.x[d] - offset_6[d];
                        sum1 += xd * xd;
                        sum2 += Math.Cos(2.0 * Math.PI * xd);
                    }
                    sum1 = -0.2 * Math.Sqrt(sum1 / xs.size);
                    sum2 /= xs.size;
                    f = f + 20.0 + Math.E - 20.0 * Math.Exp(sum1) - Math.Exp(sum2);
                    break;


                case 0:		// Parabola (Sphere)
                    f = 0.0;
                    for (d = 0; d < xs.size; d++)
                    {
                        xd = xs.x[d];
                        //xd = xs.x[d]-50; // Shifted

                        f = f + xd * xd;
                    }
                    //printf("\n2442  xs.size %i f %f x0 %f",xs.size,f,xs.x[0]);
                    break;

                case 1:		// Griewank
                    f = 0;
                    p = 1;

                    for (d = 0; d < xs.size; d++)
                    {
                        xd = xs.x[d];
                        //xd=xd-150;
                        f = f + xd * xd;
                        p = p * Math.Cos(xd / Math.Sqrt((double)(d + 1)));
                    }
                    f = f / 4000 - p + 1;
                    break;

                case 2:		// Rosenbrock
                    f = 0;
                    t0 = xs.x[0] + 1;	// Solution on (0,...0)
                    for (d = 1; d < xs.size; d++)
                    {

                        t1 = xs.x[d] + 1;
                        tt = 1 - t0;
                        f += tt * tt;
                        tt = t1 - t0 * t0;
                        f += 100 * tt * tt;
                        t0 = t1;
                    }
                    break;

                case 3:		// Rastrigin
                    f = 0;
                    k = 10;

                    for (d = 0; d < xs.size; d++)
                    {
                        xd = xs.x[d];
                        f = f + xd * xd - k * Math.Cos(2 * Math.PI * xd);
                    }
                    f = f + xs.size * k;
                    break;

                case 4:		// 2D Tripod function
                    // Note that there is a big discontinuity right on the solution
                    // point.
                    x1 = xs.x[0];
                    x2 = xs.x[1];
                    s11 = (1.0 - Tools.sign(x1)) / 2;
                    s12 = (1.0 + Tools.sign(x1)) / 2;
                    s21 = (1.0 - Tools.sign(x2)) / 2;
                    s22 = (1.0 + Tools.sign(x2)) / 2;

                    //f = s21 * (Math.Abs(x1) - x2); // Solution on (0,0)
                    f = s21 * (Math.Abs(x1) + Math.Abs(x2 + 50)); // Solution on (0,-50)
                    f = f + s22 * (s11 * (1 + Math.Abs(x1 + 50) + Math.Abs(x2 - 50))
                                   + s12 * (2 + Math.Abs(x1 - 50) + Math.Abs(x2 - 50)));
                    break;

                case 5:  // Ackley
                    sum1 = 0;
                    sum2 = 0;
                    DD = x.size;
                    for (d = 0; d < x.size; d++)
                    {
                        xd = xs.x[d];
                        sum1 = sum1 + xd * xd;
                        sum2 = sum2 + Math.Cos(2 * Math.PI * xd);
                    }
                    f = -20 * Math.Exp(-0.2 * Math.Sqrt(sum1 / DD)) - Math.Exp(sum2 / DD) + 20 + Math.Exp(1);
                    break;

                case 6:  // Center-bias test function
                    // on [-100,100]^2
                    for (d = 0; d < x.size; d++) summit.x[d] = shift * center[d];

                    z = Tools.distanceL(x, summit, 2);
                    if (z > radius) f = radius; else f = z;
                    break;

                case 7: // Pressure vessel (penalty method)
                case 1007: // confinement method
                    // Ref New Optim. Tech. in Eng. p 638
                    // D = 4
                    x1 = xs.x[0]; // [1.1,12.5] granularity 0.0625
                    x2 = xs.x[1];// [0.6,12.5] granularity 0.0625
                    x3 = xs.x[2]; // [0,240]
                    x4 = xs.x[3];// [0,240]

                    f = 0.6224 * x1 * x3 * x4 + 1.7781 * x2 * x3 * x3 + 3.1611 * x1 * x1 * x4 + 19.84 * x1 * x1 * x3;

                    ff = Position.Constraint(xs, pb.function, pb.EpsilonConstraint);

                    if (pb.Constraint == 0)// Constraints, by penalty method
                    {
                        if (ff.f[1] > 0) { c = 1 + Math.Pow(10, 10) * ff.f[1]; f = f * c * c; }
                        if (ff.f[2] > 0) { c = 1 + ff.f[2]; f = f * c * c; }
                        if (ff.f[3] > 0) { c = 1 + ff.f[3]; f = f * c * c; }
                    }
                    break;

                case 8: // Coil compression spring  (penalty method)
                    // Ref New Optim. Tech. in Eng. p 644
                case 1008: // confinement method

                    x1 = xs.x[0]; // {1,2, ... 70}
                    x2 = xs.x[1];//[0.6, 3]
                    x3 = xs.x[2];// relaxed form [0.207,0.5]  dx=0.001
                    // In the original problem, it is a list of
                    // acceptable values
                    // {0.207,0.225,0.244,0.263,0.283,0.307,0.331,0.362,0.394,0.4375,0.5}

                    f = Math.PI * Math.PI * x2 * x3 * x3 * (x1 + 2) * 0.25;
                    // Constraints
                    ff = Position.Constraint(xs, pb.function, pb.EpsilonConstraint);
                    if (pb.Constraint == 0)
                    {
                        if (ff.f[1] > 0) { c = 1 + ff.f[1]; f = f * c * c * c; }
                        if (ff.f[2] > 0) { c = 1 + ff.f[1]; f = f * c * c * c; }
                        if (ff.f[3] > 0) { c = 1 + ff.f[3]; f = f * c * c * c; }
                        if (ff.f[4] > 0) { c = 1 + Math.Pow(10, 10) * ff.f[4]; f = f * c * c * c; }
                        if (ff.f[5] > 0) { c = 1 + Math.Pow(10, 10) * ff.f[5]; f = f * c * c * c; }
                    }
                    break;

                case 9: // Gear train
                    x1 = xs.x[0]; // {12,13, ... 60}
                    x2 = xs.x[1];// {12,13, ... 60}
                    x3 = xs.x[2];// {12,13, ... 60}
                    x4 = xs.x[3];// {12,13, ... 60}

                    f = 1.0 / 6.931 - x1 * x2 / (x3 * x4);
                    f = f * f;
                    break;

                    // case 10
                case 10: // Cellular phone
                    // Grid 100 x 100. You may modify it for more or less precision
                    // (which implies, of course, more or less computation time)
                    grid = 100;
                    f = 0;

                    dx1 = (pb.SwarmSize.max[0] - pb.SwarmSize.min[0]) / grid;

                    for (i = 0; i <= grid; i++)
                    {
                        y1 = pb.SwarmSize.min[0] + i * dx1;
                        dx2 = (pb.SwarmSize.max[1] - pb.SwarmSize.min[1]) / grid;
                        for (j = 0; j <= grid; j++)
                        {
                            y2 = pb.SwarmSize.min[1] + j * dx2;
                            z = 0; // Max field
                            for (d = 0; d < xs.size; d = d + 2) // For each station
                            {
                                x1 = xs.x[d];
                                x2 = xs.x[d + 1];

                                z2 = 1.0 / ((x1 - i) * (x1 - y1) + (x2 - j) * (x2 - y2) + 1);
                                if (z2 > z) z = z2;
                            }
                            f = f + z;
                        }
                    }
                    f = 1.0 / f; // In order to have something to minimise
                    break;


                case 12: // Schwefel
                    f = 418.98288727243369 * pb.SwarmSize.D;
                    for (d = 0; d < pb.SwarmSize.D; d++)
                    {
                        xd = xs.x[d];
                        f = f - xd * Math.Sin(Math.Sqrt(Math.Abs(xd)));
                    }
                    break;
                case 13: // Cutting problem
                    sum1 = 0;
                    f = 0;
                    for (d = 0; d < pb.SwarmSize.D; d++) // Waste
                    {
                        sum1 = sum1 + xs.x[d];
                        if (sum1 < Constants.toCut && d < pb.SwarmSize.D - 1) continue;
                        if (sum1 < Constants.toCut && d == pb.SwarmSize.D - 1)
                        {
                            f = f + Constants.toCut - sum1;
                            continue;
                        }
                        f = f + Constants.toCut - sum1 + xs.x[d];
                        sum1 = xs.x[d];
                    }

                    //Constraints
                    for (i = 0; i < pb.SwarmSize.valueNb; i++) // For each value ...
                    {
                        n = 0;
                        for (d = 0; d < pb.SwarmSize.D; d++) // ... count how many times it appears
                        {
                            if (xs.x[d] >= Constants.valueList[i] - Constants.Zero && xs.x[d] <= Constants.valueList[i] + Constants.Zero)
                                n = n + 1;
                            //printf("\n 2732  %f %f %f",valueList[i]-zero,xs.x[d], valueList[i]+zero);
                        }
                        if (n != Constants.pieceNb[i]) // If not exactly the right number of times => penalty
                        {
                            f = f + Constants.toCut; // f+100000
                            //printf("\n 2716, %i, %i, should be %i",i,n,pieceNb[i]);
                        }
                    }
                    //printf("\n2740 f %f",f);
                    break;

                case 14: // Polygon. longest perimeter
                    // Cartesian coordinates
                    xy[0, 0] = 0; xy[0, 1] = 0;
                    xy[1, 0] = xs.x[0]; xy[1, 1] = 0;
                    i = 2;
                    for (d = 1; d < pb.SwarmSize.D; d = d + 2)
                    {
                        xy[i, 0] = xs.x[d + 1] * Math.Cos(xs.x[d]); // rho*Math.Cos(theta)
                        xy[i, 1] = xs.x[d + 1] * Math.Sin(xs.x[d]); // rho*sin(theta)
                        i = i + 1;
                    }

                    // Distances
                    n = i;
                    z2 = 0; // Distance max
                    y = 0; // Sum of the distances
                    for (i = 1; i < n; i++)
                    {
                        for (j = 0; j < i; j++)
                        {
                            zz[i, j] = Math.Sqrt(Math.Pow(xy[i, 0] - xy[j, 0], 2) + Math.Pow(xy[i, 1] - xy[j, 1], 2));
                            y = y + zz[i, j];
                            if (zz[i, j] > z2) z2 = zz[i, j];
                        }
                    }

                    // Perimeter
                    //	f=zz[n-1][0];
                    //	for(i=1;i<n;i++) f=f+zz[i][i-1];

                    // To minimise
                    //f=2+4*sin(pi/12)-f;
                    f = Math.PI - f;

                    // Max sum of the distances y
                    f = n * (n - 1) - y; //to minimise

                    // Constraints -----
                    p = (n - 2) * Math.PI; // penalty
                    //dTheta=(n-2)*pi/(2*n);
                    dTheta = Constants.Zero;
                    // Diameter = 1
                    if (z2 > 1 + Constants.Zero) f = f + p;
                    if (z2 < 1 - Constants.Zero) f = f + p;

                    // Convex
                    for (d = 1; d < pb.SwarmSize.D - 2; d = d + 2) // Increasing theta
                    {
                        //if(xs.x[d]>xs.x[d+2]) f=f+p;
                        if (xs.x[d + 2] < xs.x[d] + dTheta) f = f + p;
                    }

                    // Angle max
                    if (xs.x[pb.SwarmSize.D - 2] > Math.PI) f = f + p;

                    // Sides
                    for (d = 2; d < pb.SwarmSize.D - 2; d = d + 2)
                    {
                        if (xs.x[d] < xs.x[d - 2] && xs.x[d] < xs.x[d + 2]) f = f + p;
                    }

                    //if(xs.x[0]<0.2) f=f+p;
                    /*
                     for(d=2;d<pb.SwarmSize.D;d=d+2) //
                     {
                         if(xs.x[d]<xs.x[0]) f=f+p;
                     }
                     */
                    break;

                case 15: // Constrained g13, Penalty method
                case 1015: // 			Confinement method
                    // D = 5
                    // x1, x2 in [-2.3, 2.3]
                    // x3, x4, x5 [-3.2, 3.2]

                    f = Math.Exp(xs.x[0] * xs.x[1] * xs.x[2] * xs.x[3] * xs.x[4]);
                    ff = Position.Constraint(xs, pb.function, pb.EpsilonConstraint);

                    if (pb.Constraint == 0) // Penalty method
                    {
                        if (ff.f[1] > 0) f = f + ff.f[1];
                        if (ff.f[2] > 0) f = f + ff.f[2];
                        if (ff.f[3] > 0) f = f + ff.f[3];
                    }

                    break;

                case 16: // G3 (constrained)
                    // min =0 on (1/Math.Sqrt(D), ...)
                    f = 1;
                    sum1 = 0;
                    for (d = 0; d < x.size; d++)
                    {
                        xd = xs.x[d];
                        f = f * xd;
                        sum1 = sum1 + xd * xd;
                    }
                    f = Math.Abs(1 - Math.Pow(x.size, x.size / 2.0) * f) + x.size * Math.Abs(sum1 - 1);
                    break;

                case 17: // Lennard-Jones
                    f = Program.lennard_jones(xs);
                    break;

                case 99: // Test
                    x1 = xs.x[0];
                    x2 = xs.x[1];
                    xd = x1 * x1 + x2 * x2;
                    f = xd - Math.Cos(18 * x1) - Math.Cos(18 * x2);
                    f = f + Math.Abs(xd - 0.6); // One constraint (penalty method);

                    break;

                    // 1D Test "deceptive" on [0 1]
                    xd = xs.x[0];
                    if (xd < 0.6) f = xd;
                    else f = 0.6 - (0.5 / 0.4) * (xd - 0.6);
                    break;

                    // 1D test  [-50, 50]
                    xd = xs.x[0];
                    f = xd * (xd + 1) * Math.Cos(xd);
                    break;

                    // KrishnaKumar
                    f = 0;
                    for (d = 0; d < xs.size - 1; d++)
                    {
                        f = f + Math.Sin(xs.x[d] + xs.x[d + 1]) + Math.Sin(2 * xs.x[d] * xs.x[d + 1] / 3);
                    }
                    break;

                    xd = xs.x[0];
                    //if (xd<9) f=10-xd; else f=10*xd-90;
                    if (xd <= 1) f = 10 * xd; else f = 11 - xd;
                    break;

                case 1000: //1-D  Landscape on a file
                    // Find the nearest x
                    xd = xs.x[0];
                    z = Math.Abs(xd - Program.funct.x[0]);
                    f = Program.funct.fx[0];
                    for (d = 1; d < Program.funct.N; d++)
                    {
                        z2 = Math.Abs(xd - Program.funct.x[d]);
                        if (z2 < z)
                        {
                            z = z2;
                            f = Program.funct.fx[d];
                        }
                    }
                    break;

            }

            ff.f[0] = Math.Abs(f - pb.ObjectiveValue);
            return ff;

        }
    }
}
