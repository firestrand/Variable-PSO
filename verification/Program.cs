extern alias baseline;
using System;
using System.Collections.Generic;
using System.IO;
using Old = baseline::VPSO;
using Current = VPSO;

class Runner
{
    static void Main(string[] args)
    {
        Study.PublicApi(typeof(Old.Problem.Problem).Assembly, typeof(Current.Problem.Problem).Assembly);
        var random = new Random(1937);
        int checks = 0;
        var measurements = new List<object>();
        foreach (int function in new[] { 100, 102, 103, 104, 105, 106 })
        {
            var a = new Old.Position(30) { size = 30 };
            var b = new Current.Position(30) { size = 30 };
            var oldProblem = new Old.Problem.Problem { function = function };
            var problem = new Current.Problem.Problem { function = function };
            var input = new double[30];
            for (int sample = 0; sample < 1000; sample++)
            {
                for (int d = 0; d < 30; d++) input[d] = sample == 0 ? 0 : sample == 1 ? 100 : sample == 2 ? -100 : (random.NextDouble() * 200 - 100);
                Array.Copy(input, a.x, 30); Array.Copy(input, b.x, 30);
                var expected = Old.Problem.Problem.perf(a, oldProblem);
                var actual = Current.Problem.Problem.perf(b, problem);
                if (expected.size != actual.size) throw new Exception("Fitness size changed");
                for (int d = 0; d < expected.f.Length; d++) { Study.Equal(expected.f[d], actual.f[d]); checks++; }
                for (int d = 0; d < 30; d++) { Study.Equal(a.x[d], b.x[d]); checks++; }
            }
            measurements.Add(Study.Measure($"Function {function}: constant-table reuse and evaluator extraction",
                () => { Array.Copy(input, a.x, 30); return Old.Problem.Problem.perf(a, oldProblem).f[0]; },
                () => { Array.Copy(input, b.x, 30); return Current.Problem.Problem.perf(b, problem).f[0]; }));
        }
        Old.Program.funct = new Old.Landscape(3) { N = 3, x = new[] { -1.0, 0.0, 1.0 }, fx = new[] { 2.0, 0.0, 3.0 } };
        Current.Program.funct = new Current.Landscape(3) { N = 3, x = new[] { -1.0, 0.0, 1.0 }, fx = new[] { 2.0, 0.0, 3.0 } };
        foreach (int function in new[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 99, 1000, 1007, 1008, 1015 })
        {
            int definition = function == 1007 ? 7 : function == 1008 ? 8 : function == 1015 ? 15 : function;
            Old.Problem.Problem oldProblem;
            Current.Problem.Problem problem;
            if (function == 13)
            {
                // The legacy factory allocates 32 coordinates for a 51-piece problem.
                // Verify its existing exception, then test the evaluator with valid-sized state.
                try { Old.Problem.Problem.problemDef(13, null); throw new Exception("Baseline cutting-stock factory unexpectedly succeeded"); }
                catch (IndexOutOfRangeException) { }
                try { Current.Problem.Problem.problemDef(13, null); throw new Exception("Cutting-stock factory behavior changed"); }
                catch (IndexOutOfRangeException) { }
                checks++;
                oldProblem = new Old.Problem.Problem { function = 13, SwarmSize = new Old.SwarmSize(51) { D = 51, valueNb = 5 }, ObjectiveValue = 940 };
                problem = new Current.Problem.Problem { function = 13, SwarmSize = new Current.SwarmSize(51) { D = 51, valueNb = 5 }, ObjectiveValue = 940 };
                for (int d = 0; d < 51; d++)
                {
                    oldProblem.SwarmSize.min[d] = problem.SwarmSize.min[d] = 1380;
                    oldProblem.SwarmSize.max[d] = problem.SwarmSize.max[d] = 2150;
                }
            }
            else
            {
                oldProblem = Old.Problem.Problem.problemDef(definition, null);
                problem = Current.Problem.Problem.problemDef(definition, null);
            }
            oldProblem.function = function; problem.function = function;
            if (oldProblem.SwarmSize.D != problem.SwarmSize.D) throw new Exception("Problem dimensions changed");
            Study.Equal(oldProblem.ObjectiveValue, problem.ObjectiveValue); checks++;
            var a = new Old.Position(Math.Max(Old.Constants.DMax, oldProblem.SwarmSize.D)) { size = oldProblem.SwarmSize.D };
            var b = new Current.Position(Math.Max(Current.Constants.DMax, problem.SwarmSize.D)) { size = problem.SwarmSize.D };
            for (int sample = 0; sample < 100; sample++)
            {
                for (int d = 0; d < a.size; d++)
                {
                    Study.Equal(oldProblem.SwarmSize.min[d], problem.SwarmSize.min[d]); checks++;
                    Study.Equal(oldProblem.SwarmSize.max[d], problem.SwarmSize.max[d]); checks++;
                    double value = function == 1000 ? random.NextDouble() * 4 - 2 : oldProblem.SwarmSize.min[d] + (oldProblem.SwarmSize.max[d] - oldProblem.SwarmSize.min[d]) * random.NextDouble();
                    a.x[d] = value; b.x[d] = value;
                }
                double oldCount = Old.Program.nEval, newCount = Current.Program.nEval;
                var expected = Old.Problem.Problem.perf(a, oldProblem);
                var actual = Current.Problem.Problem.perf(b, problem);
                if (expected.size != actual.size) throw new Exception("Fitness shape changed");
                Study.Equal(Old.Program.nEval - oldCount, Current.Program.nEval - newCount); checks++;
                for (int d = 0; d < expected.f.Length; d++) { Study.Equal(expected.f[d], actual.f[d]); checks++; }
                for (int d = 0; d < a.size; d++) { Study.Equal(a.x[d], b.x[d]); checks++; }
            }
        }
        var oldSphere = new Old.Problem.CEC2005F1Circle();
        var sphere = new Current.Problem.CEC2005F1Circle();
        for (int sample = 0; sample < 1000; sample++)
        {
            var a = new Old.Position(30) { size = 30 };
            var b = new Current.Position(30) { size = 30 };
            for (int d = 0; d < 30; d++) a.x[d] = b.x[d] = random.NextDouble() * 200 - 100;
            double oldCount = Old.Program.nEval, newCount = Current.Program.nEval;
            var expected = oldSphere.Evaluate(a); var actual = sphere.Evaluate(b);
            if (expected.size != actual.size) throw new Exception("Sphere fitness shape changed");
            Study.Equal(Old.Program.nEval - oldCount, Current.Program.nEval - newCount); checks++;
            for (int d = 0; d < expected.f.Length; d++) { Study.Equal(expected.f[d], actual.f[d]); checks++; }
            for (int d = 0; d < 30; d++) { Study.Equal(a.x[d], b.x[d]); checks++; }
        }
        Study.Save(args[0], checks, measurements, "Pinned Git baseline; constant-table reuse, responsibility extraction and shared sphere computation; nine alternating trials after warmup. Compare fitness and input mutation exactly. Existing shifted-sphere formula behavior is preserved.");
    }
}
static class Study
{
    public static void PublicApi(System.Reflection.Assembly baseline, System.Reflection.Assembly current)
    {
        string[] Describe(System.Reflection.Assembly assembly)
        {
            var members = new List<string>();
            foreach (var type in assembly.GetExportedTypes())
            {
                members.Add(type.FullName);
                foreach (var member in type.GetMembers(System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.Static | System.Reflection.BindingFlags.DeclaredOnly))
                    members.Add(type.FullName + ":" + member.MemberType + ":" + member);
            }
            members.Sort(StringComparer.Ordinal);
            return members.ToArray();
        }
        var expected = Describe(baseline); var actual = Describe(current);
        if (expected.Length != actual.Length) throw new Exception("Public API member count changed");
        for (int i = 0; i < expected.Length; i++)
            if (expected[i] != actual[i]) throw new Exception($"Public API changed: {expected[i]} versus {actual[i]}");
        Console.WriteLine($"PASS: {expected.Length} public type/member signatures unchanged");
    }

    public static void Equal(double expected, double actual)
    {
        if (BitConverter.DoubleToInt64Bits(expected) != BitConverter.DoubleToInt64Bits(actual))
            throw new Exception($"Numerical mismatch: {expected:R} versus {actual:R}");
    }
    public static object Measure(string name, Func<double> baseline, Func<double> candidate, int iterations = 20000)
    {
        for (int i = 0; i < 10000; i++) { baseline(); candidate(); }
        var oldTimes = new double[9]; var newTimes = new double[9];
        var oldBytes = new long[9]; var newBytes = new long[9];
        double sink = 0;
        void Batch(Func<double> call, double[] times, long[] bytes, int trial)
        {
            long before = GC.GetAllocatedBytesForCurrentThread();
            long start = System.Diagnostics.Stopwatch.GetTimestamp();
            for (int i = 0; i < iterations; i++) sink += call();
            times[trial] = System.Diagnostics.Stopwatch.GetElapsedTime(start).TotalMilliseconds;
            bytes[trial] = GC.GetAllocatedBytesForCurrentThread() - before;
        }
        for (int trial = 0; trial < 9; trial++)
        {
            if (trial % 2 == 0) { Batch(baseline, oldTimes, oldBytes, trial); Batch(candidate, newTimes, newBytes, trial); }
            else { Batch(candidate, newTimes, newBytes, trial); Batch(baseline, oldTimes, oldBytes, trial); }
        }
        GC.KeepAlive(sink);
        double[] sortedOld = (double[])oldTimes.Clone(), sortedNew = (double[])newTimes.Clone();
        Array.Sort(sortedOld); Array.Sort(sortedNew);
        return new { name, iterations, baselineMilliseconds = oldTimes, candidateMilliseconds = newTimes,
            baselineMedianMilliseconds = sortedOld[4], candidateMedianMilliseconds = sortedNew[4],
            speedup = sortedOld[4] / sortedNew[4], baselineBytesPerCall = oldBytes[4] / (double)iterations,
            candidateBytesPerCall = newBytes[4] / (double)iterations };
    }
    public static void Save(string path, int checks, List<object> measurements, string notes)
    {
        File.WriteAllText(path, System.Text.Json.JsonSerializer.Serialize(new {
            runtime = System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription,
            architecture = System.Runtime.InteropServices.RuntimeInformation.ProcessArchitecture.ToString(),
            numericalChecks = checks, maximumAbsoluteError = 0, maximumRelativeError = 0,
            measurements, notes
        }, new System.Text.Json.JsonSerializerOptions { WriteIndented = true }) + "\n");
        Console.WriteLine($"PASS: {checks} exact numerical checks; results: {path}");
    }
}
