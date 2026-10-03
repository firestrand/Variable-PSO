using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace VPSOTest
{
    [TestClass]
    public class ProblemTests
    {
        [TestMethod]
        public void CanInitializeCEC2005F1Sphere()
        {
            var problem = new VPSO.Problem.CEC2005F1Circle();
            Assert.AreEqual(30, problem.SwarmSize.D);
            Assert.AreEqual(-450.0, problem.ObjectiveValue);
            for (int d = 0; d < problem.SwarmSize.D; d++)
            {
                Assert.AreEqual(-100.0, problem.SwarmSize.min[d]);
                Assert.AreEqual(100.0, problem.SwarmSize.max[d]);
                Assert.AreEqual(0.0, problem.SwarmSize.q.Q[d]);
            }
        }
        [TestMethod]
        public void SphereEvaluationPreservesLegacyFitnessShapeAndCountsOnce()
        {
            var problem = new VPSO.Problem.CEC2005F1Circle();
            var position = new VPSO.Position(30) { size = 1 };
            position.x[0] = 2.0;
            double previousCount = VPSO.Program.nEval;
            try
            {
                var fitness = problem.Evaluate(position);
                Assert.AreEqual(4.0, fitness.f[0]);
                Assert.AreEqual(0, fitness.size);
                Assert.AreEqual(41.3119, position.x[0], 1e-12);
                Assert.AreEqual(previousCount + 1, VPSO.Program.nEval);
            }
            finally
            {
                VPSO.Program.nEval = previousCount;
            }
        }
    }
}
