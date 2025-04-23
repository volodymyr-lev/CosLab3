using ScottPlot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Lab2DPF
{
    internal static class CalcUnit
    {
        public static List<double> GivenValues = new List<double>() { 2.81, 2.95, 3.21, 3.39, 3.54, 3.76 };

        public static double[] LeastSquaresQuadratic(double[] x, double[] y)
        {
            int n = x.Length;
            double sumX = x.Sum();
            double sumX2 = x.Sum(val => val * val);
            double sumX3 = x.Sum(val => val * val * val);
            double sumX4 = x.Sum(val => val * val * val * val);
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumX2Y = x.Zip(y, (xi, yi) => xi * xi * yi).Sum();

            double[,] A =
            {
                { n, sumX, sumX2 },
                { sumX, sumX2, sumX3 },
                { sumX2, sumX3, sumX4 }
            };

            double[] B = { sumY, sumXY, sumX2Y };
            return SolveLinearSystem(A, B);
        }
        private static double[] SolveLinearSystem(double[,] A, double[] B)
        {
            int n = B.Length;
            double[] X = new double[n];

            for (int i = 0; i < n; i++)
            {
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                    if (Math.Abs(A[k, i]) > Math.Abs(A[maxRow, i]))
                        maxRow = k;

                for (int k = i; k < n; k++)
                    (A[maxRow, k], A[i, k]) = (A[i, k], A[maxRow, k]);
                (B[maxRow], B[i]) = (B[i], B[maxRow]);

                for (int k = i + 1; k < n; k++)
                {
                    double factor = A[k, i] / A[i, i];
                    for (int j = i; j < n; j++)
                        A[k, j] -= factor * A[i, j];
                    B[k] -= factor * B[i];
                }
            }

            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i + 1; j < n; j++)
                    sum += A[i, j] * X[j];
                X[i] = (B[i] - sum) / A[i, i];
            }
            return X;
        }


        private static List<double> Ans { get; set; } = new List<double>();
        private static List<double> Bns { get; set; } = new List<double>();

        public static Dictionary<double, double> MNKcoords = new Dictionary<double, double>();

        public static double Harmonics { get; private set; }

        public static double F(double x)
        {
            return MNKcoords[x];
        }
        public static double An(int n, double I1, double I2)
        {
            if (n == 0)
            {
                // a0 = (1 / π) * ∫ f(x) dx
                return (1.0 / Math.PI) * Integral(x => InterpolatedF(x), I1, I2, 4000);
            }

            // an = (1 / π) * ∫ f(x) * cos(n * x) dx
            return (1.0 / Math.PI) * Integral(x => InterpolatedF(x) * Math.Cos(n * x), I1, I2, 4000);
        }

        public static double Bn(int n, double I1, double I2)
        {
            if (n == 0)
            {
                // b0 = 0, оскільки sin(0 * x) = 0
                return 0;
            }

            // bn = (1 / π) * ∫ f(x) * sin(n * x) dx
            return (1.0 / Math.PI) * Integral(x => InterpolatedF(x) * Math.Sin(n * x), I1, I2, 4000);
        }

        private static double InterpolatedF(double x)
        {
            // Інтерполяція значень з MNKcoords
            var keys = MNKcoords.Keys.ToArray();
            var values = MNKcoords.Values.ToArray();

            for (int i = 0; i < keys.Length - 1; i++)
            {
                if (x >= keys[i] && x <= keys[i + 1])
                {
                    // Лінійна інтерполяція
                    double t = (x - keys[i]) / (keys[i + 1] - keys[i]);
                    return values[i] * (1 - t) + values[i + 1] * t;
                }
            }

            // Якщо x поза межами то повертаємо крайні значення
            return x < keys[0] ? values[0] : values[^1];
        }

        private static double Integral(Func<double, double> f, double a, double b, int steps)
        {
            double h = (b - a) / steps;
            double sum = 0.5 * (f(a) + f(b));

            for (int i = 1; i < steps; i++)
            {
                sum += f(a + i * h);
            }

            return sum * h;
        }
        

        public static double CountFourier(int harmonics, double x, double I1, double I2)
        {
            Harmonics = harmonics;

            double result = Ans[0] / 2.0;

            for (int n = 1; n <= harmonics; n++)
            {
                double an = Ans[n];
                double bn = Bns[n];

                result += an * Math.Cos(n * x) + bn * Math.Sin(n * x);
            }

            return result;
        }




        public static void countCoefficients(int N, double I1, double I2)
        {
            Ans.Clear(); Bns.Clear();

            double a0 = An(0, I1, I2); Ans.Add(a0);
            double b0 = Bn(0, I1, I2); Bns.Add(b0);

            for (int i = 1; i <= N; i++)
            {
                double an = An(i, I1, I2);
                double bn = Bn(i, I1, I2);

                Ans.Add(an); Bns.Add(bn);

            }
        }

        //public static Complex[] DFT(List<double> values)
        //{
        //    int N = values.Count;
        //    Complex[] result = new Complex[N];

        //    for (int k = 0; k < N; k++)
        //    {
        //        Complex sum = Complex.Zero;
        //        for (int n = 0; n < N; n++)
        //        {
        //            double angle = -2.0 * Math.PI * k * n / N;
        //            sum += values[n] * Complex.Exp(new Complex(0, angle));
        //        }
        //        result[k] = sum;
        //    }

        //    return result;
        //}

        //public static List<DataPoint> ReconstructSignal(Complex[] dftValues, double I1, double I2, int pointsCount)
        //{
        //    int N = dftValues.Length;
        //    List<DataPoint> result = new List<DataPoint>();
        //    double step = (I2 - I1) / (pointsCount - 1);

        //    for (int i = 0; i < pointsCount; i++)
        //    {
        //        double t = I1 + i * step;
        //        Complex sum = Complex.Zero;
        //        double normalizedT = (t - I1) / (I2 - I1) * N;

        //        for (int k = 0; k < N; k++)
        //        {
        //            double angle = 2.0 * Math.PI * k * normalizedT / N;
        //            sum += dftValues[k] * Complex.Exp(new Complex(0, angle));
        //        }
        //        result.Add(new DataPoint(t, sum.Real/N, 0));
        //    }
        //    return result;
        //}
    }
}
