﻿using ScottPlot;
using ScottPlot.WPF;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Lab2DPF
{
    internal class Plotter
    {
        public WpfPlot PlotWindow { get; init; }

        public double Step { get; set; }
        List<DataPoint> dataXYMNK = new List<DataPoint>();

        public Plotter(in WpfPlot PlotWindow)
        {
            this.PlotWindow = PlotWindow;
            PlotWindow.Plot.Add.HorizontalLine(0, color: Color.FromSDColor(System.Drawing.Color.Black));
            PlotWindow.Plot.Add.VerticalLine(0, color: Color.FromSDColor(System.Drawing.Color.Black));
        }

        public void PlotGivenDiscrete(double I1, double I2)
        {
            Step = Math.Abs(I2 - I1) / (CalcUnit.GivenValues.Count()-1);

            List<double> dataX = new List<double>();
            List<double> dataY = new List<double>();

            int n = 0;
            
            for (double i = I1; i <= I2 && n < CalcUnit.GivenValues.Count(); i += Step)
            {
                dataX.Add(i);
                dataY.Add(CalcUnit.GivenValues[n]);
                n++;
            }

            var scatter = PlotWindow.Plot.Add.ScatterPoints(dataX.ToArray(), dataY.ToArray());
            PlotWindow.Refresh();
        }

        public (double, double, double, double) PlotMNK(int N, double i1, double i2)
        {

            List<double> dataX = new List<double>();
            List<double> dataY = new List<double>();

            CalcUnit.MNKcoords.Clear();
            dataXYMNK.Clear();

            int n = 0;
            for (double i = i1; i <= i2 && n < CalcUnit.GivenValues.Count(); i += Step)
            {
                dataX.Add(i);
                dataY.Add(CalcUnit.GivenValues[n]);
                n++;
            }

            var coeffs = CalcUnit.LeastSquaresQuadratic(dataX.ToArray(), dataY.ToArray());
            double a0 = coeffs[0], a1 = coeffs[1], a2 = coeffs[2];

            List<double> dataXMNK = new List<double>();
            List<double> dataYMNK = new List<double>();

            for (double xi = dataX.Min(); xi <= dataX.Max(); xi += 0.01)
            {
                double yi = a0 + a1 * xi + a2 * xi * xi;
                dataXMNK.Add(xi);
                dataYMNK.Add(yi);

                CalcUnit.MNKcoords.Add(xi, yi);

                //for error
                dataXYMNK.Add(new DataPoint(xi, yi, 0));
            }

            var scatter = PlotWindow.Plot.Add.Scatter(dataXMNK.ToArray(), dataYMNK.ToArray());

            return (CalculateAbsoluteError(i1,i2,dataXYMNK,CalcUnit.GivenValues), a0, a1, a2);
        }

        public double PlotFourier(int N, double I1, double I2)
        {
            List<double> dataX = new List<double>();
            List<double> dataY = new List<double>();
            List<DataPoint> dataXY = new List<DataPoint>();

            double step = 0.01;

            CalcUnit.countCoefficients(N, I1, I2);

            for (double i = I1; i < I2; i += step)
            {
                double fourier = CalcUnit.CountFourier(N, i, I1, I2);
                dataX.Add(i);
                dataY.Add(fourier);


                //for error
                dataXY.Add(new DataPoint(i,fourier,0));
            }

            var scatter = PlotWindow.Plot.Add.Scatter(dataX.ToArray(), dataY.ToArray());



            double deltaAprox = 0.0;
            for (int i = 0; i < dataXY.Count; i++)
            {
                deltaAprox += Math.Abs(dataXY[i].Y - dataXYMNK[i].Y);
            }

            return deltaAprox/dataXY.Count;
        }

        private double CalculateAbsoluteError(double I1, double I2, List<DataPoint> restored, List<double> originalValues)
        {
            int N = originalValues.Count;
            double errorSum = 0.0;

            Step = Math.Abs(I2 - I1) / N;

            double epsilon = 1e-2;

            int n = 0;
            for (double i = I1; i <= I2 && n < N; i+=Step)
            {
                double fRestoredValue = restored.FirstOrDefault(p=>Math.Abs(p.X - i)<epsilon).Y;
                double fOriginalValue = originalValues[n];
            
                errorSum += Math.Abs(fOriginalValue - fRestoredValue);
            
                n++;
            }

            return errorSum / N; 
        }

        public void ClearPlot()
        {
            PlotWindow.Plot.Clear();

            PlotWindow.Plot.Add.HorizontalLine(0, color: Color.FromSDColor(System.Drawing.Color.Black));
            PlotWindow.Plot.Add.VerticalLine(0, color: Color.FromSDColor(System.Drawing.Color.Black));

            PlotWindow.Refresh();
        }
    }
}