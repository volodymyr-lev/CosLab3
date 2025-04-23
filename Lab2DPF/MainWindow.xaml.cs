using System.Numerics;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Lab2DPF
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private Plotter plotter;
        private int N;
        private double I1, I2;

        public MainWindow()
        {
            InitializeComponent();

            plotter = new Plotter(WpfPlot1);
        }

        private void BtnDrawGiven_Click(object sender, RoutedEventArgs e)
        {
            GetData();
            plotter.PlotGivenDiscrete(I1, I2);
        }

        private void BtnDrawAproxDots_Click(object sender, RoutedEventArgs e)
        {
            if (!GetData()) return;

            (double deltaAprox, double a0, double a1, double a2) = plotter.PlotMNK(N, I1, I2);

            LblDeltaAprox.Content = $"Delta Aprox: {deltaAprox}";
            LblA0.Content = $"a0: {a0}";
            LblA1.Content = $"a1: {a1}";
            LblA2.Content = $"a2: {a2}";

        }

        private void BtnDrawAprox_Click(object sender, RoutedEventArgs e)
        {
            DataGridComplex.Items.Clear();
            if (!GetData()) return;
            
            (Complex[] complexVals, double error) = plotter.PlotFourier(N, I1, I2);

            //LblDelta.Content = $"Delta Fourier: {error}";
        }

        private void BtnClear_Click(object sender, RoutedEventArgs e)
        {
            plotter.ClearPlot();
        }


        private bool GetData()
        {
            if (string.IsNullOrEmpty(TextBoxI1.Text) || string.IsNullOrEmpty(TextBoxI2.Text)) return false;

            N = InputConverter.GetN(TextBoxN.Text);
            (I1, I2) = InputConverter.GetIntervals(TextBoxI1.Text, TextBoxI2.Text);
            return true;
        }
    }
}