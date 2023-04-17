internal class Program
{
    private static void Main(string[] args)
    {
        int n = 7;
        int p = 3;
        double[,] A = new double[n, n];
        int[,] B = new int[n, 1];

        for (int i = 0; i < n; i++)
        {
            B[i, 0] = n + 1 - i;

            for (int j = 0; j < n; j++)
            {
                double ij = (i+1) * (j+1);
                A[i, j] = 1500 / Math.Pow(3 + 0.5 * ij, p);
            }
        }

        Console.WriteLine("Matrix A:");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Console.Write("{0, -18:F8}", A[i, j]);
            }
            Console.WriteLine();
        }

        Console.WriteLine("Matrix B:");
        for (int i = 0; i < n; i++)
        {
            Console.WriteLine("{0, -1}", B[i, 0]);
        }

        Console.ReadLine();
    }
}