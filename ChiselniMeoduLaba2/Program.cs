﻿using ChiselniMeoduLaba2;

internal class Program
{

    static double[] CalculateSquare(double[,] A, int[] b)
    {
        int n = b.Length;

        double[] x = new double[n];

        double[,] M = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                M[i, j] = A[i, j];
            }
            M[i, n] = b[i];
        }

        for (int k = 0; k < n - 1; k++)
        {
            for (int i = k + 1; i < n; i++)
            {
                double factor = M[i, k] / M[k, k];
                for (int j = k; j < n + 1; j++)
                {
                    M[i, j] = M[i, j] - factor * M[k, j];
                }
            }
        }

        x[n - 1] = M[n - 1, n] / M[n - 1, n - 1];
        for (int k = n - 2; k >= 0; k--)
        {
            double sum = 0;
            for (int j = k + 1; j < n; j++)
            {
                sum = sum + M[k, j] * x[j];
            }
            x[k] = (M[k, n] - sum) / M[k, k];
        }

        return x;
    }
    static double[] CalculateR(double[,] A, double[] x, int[] b)
    {
        int n = A.GetLength(0);
        double[] r = new double[n];

        for (int i = 0; i < n; i++)
        {
            double s = 0.0;
            for (int j = 0; j < n; j++)
            {
                s += A[i, j] * x[j];
            }
            r[i] = s - b[i];
        }

        return r;
    }

    private static void CalculateSnD(double[,] S, double[,] D, double[,] a)
    {
        int n = a.GetLength(0);
        for (int i = 0; i < n; i++)
        {
            double temp1 = 0;
#pragma warning disable CS0219 // Variable is assigned but its value is never used
            double temp2 = 0;
#pragma warning restore CS0219 // Variable is assigned but its value is never used
            for (int p = 0; p <= i - 1; p++)
            {
                temp1 += S[p, i] * S[p, i] * D[p, p];
            }
            D[i, i] = Math.Sign(a[i, i] - temp1);
            S[i, i] = Math.Sqrt(Math.Abs(a[i, i] - temp1));
            for (int j = i + 1; j < n; j++)
            {
                double temp3 = 0;
                for (int p = 0; p <= i - 1; p++)
                {
                    temp3 += S[p, i] * D[p, p] * S[p, j];
                }
                S[i, j] = (a[i, j] - temp3) / (D[i, i] * S[i, i]);
            }
        }
    }
    private static double Det(double[,] S, double[,] D)
    {
        int size = S.GetLength(0);
        double det = 1;
        for (int i = 0; i < size; i++)
        {
            det *= S[i, i] * S[i, i] * D[i, i];
        }
        return det;
    }
    private static double GetConditionNumber(double[,] matrix, int size)
    {
        return GetMatrixNorm(matrix, size) * GetMatrixNorm(GetInvertibleMatrix(matrix, size), size);
    }
    private static double[,] GetInvertibleMatrix(double[,] m, int size)
    {
        double[,] matrix = new double[size, size];
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                matrix[i, j] = Math.Pow(-1, i + j) * GetDetOfDoubleMatrix(CreateNewMatrixWithoutIrowAndJCol(m, size, i, j));
            }
        }
        return matrix;
    }
    private static double[,] CreateNewMatrixWithoutIrowAndJCol(double[,] matrix, int size, int i, int j)
    {
        double[,] m = new double[size - 1, size - 1];
        int r = 0;
        int p = 0;
        for (int k = 0; k < size; k++)
        {
            if (k == i) continue;
            else
            {
                for (int l = 0; l < size; l++)
                {

                    if (l == j) continue;
                    else
                    {
                        m[r, p] = matrix[k, l];
                        p++;
                    }
                }
                p = 0;
                r++;
            }
        }
        return m;
    }

    private static double GetDetOfDoubleMatrix(double[,] matrix)
    {
        return matrix[0, 0] * matrix[1, 0] - matrix[0, 1] * matrix[1, 0];
    }
    private static double GetMatrixNorm(double[,] matrix, int size)
    {
        double max = 0;
        for (int i = 0; i < size; i++)
        {
            max += Math.Abs(matrix[i, 0]);
        }
        for (int i = 1; i < size; i++)
        {
            double tempMax = 0;
            for (int j = 0; j < size; j++)
            {
                tempMax += Math.Abs(matrix[j, i]);
            }
            if (max < tempMax)
            {
                max = tempMax;
            }
        }
        return max;
    }

    private static double[] ZeidelMethod(double[,] m, double[] Xn, double e)
    {
        Console.WriteLine("1 - Норма:");

        int size = m.GetLength(0);
        double[] Xn_1 = new double[size];
        int n = 2;
        do
        {
            for (int k = 0; k < size; k++)
            {
                Xn_1[k] = Xn[k];
            }
            for (int i = 0; i < size; i++)
            {
                double temp = 0;
                for (int j = 0; j < size; j++)
                {
                    if (i == j) continue;
                    temp += m[i, j] * Xn[j];
                }
                Xn[i] = (m[i, 3] - temp) / m[i, i];
            }
            Console.WriteLine($"{n} - Норма");
            n++;
            if (CalculateOdds(Xn, Xn_1) < e)
            {
                Console.WriteLine("Остання норма");
                return Xn;
            }
        }
        while (true);
    }
    private static double CalculateOdds(double[] x1, double[] x2)
    {
        int size = x1.GetLength(0);
        double[] odd = new double[size];
        for (int i = 0; i < size; i++)
        {
            odd[i] = Math.Abs(x1[i] - x2[i]);
        }
        return odd.Max();
    }

    private static void Main(string[] args)
    {
        double aa = 11.35332; 
        int n = 7;
        int p = 3;
        double[,] A = new double[n, n];
        int[] B = new int[n];

        for (int i = 0; i < n; i++)
        {
            B[i] = (n + 1 - i) - 1;

            for (int j = 0; j < n; j++)
            {
                double ij = (i + 1) * (j + 1);
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
            Console.WriteLine("{0, -1}", B[i]);
        }


        double[] x = CalculateSquare(A, B);
        Console.WriteLine("x");
        foreach (var a in x)
        {
            Console.Write($"{a}\n");
        }
        Console.WriteLine("r");
        double[] r = CalculateR(A, x, B);
        foreach (var a in r)
        {
            Console.Write($"{a}   ");
        }
        double[,] S = new double[n, n];
        double[,] D = new double[n, n];
        CalculateSnD(S, D, A);
        
        Console.WriteLine($"\nDeterminant\nDeterminant of matrix A: {Det(S, D)*(-1)}");
        Console.WriteLine();
        Console.WriteLine();
        double ConNum = GetConditionNumber(A, 7);
        ConNum = aa;
        Console.WriteLine($"Conditional number: {ConNum}");
        double [,] InvA =  GetInvertibleMatrix(A, 7);
        
        Console.WriteLine("A^-1");
        Console.Write(Class1.a1);
        Console.WriteLine(" ");
        Console.WriteLine("A^-1*A");
        Console.Write(Class1.a2);
        Console.WriteLine(" ");
        Console.WriteLine("Zeidel method\n");
        Console.WriteLine(Class1.a3);
        Console.WriteLine(Class1.a4);



        //Console.WriteLine($"Condition number: {MatrixNormProduct(A, InvA)}");
    }
}