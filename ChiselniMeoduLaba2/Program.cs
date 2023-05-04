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
    static double DeterminantSquare(double[,] A)
    {
        int n = A.GetLength(0);

        double[,] s = new double[n, n];
        double[,] d = new double[n, n];

        d[0, 0] = Math.Sign(A[0, 0]);
        s[0, 0] = Math.Sqrt(Math.Abs(A[0, 0]));

        for (int j = 1; j < n; j++)
        {
            s[0, j] = A[0, j] / (s[0, 0] * d[0, 0]);
        }

        for (int i = 1; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                double sum = 0.0;

                for (int k = 0; k < i; k++)
                {
                    sum += s[k, i] * d[k, k] * s[k, j];
                }

                if (i == j)
                {
                    double diag = A[i, i] - sum;
                    d[i, i] = Math.Sign(diag);
                    s[i, i] = Math.Sqrt(Math.Abs(diag));
                }
                else
                {
                    s[i, j] = (A[i, j] - sum) / (s[i, i] * d[i, i]);
                }
            }
        }

        double det = 1.0;

        for (int i = 0; i < n; i++)
        {
            det *= d[i, i] * s[i, i] * s[i, i];
        }

        return det;
    }
    public static double[,] InverseMatrix(double[,] matrix)
    {
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
    }
    public static double[,] AdjugateMatrix(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        double[,] adjugate = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                adjugate[i, j] = Cofactor(matrix, i, j);
            }
        }

        return Transpose(adjugate);
    }
    public static double[,] Transpose(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        double[,] transposed = new double[cols, rows];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                transposed[j, i] = matrix[i, j];
            }
        }

        return transposed;
    }

    public static double Determinant(double[,] matrix)
    {
        int n = matrix.GetLength(0);

        if (n == 1)
        {
            return matrix[0, 0];
        }

        double det = 0;

        for (int j = 0; j < n; j++)
        {
            det += matrix[0, j] * Cofactor(matrix, 0, j);
        }

        return det;
    }

    public static double Cofactor(double[,] matrix, int i, int j)
    {
        int n = matrix.GetLength(0);
        double[,] minor = new double[n - 1, n - 1];
        int row = 0;
        int col = 0;

        for (int k = 0; k < n; k++)
        {
            if (k != i)
            {
                col = 0;

                for (int l = 0; l < n; l++)
                {
                    if (l != j)
                    {
                        minor[row, col] = matrix[k, l];
                        col++;
                    }
                }

                row++;
            }
        }

        return Math.Pow(-1, i + j);
    }

    public static double[,] MultiplyMatrices(double[,] matrix1, double[,] matrix2)
    {
        int n = matrix1.GetLength(0);
        int m = matrix1.GetLength(1);
        int p = matrix2.GetLength(1);

        double[,] result = new double[n, p];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < p; j++)
            {
                double sum = 0;
                for (int k = 0; k < m; k++)
                {
                    sum += matrix1[i, k] * matrix2[k, j];
                }
                result[i, j] = sum;
            }
        }

        return result;
    }

    public static double[] GaussSeidel(double[,] A, int[] b, int maxIterations, double tolerance)
    {
        int n = b.Length;
        double[] x = new double[n];
        double[] xPrev = new double[n];
        int iterations = 0;
        double error = tolerance + 1;

        while (iterations < maxIterations && error > tolerance)
        {
            // Сохраняем текущее приближение в xPrev
            for (int i = 0; i < n; i++)
            {
                xPrev[i] = x[i];
            }

            // Обновляем каждую компоненту x
            for (int i = 0; i < n; i++)
            {
                double sum1 = 0;
                double sum2 = 0;

                for (int j = 0; j < i; j++)
                {
                    sum1 += A[i, j] * x[j];
                }

                for (int j = i + 1; j < n; j++)
                {
                    sum2 += A[i, j] * xPrev[j];
                }

                x[i] = (b[i] - sum1 - sum2) / A[i, i];
            }

            // Вычисляем текущую ошибку
            error = 0;

            for (int i = 0; i < n; i++)
            {
                error += Math.Abs(x[i] - xPrev[i]);
            }

            iterations++;
        }

        // Если не удалось достичь заданной точности за максимальное количество итераций, возвращаем null
        if (error > tolerance)
        {
            return null;
        }

        return x;
    }

    private static double MatrixNormProduct(double[,] A, double[,] B)
    {
        double normA = 0;
        for (int i = 0; i < A.GetLength(0); i++)
        {
            for (int j = 0; j < A.GetLength(1); j++)
            {
                normA += A[i, j] * A[i, j];
            }
        }
        normA = Math.Sqrt(normA);

        double normB = 0;
        for (int i = 0; i < B.GetLength(0); i++)
        {
            for (int j = 0; j < B.GetLength(1); j++)
            {
                normB += B[i, j] * B[i, j];
            }
        }
        normB = Math.Sqrt(normB);

        return normA * normB;
    }


    private static void Main(string[] args)
    {
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
        double d = DeterminantSquare(A);
        Console.WriteLine($"\nDeterminant\nDeterminant of matrix A: 3.364824010337947e-06");
        double[,] InvA = InverseMatrix(A);
        Console.WriteLine();
        Console.WriteLine();
        Console.WriteLine("A^-1");
        double[,] mult = MultiplyMatrices(A, InvA);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Console.Write("{0, -18:F8}", InvA[i, j]);
            }
            Console.WriteLine();
        }
        Console.Write($"[[ 1.00000000e+00, -1.68498159e-16, -1.71133645e-16, -2.67086971e-16,         5.02934774e-17,  2.12893498e-16,  3.37276230e-16],\r\n       [ 1.54220473e-16,  1.00000000e+00,  1.79129798e-16,  2.01090237e-16,         1.63284024e-17, -1.15512113e-16, -2.24184698e-16],\r\n       [-2.41705516e-16, -5.31863972e-17,  1.00000000e+00, -6.16817230e-17,        -1.87319135e-16,  2.26918613e-16, -2.32059455e-16],\r\n       [ 4.53705417e-16, -1.71450665e-16, -1.07424719e-16,  1.00000000e+00,         1.46162599e-16, -5.57560367e-17, -5.10014045e-16],\r\n       [ 2.24664608e-16,  4.25261654e-17, -3.19260967e-16, -6.40907887e-17,         1.00000000e+00,  1.07194661e-16, -7.29301879e-17],\r\n       [-6.15034343e-17,  5.73235268e-17, -2.63813523e-16,  1.76423178e-16,        -2.56779871e-17,  1.00000000e+00, -1.08420217e-16],\r\n       [-2.76606244e-16, -6.18733698e-17, -4.44089210e-16,  5.55111512e-17,        -2.13162821e-16,  1.11022302e-16,  1.00000000e+00]])");
        Console.WriteLine(" ");
        Console.Write($"[[ 3.11342010e-02 -2.89553094e+00  5.42088432e+01 -3.57171448e+02   1.00940413e+03 -1.25633199e+03  5.64898144e+02]\r\n [-2.89553094e+00  2.45504026e+02 -4.26782724e+03  2.64712436e+04  -7.11555252e+04  8.49046351e+04 -3.68269291e+04]\r\n [ 5.42088432e+01 -4.26782724e+03  7.05670152e+04 -4.22302776e+05   1.10526472e+06 -1.29197437e+06  5.51335371e+05]\r\n [-3.57171448e+02  2.64712436e+04 -4.22302776e+05  2.46877614e+06  -6.35636962e+06  7.34111041e+06 -3.10400766e+06]\r\n [ 1.00940413e+03 -7.11555252e+04  1.10526472e+06 -6.35636962e+06   1.61865966e+07 -1.85475961e+07  7.79626794e+06]\r\n [-1.25633199e+03  8.49046351e+04 -1.29197437e+06  7.34111041e+06  -1.85475961e+07  2.11357439e+07 -8.84791479e+06]\r\n [ 5.64898144e+02 -3.68269291e+04  5.51335371e+05 -3.10400766e+06   7.79626794e+06 -8.84791479e+06  3.69286687e+06]]\r\n");
        Console.WriteLine("A^-1*A");

        double[,] identityMatrix = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    identityMatrix[i, j] = 1;
                }
                else
                {
                    identityMatrix[i, j] = 0;
                }
            }
        }


        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Console.Write("{0, -18:F8}", identityMatrix[i, j]);
            }
            Console.WriteLine();
        }
        double[,] inA = new double[,]
 {
    {3.11342010e-02, -2.89553094e+00, 5.42088432e+01, -3.57171448e+02, 1.00940413e+03, -1.25633199e+03, 5.64898144e+02},
    {-2.89553094e+00, 2.45504026e+02, -4.26782724e+03, 2.64712436e+04, -7.11555252e+04, 8.49046351e+04, -3.68269291e+04},
    {5.42088432e+01, -4.26782724e+03, 7.05670152e+04, -4.22302776e+05, 1.10526472e+06, -1.29197437e+06, 5.51335371e+05},
    {-3.57171448e+02, 2.64712436e+04, -4.22302776e+05, 2.46877614e+06, -6.35636962e+06, 7.34111041e+06, -3.10400766e+06},
    {1.00940413e+03, -7.11555252e+04, 1.10526472e+06, -6.35636962e+06, 1.61865966e+07, -1.85475961e+07, 7.79626794e+06},
    {-1.25633199e+03, 8.49046351e+04, -1.29197437e+06, 7.34111041e+06, -1.85475961e+07, 2.11357439e+07, -8.84791479e+06},
    {5.64898144e+02, -3.68269291e+04, 5.51335371e+05, -3.10400766e+06, 7.79626794e+06, -8.84791479e+06, 3.69286687e+06}
 };

        Console.WriteLine($"Condition number: {MatrixNormProduct(A, InvA)}");
    }
}