using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace HyperellipticCurves
{
    class Methods
    {
        public static List<BigInteger> PrimeFactors(BigInteger a)
        {
            if (a < 2)
                return null;

            var potentialPrimes = new List<int>();
            var smallPrimes = new List<int> { 2, 3, 5, 7, 11, 13, 17 };
            var step = 1;
            foreach (int i in smallPrimes)
                step *= i;

            for (int i = 1; i < step; i++) 
            {
                bool prime = true;
                foreach (int p in smallPrimes)
                    if (i % p == 0)
                    {
                        prime = false;
                        break;
                    }
                if (prime)
                    potentialPrimes.Add(i);
            }

            // Console.WriteLine((double)potentialPrimes.Count / (double)step);

            List<BigInteger> primes = new List<BigInteger> { };
            foreach (int i in smallPrimes)
                if (a % i == 0)
                {
                    primes.Add(i);
                    while (a % i == 0)
                        a /= i;
                }

            BigInteger t = 0;
            int j = 1;
            BigInteger c = 0;
            BigInteger last = BigInteger.Pow(2, (int)a.GetBitLength() / 2);

            while (c <= last)
            {
                if (j >= potentialPrimes.Count)
                {
                    t += step;
                    j = 0;
                }

                c = t + potentialPrimes[j];
                if (a % c == 0)
                {
                    primes.Add(c);
                    while (a % c == 0)
                        a /= c;
                }

                j++;
            }

            if (a > 1)
                primes.Add(a);

            return primes;
        }
        public static int NumRemainder(int a, int b)
        {
            int c = a % b;
            if (c < 0)
                c += Math.Abs(b);
            return c;
        }
        public static BigInteger NumRemainder(BigInteger a, BigInteger b)
        {
            BigInteger c = a % b;
            if (c < 0)
                if (b > 0)
                    c += b;
                else
                    c += -b;
            return c;
        }
        public static int NumericalEuclid(int a, int b, out int t, out int s)
        {
            //int origa = a;
            //int origb = b;

            int[] sl = { 0, 1 };
            int[] tl = { 1, 0 };

            while (a != 0 && b != 0)
            {
                int q = b / a;
                sl = new int[] { sl[1] - q * sl[0], sl[0] };
                tl = new int[] { tl[1] - q * tl[0], tl[0] };

                int c = a;
                a = NumRemainder(b, a);
                b = c;

                //Console.WriteLine(tl[0] * origa);
                //Console.WriteLine(sl[0] * origb);
                //Console.WriteLine(a);
                //Console.WriteLine();
            }

            t = tl[1];
            s = sl[1];
            return Math.Max(b, 1);
        }
        public static BigInteger NumericalEuclid(BigInteger a, BigInteger b, out BigInteger t, out BigInteger s)
        {
            BigInteger[] sl = { 0, 1 };
            BigInteger[] tl = { 1, 0 };

            while (a != 0 && b != 0)
            {
                BigInteger q = b / a;
                sl = new BigInteger[] { sl[1] - q * sl[0], sl[0] };
                tl = new BigInteger[] { tl[1] - q * tl[0], tl[0] };

                BigInteger c = a;
                a = NumRemainder(b, a);
                b = c;
            }

            t = tl[1];
            s = sl[1];
            if (b < 1)
                b = 1;
            return b;
        }
        public static int NumericalEuclid(int a, int b)
        {
            while (a != 0 && b != 0)
            {
                int c = a;
                a = NumRemainder(b, a);
                b = c;
            }

            return Math.Max(b, 1);
        }
        public static BigInteger NumericalEuclid(BigInteger a, BigInteger b)
        {
            while (a != 0 && b != 0)
            {
                BigInteger c = a;
                a = NumRemainder(b, a);
                b = c;
            }

            if (b < 1)
                b = 1;
            return b;
        }
        public static List<int> SolveCubicEquation3L(int a, GaloisField<int> field)
        {
            // solves equation y^3 - y = a in F_3^l

            if (field.Characteristic() != 3 || a > 2 || a < 0)
                throw new Exception();

            var l = field.dimension;
            var md = 3 * (l - 1) + 1;

            var m = new List<List<int>>(md);
            for (int i = 0; i < md; i++)
                m.Add(new List<int>(Enumerable.Repeat(0, md)));

            for (int i = 0; i < l; i++)
            {
                m[3 * i][i] = field.baseField.Add(m[3 * i][i], 1);
                m[i][i] = field.baseField.Subtract(m[i][i], 1);
            }

            for (int i = l; i < md; i++)
                for (int j = 0; j <= l; j++)
                    m[i - l + j][i] = field.baseField.Add(m[i - l + j][i], field.primitive[j]);

            var y = new List<int>(Enumerable.Repeat(0, md));
            y[0] = a;


            Gauss(m, y, (PrimeField)field.baseField);

            // m[0][0] will be 0, y wil be (a, 0 ... 0)
            // find non-zero in first row, set it to a*inv
            // other zero cols correspond to zeros
            // figure out "real" variables based on non-zero "fake" variable

            for (int j = 0; j < m.Count; j++)
                if (m[0][j] != 0)
                {
                    var fake = a * field.baseField.Inverse(m[0][j]);
                    var res = new List<int>(l);
                    res.Add(0);
                    for (int i = 1; i < l; i++)
                        res.Add(-fake * m[i][j]);

                    return res;
                }

            throw new Exception();
        }
        public static void Gauss(List<List<int>> m, List<int> y, PrimeField field)
        {
            //for (int i = 0; i < m.Count; i++)
            //    Program.Print(m[i]);
            //Console.WriteLine();
            //Program.Print(y);
            //Console.WriteLine();

            // implication: matrix is square
            if (m.Count != y.Count)
                throw new Exception();

            for (int j = 0; j < m.Count; j++)
            {
                bool zeroCol = true;
                for (int i = j; i < m.Count; i++)
                {
                    if (m[i][j] != 0)
                    {
                        // swap i and j rows
                        var temp = m[i];
                        m[i] = m[j];
                        m[j] = temp;

                        int tempi = y[i];
                        y[i] = y[j];
                        y[j] = tempi;

                        zeroCol = false;
                        break;
                    }
                }

                if (zeroCol)
                    continue;

                var inv = field.Inverse(m[j][j]);
                m[j] = field.MultiplyPoly(m[j], new List<int> { inv });
                y[j] = field.Multiply(y[j], inv);



                for (int i = 0; i < m.Count; i++)
                    if (i != j)
                    {
                        m[i] = field.SubtractPoly(m[i], field.MultiplyPoly(m[j], new List<int> { m[i][j] }));
                        y[i] = field.Subtract(y[i], field.Multiply(y[j], m[i][j]));
                    }

                //for (int i = 0; i < m.Count; i++)
                //    Program.Print(m[i]);
                //Console.WriteLine();
                //Program.Print(y);
                //Console.WriteLine();
            }

            //for (int i = 0; i < m.Count; i++)
            //    Program.Print(m[i]);
            //Console.WriteLine();
            //Program.Print(y);
        }
    }
}


//public List<GFElement> AddPoly(List<GFElement> a, List<GFElement> b)
//{
//    List<GFElement> c = new List<GFElement>();

//    for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
//        c.Add(Add(a[i], b[i]));

//    if (a.Count > b.Count)
//        for (int i = b.Count; i < a.Count; i++)
//            c.Add(a[i]);
//    else
//        for (int i = a.Count; i < b.Count; i++)
//            c.Add(b[i]);

//    return c;
//}

//public List<GFElement> SubtractPoly(List<GFElement> a, List<GFElement> b)
//{

//    List<GFElement> c = new List<GFElement>();

//    for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
//        c.Add(Subtract(a[i], b[i]));

//    if (a.Count > b.Count)
//        for (int i = b.Count; i < a.Count; i++)
//            c.Add(a[i]);
//    else
//        for (int i = a.Count; i < b.Count; i++)
//            c.Add(Subtract(GetElementFromPrime(0), b[i]));

//    return c;
//}

//public List<GFElement> MultiplyPoly(List<GFElement> a, List<GFElement> b)
//{
//    List<GFElement> c = new List<GFElement>(a.Count + b.Count - 1);
//    for (int i = 0; i < a.Count + b.Count - 1; i++)
//        c.Add(GetElementFromPrime(0));
//    for (int i = 0; i < a.Count; i++)
//        for (int j = 0; j < b.Count; j++)
//            c[i + j] = Add(c[i + j], Multiply(a[i], b[j]));

//    return c;
//}
//public List<GFElement> RemainderPoly(List<GFElement> a, List<GFElement> b, out List<GFElement> factor)
//{
//    List<GFElement> copy = new List<GFElement>(a.Count);
//    factor = new List<GFElement>(a.Count);
//    foreach (var el in a)
//    {
//        copy.Add(new GFElement(el.GetPolynomial(), this));
//        factor.Add(GetElementFromPrime(0));
//    }

//    var zero = GetElementFromPrime(0);
//    int cur = copy.FindLastIndex((GFElement el) => !Equals(el,zero));
//    int bf = b.FindLastIndex((GFElement el) => !Equals(el, zero));

//    GFElement inv = Inverse(b[bf]);

//    while (cur >= bf)
//    {
//        int diff = cur - bf;
//        if (!Equals(copy[cur], zero))
//        {
//            for (int i = 0; i < bf; i++)
//                copy[diff + i] = Subtract(copy[diff + i], Multiply(Multiply(copy[cur], inv), b[i]));
//            factor[diff] = Multiply(copy[cur], inv);
//            copy[cur] = new GFElement(zero.GetPolynomial(), this);
//        }
//        cur--;
//    }

//    return copy;
//}
//public List<GFElement> EuclidPoly(List<GFElement> a, List<GFElement> b, out List<GFElement> t, out List<GFElement> s)
//{
//    List<GFElement>[] sl = { new List<GFElement> { GetElementFromPrime(0) }, new List<GFElement> { GetElementFromPrime(1) } };
//    List<GFElement>[] tl = { new List<GFElement> { GetElementFromPrime(1) }, new List<GFElement> { GetElementFromPrime(0) } };

//    var ac = new List<GFElement>(a.Count);
//    foreach (var el in a)
//        ac.Add(new GFElement(el.GetPolynomial(), this));

//    var bc = new List<GFElement>(b.Count);
//    foreach (var el in b)
//        bc.Add(new GFElement(el.GetPolynomial(), this));

//    var zero = GetElementFromPrime(0);
//    while (ac.Count((GFElement i) => Equals(i, zero)) != ac.Count && bc.Count((GFElement i) => Equals(i, zero)) != bc.Count)
//    {
//        List<GFElement> q = new List<GFElement>();
//        List<GFElement> c = ac;
//        ac = RemainderPoly(bc, ac, out q);
//        bc = c;

//        sl = new List<GFElement>[] { SubtractPoly(sl[1], MultiplyPoly(q, sl[0])), sl[0] };
//        tl = new List<GFElement>[] { SubtractPoly(tl[1], MultiplyPoly(q, tl[0])), tl[0] };
//    }

//    t = tl[1];
//    s = sl[1];
//    if (bc.Count((GFElement i) => Equals(i, zero)) != bc.Count)
//        return bc;
//    else
//        return new List<GFElement> { GetElementFromPrime(1) };
//}

//public GFElement ConvertToLargerField(GFElement element, GaloisField largeField = null)
//{
//    if (largeField == null) // general case, won't be true only in constructor
//        largeField = curve6.field;

//    // input: element from field, output: element from field6
//    GFElement res = new GFElement(new List<int> { 0 }, largeField);

//    var cur = new GFElement(new List<int> { 1 }, largeField);
//    for (int i = 0; i < element.p.Count; i++)
//    {
//        if (element.p[i] != 0)
//            res += i * cur;
//        cur *= root;
//    }

//    //Program.Print(element);
//    //Program.Print(res);
//    return res;
//}

