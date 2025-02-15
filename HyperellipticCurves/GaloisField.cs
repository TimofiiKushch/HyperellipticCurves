using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using System.IO;
using System.Threading.Tasks;
using System.Numerics;
using System.Security.Cryptography;

namespace HyperellipticCurves
{
    public class GFElement
    {
        public List<int> p { get; }
        public GaloisField field { get; }

        public GFElement(List<int> polynomial, GaloisField field)
        {
            p = field.FixRepresentation(polynomial);
            this.field = field;
        }

        public static GFElement operator -(GFElement a)
        {
            return a.field.Subtract(new GFElement(new List<int> { 0 }, a.field), a);
        }

        public static GFElement operator +(GFElement a, GFElement b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return a.field.Add(a, b);
            else
                return null;
        }

        public static GFElement operator -(GFElement a, GFElement b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return a.field.Subtract(a, b);
            else
                return null;
        }

        public static GFElement operator *(GFElement a, GFElement b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return a.field.Multiply(a, b);
            else
                return null;
        }
        public static GFElement operator *(int a, GFElement b)
        {
            return b.field.Multiply(new GFElement(new List<int> { a }, b.field), b);
        }

        public static GFElement operator /(GFElement a, GFElement b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                return a * a.field.Inverse(b);
            }
            else
                return null;
        }

        public static bool operator ==(GFElement a, GFElement b)
        {
            if (a.field != b.field)
                return false;

            return a.field.IsEqual(a, b);
        }

        public static bool operator !=(GFElement a, GFElement b)
        {
            return !(a == b);
        }
    }
    public class ExtensionElement 
    {
        public List<GFElement> p { get; }
        public GaloisFieldExtension field { get; }

        public ExtensionElement(List<GFElement> polynomial, GaloisFieldExtension field)
        {
            p = field.FixRepresentation(polynomial);
            this.field = field;
        }

        public static ExtensionElement operator -(ExtensionElement a)
        {
            return new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(0) }, a.field) - a;
        }

        public static ExtensionElement operator +(ExtensionElement a, ExtensionElement b)
        {
            if (a.field == b.field)
                return new(a.field.baseField.AddPoly(a.p, b.p), a.field);
            else
                return null;
        }

        public static ExtensionElement operator -(ExtensionElement a, ExtensionElement b)
        {
            if (a.field == b.field)
                return new(a.field.baseField.SubtractPoly(a.p, b.p), a.field);
            else
                return null;
        }

        public static ExtensionElement operator *(ExtensionElement a, ExtensionElement b)
        {
            if (a.field == b.field)
                return new(a.field.baseField.MultiplyPoly(a.p, b.p), a.field);
            else
                return null;
        }
        public static ExtensionElement operator *(GFElement a, ExtensionElement b)
        {
            return new ExtensionElement(new List<GFElement> { a }, b.field) * b;
        }

        public static ExtensionElement operator /(ExtensionElement a, ExtensionElement b)
        {
            if (a.field == b.field)
            {
                return a * a.field.Inverse(b);
            }
            else
                return null;
        }

        public static bool operator ==(ExtensionElement a, ExtensionElement b)
        {
            if (a.field != b.field)
                return false;

            return a.field.IsEqual(a, b);
        }

        public static bool operator !=(ExtensionElement a, ExtensionElement b)
        {
            return !(a == b);
        }
    }

    public class ECP
    {
        public GFElement x { get; }
        public GFElement y { get; }

        public GaloisField field { get; }

        public ECP(GFElement x, GFElement y)
        {
            if (x.field == y.field)
            {
                this.x = x;
                this.y = y;
                field = x.field;
            }
        }

        public static bool operator ==(ECP a, ECP b)
        {
            return a.x == b.x && a.y == b.y;
        }

        public static bool operator !=(ECP a, ECP b)
        {
            return !(a == b);
        }
    }
    public class ECPExtension
    {
        public ExtensionElement x { get; }
        public ExtensionElement y { get; }

        public GaloisFieldExtension field { get; }

        public ECPExtension(ExtensionElement x, ExtensionElement y)
        {
            if (x.field == y.field)
            {
                this.x = x;
                this.y = y;
                field = x.field;
            }
        }

        public static bool operator ==(ECPExtension a, ECPExtension b)
        {
            return a.x == b.x && a.y == b.y;
        }

        public static bool operator !=(ECPExtension a, ECPExtension b)
        {
            return !(a == b);
        }
    }


    public class GaloisField
    {
        public int characteristic { get; }
        public int dimension { get; }
        private List<int> primitive = null;
        readonly PrimeField primeField;

        public BigInteger size { get; }

        public GaloisField(int characteristic, List<int> primitive)
        {
            this.characteristic = characteristic;
            dimension = primitive.Count - 1;
            this.primitive = primitive;

            primeField = new PrimeField(characteristic);

            size = BigInteger.Pow(characteristic, dimension);
        }

        public GFElement Scalar(int a)
        {
            return new(new List<int> { a }, this);
        }
        public GFElement Add(GFElement a, GFElement b)
        {
            return new GFElement(primeField.AddPoly(a.p, b.p), this);
        }
        public GFElement Subtract(GFElement a, GFElement b)
        {
            return new GFElement(primeField.SubtractPoly(a.p, b.p), this);
        }
        public GFElement Multiply(GFElement a, GFElement b)
        {
            return new GFElement(primeField.MultiplyPoly(a.p, b.p), this);
        }
        public GFElement Inverse(GFElement a)
        {
            if (primeField.Degree(a.p) == 0)
                return new GFElement(new List<int> { primeField.Inverse(a.p[0]) }, this);

            PrimePolynomial inv, s;
            var gcd = primeField.EuclidPoly(new(a.p, primeField), new(primitive, primeField), out inv, out s);
            var res = new GFElement(inv.poly, this);
            return primeField.Inverse(gcd.poly[0]) * res;
        }
        public GFElement Pow(GFElement a, BigInteger k)
        {
            GFElement res = null;
            GFElement last = a;

            while (true)
            {
                if (k % 2 == 1)
                {
                    if (!(res is null))
                        res = res * last;
                    else
                        res = last;
                }

                k /= 2;
                if (k == 0)
                    break;
                last = last * last;
            }

            if (res is null)
                return new GFElement(new List<int> { 1 }, this);

            return res;
        }
        public static List<GFElement> Root(GFElement a)
        {
            if (a.field.size % 2 == 1)
            {
                // quadratic residue check
                if (a.field.Pow(a, (a.field.size - 1) / 2) == new GFElement(new List<int> { 1 }, a.field))
                {
                    if (a.field.size % 4 == 3)
                    {
                        var root = a.field.Pow(a, (a.field.size + 1) / 4);
                        return new List<GFElement> { root, -root };
                    }
                    else if (a.field.size % 4 == 1)
                    {
                        // tonelli - shanks
                        BigInteger p1 = a.field.size - 1;
                        BigInteger odd = p1;

                        int degree2 = 0;
                        while (odd % 2 == 0)
                        {
                            odd /= 2;
                            degree2++;
                        }

                        GFElement re;
                        while (true)
                        {
                            re = a.field.RandomNonZero();
                            if (a.field.Pow(re, p1 / 2) == a.field.Scalar(-1))
                                break;
                        }

                        // precalculate all necessary powers
                        var apow = new List<GFElement>();
                        var repow = new List<GFElement>();
                        for (int i = 0; i < degree2; i++)
                        {
                            if (i == 0)
                            {
                                apow.Add(a.field.Pow(a, (odd + 1) / 2));
                                repow.Add(a.field.Pow(re, odd));
                            }
                            else if (i == 1)
                            {
                                apow.Add(a.field.Pow(a, odd));
                                repow.Add(repow[i - 1] * repow[i - 1]);
                            }
                            else
                            {
                                apow.Add(apow[i - 1] * apow[i - 1]);
                                repow.Add(repow[i - 1] * repow[i - 1]);
                            }
                        }

                        int n = 0;
                        BigInteger k = 1;
                        for (int i = degree2 - 1; i >= 1; i--)
                        {
                            n++;
                            if (apow[i] * a.field.Pow(repow[i], k) == a.field.Scalar(-1))
                                k += BigInteger.Pow(2, n - 1);
                        }

                        var root = apow[0] * a.field.Pow(repow[0], k);
                        return new List<GFElement> { root, -root };
                    }
                }
            }

            return null;
        }
        public static GFElement SolveCubicEquation3L(int a, GaloisField field)
        {
            // solves equation y^3 - y = a in F_3^l

            if (field.characteristic != 3 || a > 2 || a < 0)
                throw new Exception();

            var l = field.dimension;
            var md = 3 * (l - 1) + 1;

            var m = new List<List<int>>(md);
            for (int i = 0; i < md; i++)
                m.Add(new List<int>(Enumerable.Repeat(0, md)));

            for (int i = 0; i < l; i++)
            {
                m[3 * i][i] = field.primeField.Add(m[3 * i][i], 1);
                m[i][i] = field.primeField.Subtract(m[i][i], 1);
            }

            for (int i = l; i < md; i++)
                for (int j = 0; j <= l; j++)
                    m[i - l + j][i] = field.primeField.Add(m[i - l + j][i], field.primitive[j]);

            var y = new List<int>(Enumerable.Repeat(0, md));
            y[0] = a;


            field.primeField.Gauss(m, y);

            // m[0][0] will be 0, y wil be (a, 0 ... 0)
            // find non-zero in first row, set it to a*inv
            // other zero cols correspond to zeros
            // figure out "real" variables based on non-zero "fake" variable

            for (int j = 0; j < m.Count; j++)
                if (m[0][j] != 0)
                {
                    var fake = a * field.primeField.Inverse(m[0][j]);
                    var res = new List<int>(l);
                    res.Add(0);
                    for (int i = 1; i < l; i++)
                        res.Add(-fake * m[i][j]);

                    return new GFElement(res, field);
                }

            throw new Exception();
        }
        public GFElement RandomNonZero(int seed = -1)
        {
            var rand = new Random();
            if (seed >= 0)
                rand = new Random(seed);

            var poly = new List<int>(dimension);
            bool zero = true;

            while (zero)
            {
                poly.Clear();
                for (int i = 0; i < dimension; i++)
                {
                    poly.Add(rand.Next(0, characteristic));
                    if (poly[i] != 0)
                        zero = false;
                }
            }

            return new GFElement(poly, this);
        }
        public List<int> FixRepresentation(List<int> polynomial)
        {
            List<int> copy = new List<int>(polynomial);
            while (copy.Count < dimension)
                copy.Add(0);
            for (int i = 0; i < copy.Count; i++)
                copy[i] = Methods.NumRemainder(copy[i], characteristic);
            if (copy.Count > dimension)
            {
                var temp = new List<int>();
                copy = primeField.RemainderPoly(copy, primitive, out temp);
                copy.RemoveRange(dimension, copy.Count - dimension);
            }

            return copy;
        }
        public bool IsEqual(GFElement a, GFElement b)
        {
            // implication: elements from the same field
            List<int> ap = a.p, bp = b.p;
            for (int i = 0; i < ap.Count; i++)
                if (ap[i] != bp[i])
                    return false;
            return true;
        }
        public int Degree(List<GFElement> polynomial)
        {
            var zero = Scalar(0);
            var res = polynomial.FindLastIndex((GFElement el) => el != zero);
            return res < 0 ? 0 : res;
        }
        public GFElement LeadingCoeff(List<GFElement> polynomial)
        {
            return polynomial[Degree(polynomial)];
        }

        public List<GFElement> AddPoly(List<GFElement> a, List<GFElement> b)
        {
            List<GFElement> c = new List<GFElement>();

            for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
                c.Add(a[i] + b[i]);

            if (a.Count > b.Count)
                for (int i = b.Count; i < a.Count; i++)
                    c.Add(a[i]);
            else
                for (int i = a.Count; i < b.Count; i++)
                    c.Add(b[i]);

            return c;
        }

        public List<GFElement> SubtractPoly(List<GFElement> a, List<GFElement> b)
        {

            List<GFElement> c = new List<GFElement>();

            for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
                c.Add(a[i] - b[i]);

            if (a.Count > b.Count)
                for (int i = b.Count; i < a.Count; i++)
                    c.Add(a[i]);
            else
                for (int i = a.Count; i < b.Count; i++)
                    c.Add(-b[i]);

            return c;
        }

        public List<GFElement> MultiplyPoly(List<GFElement> a, List<GFElement> b)
        {
            List<GFElement> c = new List<GFElement>(a.Count + b.Count - 1);
            for (int i = 0; i < a.Count + b.Count - 1; i++)
                c.Add(Scalar(0));
            for (int i = 0; i < a.Count; i++)
                for (int j = 0; j < b.Count; j++)
                    c[i + j] += a[i] * b[j];

            return c;
        }
        public List<GFElement> RemainderPoly(List<GFElement> a, List<GFElement> b, out List<GFElement> factor)
        {
            List<GFElement> copy = new List<GFElement>(a.Count);
            factor = new List<GFElement>(a.Count);
            foreach (var el in a)
            {
                copy.Add(new GFElement(el.p, this));
                factor.Add(Scalar(0));
            }

            var zero = Scalar(0);
            int cur = copy.FindLastIndex((GFElement el) => el != zero);
            int bf = b.FindLastIndex((GFElement el) => el != zero);

            GFElement inv = Inverse(b[bf]);

            while (cur >= bf)
            {
                int diff = cur - bf;
                if (copy[cur] != zero)
                {
                    for (int i = 0; i < bf; i++)
                        copy[diff + i] -= copy[cur] * inv * b[i];
                    factor[diff] = copy[cur] * inv;
                    copy[cur] = new GFElement(zero.p, this);
                }
                cur--;
            }

            return copy;
        }
        public List<GFElement> EuclidPoly(List<GFElement> a, List<GFElement> b, out List<GFElement> t, out List<GFElement> s)
        {
            List<GFElement>[] sl = { new List<GFElement> { Scalar(0) }, new List<GFElement> { Scalar(1) } };
            List<GFElement>[] tl = { new List<GFElement> { Scalar(1) }, new List<GFElement> { Scalar(0) } };

            var ac = new List<GFElement>(a.Count);
            foreach (var el in a)
                ac.Add(new GFElement(el.p, this));

            var bc = new List<GFElement>(b.Count);
            foreach (var el in b)
                bc.Add(new GFElement(el.p, this));

            var zero = Scalar(0);
            while (LeadingCoeff(ac) != zero && LeadingCoeff(bc) != zero)
            {
                List<GFElement> q = new List<GFElement>();
                List<GFElement> c = ac;
                ac = RemainderPoly(bc, ac, out q);
                bc = c;

                sl = new List<GFElement>[] { SubtractPoly(sl[1], MultiplyPoly(q, sl[0])), sl[0] };
                tl = new List<GFElement>[] { SubtractPoly(tl[1], MultiplyPoly(q, tl[0])), tl[0] };
            }

            t = tl[1];
            s = sl[1];
            if (LeadingCoeff(bc) != zero)
                return bc;
            else
                return new List<GFElement> { Scalar(1) };
        }
        public void Gauss(List<List<GFElement>> m, List<GFElement> y)
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
                    if (m[i][j] != Scalar(0))
                    {
                        // swap i and j rows
                        var temp = m[i];
                        m[i] = m[j];
                        m[j] = temp;

                        GFElement tempi = y[i];
                        y[i] = y[j];
                        y[j] = tempi;

                        zeroCol = false;
                        break;
                    }
                }

                if (zeroCol)
                    continue;

                var inv = Inverse(m[j][j]);
                m[j] = MultiplyPoly(m[j], new List<GFElement> { inv });
                y[j] = Multiply(y[j], inv);



                for (int i = 0; i < m.Count; i++)
                    if (i != j)
                    {
                        m[i] = SubtractPoly(m[i], MultiplyPoly(m[j], new List<GFElement> { m[i][j] }));
                        y[i] = Subtract(y[i], Multiply(y[j], m[i][j]));
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

        public static bool operator ==(GaloisField a, GaloisField b)
        {
            if (a.dimension != b.dimension || a.characteristic != b.characteristic)
                return false;
            for (int i = 0; i < a.dimension; i++)
                if (a.primitive[i] != b.primitive[i])
                    return false;

            return true;
        }
        public static bool operator !=(GaloisField a, GaloisField b)
        {
            return !(a == b);
        }
    }
    public class GaloisFieldExtension
    {
        public int dimension { get; }
        private List<GFElement> primitive = null;
        readonly public GaloisField baseField;
        public BigInteger size { get; }
        public static ExtensionElement nonSquare = null;

        public GaloisFieldExtension(GaloisField baseField, List<int> primitive)
        {
            this.baseField = baseField;
            dimension = primitive.Count - 1;
            this.primitive = (from el in primitive select baseField.Scalar(el)).ToList();
            size = BigInteger.Pow(baseField.size, dimension);
        }
        public ExtensionElement Inverse(ExtensionElement a)
        {
            if (baseField.Degree(a.p) == 0)
                return new ExtensionElement(new List<GFElement> { baseField.Inverse(a.p[0]) }, this);

            List<GFElement> inv, s;
            var gcd = baseField.EuclidPoly(a.p, primitive, out inv, out s);
            var res = new ExtensionElement(inv, this);
            return baseField.Inverse(gcd[0]) * res;
        }

        public List<GFElement> FixRepresentation(List<GFElement> polynomial)
        {
            List<GFElement> copy = new List<GFElement>(polynomial);
            while (copy.Count < dimension)
                copy.Add(baseField.Scalar(0));
            if (copy.Count > dimension)
            {
                var temp = new List<GFElement>();
                copy = baseField.RemainderPoly(copy, primitive, out temp);
                copy.RemoveRange(dimension, copy.Count - dimension);
            }

            return copy;
        }
        public bool IsEqual(ExtensionElement a, ExtensionElement b)
        {
            // implication: elements from the same field
            List<GFElement> ap = a.p, bp = b.p;
            for (int i = 0; i < ap.Count; i++)
                if (ap[i] != bp[i])
                    return false;
            return true;
        }
        public ExtensionElement Pow(ExtensionElement a, BigInteger k)
        {
            ExtensionElement res = null;
            ExtensionElement last = a;

            while (true)
            {
                if (k % 2 == 1)
                {
                    if (!(res is null))
                        res = res * last;
                    else
                        res = last;
                }

                k /= 2;
                if (k == 0)
                    break;
                last = last * last;
            }

            if (res is null)
                return new ExtensionElement(new List<GFElement> { baseField.Scalar(1) }, this);

            return res;
        }
        public ExtensionElement RandomNonZero(int seed = -1)
        {
            var rand = new Random();
            if (seed >= 0)
                rand = new Random(seed);

            var poly = new List<GFElement>(dimension);
            bool zero = true;

            while (zero)
            {
                poly.Clear();
                for (int i = 0; i < dimension; i++)
                {
                    if (rand.Next() % baseField.size == 0)
                        poly.Add(baseField.Scalar(0));
                    else
                    {
                        poly.Add(baseField.RandomNonZero(seed >= 0 ? seed + i : -1));
                        zero = false;
                    }
                }
            }

            return new ExtensionElement(poly, this);
        }

        public static List<ExtensionElement> Root(ExtensionElement a)
        {
            if (a.field.size % 2 == 1)
            {
                // quadratic residue check
                if (a.field.Pow(a, (a.field.size - 1) / 2) == new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(1) }, a.field))
                {
                    if (a.field.size % 4 == 3)
                    {
                        var root = a.field.Pow(a, (a.field.size + 1) / 4);
                        return new List<ExtensionElement> { root, -root };
                    }
                    else if (a.field.size % 4 == 1)
                    {
                        // tonelli - shanks
                        BigInteger p1 = a.field.size - 1;
                        BigInteger odd = p1;

                        int degree2 = 0;
                        while (odd % 2 == 0)
                        {
                            odd /= 2;
                            degree2++;
                        }

                        if (nonSquare is not null && a.field != nonSquare.field)
                            nonSquare = null;

                        while (nonSquare is null)
                        {
                            var re = a.field.RandomNonZero();
                            if (a.field.Pow(re, p1 / 2) == new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(-1) }, a.field))
                                nonSquare = re;
                        }

                        // precalculate all necessary powers
                        var apow = new List<ExtensionElement>();
                        var repow = new List<ExtensionElement>();
                        for (int i = 0; i < degree2; i++)
                        {
                            if (i == 0)
                            {
                                apow.Add(a.field.Pow(a, (odd + 1) / 2));
                                repow.Add(a.field.Pow(nonSquare, odd));
                            }
                            else if (i == 1)
                            {
                                apow.Add(a.field.Pow(a, odd));
                                repow.Add(repow[i - 1] * repow[i - 1]);
                            }
                            else
                            {
                                apow.Add(apow[i - 1] * apow[i - 1]);
                                repow.Add(repow[i - 1] * repow[i - 1]);
                            }
                        }

                        int n = 0;
                        BigInteger k = 1;
                        for (int i = degree2 - 1; i >= 1; i--)
                        {
                            n++;
                            if (apow[i] * a.field.Pow(repow[i], k) == new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(-1) }, a.field))
                                k += BigInteger.Pow(2, n - 1);
                        }

                        var root = apow[0] * a.field.Pow(repow[0], k);
                        return new List<ExtensionElement> { root, -root };
                    }
                }
            }

            return null;
        }
        public static List<ExtensionElement> RootBasic(ExtensionElement a)
        {
            if (a.field.size % 2 == 1)
            {
                // quadratic residue check
                if (a.field.Pow(a, (a.field.size - 1) / 2) == new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(1) }, a.field))
                {
                    if (a.field.size % 4 == 3)
                    {
                        var root = a.field.Pow(a, (a.field.size + 1) / 4);
                        return new List<ExtensionElement> { root, -root };
                    }
                    else if (a.field.size % 4 == 1)
                    {
                        // tonelli - shanks
                        BigInteger p1 = a.field.size - 1;
                        BigInteger p2 = p1 / 2;

                        ExtensionElement neg = new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(-1) }, a.field);
                        ExtensionElement re;
                        while (true)
                        {
                            re = a.field.RandomNonZero();
                            if (a.field.Pow(re, p2) == neg)
                                break;
                        }

                        // precalculate all necessary powers
                        int n = 0;
                        BigInteger k = 1;
                        BigInteger adeg, redeg;


                        while ((p2 / BigInteger.Pow(2, n)) % 2 == 0)
                        {
                            n++;
                            adeg = p2 / BigInteger.Pow(2, n);
                            redeg = p1 * k / BigInteger.Pow(2, n);

                            if (a.field.Pow(a, adeg) * a.field.Pow(re, redeg) == neg)
                                k += BigInteger.Pow(2, n - 1);
                        }

                        adeg = (p2 / BigInteger.Pow(2, n) + 1) / 2;
                        redeg = p1 * k / BigInteger.Pow(2, n + 1);

                        var root = a.field.Pow(a, adeg) * a.field.Pow(re, redeg);
                        return new List<ExtensionElement> { root, -root };
                    }
                }
            }

            return null;
        }
        public static ExtensionElement SolveCubicEquation3L(int a, GaloisFieldExtension field)
        {
            // solves equation y^3 - y = a in F_3^l

            if (field.baseField.characteristic != 3 || a > 2 || a < 0)
                throw new Exception();

            var l = field.dimension;
            var md = 3 * (l - 1) + 1;

            var m = new List<List<GFElement>>(md);
            for (int i = 0; i < md; i++)
                m.Add(new List<GFElement>(Enumerable.Repeat(field.baseField.Scalar(0), md)));

            for (int i = 0; i < l; i++)
            {
                m[3 * i][i] += field.baseField.Scalar(1);
                m[i][i] -= field.baseField.Scalar(1);
            }

            for (int i = l; i < md; i++)
                for (int j = 0; j <= l; j++)
                    m[i - l + j][i] += field.primitive[j];

            var y = new List<GFElement>(Enumerable.Repeat(field.baseField.Scalar(0), md));
            y[0] = field.baseField.Scalar(a);


            field.baseField.Gauss(m, y);

            // m[0][0] will be 0, y wil be (a, 0 ... 0)
            // find non-zero in first row, set it to a*inv
            // other zero cols correspond to zeros
            // figure out "real" variables based on non-zero "fake" variable

            for (int j = 0; j < m.Count; j++)
                if (m[0][j] != field.baseField.Scalar(0))
                {
                    var fake = a * field.baseField.Inverse(m[0][j]);
                    var res = new List<GFElement>(l);
                    res.Add(field.baseField.Scalar(0));
                    for (int i = 1; i < l; i++)
                        res.Add(-fake * m[i][j]);

                    return new ExtensionElement(res, field);
                }

            throw new Exception();
        }
        public static bool operator ==(GaloisFieldExtension a, GaloisFieldExtension b)
        {
            if (a.dimension != b.dimension || a.baseField != b.baseField)
                return false;
            for (int i = 0; i < a.dimension; i++)
                if (a.primitive[i] != b.primitive[i])
                    return false;

            return true;
        }
        public static bool operator !=(GaloisFieldExtension a, GaloisFieldExtension b)
        {
            return !(a == b);
        }
    }

    //public interface Curve<T1,T2>
    //{
    //    public T1 RightSide(T1 x);
    //    public bool IsOnTheCurve(T2 a);
    //    public T2 Action(T2 a1, T2 a2);
    //    public T2 Mult(T2 a, BigInteger k);
    //}
    public class EllipticCurve // curve of the form y^2 = x^3 + ax + b
    {
        public GaloisField field { get; private set; }
        public GFElement a { get; private set; }
        public GFElement b { get; private set; }
        public EllipticCurve(GFElement a, GFElement b)
        {
            field = a.field;
            this.a = a;
            this.b = b;
        }

        public GFElement RightSide(GFElement x)
        {
            return x * x * x + a * x + b;
        }
        public bool IsOnTheCurve(ECP a)
        {
            return a.field.IsEqual(a.y * a.y, RightSide(a.x));
        }

        private ECP Double(ECP a1)
        {
            var x1 = a1.x;
            var y1 = a1.y;

            if (field.characteristic == 2)
            {
                throw new Exception($"Action not defined for field with characteristic {field.characteristic}");
            }
            else if (field.characteristic == 3)
            {
                //var y2 = 2 * y1;
                //var y2inv = field.Inverse(y2);
                //var check = y2 * y2inv;
                //Console.WriteLine("check");
                //Program.Print(check);
                //Console.WriteLine();
                //var test = a * y2inv;

                var s = a / (2 * y1);
                var x3 = s * s - 2 * x1;
                var y3 = s * (x1 - x3) - y1;
                return new ECP(x3, y3);
            }
            else
            {
                var s = (3 * x1 * x1 + a) / (2 * y1);
                var x3 = s * s - 2 * x1;
                var y3 = s * (x1 - x3) - y1;
                return new ECP(x3, y3);
            }

        }
        private ECP Add(ECP a1, ECP a2)
        {
            var x1 = a1.x;
            var x2 = a2.x;
            var y1 = a1.y;
            var y2 = a2.y;

            switch (field.characteristic)
            {
                case 2:
                    throw new Exception($"Action not defined for field with characteristic {field.characteristic}");
                default: // the same for characteristic >= 3
                    var s = (y2 - y1) / (x2 - x1);
                    var x3 = s * s - x2 - x1;
                    var y3 = s * (x1 - x3) - y1;
                    return new ECP(x3, y3);
            }
        }

        public ECP Action(ECP a1, ECP a2)
        {
            if (a1.field != a2.field)
                throw new Exception();

            var field = a1.field;

            var x1 = a1.x;
            var x2 = a2.x;
            var y1 = a1.y;
            var y2 = a2.y;

            if (field.IsEqual(x1, x2))
            {
                if (field.IsEqual(y1, -y2))
                {
                    return null;
                }
                else
                {
                    // double
                    return Double(a1);
                }
            }
            else
            {
                // add
                return Add(a1, a2);
            }
        }

        public ECP Mult(ECP a, BigInteger k)
        {
            ECP res = null;
            ECP last = a;

            while (last is not null)
            {
                if (k % 2 == 1)
                {
                    if (res is not null)
                        res = Action(res, last);
                    else
                        res = last;
                }

                k /= 2;
                if (k == 0)
                    break;
                last = Action(last, last);
                //Program.Print(last.x);
            }

            return res;
        }
    }
    public class EllipticCurveExtension
    {
        public GaloisFieldExtension field { get; private set; }
        public ExtensionElement a { get; private set; }
        public ExtensionElement b { get; private set; }
        public EllipticCurveExtension(ExtensionElement a, ExtensionElement b)
        {
            field = a.field;
            this.a = a;
            this.b = b;
        }
        private ECPExtension Double(ECPExtension a1)
        {
            var x1 = a1.x;
            var y1 = a1.y;

            if (field.baseField.characteristic == 2)
            {
                throw new Exception($"Action not defined for field with characteristic {field.baseField.characteristic}");
            }
            else if (field.baseField.characteristic == 3)
            {
                //var y2 = 2 * y1;
                //var y2inv = field.Inverse(y2);
                //var check = y2 * y2inv;
                //Console.WriteLine("check");
                //Program.Print(check);
                //Console.WriteLine();
                //var test = a * y2inv;

                var s = a / (field.baseField.Scalar(2) * y1);
                var x3 = s * s - field.baseField.Scalar(2) * x1;
                var y3 = s * (x1 - x3) - y1;
                return new ECPExtension(x3, y3);
            }
            else
            {
                var s = (field.baseField.Scalar(3) * x1 * x1 + a) / (field.baseField.Scalar(2) * y1);
                var x3 = s * s - field.baseField.Scalar(2) * x1;
                var y3 = s * (x1 - x3) - y1;
                return new ECPExtension(x3, y3);
            }

        }
        private ECPExtension Add(ECPExtension a1, ECPExtension a2)
        {
            var x1 = a1.x;
            var x2 = a2.x;
            var y1 = a1.y;
            var y2 = a2.y;

            switch (field.baseField.characteristic)
            {
                case 2:
                    throw new Exception($"Action not defined for field with characteristic {field.baseField.characteristic}");
                default: // the same for characteristic >= 3
                    var s = (y2 - y1) / (x2 - x1);
                    var x3 = s * s - x2 - x1;
                    var y3 = s * (x1 - x3) - y1;
                    return new ECPExtension(x3, y3);
            }
        }

        public ECPExtension Action(ECPExtension a1, ECPExtension a2)
        {
            if (a1.field != a2.field)
                throw new Exception();

            var x1 = a1.x;
            var x2 = a2.x;
            var y1 = a1.y;
            var y2 = a2.y;

            if (x1 == x2)
            {
                if (y1 == -y2)
                {
                    return null;
                }
                else
                {
                    // double
                    return Double(a1);
                }
            }
            else
            {
                // add
                return Add(a1, a2);
            }
        }

        public bool IsOnTheCurve(ECPExtension a)
        {
            return a.y * a.y == RightSide(a.x);
        }

        public ECPExtension Mult(ECPExtension a, BigInteger k)
        {
            ECPExtension res = null;
            ECPExtension last = a;

            while (last is not null)
            {
                Console.WriteLine(k);
                if (k % 2 == 1)
                {
                    if (res is not null)
                        res = Action(res, last);
                    else
                        res = last;
                }

                k /= 2;
                if (k == 0)
                    break;
                last = Action(last, last);
                //Program.Print(last.x);
            }

            return res;
        }

        public ExtensionElement RightSide(ExtensionElement x)
        {
            return x * x * x + a * x + b;
        }
    }

    public class EllipticCurveManager
    {
        public EllipticCurve curve { get; private set; }
        public EllipticCurveExtension curve6 { get; private set; }

        //GFElement root;

        BigInteger? m;
        public BigInteger? q { get; private set; }

        ExtensionElement u;
        ExtensionElement rp;
        ExtensionElement rm;

        static EllipticCurveExtension lastCurve;
        static BigInteger order;
        static ECPExtension s;
        private class FieldData
        {
            public FieldData(int dimension, Dictionary<int, int> p1)
            {
                this.dimension = dimension;
                this.p1 = p1;
            }

            public int dimension;
            public Dictionary<int, int> p1;
        }
        public EllipticCurveManager(List<int> a, List<int> b, int l, bool positive)
        {
            // construct necessary fields
            // smaller field requires primitive polynomial to convert elements
            GaloisField field;
            GaloisFieldExtension field6;
            (field, field6) = InitFields(l, positive);
            InitMQ(field, positive);

            curve = new EllipticCurve(new GFElement(a, field), new GFElement(b, field));
            curve6 = new EllipticCurveExtension(new(new List<GFElement> { curve.a }, field6), new(new List<GFElement> { curve.b }, field6));
        }
        private BigInteger BonehCurveSize(GaloisField field, bool positive)
        {
            if (positive)
            {
                if (field.dimension % 12 == 1 || field.dimension % 12 == 11)
                    return BigInteger.Pow(3, field.dimension) + BigInteger.Pow(3, (field.dimension + 1) / 2) + 1;
                else if (field.dimension % 12 == 5 || field.dimension % 12 == 7)
                    return BigInteger.Pow(3, field.dimension) - BigInteger.Pow(3, (field.dimension + 1) / 2) + 1;
                else
                    throw new Exception("Wrong field dimension");
            }
            else
            {
                if (field.dimension % 12 == 1 || field.dimension % 12 == 11)
                    return BigInteger.Pow(3, field.dimension) - BigInteger.Pow(3, (field.dimension + 1) / 2) + 1;
                else if (field.dimension % 12 == 5 || field.dimension % 12 == 7)
                    return BigInteger.Pow(3, field.dimension) + BigInteger.Pow(3, (field.dimension + 1) / 2) + 1;
                else
                    throw new Exception("Wrong field dimension");
            }
        }
        private void InitMQ(GaloisField field, bool positive)
        {
            m = BonehCurveSize(field, positive);
            var factors = PrimeFactors(m ?? -1);
            q = factors?.Max<BigInteger>();

            if (q == null)
                Console.WriteLine($"Factorization for m = {m} missing, q is set to null");
        }
        private (GaloisField field, GaloisFieldExtension field6) InitFields(int l, bool positive)
        {
            GaloisField field = null;
            GaloisFieldExtension field6 = null;

            string fileName = "fields.json";
            List<FieldData> source = null;

            JsonSerializerOptions options = new();
            options.IncludeFields = true;

            using (StreamReader r = new StreamReader(fileName))
            {
                string json = r.ReadToEnd();
                source = JsonSerializer.Deserialize<List<FieldData>>(json, options);
                source.RemoveAll((FieldData data) => data.dimension == 0);

                for (int i = 0; i < source.Count; i++)
                {
                    var data = source[i];
                    if (data.dimension == l)
                    {
                        var p = new List<int>(Enumerable.Repeat(0, data.dimension + 1));
                        var p6 = new List<int>(Enumerable.Repeat(0, 7));

                        foreach (var key in data.p1.Keys)
                            p[key] = data.p1[key];
                        p6[0] = 2;
                        p6[1] = 1;
                        p6[6] = 1;

                        field = new GaloisField(3, p);
                        field6 = new GaloisFieldExtension(field, p6);

                        // need to define automorphism
                        var tempField = new GaloisField(3, p6);

                        var y2 = GaloisField.Root(new GFElement(new List<int> { -1 }, tempField))[0];
                        var y3 = GaloisField.SolveCubicEquation3L(positive ? 1 : 2, tempField);

                        u = new((from e in y2.p select field.Scalar(e)).ToList(), field6);
                        rm = new((from e in y3.p select field.Scalar(e)).ToList(), field6);
                        rp = -rm;

                        break;
                    }
                }
            }

            return (field, field6);
        }
        private List<BigInteger> PrimeFactors(BigInteger m)
        {
            var res = new List<BigInteger>();

            // possible group sizes for elliptic curves
            if (m == BigInteger.Parse("2269"))
            {
                res.Add(BigInteger.Parse("2269"));
                return res;
            }
            else if (m == BigInteger.Parse("49269609804781974450852068861184694669"))
            {
                res.Add(BigInteger.Parse("49269609804781974450852068861184694669"));
                return res;
            }
            // possible non-zero elements of a field
            else if (m == BigInteger.Pow(3, 7) - 1)
            {
                res.Add(BigInteger.Parse("2"));
                res.Add(BigInteger.Parse("1093"));
                return res;
            }
            else if (m == BigInteger.Pow(3, 42) - 1)
            {
                res.Add(BigInteger.Parse("2"));
                res.Add(BigInteger.Parse("7"));
                res.Add(BigInteger.Parse("13"));
                res.Add(BigInteger.Parse("547"));
                res.Add(BigInteger.Parse("1093"));
                res.Add(BigInteger.Parse("2269"));
                res.Add(BigInteger.Parse("368089"));
                return res;
            }
            else if (m == BigInteger.Pow(3, 79) - 1)
            {
                res.Add(BigInteger.Parse("2"));
                res.Add(BigInteger.Parse("432853009"));
                res.Add(BigInteger.Parse("392038110671"));
                res.Add(BigInteger.Parse("145171177264407947"));
                return res;
            }
            else if (m == BigInteger.Pow(3, 474) - 1)
            {
                res.Add(BigInteger.Parse("2"));
                res.Add(BigInteger.Parse("7"));
                res.Add(BigInteger.Parse("13"));
                res.Add(BigInteger.Parse("14221"));
                res.Add(BigInteger.Parse("66361"));
                res.Add(BigInteger.Parse("2275201"));
                res.Add(BigInteger.Parse("37520893"));
                res.Add(BigInteger.Parse("432853009"));
                res.Add(BigInteger.Parse("327220181191"));
                res.Add(BigInteger.Parse("392038110671"));
                res.Add(BigInteger.Parse("145171177264407947"));
                res.Add(BigInteger.Parse("567239060331150635317"));
                res.Add(BigInteger.Parse("217536018992883276988362961"));
                res.Add(BigInteger.Parse("52802131282482966001751707837"));
                res.Add(BigInteger.Parse("49269609804781974450852068861184694669"));
                res.Add(BigInteger.Parse("94251934933904430405353683133566414129"));
                return res;
            }

            return null;
        }
        private List<int> FindIrreducible(int fp, int n)
        {
            var primeField = new PrimeField(fp);

            var nfactors = Methods.PrimeFactors(n);
            var np = new List<int>();
            foreach (int p in nfactors)
                np.Add(n / p);

            var irreducible = new List<int>(Enumerable.Repeat(0, n + 1));
            irreducible[n] = 1;

            int iteration = 1;
            while (true)
            {
                Console.WriteLine($"Iteration: {iteration}");
                iteration++;

                for (int i = 0; i < n; i++)
                {
                    irreducible[i] = primeField.Add(irreducible[i], 1);
                    if (irreducible[0] == 0)
                    {
                        irreducible[0] = 1;
                        continue;
                    }
                    if (irreducible[i] != 0)
                        break;
                }

                //irreducible = new List<int> { 1, 2, 0, 0, 0, 1 };

                // need to compute gcd(x^p^n - x, f)
                // first step can't be computed right away, need to get rid of the large degree
                // solution: compute all remainders of x^p^k, 1 <= k <= np step by step
                // check if gcd(x^p^n/p - x, f) == 1 at the same time
                bool reducible = false;
                var q = new List<int>();
                var f = new PrimePolynomial(irreducible, primeField);
                var x = new PrimePolynomial(new List<int> { 0, 1 }, primeField);
                var k = x % f; // start
                for (int i = 0; i < n; i++)
                {
                    if (np.Contains(i))
                    {
                        var kx = k - x % f; // this is now (x^p^np - x) % f
                                            // run euclid
                        PrimePolynomial t, s;
                        var gcd = primeField.EuclidPoly(kx, f, out t, out s);
                        if (primeField.Degree(gcd.poly) > 1)
                        {
                            reducible = true;
                            break;
                        }
                    }

                    k = PrimePolynomial.Pow(k, fp);
                    k %= f;
                }

                //Console.Write(reducible.ToString() + " ");
                //Print(irreducible);

                if (!reducible) // if (x^p^n - x) % f == 0 => f irreducible
                {
                    var xpn = k - x % f; // has to be zero if f is irreducible
                    if (primeField.Degree(xpn.poly) == 0 && xpn.poly[0] == 0)
                        return irreducible;
                }
            }
        }
        public ECP RandomPointOfPrimeOrder(BigInteger order)
        {
            EllipticCurve c = curve;

            ECP p = null;
            while (true)
            {
                GFElement px;
                GFElement py;

                px = c.field.RandomNonZero();
                var roots = GaloisField.Root(c.RightSide(px));
                if (roots != null)
                    py = roots[0];
                else
                    continue;

                p = new ECP(px, py);
                var pq = c.Mult(p, order);

                if (pq is null)
                    break;
            }

            return p;
        }
        class HPQ
        {
            ECPExtension p;
            ECPExtension q;
            ExtensionElement lambda;

            public HPQ(ECPExtension p, ECPExtension q, ExtensionElement a)
            {
                this.p = p;
                this.q = q;

                if (p.field.IsEqual(p.x, q.x))
                {
                    if (p.field.IsEqual(p.y, q.y))
                    {
                        // tangent
                        if (p.y == new ExtensionElement(new List<GFElement> { p.field.baseField.Scalar(0) }, p.field))
                            lambda = null;
                        else if (p.field.baseField.characteristic == 3)
                            lambda = a / (p.field.baseField.Scalar(2) * p.y);
                        else if (p.field.baseField.characteristic > 3)
                            lambda = (p.field.baseField.Scalar(3) * p.x * p.x + a) / (p.field.baseField.Scalar(2) * p.y);
                        else
                            throw new Exception("Slope not defined");
                    }
                    else
                    {
                        // vertical
                        lambda = null;
                    }
                }
                else
                {
                    // slope
                    lambda = (q.y - p.y) / (q.x - p.x);
                }
            }

            public ExtensionElement Value(ECPExtension a)
            {
                if (lambda is not null)
                {
                    var num = a.y - p.y - lambda * (a.x - p.x);
                    var den = a.x + p.x + q.x - lambda * lambda; // some terms skipped cause curve is simple
                    return num / den;
                }
                else
                    return a.x - p.x;
            }
        }

        static ExtensionElement EvaluateF(List<HPQ> f, ECPExtension a)
        {
            var res = new ExtensionElement(new List<GFElement> { a.field.baseField.Scalar(1) }, a.field);
            foreach (HPQ h in f)
            {
                if (h is null)
                    res = res * res;
                else
                    res = res * h.Value(a);
            }

            return res;
        }
        static List<HPQ> Miller(ECPExtension p, BigInteger n, EllipticCurveExtension curve)
        {
            var binary = new List<int>();
            while (n != 1)
            {
                binary.Add((int)(n % 2));
                n /= 2;
            }
            binary.Add(1);

            var t = binary.Count - 2;

            var f = new List<HPQ>();
            var T = p;
            for (int i = t; i >= 0; i--)
            {
                f.Add(null);
                f.Add(new HPQ(T, T, curve.a));

                T = curve.Action(T, T);
                if (binary[i] == 1)
                {
                    f.Add(new HPQ(T, p, curve.a));
                    T = curve.Action(T, p);
                }
            }

            return f;
        }
        static ExtensionElement MillerUpdated(ECPExtension p, ECPExtension r, BigInteger n, EllipticCurveExtension curve)
        {
            var binary = new List<int>();
            while (n != 1)
            {
                binary.Add((int)(n % 2));
                n /= 2;
            }
            binary.Add(1);
            //binary.Reverse();

            var t = binary.Count - 2;

            var f = new ExtensionElement(new List<GFElement> { curve.field.baseField.Scalar(1) }, curve.field);
            var T = p;
            //int degree = 1;
            for (int i = t; i >= 0; i--)
            {
                Console.WriteLine(i);
                f *= f;
                //degree *= 2;
                f *= new HPQ(T, T, curve.a).Value(r);

                T = curve.Action(T, T);
                

                if (binary[i] == 1)
                {
                    f *= new HPQ(T, p, curve.a).Value(r);
                    //degree += 1;
                    T = curve.Action(T, p);
                }
            }
            //Console.WriteLine(degree);

            return f;
        }
        static public ExtensionElement WeilPairing(ECPExtension p, ECPExtension q, BigInteger n, EllipticCurveExtension curve)
        {
            if (n != order || curve != lastCurve)
            {
                lastCurve = curve;
                order = n;
                s = null;
            }

            while (s is null)
            {
                Console.WriteLine("trying to find s");
                var sx = p.field.RandomNonZero();
                var sy2 = curve.RightSide(sx);
                var sy = GaloisFieldExtension.Root(sy2);

                if (sy == null)
                    continue;
                else
                    Console.WriteLine("square");

                s = new ECPExtension(sx, sy[0]);

                if (curve.Mult(s, n) is not null)
                    break;
            }
            Console.WriteLine("Found s");
            //var ns = curve.Mult(s, n);
            //Program.Print(MillerUpdated(s, ns, n, curve));


            //Program.Print(s.x);
            //Program.Print(s.y);
            //s = new(new(new List<GFElement> { curve.field.baseField.Scalar(0) }, curve.field), new(new List<GFElement> { curve.field.baseField.Scalar(36) }, curve.field));

            //var fp = Miller(p, n, curve);
            //var fq = Miller(q, n, curve);
            var sn = new ECPExtension(s.x, -s.y);
            //Console.WriteLine(curve.IsOnTheCurve(sn));
            //Console.WriteLine(curve.IsOnTheCurve(s));

            //Program.Print(EvaluateF(fp, p));
            //Program.Print(EvaluateF(fp, q));
            //var num = EvaluateF(fp, curve.Action(q, s)) / EvaluateF(fp, s);
            //var den = EvaluateF(fq, curve.Action(p, sn)) / EvaluateF(fq, sn);

            var num = MillerUpdated(p, curve.Action(q, s), n, curve) / MillerUpdated(p, s, n, curve);
            var den = MillerUpdated(q, curve.Action(p, sn), n, curve) / MillerUpdated(q, sn, n, curve);
            return num / den;
        }
        public ECP MapToGroup(byte[] m, float delta)
        {
            int I = (int)Math.Ceiling(Math.Log2(Math.Log2(1 / delta)));
            int maxIter = (int)Math.Pow(2, I);

            var hashAlg = SHA256.Create();
            for (int iter = 0; iter < maxIter; iter++)
            {
                BigInteger h = new BigInteger(m);
                h += iter;
                var hash = hashAlg.ComputeHash(h.ToByteArray());

                h = new BigInteger(hash);
                //Console.WriteLine(h);

                BigInteger mod = curve.field.size * 2;
                h %= mod;

                bool b = h < curve.field.size;
                if (!b)
                    h -= curve.field.size;

                List<int> poly = new List<int>(curve.field.dimension);
                for (int i = 0; i < curve.field.dimension; i++)
                {
                    poly.Add((int)(h % 3));
                    h /= 3;
                }

                var x = new GFElement(poly, curve.field);

                var roots = GaloisField.Root(curve.RightSide(x));
                if (roots == null)
                    continue;

                if (roots[0].p[0] > roots[1].p[0])
                {
                    // swap
                    var temp = roots[0];
                    roots[0] = roots[1];
                    roots[1] = temp;
                }

                int index = 0;
                if (b)
                    index = 1;
                var pm = new ECP(x, roots[index]);

                var mult = this.m / this.q; // must know m and q
                pm = curve.Mult(pm, mult ?? throw new Exception("m or q is null"));

                if (pm is not null)
                    return pm;
            }

            throw new Exception();
        }
        public ECPExtension Auto6L(ECPExtension a)
        {
            // check if curve is positive for small curve, that's ok
            if (curve.b == new GFElement(new List<int> { 1 }, curve.field))
                return new ECPExtension(-a.x + rp, u * a.y); // +++
            else
                return new ECPExtension(-a.x + rm, u * a.y); // ---
        }
    }
}

