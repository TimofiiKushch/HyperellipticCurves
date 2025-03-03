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
    public class GFElement<T>
    {
        public List<T> p { get; }
        public GaloisField<T> field { get; }

        public GFElement(List<T> polynomial, GaloisField<T> field)
        {
            p = field.FixRepresentation(polynomial);
            this.field = field;
        }

        public static GFElement<T> operator -(GFElement<T> a)
        {
            return a.field.Subtract(a.field.Scalar(0), a);
        }

        public static GFElement<T> operator +(GFElement<T> a, GFElement<T> b)
        {
            if (a.field == b.field)
                return a.field.Add(a, b);
            else
                return null;
        }

        public static GFElement<T> operator -(GFElement<T> a, GFElement<T> b)
        {
            if (a.field == b.field)
                return a.field.Subtract(a, b);
            else
                return null;
        }

        public static GFElement<T> operator *(GFElement<T> a, GFElement<T> b)
        {
            if (a.field == b.field)
                return a.field.Multiply(a, b);
            else
                return null;
        }
        public static GFElement<T> operator *(int a, GFElement<T> b)
        {
            return b.field.Multiply(b.field.Scalar(a), b);
        }

        public static GFElement<T> operator /(GFElement<T> a, GFElement<T> b)
        {
            if (a.field == b.field)
            {
                return a * a.field.Inverse(b);
            }
            else
                return null;
        }

        public static bool operator ==(GFElement<T> a, GFElement<T> b)
        {
            if (a.field != b.field)
                return false;

            return a.field.IsEqual(a, b);
        }

        public static bool operator !=(GFElement<T> a, GFElement<T> b)
        {
            return !(a == b);
        }
    }
    public class GaloisField<T> : IField<GFElement<T>>
    {
        //public int characteristic { get; }
        public int dimension { get; }
        readonly public List<T> primitive = null;
        readonly public IField<T> baseField;

        public BigInteger size { get; }

        public GaloisField(IField<T> baseField, List<int> primitive)
        {
            this.baseField = baseField;

            this.primitive = new List<T>();
            foreach (int i in primitive)
                this.primitive.Add(baseField.Scalar(i));

            dimension = primitive.Count - 1;
            size = BigInteger.Pow(baseField.Size(), dimension);
        }
        public GaloisField(IField<T> baseField, List<T> primitive)
        {
            this.baseField = baseField;
            this.primitive = primitive;
            dimension = primitive.Count - 1;
            size = BigInteger.Pow(baseField.Size(), dimension);
        }

        public GFElement<T> Scalar(int a)
        {
            return new(new List<T> { baseField.Scalar(a) }, this);
        }
        public GFElement<T> Inverse(GFElement<T> a)
        {
            if (Degree(a.p) == 0)
                return new GFElement<T>(new List<T> { baseField.Inverse(a.p[0]) }, this);

            List<T> inv, s;
            var gcd = EuclidPoly(a.p, primitive, out inv, out s);
            inv = MultiplyPoly(new List<T> { baseField.Inverse(gcd[0]) }, inv);
            return new GFElement<T>(inv, this);
        }
        public GFElement<T> Pow(GFElement<T> a, BigInteger k)
        {
            GFElement<T> res = null;
            GFElement<T> last = a;

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
                return Scalar(1);

            return res;
        }
        public static GFElement<T> Root(GFElement<T> a)
        {
            if (a.field.size % 2 == 1)
            {
                // quadratic residue check
                if (a.field.Pow(a, (a.field.size - 1) / 2) == a.field.Scalar(1))
                {
                    if (a.field.size % 4 == 3)
                    {
                        var root = a.field.Pow(a, (a.field.size + 1) / 4);
                        return root;
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

                        GFElement<T> re;
                        while (true)
                        {
                            re = a.field.RandomNonZero();
                            if (a.field.Pow(re, p1 / 2) == a.field.Scalar(-1))
                                break;
                        }

                        // precalculate all necessary powers
                        var apow = new List<GFElement<T>>();
                        var repow = new List<GFElement<T>>();
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
                        return root;
                    }
                }
            }

            return null;
        }
        public GFElement<T> RandomNonZero(int seed = -1)
        {
            var rand = new Random();
            if (seed >= 0)
                rand = new Random(seed);

            var poly = new List<T>(dimension);
            bool zero = true;

            while (zero)
            {
                poly.Clear();
                for (int i = 0; i < dimension; i++)
                {
                    poly.Add(baseField.RandomNonZero());
                    if (!baseField.IsEqual(poly[i], baseField.Scalar(0)))
                        zero = false;
                }
            }

            return new GFElement<T>(poly, this);
        }
        public List<T> FixRepresentation(List<T> polynomial)
        {
            List<T> copy = new List<T>(polynomial);
            if (typeof(T) == typeof(int))
                for (int i = 0; i < copy.Count; i++)
                    copy[i] = (T)(object)Methods.NumRemainder((int)(object)copy[i], Characteristic());

            if (copy.Count > dimension)
            {
                var temp = new List<T>();
                copy = RemainderPoly(copy, primitive, out temp);
            }

            while (copy.Count < dimension)
                copy.Add(baseField.Scalar(0));

            return copy;
        }
        public bool IsEqual(GFElement<T> a, GFElement<T> b)
        {
            // implication: elements from the same field
            List<T> ap = a.p, bp = b.p;
            for (int i = 0; i < ap.Count; i++)
                if (!baseField.IsEqual(ap[i], bp[i]))
                    return false;
            return true;
        }
        public int Degree(List<T> polynomial)
        {
            var zero = baseField.Scalar(0);
            var res = polynomial.FindLastIndex((T el) => !baseField.IsEqual(el, zero));
            return res < 0 ? 0 : res;
        }
        public T LeadingCoeff(List<T> polynomial)
        {
            if (polynomial.Count == 0)
                return baseField.Scalar(0);
            return polynomial[Degree(polynomial)];
        }
        public GFElement<T> Add(GFElement<T> a, GFElement<T> b)
        {
            return new GFElement<T>(AddPoly(a.p, b.p), this);
        }
        public GFElement<T> Subtract(GFElement<T> a, GFElement<T> b)
        {
            return new GFElement<T>(SubtractPoly(a.p, b.p), this);
        }
        public GFElement<T> Multiply(GFElement<T> a, GFElement<T> b)
        {
            return new GFElement<T>(MultiplyPoly(a.p, b.p), this);
        }

        List<T> AddPoly(List<T> a, List<T> b)
        {
            List<T> c = new List<T>();

            for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
                c.Add(baseField.Add(a[i], b[i]));

            if (a.Count > b.Count)
                for (int i = b.Count; i < a.Count; i++)
                    c.Add(a[i]);
            else
                for (int i = a.Count; i < b.Count; i++)
                    c.Add(b[i]);

            return c;
        }
        List<T> SubtractPoly(List<T> a, List<T> b)
        {

            List<T> c = new List<T>();

            for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
                c.Add(baseField.Subtract(a[i], b[i]));

            if (a.Count > b.Count)
                for (int i = b.Count; i < a.Count; i++)
                    c.Add(a[i]);
            else
                for (int i = a.Count; i < b.Count; i++)
                    c.Add(-(dynamic)b[i]);

            return c;
        }
        List<T> MultiplyPoly(List<T> a, List<T> b)
        {
            List<T> c = new List<T>(a.Count + b.Count - 1);
            for (int i = 0; i < a.Count + b.Count - 1; i++)
                c.Add(baseField.Scalar(0));
            for (int i = 0; i < a.Count; i++)
                for (int j = 0; j < b.Count; j++)
                    c[i + j] = baseField.Add(c[i + j], baseField.Multiply(a[i], b[j]));

            return c;
        }
        public List<T> RemainderPoly(List<T> a, List<T> b, out List<T> factor)
        {
            List<T> copy = new List<T>(a.Count);
            factor = new List<T>(a.Count);
            foreach (var el in a)
            {
                copy.Add(baseField.Clone(el));
                factor.Add(baseField.Scalar(0));
            }

            var zero = baseField.Scalar(0);
            int cur = copy.FindLastIndex((T el) => !baseField.IsEqual(el, zero));
            int bf = b.FindLastIndex((T el) => !baseField.IsEqual(el, zero));

            T inv = baseField.Inverse(b[bf]);

            while (cur >= bf)
            {
                int diff = cur - bf;
                if (!baseField.IsEqual(copy[cur], zero))
                {
                    for (int i = 0; i < bf; i++)
                        copy[diff + i] = baseField.Subtract(copy[diff + i], baseField.Multiply(baseField.Multiply(copy[cur], inv), b[i]));
                    factor[diff] = baseField.Multiply(copy[cur], inv);
                    copy[cur] = baseField.Clone(zero);
                }
                cur--;
            }

            int last = copy.FindLastIndex((T el) => !baseField.IsEqual(el, zero));
            if (last < copy.Count - 1)
                copy.RemoveRange(last + 1, copy.Count - last - 1);

            last = factor.FindLastIndex((T el) => !baseField.IsEqual(el, zero));
            if (last < factor.Count - 1)
                factor.RemoveRange(last + 1, factor.Count - last - 1);

            return copy;
        }
        public List<T> EuclidPoly(List<T> a, List<T> b, out List<T> t, out List<T> s)
        {
            List<T>[] sl = new List<T>[]{ new List<T> { baseField.Scalar(0) }, new List<T> { baseField.Scalar(1) } };
            List<T>[] tl = new List<T>[]{ new List<T> { baseField.Scalar(1) }, new List<T> { baseField.Scalar(0) } };

            var ac = new List<T>(a.Count);
            foreach (var el in a)
                ac.Add(baseField.Clone(el));

            var bc = new List<T>(b.Count);
            foreach (var el in b)
                bc.Add(baseField.Clone(el));

            var zero = baseField.Scalar(0);
            while (!baseField.IsEqual(LeadingCoeff(ac), zero) && !baseField.IsEqual(LeadingCoeff(bc), zero))
            {
                List<T> q = new List<T>();
                List<T> c = ac;
                ac = RemainderPoly(bc, ac, out q);
                bc = c;

                sl = new List<T>[] { SubtractPoly(sl[1], MultiplyPoly(q, sl[0])), sl[0] };
                tl = new List<T>[] { SubtractPoly(tl[1], MultiplyPoly(q, tl[0])), tl[0] };
            }

            t = tl[1];
            s = sl[1];
            if (!baseField.IsEqual(LeadingCoeff(bc), zero))
                return bc;
            else
                return new List<T> { baseField.Scalar(1) };
        }

        public GFElement<T> Clone(GFElement<T> a)
        {
            var p = new List<T>();
            foreach (var el in a.p)
                p.Add(baseField.Clone(el));
            return new GFElement<T>(p, this);
        }

        public BigInteger Size()
        {
            return size;
        }

        public int Characteristic()
        {
            return baseField.Characteristic();
        }

        public bool IsEqual(IField<GFElement<T>> b)
        {
            return this == (GaloisField<T>)b;
        }

        public static bool operator ==(GaloisField<T> a, GaloisField<T> b)
        {
            if (a.dimension != b.dimension || !a.baseField.IsEqual(b.baseField))
                return false;
            for (int i = 0; i < a.dimension; i++)
                if (!a.baseField.IsEqual(a.primitive[i], b.primitive[i]))
                    return false;

            return true;
        }
        public static bool operator !=(GaloisField<T> a, GaloisField<T> b)
        {
            return !(a == b);
        }
    }
    public class ECP<T>
    {
        public GFElement<T> x { get; }
        public GFElement<T> y { get; }

        public GaloisField<T> field { get; }

        public ECP(GFElement<T> x, GFElement<T> y)
        {
            if (x.field == y.field)
            {
                this.x = x;
                this.y = y;
                field = x.field;
            }
        }

        public static bool operator ==(ECP<T> a, ECP<T> b)
        {
            return a.x == b.x && a.y == b.y;
        }

        public static bool operator !=(ECP<T> a, ECP<T> b)
        {
            return !(a == b);
        }
    }
    public class EllipticCurve<T> // curve of the form y^2 = x^3 + ax + b
    {
        public GaloisField<T> field { get; private set; }
        public GFElement<T> a { get; private set; }
        public GFElement<T> b { get; private set; }
        public EllipticCurve(GFElement<T> a, GFElement<T> b)
        {
            field = a.field;
            this.a = a;
            this.b = b;
        }

        public GFElement<T> RightSide(GFElement<T> x)
        {
            return x * x * x + a * x + b;
        }
        public bool IsOnTheCurve(ECP<T> a)
        {
            if (a is null)
                return true;
            return a.field.IsEqual(a.y * a.y, RightSide(a.x));
        }
        private ECP<T> Double(ECP<T> a1)
        {
            var x1 = a1.x;
            var y1 = a1.y;

            if (field.Characteristic() == 2)
            {
                throw new Exception($"Action not defined for field with characteristic {field.Characteristic()}");
            }
            else if (field.Characteristic() == 3)
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
                return new ECP<T>(x3, y3);
            }
            else
            {
                var s = (3 * x1 * x1 + a) / (2 * y1);
                var x3 = s * s - 2 * x1;
                var y3 = s * (x1 - x3) - y1;
                return new ECP<T>(x3, y3);
            }

        }
        private ECP<T> Add(ECP<T> a1, ECP<T> a2)
        {
            var x1 = a1.x;
            var x2 = a2.x;
            var y1 = a1.y;
            var y2 = a2.y;

            switch (field.Characteristic())
            {
                case 2:
                    throw new Exception($"Action not defined for field with characteristic {field.Characteristic()}");
                default: // the same for characteristic >= 3
                    var s = (y2 - y1) / (x2 - x1);
                    var x3 = s * s - x2 - x1;
                    var y3 = s * (x1 - x3) - y1;
                    return new ECP<T>(x3, y3);
            }
        }
        public ECP<T> Action(ECP<T> a1, ECP<T> a2)
        {
            if (a1.field != a2.field)
                throw new Exception();

            var field = a1.field;

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
        public ECP<T> Mult(ECP<T> a, BigInteger k)
        {
            ECP<T> res = null;
            ECP<T> last = a;

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
                //Console.WriteLine(k);
                if (k == 0)
                    break;
                last = Action(last, last);
                //Program.Print(last.x);
            }

            return res;
        }
        public static bool operator ==(EllipticCurve<T> a, EllipticCurve<T> b)
        {
            return a.field == b.field && a.a == b.a && a.b == b.b;
        }
        public static bool operator !=(EllipticCurve<T> a, EllipticCurve<T> b)
        {
            return !(a == b);
        }
    }

    public class EllipticCurveManager
    {
        public EllipticCurve<int> curve { get; private set; }
        public EllipticCurve<GFElement<int>> curve6 { get; private set; }

        BigInteger m;
        public BigInteger q { get; private set; }
        public ECP<int> pq { get; private set; }
        public ECP<GFElement<int>> s { get; private set; }

        GFElement<GFElement<int>> u;
        GFElement<GFElement<int>> rp;
        GFElement<GFElement<int>> rm;

        static EllipticCurve<GFElement<int>> lastCurve;
        static BigInteger order;
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
        private class CurveData
        {
            public CurveData(int l, bool positive, string q, List<int> pqx, List<int> pqy, List<List<int>> sx, List<List<int>> sy)
            {
                this.l = l;
                this.positive = positive;
                this.q = q;
                this.pqx = pqx;
                this.pqy = pqy;
                this.sx = sx;
                this.sy = sy;
            }

            public int l;
            public bool positive;
            public string q;
            public List<int> pqx;
            public List<int> pqy;
            public List<List<int>> sx;
            public List<List<int>> sy;
        }
        public EllipticCurveManager(List<int> a, List<int> b, int l, bool positive, bool initCurveData = false)
        {
            // construct necessary fields
            // smaller field requires primitive polynomial to convert elements
            GaloisField<int> field;
            GaloisField<GFElement<int>> field6;
            (field, field6) = InitFields(l, positive);
            curve = new(new GFElement<int>(a, field), new GFElement<int>(b, field));
            curve6 = new(new(new List<GFElement<int>> { curve.a }, field6), new(new List<GFElement<int>> { curve.b }, field6));

            if (initCurveData)
                InitCurveData(field, field6, positive);
        }
        private BigInteger BonehCurveSize(GaloisField<int> field, bool positive)
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
        private void InitCurveData(GaloisField<int> field, GaloisField<GFElement<int>> field6, bool positive)
        {
            m = BonehCurveSize(field, positive);
            
            string fileName = "curves.json";
            List<CurveData> source = null;

            JsonSerializerOptions options = new();
            options.IncludeFields = true;

            using (StreamReader r = new StreamReader(fileName))
            {
                string json = r.ReadToEnd();
                source = JsonSerializer.Deserialize<List<CurveData>>(json, options);
                source.RemoveAll((CurveData data) => data.l == 0);

                for (int i = 0; i < source.Count; i++)
                {
                    var data = source[i];
                    if (data.l == field.dimension && data.positive == positive)
                    {
                        q = BigInteger.Parse(data.q);
                        pq = new(new(data.pqx, field), new(data.pqy, field));
                        s = new(new((from list in data.sx select new GFElement<int>(list, field)).ToList(), 
                            field6), new((from list in data.sy select new GFElement<int>(list, field)).ToList(), 
                            field6));
                        break;
                    }
                }
            }

            if (pq is null)
            {
                var factors = PrimeFactors(m);
                q = factors.Max<BigInteger>();

                pq = RandomPointOfPrimeOrder(q, curve, true);
                s = RandomPointOfPrimeOrder(q, curve6, false);
                source.Add(new(field.dimension, positive, q.ToString(), pq.x.p, pq.y.p,
                    (from el in s.x.p select el.p).ToList(),
                    (from el in s.y.p select el.p).ToList()));
                source.Sort((data1, data2) => data1.l - data2.l);

                string json = JsonSerializer.Serialize(source, options);
                File.WriteAllText(fileName, json);
            }

        }
        private (GaloisField<int> field, GaloisField<GFElement<int>> field6) InitFields(int l, bool positive)
        {
            GaloisField<int> field = null;
            GaloisField<GFElement<int>> field6 = null;

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

                        field = new(new PrimeField(3), p);
                        field6 = new(field, p6);

                        // need to define automorphism
                        var tempField = new GaloisField<int>(new PrimeField(3), p6);

                        var y2 = GaloisField<int>.Root(new GFElement<int>(new List<int> { -1 }, tempField));
                        var y3 = Methods.SolveCubicEquation3L(positive ? 1 : 2, tempField);

                        u = new((from e in y2.p select field.Scalar(e)).ToList(), field6);
                        rm = new((from e in y3 select field.Scalar(e)).ToList(), field6);
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
        public static ECP<T> RandomPointOfPrimeOrder<T>(BigInteger order, EllipticCurve<T> c, bool mustBeOfThisOrder)
        {
            ECP<T> p = null;
            while (true)
            {
                GFElement<T> px, py;

                px = c.field.RandomNonZero();
                var root = GaloisField<T>.Root(c.RightSide(px));
                if (root is not null)
                    py = root;
                else
                    continue;

                p = new ECP<T>(px, py);
                var pq = c.Mult(p, order);

                if ((pq is null && mustBeOfThisOrder) || (pq is not null && !mustBeOfThisOrder))
                    break;
            }

            return p;
        }
        class HPQ
        {
            ECP<GFElement<int>> p, q;
            GFElement<GFElement<int>> lambda;

            public HPQ(ECP<GFElement<int>> p, ECP<GFElement<int>> q, GFElement<GFElement<int>> a)
            {
                this.p = p;
                this.q = q;

                if (p.field.IsEqual(p.x, q.x))
                {
                    if (p.field.IsEqual(p.y, q.y))
                    {
                        // tangent
                        if (p.y == new GFElement<GFElement<int>>(new List<GFElement<int>> { p.field.baseField.Scalar(0) }, p.field))
                            lambda = null;
                        else if (p.field.Characteristic() == 3)
                            lambda = a / (p.field.Scalar(2) * p.y);
                        else if (p.field.Characteristic() > 3)
                            lambda = (p.field.Scalar(3) * p.x * p.x + a) / (p.field.Scalar(2) * p.y);
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

            public GFElement<GFElement<int>> Value(ECP<GFElement<int>> a)
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
        static List<GFElement<GFElement<int>>> MillerUpdated(ECP<GFElement<int>> p, List<ECP<GFElement<int>>> r, BigInteger n, EllipticCurve<GFElement<int>> curve)
        {
            var binary = new List<int>();
            while (n != 1)
            {
                binary.Add((int)(n % 2));
                n /= 2;
            }
            binary.Add(1);

            var t = binary.Count - 2;

            var res = new List<GFElement<GFElement<int>>>();
            for (int i = 0; i < r.Count; i++)
                res.Add(new GFElement<GFElement<int>>(new List<GFElement<int>> { curve.field.baseField.Scalar(1) }, curve.field));
            var T = p;

            for (int i = t; i >= 0; i--)
            {
                Console.WriteLine(i);
                for (int j = 0; j < r.Count; j++)
                    res[j] *= res[j];

                var hpq = new HPQ(T, T, curve.a);
                for (int j = 0; j < r.Count; j++)
                    res[j] *= hpq.Value(r[j]);

                T = curve.Action(T, T);

                if (binary[i] == 1)
                {
                    hpq = new HPQ(T, p, curve.a);
                    for (int j = 0; j < r.Count; j++)
                        res[j] *= hpq.Value(r[j]);
                    T = curve.Action(T, p);
                }
            }

            return res;
        }
        static public GFElement<GFElement<int>> WeilPairing(ECP<GFElement<int>> p, ECP<GFElement<int>> q, BigInteger n, EllipticCurve<GFElement<int>> curve, ECP<GFElement<int>> s = null)
        {
            if (s is null)
                s = RandomPointOfPrimeOrder(n, curve, false);

            //s = new(new(new List<GFElement> { curve.field.baseField.Scalar(0) }, curve.field), new(new List<GFElement> { curve.field.baseField.Scalar(36) }, curve.field));

            var sn = new ECP<GFElement<int>>(s.x, -s.y);

            var args1 = new List<ECP<GFElement<int>>> { curve.Action(q, s), s };
            var args2 = new List<ECP<GFElement<int>>> { curve.Action(p, sn), sn };

            var res1 = MillerUpdated(p, args1, n, curve);
            var res2 = MillerUpdated(q, args2, n, curve);

            var num = res1[0] / res1[1];
            var den = res2[0] / res2[1];
            return num / den;
        }
        public ECP<int> MapToGroup(byte[] m, float delta)
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

                var x = new GFElement<int>(poly, curve.field);

                var root = GaloisField<int>.Root(curve.RightSide(x));
                if (root is null)
                    continue;

                if (root.p[0] > -root.p[0])
                    root = -root;
                if (b)
                    root = -root;
                var pm = new ECP<int>(x, root);

                var mult = this.m / this.q; // must know m and q
                pm = curve.Mult(pm, mult);

                if (pm is not null)
                    return pm;
            }

            throw new Exception();
        }
        public ECP<GFElement<int>> Auto6L(ECP<GFElement<int>> a)
        {
            // check if curve is positive for the small curve
            if (curve.b == new GFElement<int>(new List<int> { 1 }, curve.field))
                return new ECP<GFElement<int>>(-a.x + rp, u * a.y); // +++
            else
                return new ECP<GFElement<int>>(-a.x + rm, u * a.y); // ---
        }
    }
}

