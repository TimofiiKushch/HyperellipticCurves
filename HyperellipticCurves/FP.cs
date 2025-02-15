using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace HyperellipticCurves
{
    public class FPInt
    {
        public static int multiplications;
        public static int inversions;
        public BigInteger value { get; }
        public FP field;
        public FPInt(BigInteger value, FP field)
        {
            this.value = Methods.NumRemainder(value, field.characteristic);
            this.field = field;
        }
        public static FPInt Zero(FP field)
        {
            return new FPInt(new(0), field);
        }
        public static FPInt One(FP field)
        {
            return new FPInt(new(1), field);
        }
        public static FPInt operator -(FPInt a)
        {
            return new FPInt(-a.value, a.field);
        }

        public static FPInt operator +(FPInt a, FPInt b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new FPInt(a.value + b.value, a.field);
            else
                return null;
        }

        public static FPInt operator -(FPInt a, FPInt b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new FPInt(a.value - b.value, a.field);
            else
                return null;
        }

        public static FPInt operator *(FPInt a, FPInt b)
        {
            multiplications += 1;
            if (a.field.characteristic == b.field.characteristic)
                return new FPInt(a.value * b.value, a.field);
            else
                return null;
        }
        public static FPInt operator *(BigInteger a, FPInt b)
        {
            multiplications += 1;
            return new FPInt(a * b.value, b.field);
        }
        public static FPInt operator *(int a, FPInt b)
        {
            return new FPInt(a * b.value, b.field);
        }

        public static bool operator ==(FPInt a, FPInt b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                return a.value == b.value;
            }
            else
                return false;
        }
        public static bool operator !=(FPInt a, FPInt b)
        {
            return !(a == b);
        }
        public virtual FPInt Inverse()
        {
            inversions += 1;
            BigInteger inv, s;
            Methods.NumericalEuclid(value, field.characteristic, out inv, out s);
            return new(inv, field);
        }

    }
    public class Point
    {
        public FPInt x { get; }
        public FPInt y { get; }

        public Point(FPInt x, FPInt y)
        {
            this.x = x;
            this.y = y;
        }

        public override string ToString() 
        {
            return x.value + " " + y.value;
        }
    }

    public class PrimeDivisorFP
    {
        public PrimePolynomialFP u;
        public PrimePolynomialFP v;

        public PrimeDivisorFP(PrimePolynomialFP u, PrimePolynomialFP v)
        {
            this.u = u;
            this.v = v;
        }
        public PrimeDivisorFP(ProjectiveDivisor divisor)
        {
            if (!divisor.u2IsZero)
                this.u = new(new List<FPInt> { divisor.u0, divisor.u1, divisor.z }, divisor.u0.field);
            else
                this.u = new(new List<FPInt> { divisor.u0, divisor.u1 }, divisor.u0.field);

            this.v = new(new List<FPInt> { divisor.v0, divisor.v1 }, divisor.v0.field);

            var zinv = divisor.z.Inverse();
            this.u = zinv * this.u;
            this.v = zinv * this.v;
        }
        public PrimeDivisorFP(NewDivisor divisor)
        {
            if (!divisor.u2IsZero)
                this.u = new(new List<FPInt> { divisor.u0, divisor.u1, divisor.z1 }, divisor.u0.field);
            else
                this.u = new(new List<FPInt> { divisor.u0, divisor.u1 }, divisor.u0.field);

            this.v = new(new List<FPInt> { divisor.v0, divisor.v1 }, divisor.v0.field);

            this.u = divisor.z1.Inverse() * this.u;
            this.v = (divisor.z1 * divisor.Z1 * divisor.Z2).Inverse() * this.v;
        }
        public void Print()
        {
            u.Print();
            v.Print();
        }
    }
    public class ProjectiveDivisor
    {
        public FPInt u1 { get; }
        public FPInt u0 { get; }
        public FPInt v1 { get; }
        public FPInt v0 { get; }
        public FPInt z { get; }
        public bool u2IsZero { get; }
        public ProjectiveDivisor(FPInt u1, FPInt u0, FPInt v1, FPInt v0, FPInt z, bool u2IsZero)
        {
            this.u1 = u1;
            this.u0 = u0;
            this.v1 = v1;
            this.v0 = v0;
            this.z = z;
            this.u2IsZero = u2IsZero;
        }
        public ProjectiveDivisor(PrimeDivisorFP divisor)
        {
            this.u1 = divisor.u[1];
            this.u0 = divisor.u[0];
            this.v1 = divisor.v[1];
            this.v0 = divisor.v[0];
            this.z = FPInt.One(u1.field);
            this.u2IsZero = divisor.u[2] == FPInt.Zero(u1.field);
        }
    }
    public class NewDivisor
    {
        public FPInt u1 { get; }
        public FPInt u0 { get; }
        public FPInt v1 { get; }
        public FPInt v0 { get; }
        public FPInt Z1 { get; }
        public FPInt Z2 { get; }
        public FPInt z1 { get; }
        public FPInt z2 { get; }
        public bool u2IsZero { get; }
        public NewDivisor(FPInt u1, FPInt u0, FPInt v1, FPInt v0, FPInt Z1, FPInt Z2, FPInt z1, FPInt z2, bool u2IsZero)
        {
            this.u1 = u1;
            this.u0 = u0;
            this.v1 = v1;
            this.v0 = v0;
            this.Z1 = Z1;
            this.Z2 = Z2;
            this.z1 = z1;
            this.z2 = z2;
            this.u2IsZero = u2IsZero;
        }
        public NewDivisor(PrimeDivisorFP divisor)
        {
            this.u1 = divisor.u[1];
            this.u0 = divisor.u[0];
            this.v1 = divisor.v[1];
            this.v0 = divisor.v[0];
            this.Z1 = FPInt.One(u1.field);
            this.Z2 = FPInt.One(u1.field);
            this.z1 = FPInt.One(u1.field);
            this.z2 = FPInt.One(u1.field);
            this.u2IsZero = divisor.u[2] == FPInt.Zero(u1.field);
        }
    }

    public class PrimePolynomialFP
    {
        private List<FPInt> poly;
        public FP field;
        public PrimePolynomialFP(List<FPInt> poly, FP field)
        {
            this.poly = new(poly);

            int degree = field.Degree(this.poly);
            if (this.poly.Count - degree - 1 >= 0)
                this.poly.RemoveRange(degree + 1, this.poly.Count - degree - 1);

            this.field = field;
        }
        public PrimePolynomialFP(List<BigInteger> poly, FP field)
        {
            this.poly = new(poly.Count);
            for (int i = 0; i < poly.Count; i++)
                this.poly.Add(new(poly[i], field));

            int degree = field.Degree(this.poly);
            this.poly.RemoveRange(degree + 1, this.poly.Count - degree - 1);

            this.field = field;
        }
        public FPInt this[int index]
        {
            get {
                if (index >= poly.Count)
                    return FPInt.Zero(field);
                else
                    return poly[index]; 
            }
        }
        public PrimePolynomialFP Copy()
        {
            return new(new List<FPInt>(poly), field);
        }
        public int Degree()
        {
            return field.Degree(poly);
        }
        public FPInt LeadingCoeff()
        {
            return field.LeadingCoeff(poly);
        }
        public void Print()
        {
            for (int i = 0; i < poly.Count; i++)
                Console.Write(poly[i].value + " ");
            Console.WriteLine();
        }
        public static PrimePolynomialFP RemainderPoly(PrimePolynomialFP a, PrimePolynomialFP b, out PrimePolynomialFP factor)
        {
            List<FPInt> copy = new List<FPInt>(a.poly);
            var factorPoly = new List<FPInt>(a.poly.Count);
            for (int i = 0; i < a.poly.Count; i++)
                factorPoly.Add(FPInt.Zero(a.field));

            int cur = a.Degree();
            int bf = b.Degree();

            FPInt inv = b[bf].Inverse();

            while (cur >= bf)
            {
                if (copy[cur].value != 0)
                {
                    int diff = cur - bf;
                    for (int i = 0; i < bf; i++)
                        copy[diff + i] = copy[diff + i] - copy[cur] * inv * b[i];
                    factorPoly[diff] = copy[cur] * inv;
                    copy[cur] = new(new(0), a.field);
                }
                cur--;
            }
            factor = new(factorPoly, a.field);

            return new(copy, a.field);
        }
        public PrimePolynomialFP Monic()
        {
            var scalar = LeadingCoeff();
            var inv = scalar.Inverse();
            int degree = Degree();

            List<FPInt> newPoly = new();
            for (int i = 0; i < degree; i++)
                newPoly.Add(poly[i] * inv);
            newPoly.Add(new(new(1), field));

            return new(newPoly, field);
        }
        public PrimePolynomialFP Derivative()
        {
            List<FPInt> derivative = new();
            for (int i = 1; i < poly.Count; i++)
                derivative.Add(poly[i] * new FPInt(new(i), field));
            return new(derivative, field);
        }
        public FPInt Value(FPInt point)
        {
            FPInt sum = new(new(0), field);
            FPInt exponent = new(new(1), field);
            for (int i = 0; i < poly.Count; i++)
            {
                sum += exponent * poly[i];
                exponent *= point;
            }
            return sum;
        }
        public static bool ResultantIsZero(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            var ap = a.poly;
            var bp = b.poly;

            int an = a.Degree();
            int bn = b.Degree();
            int n = an + bn;
            List<List<FPInt>> m = new(n);

            // 1: build the matrix
            for (int i = 0; i < n; i++)
            {
                if (i < an)
                {
                    var temp = new List<FPInt>(n);
                    for (int j = 0; j < n; j++)
                    {
                        if (j < i || j > i + an)
                            temp.Add(new(new(0), a.field));
                        else
                            temp.Add(ap[j - i]);
                    }
                    m.Add(temp);
                }
                else
                {
                    var temp = new List<FPInt>(n);
                    for (int j = 0; j < n; j++)
                    {
                        if (j < i - an || j > i - an + bn)
                            temp.Add(new(new(0), a.field));
                        else
                            temp.Add(bp[j - i + an]);
                    }
                    m.Add(temp);
                }
            }

            // 2: gauss
            for (int i = 0; i < n; i++)
            {
                // find non zero
                if (m[i][i].value == 0)
                {
                    int nonZeroI = -1;
                    for (int j = i + 1; j < n; j++)
                        if (m[j][i].value != 0)
                        {
                            nonZeroI = j;
                            break;
                        }
                    if (nonZeroI != -1)
                    {
                        var temp = m[i];
                        m[i] = m[nonZeroI];
                        m[nonZeroI] = temp;
                    }
                    else
                        return true;
                }

                for (int j = i + 1; j < n; j++)
                    m[j] = (m[i][i] * (new PrimePolynomialFP(m[j], a.field)) - m[j][i] * (new PrimePolynomialFP(m[i], a.field))).poly;
            }

            return false;
        }
        public static FPInt ResultanSmall(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            // deg a, b <= 2

            if (a[2] == FPInt.Zero(a.field) && b[2] == FPInt.Zero(b.field))
            {
                return a[0] * b[1] - a[1] * b[0];
            }
            else if (a[2] != FPInt.Zero(a.field) && b[2] != FPInt.Zero(b.field))
            {
                FPInt v0 = a[0];
                FPInt v1 = a[1];
                FPInt v2 = a[2];

                FPInt u0 = b[0];
                FPInt u1 = b[1];
                FPInt u2 = b[2];

                var w1 = v0 * u2;
                var w2 = v2 * u0;
                var w3 = v1 * u1;
                var w4 = w2 - w1;
                var res2 = w4 * w4 - w3 * (w1 + w2) + v0 * v2 * u1 * u1 + v1 * v1 * u0 * u2;

                return res2;
            }
            else
            {
                if (a[2] != FPInt.Zero(a.field))
                {
                    FPInt v0 = a[0];
                    FPInt v1 = a[1];
                    FPInt v2 = a[2];

                    FPInt u0 = b[0];
                    FPInt u1 = b[1];

                    return v0 * u1 * u1 - v1 * u0 * u1 + v2 * u0 * u0;
                }
                else
                {
                    FPInt v0 = a[0];
                    FPInt v1 = a[1];

                    FPInt u0 = b[0];
                    FPInt u1 = b[1];
                    FPInt u2 = b[2];

                    return v0 * v0 * u2 - v0 * v1 * u1 + v1 * v1 * u0;
                }
            }
        }

        public static PrimePolynomialFP operator -(PrimePolynomialFP a)
        {
            return new PrimePolynomialFP(a.field.SubtractPoly(new List<FPInt> { new(new(0), a.field) }, a.poly), a.field);
        }

        public static PrimePolynomialFP operator +(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new PrimePolynomialFP(a.field.AddPoly(a.poly, b.poly), a.field);
            else
                return null;
        }

        public static PrimePolynomialFP operator -(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new PrimePolynomialFP(a.field.SubtractPoly(a.poly, b.poly), a.field);
            else
                return null;
        }

        public static PrimePolynomialFP operator *(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new PrimePolynomialFP(a.field.MultiplyPoly(a.poly, b.poly), a.field);
            else
                return null;
        }
        public static PrimePolynomialFP operator *(FPInt a, PrimePolynomialFP b)
        {
            return new PrimePolynomialFP(b.field.MultiplyPoly(new() { a }, b.poly), b.field);
        }
        public static PrimePolynomialFP operator *(PrimePolynomialFP b, FPInt a)
        {
            return a * b;
        }

        public static PrimePolynomialFP operator %(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                PrimePolynomialFP factor;
                return RemainderPoly(a, b, out factor);
            }
            else
                return null;
        }

        public static PrimePolynomialFP operator /(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                PrimePolynomialFP factor;
                RemainderPoly(a, b, out factor);
                return factor;
            }
            else
                return null;
        }
        public static bool operator ==(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                if (a.poly.Count != b.poly.Count)
                    return false;
                else
                {
                    for (int i = 0; i < a.poly.Count; i++)
                        if (a[i] != b[i])
                            return false;
                    return true;
                }
            }
            else
                return false;
        }
        public static bool operator !=(PrimePolynomialFP a, PrimePolynomialFP b)
        {
            return !(a == b);
        }

    }

    public class FP
    {
        public static int notCommonProj = 0;
        public enum AdditionMethod
        {
            cantor,
            cantorOptimized,
            cantorProjective,
            cantorNew
        }
        public BigInteger characteristic { get; }
        public AdditionMethod method { get; set; }

        public FP(BigInteger characteristic)
        {
            this.characteristic = characteristic;
        }
        public FPInt n(BigInteger value)
        {
            return new FPInt(value, this);
        }
        public Point p(BigInteger x, BigInteger y)
        {
            return new Point(n(x), n(y));
        }
        public Point RandomPoint(PrimePolynomialFP f, Random rand)
        {
            // h assumed 0
            FPInt x0, a, legendre;
            while (true)
            {
                x0 = new FPInt(rand.Next(), this);
                a = f.Value(x0);
                legendre = Power(a, (characteristic - 1) / 2);
                if (legendre.value == 1)
                    break; 
            }

            // Tonelli and Shanks

            // 1
            FPInt e = FPInt.Zero(this);
            var r = characteristic - 1;
            while (r % 2 == 0)
            {
                e += FPInt.One(this);
                r /= 2;
            }

            // 2
            FPInt n;
            while (true)
            {
                n = new(rand.Next(), this);
                legendre = Power(n, (characteristic - 1) / 2);
                if (legendre.value == characteristic - 1)
                    break;
            }

            // 3
            var z = Power(n, r);
            var y = z;
            var s = e;
            var x = Power(a, (r - 1) / 2);

            // 4
            var b = a * x * x;
            x = a * x;

            // 5
            while (b.value != 1)
            {
                // 6
                var m = new FPInt(1, this);
                var b2m = b * b;

                // 7
                while (b2m.value != 1)
                {
                    m = m + new FPInt(1, this);
                    b2m = b2m * b2m;
                }

                // 8
                var sm1 = Power(new FPInt(2, this), (s - m - new FPInt(1, this)).value);
                var t = Power(y, sm1.value);
                y = t * t;
                s = m;

                // 9
                x = t * x;
                b = y * b;
            }

            return new(x0, x);
        }


        public List<FPInt> AddPoly(List<FPInt> a, List<FPInt> b)
        {
            List<FPInt> c = new List<FPInt>();

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

        public List<FPInt> SubtractPoly(List<FPInt> a, List<FPInt> b)
        {

            List<FPInt> c = new List<FPInt>();

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

        public List<FPInt> MultiplyPoly(List<FPInt> a, List<FPInt> b)
        {
            List<FPInt> c = new List<FPInt>(a.Count + b.Count - 1);
            for (int i = 0; i < a.Count + b.Count - 1; i++)
                c.Add(new(new(0), this));
            for (int i = 0; i < a.Count; i++)
                for (int j = 0; j < b.Count; j++)
                    c[i + j] = c[i + j] + a[i] * b[j];

            return c;
        }
        
        public PrimePolynomialFP EuclidPoly(PrimePolynomialFP a, PrimePolynomialFP b, out PrimePolynomialFP t, out PrimePolynomialFP s)
        {
            if (a.LeadingCoeff().value == 0)
            {
                t = new(new List<FPInt> { new(new(0), this) }, this);
                s = new(new List<FPInt> { new(new(1), this) }, this);
                return b.Copy();
            }
            else if (b.LeadingCoeff().value == 0)
            {
                t = new(new List<FPInt> { new(new(1), this) }, this);
                s = new(new List<FPInt> { new(new(0), this) }, this);
                return a.Copy();
            }

            PrimePolynomialFP[] sl = { new(new List<FPInt> { new(new(0), this) }, this), new(new List<FPInt> { new(new(1), this) }, this) };
            PrimePolynomialFP[] tl = { new(new List<FPInt> { new(new(1), this) }, this), new(new List<FPInt> { new(new(0), this) }, this) };

            var ac = a.Copy();
            var bc = b.Copy();

            while (ac.LeadingCoeff().value != 0 && bc.LeadingCoeff().value != 0)
            {
                PrimePolynomialFP q;
                PrimePolynomialFP c = ac;
                ac = PrimePolynomialFP.RemainderPoly(bc, ac, out q);
                bc = c;

                sl = new PrimePolynomialFP[] { sl[1] - q * sl[0], sl[0] };
                tl = new PrimePolynomialFP[] { tl[1] - q * tl[0], tl[0] };

                //Program.Print(ac.poly);
                //Program.Print(q);
                //Program.Print(sl[0].poly);
                //Program.Print(tl[0].poly);
                //Console.WriteLine();

                //Program.Print(MultiplyPoly(a, tl[1]));
                //Program.Print(MultiplyPoly(b, sl[1]));
                //Program.Print(ac);
                //Console.WriteLine();
            }

            t = tl[1];
            s = sl[1];
            if (bc.LeadingCoeff().value != 0)
                return bc;
            else
                return new PrimePolynomialFP(new List<FPInt> { new(new(1), this) }, this);
        }
        public int Degree(List<FPInt> polynomial)
        {
            int li = polynomial.FindLastIndex((FPInt i) => i.value != 0);
            if (li < 0)
                return 0;
            else
                return li;
        }

        public virtual FPInt LeadingCoeff(List<FPInt> polynomial)
        {
            return polynomial[Degree(polynomial)];
        }

        public PrimeDivisorFP ScalarMult(PrimeDivisorFP a, BigInteger k, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // h assumed 0 fo projective and new coordinates

            if (method == AdditionMethod.cantor) 
            {
                PrimeDivisorFP res = null;
                PrimeDivisorFP last = a;

                while (true)
                {
                    if (k % 2 == 1)
                    {
                        if (res != null)
                            res = Cantor(res, last, h, f);
                        else
                            res = last;
                    }

                    k /= 2;
                    if (k == 0)
                        break;
                    last = Cantor(last, last, h, f);

                    //if (res != null)
                    //    res.Print();
                    //last.Print();
                }
                return res;
            }
            else if (method == AdditionMethod.cantorOptimized)
            {
                PrimeDivisorFP res = null;
                PrimeDivisorFP last = a;

                while (true)
                {
                    if (k % 2 == 1)
                    {
                        if (res != null)
                            res = CantorOptimized(res, last, h, f);
                        else
                            res = last;
                    }

                    k /= 2;
                    if (k == 0)
                        break;
                    last = CantorOptimized(last, last, h, f);
                }
                return res;
            }
            else if (method == AdditionMethod.cantorProjective)
            {
                ProjectiveDivisor res = null;
                ProjectiveDivisor last = new(a);

                while (true)
                {
                    if (k % 2 == 1)
                    {
                        if (res != null)
                            res = CantorProjective(res, last, h, f);
                        else
                            res = last;
                    }

                    k /= 2;
                    if (k == 0)
                        break;
                    last = CantorProjective(last, last, h, f);
                }
                return new(res);
            }
            else if (method == AdditionMethod.cantorNew)
            {
                NewDivisor res = null;
                NewDivisor last = new(a);

                while (true)
                {
                    if (k % 2 == 1)
                    {
                        if (res != null)
                            res = CantorNew(res, last, h, f);
                        else
                            res = last;
                    }

                    k /= 2;
                    if (k == 0)
                        break;

                    //new PrimeDivisorFP(last).Print();
                    last = CantorNew(last, last, h, f);

                    //if (res != null)
                    //    new PrimeDivisorFP(res).Print();
                }
                return new(res);
            }

            return null;
        }
        public FPInt Power(FPInt a, BigInteger k)
        {
            FPInt res = null;
            FPInt last = a;

            while (true)
            {
                if (k % 2 == 1)
                {
                    if (!(res is null))
                        res = res*last;
                    else
                        res = last;
                }

                k /= 2;
                if (k == 0)
                    break;
                last = last*last;
            }

            if (res is null)
                return FPInt.One(this);

            return res;
        }
        public PrimeDivisorFP Cantor(PrimeDivisorFP a, PrimeDivisorFP b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // 0
            var u1 = a.u;
            var u2 = b.u;
            var v1 = a.v;
            var v2 = b.v;

            int g = (f.Degree() - 1) / 2;

            // 1
            PrimePolynomialFP e1, e2;
            var d1 = EuclidPoly(u1, u2, out e1, out e2);

            //Console.WriteLine("euclid");
            //(u1 * e1 + u2 * e2).Print();
            //d1.Print();
            //(u1 % d1).Print();
            //(u2 % d1).Print();
            //Console.WriteLine();
            //u1.Print();
            //u2.Print();
            //Console.WriteLine();
            //e1.Print();
            //e2.Print();
            //Console.WriteLine("euclid");

            // 2
            PrimePolynomialFP c1, c2;
            var d = EuclidPoly(d1, v1 + v2 + h, out c1, out c2);

            //Console.WriteLine("euclid");
            //(d1 * c1 + (v1 + v2 + h) * c2).Print();
            //d.Print();
            //(d1 % d).Print();
            //((v1 + v2 + h) % d).Print();
            //Console.WriteLine();
            //d1.Print();
            //(v1 + v2 + h).Print();
            //Console.WriteLine();
            //c1.Print();
            //c2.Print();
            //Console.WriteLine("euclid");

            // 3
            var s1 = c1 * e1;
            var s2 = c1 * e2;
            var s3 = c2;

            // 4
            var u = (u1 / d) * (u2 / d);

            var m1 = s1 * u1 * v2;
            var m2 = s2 * u2 * v1;
            var m3 = (v1 * v2 + f) * s3;
            var v = ((m1 + m2 + m3) / d) % u;

            //Program.Print(s1.poly);
            //Program.Print(s2.poly);
            //Program.Print(m1.poly);
            //Program.Print(m2.poly);

            // 5
            while (u.Degree() > g)
            {
                //u.Print();
                //v.Print();

                // 6
                var num = f - v * h - v * v;
                var un = num / u;

                num = -h - v;
                var vn = num % un;

                // 7
                u = un;
                v = vn;

                //Console.WriteLine();

            } // 8

            //u.Print();
            //v.Print();
            //Console.WriteLine();

            // 9
            u = u.Monic();

            // 10
            return new PrimeDivisorFP(u, v);
        }
        public PrimeDivisorFP CantorOptimized(PrimeDivisorFP a, PrimeDivisorFP b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // divisors must have mumford representation, g must be 2

            // -1
            if (a.u.Degree() > b.u.Degree())
                return CantorOptimized(b, a, h, f);

            // 0
            var u1 = a.u;
            var u2 = b.u;
            var v1 = a.v;
            var v2 = b.v;

            int g = (f.Degree() - 1) / 2;
            if (g != 2)
                throw new Exception("Sorry, not implemented for g != 2");

            if (u1.Degree() == 0)
            {
                return new PrimeDivisorFP(u2.Copy(), v2.Copy());
            }
            if (u1.Degree() == 1 && u2.Degree() == 1)
            {
                // case u1 == u2
                if (u1[0] == u2[0])
                {
                    // case v1 != v2
                    if (v1[0] == -v2[0] - h.Value(-u1[0]))
                    {
                        return new PrimeDivisorFP(new(new List<BigInteger>() { 1 }, this), new(new List<BigInteger>() { 0 }, this));
                    }
                    // case v1 == v2
                    else
                    {
                        return COD1(a, h, f);
                    }
                }
                // case u1 != u2
                {
                    var u = u1 * u2;

                    var inv = (u1[0] - u2[0]).Inverse();
                    var v = new PrimePolynomialFP(new List<FPInt>()
                    { inv * (v2[0] * u1[0] - v1[0] * u2[0]), inv*(v2[0] - v1[0]) }, this);

                    return new(u, v);
                }
            }
            else if (u1.Degree() == 1 && u2.Degree() == 2)
            {
                var u2v = u2.Value(-u1[0]);

                if (u2v.value != 0)
                {
                    return CO2B(a, b, h, f);
                }
                else
                {
                    if (u2[1] == u1[0] && u2[0] == u1[0] * u1[0])
                    {
                        var d2 = COD2(b, h, f);
                        var d1inv = InverseDeg1(a, h, f);
                        return CO2B(d2, d1inv, h, f);
                    }
                    else if (-v2.Value(-u1[0]) == v1[0] + h.Value(-u1[0])) // it seems there should be a minus sign
                    {
                        return new PrimeDivisorFP(new(new List<FPInt>() { u2[1] - u1[0], FPInt.One(this) }, this),
                            new(new List<FPInt>() { v2.Value(u1[0] - u2[1]) }, this));
                    }
                    else
                    {
                        var a2 = COD1(a, h, f);
                        var add = new PrimeDivisorFP(new(new List<FPInt>() { u2[1] - u1[0], FPInt.One(this) }, this),
                            new(new List<FPInt>() { v2.Value(u1[0] - u2[1]) }, this));
                        return CO2B(a2, add, h, f);
                    }
                }
            }
            else if (u1.Degree() == 2 && u2.Degree() == 2)
            {
                if (u1 == u2)
                {
                    if (v1 == (-v2 - h) % u1)
                    {
                        return new PrimeDivisorFP(new(new List<BigInteger>() { 1 }, this), new(new List<BigInteger>() { 0 }, this));
                    }
                    else if (v1 == v2)
                    {
                        return COD2(a, h, f);
                    }
                    else
                    {
                        var inv = (v2[1] - v1[1]).Inverse();
                        var temp = (v1[0] - v2[0]) * inv;

                        var toDouble = new PrimeDivisorFP(new(new List<FPInt>() { -temp, FPInt.One(this) }, this), new(new List<FPInt>() { v1.Value(temp) }, this));
                        return COD1(toDouble, h, f);
                    }
                }
                else
                {
                    if (PrimePolynomialFP.ResultanSmall(u1, u2).value == 0)
                    {
                        PrimePolynomialFP t, s;
                        var gcd = EuclidPoly(u1, u2, out t, out s).Monic();
                        var x1 = -gcd[0];
                        var y1 = v1.Value(x1);

                        // divisors that correspond to p1, p2, p3
                        var d1 = new PrimeDivisorFP(gcd, new(new List<FPInt>() { y1 }, this));
                        var d2 = new PrimeDivisorFP(new(new List<FPInt>() { u1[1] + x1, FPInt.One(this) }, this),
                            new(new List<FPInt>() { v1.Value(-(u1[1] + x1)) }, this));
                        var d3 = new PrimeDivisorFP(new(new List<FPInt>() { u2[1] + x1, FPInt.One(this) }, this),
                            new(new List<FPInt>() { v2.Value(-(u2[1] + x1)) }, this));

                        //Console.WriteLine(x1.value + " " + y1.value);
                        //Console.WriteLine((-d2.u[0]).value + " " + d2.v[0].value);
                        //Console.WriteLine((-d3.u[0]).value + " " + d3.v[0].value);

                        if (y1 == v2.Value(x1))
                        {
                            // double the divisor that corresponds to p1 - pinf (gcd is x - x1)
                            var co1 = CantorOptimized(d1, d1, h, f);

                            // need to compute d1 + d2 + d3
                            // since degree d2.u, d3.u == 1 the problem is reduced to a previously solved case

                            return CantorOptimized(CantorOptimized(co1, d2, h, f), d3, h, f);

                        }
                        else
                        {
                            // need to compute d2 + d3
                            // since degree d2.u, d3.u == 1 the problem is reduced to a previously solved case
                            return CantorOptimized(d2, d3, h, f);
                        }
                    }
                    else
                    {
                        // most common addition

                        return CO2A(a, b, h, f);
                    }
                }
            }

            return null;
        }
        public ProjectiveDivisor CantorProjective(ProjectiveDivisor a, ProjectiveDivisor b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // assumption: h == 0

            int g = (f.Degree() - 1) / 2;
            if (g != 2)
                throw new Exception("Sorry, not implemented for g != 2");

            // condition 1: must be of degree 2
            if (!a.u2IsZero && !b.u2IsZero)
            {
                var u1 = new PrimePolynomialFP(new List<FPInt> { a.u0, a.u1, a.z }, this) * b.z;
                var u2 = new PrimePolynomialFP(new List<FPInt> { b.u0, b.u1, b.z }, this) * a.z;

                // check if equal and not opposite
                if (u1 == u2)
                {
                    var v1 = new PrimePolynomialFP(new List<FPInt> { a.v0, a.v1}, this) * b.z;
                    var v2 = new PrimePolynomialFP(new List<FPInt> { b.v0, b.v1}, this) * a.z;

                    if (v1 != -v2 && v1 == v2)
                    {
                        // check if most common doubling
                        if (PrimePolynomialFP.ResultanSmall(new FPInt(new(2), this) * v1, u1).value != 0)
                        {
                            return ProjectiveDoubling(a, h, f);
                        }
                    }
                }
                else
                {
                    // check if most common addition
                    if (PrimePolynomialFP.ResultanSmall(u1, u2).value != 0)
                    {
                        return ProjectiveAddition(a, b, h, f);
                    }
                }
            }

            // not one of 2 most common cases => convert into affine
            notCommonProj += 1;
            return new(CantorOptimized(new(a), new(b), h, f));
        }
        public NewDivisor CantorNew(NewDivisor a, NewDivisor b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // assumption: h == 0

            int g = (f.Degree() - 1) / 2;
            if (g != 2)
                throw new Exception("Sorry, not implemented for g != 2");

            // condition 1: must be of degree 2
            if (!a.u2IsZero && !b.u2IsZero)
            {
                var u1 = new PrimePolynomialFP(new List<FPInt> { a.u0, a.u1, a.z1 }, this) * b.z1;
                var u2 = new PrimePolynomialFP(new List<FPInt> { b.u0, b.u1, b.z1 }, this) * a.z1;

                // check if equal and not opposite
                if (u1 == u2)
                {
                    var v1 = new PrimePolynomialFP(new List<FPInt> { a.v0, a.v1 }, this) * b.z1 * b.Z1 * b.Z2;
                    var v2 = new PrimePolynomialFP(new List<FPInt> { b.v0, b.v1 }, this) * a.z1 * a.Z1 * a.Z2;

                    if (v1 != -v2 && v1 == v2)
                    {
                        // check if most common doubling
                        if (PrimePolynomialFP.ResultanSmall(new FPInt(new(2), this) * v1, u1).value != 0)
                        {
                            return NewDoubling(a, h, f);
                        }
                    }
                }
                else
                {
                    // check if most common addition
                    if (PrimePolynomialFP.ResultanSmall(u1, u2).value != 0)
                    {
                        return NewAddition(a, b, h, f);
                    }
                }
            }

            // not one of 2 most common cases => convert into affine
            return new(CantorOptimized(new(a), new(b), h, f));
        }
        private PrimeDivisorFP InverseDeg1(PrimeDivisorFP a, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            var u1 = a.u;
            var v1 = a.v;

            var copy = u1.Copy();

            return new(copy, new(new List<FPInt>{ -v1[0] - h.Value(-u1[0]) }, this));
        }

        // double a, deg a == 1
        private PrimeDivisorFP COD1(PrimeDivisorFP a, PrimePolynomialFP h, PrimePolynomialFP f) 
        {
            //Console.WriteLine("Case d1");

            var u1 = a.u;
            var v1 = a.v;

            var u = u1 * u1;

            var fd = f.Derivative();
            var hd = h.Derivative();
            var fdu = fd.Value(-u1[0]);
            var hdu = hd.Value(-u1[0]);
            var hu = h.Value(-u1[0]);

            var inv = (2 * v1[0] + hu).Inverse();
            var temp = (fdu - v1[0] * hdu) * inv;
            var v = new PrimePolynomialFP(new List<FPInt>() { v1[0] + u1[0] * temp, temp }, this);

            return new(u, v);
        }
        // double a, deg a == 2
        private PrimeDivisorFP COD2(PrimeDivisorFP a, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            //Console.WriteLine("Case d2");

            var u1 = a.u;
            var v1 = a.v;

            if (PrimePolynomialFP.ResultanSmall(h + new FPInt(new(2), this) * v1, u1).value == 0)
            {
                PrimePolynomialFP t, s;
                var gcd = EuclidPoly(h + new FPInt(new(2), this) * v1, u1, out t, out s).Monic();
                var temp = u1[1] - gcd[0];
                var toDouble = new PrimeDivisorFP(new(new List<FPInt>() { temp, FPInt.One(this) }, this), new(new List<FPInt>() { v1.Value(-temp) }, this));
                return COD1(toDouble, h, f);
            }

            else
            {
                return CO2C(a, h, f);
            }
        }
        private PrimeDivisorFP CO2A(PrimeDivisorFP a, PrimeDivisorFP b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // g == 2, deg u1 == deg u2 == 2
            //Console.WriteLine("Case 2a");

            // 0
            var u1 = a.u;
            var u2 = b.u;
            var v1 = a.v;
            var v2 = b.v;

            // test
            //var t = (f - v2 * h - v2 * v2) / u2; // exact
            //var sf = ((v1 - v2) / u2) % u1;
            //var s = (v1 - v2 / u2 ) % u1;
            //var sr = (v1 - v2) % u2; // implication: sr is the numerator in .../u2
            //s.Print();
            //var l = s * u2;
            //l.Print();

            //var t1 = (l + h + new FPInt(2, this) * v2) % u2;
            //var t2 = l * (l + h + new FPInt(2, this) * v2) % u2;
            //t1.Print();
            //t2.Print();

            //var u = (t - s * (l + h + new FPInt(2, this) * v2)) / u1; // exact
            //var ud = u.Monic();
            //var vd = (-h - (l + v2)) % ud;
            //Console.WriteLine("simplified");
            //ud.Print();
            //vd.Print();

            // 1
            var z1 = u1[1] - u2[1];
            var z2 = u2[0] - u1[0];
            var z3 = u1[1] * z1 + z2;
            var r = z2 * z3 + z1 * z1 * u1[0];

            // 2
            var inv1 = z1;
            var inv0 = z3;

            // test
            //var invt = new PrimePolynomialFP(new List<FPInt>{ inv0, inv1}, this);
            //(((v1 - v2) * invt) % u1).Print();

            // 3
            var w0 = v1[0] - v2[0];
            var w1 = v1[1] - v2[1];
            var w2 = inv0 * w0;
            var w3 = inv1 * w1;

            var sd1 = (inv0 + inv1) * (w0 + w1) - w2 - w3 * (FPInt.One(this) + u1[1]);
            var sd0 = w2 - u1[0] * w3;

            if (sd1.value == 0)
            {
                // 4
                var inv = r.Inverse();
                var s0 = sd0 * inv; // s0 matches up

                // 5
                var ud0 = f[4] - u2[1] - u1[1] - s0 * s0 - s0 * h[2];

                // 6
                w1 = s0 * (u2[1] - ud0) + h[1] + v2[1] - h[2] * ud0;
                w2 = u2[0] * s0 + v2[0] + h[0];
                var vd0 = ud0 * w1 - w2;

                //PrimePolynomialFP testu = new(new List<FPInt> { ud0, FPInt.One(this) }, this);
                //PrimePolynomialFP testv = new(new List<FPInt> { new(6, this) }, this);
                //((testv * testv + testv * h - f) % testu).Print();

                // 8
                return new(new(new List<FPInt> { ud0, FPInt.One(this) }, this), new(new List<FPInt> { vd0 }, this));
            }
            else
            {
                // 4
                w1 = (r * sd1).Inverse();
                w2 = r * w1;
                w3 = sd1 * sd1 * w1;
                var w4 = r * w2;
                var w5 = w4 * w4;
                var sdd0 = sd0 * w2;

                // 5
                var ld2 = u2[1] + sdd0;
                var ld1 = u2[1] * sdd0 + u2[0];
                var ld0 = u2[0] * sdd0;

                // 6
                var ud0 = (sdd0 - u1[1]) * (sdd0 - z1 + h[2] * w4) - u1[0];
                ud0 = ud0 + ld1 + (h[1] + 2 * v2[1]) * w4 + (2 * u2[1] + z1 - f[4]) * w5;
                var ud1 = 2 * sdd0 - z1 + h[2] * w4 - w5;

                // 7
                w1 = ld2 - ud1;
                w2 = ud1 * w1 + ud0 - ld1;
                var vd1 = w2 * w3 - v2[1] - h[1] + h[2] * ud1;
                w2 = ud0 * w1 - ld0;
                var vd0 = w2 * w3 - v2[0] - h[0] + h[2] * ud0;

                // 8
                return new(new(new List<FPInt> { ud0, ud1, FPInt.One(this) }, this), new(new List<FPInt> { vd0, vd1 }, this));
            }
        }
        private PrimeDivisorFP CO2B(PrimeDivisorFP a, PrimeDivisorFP b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // g == 2, deg u1 == 1, deg u2 == 2
            //Console.WriteLine("Case 2b");

            // -1
            if (a.u.Degree() > b.u.Degree())
                return CO2B(b, a, h, f);

            // 0
            var u1 = a.u;
            var u2 = b.u;
            var v1 = a.v;
            var v2 = b.v;

            // 1
            var r = u2[0] - (u2[1] - u1[0]) * u1[0];

            // 2
            var inv = r.Inverse();

            // 3
            var s0 = inv * (v1[0] - v2[0] + v2[1] * u1[0]);

            // 4
            var l1 = s0 * u2[1];
            var l0 = s0 * u2[0];

            // 5
            var t2 = f[4] - u2[1];
            var t1 = f[3] - (f[4] - u2[1]) * u2[1] - v2[1] * h[2] - u2[0];

            // 6
            var ud1 = t2 - s0 * s0 - s0 * h[2] - u1[0];
            var ud0 = t1 - s0 * (l1 + h[1] + 2 * v2[1]) - u1[0] * ud1;

            // 7
            var vd1 = (h[2] + s0) * ud1 - (h[1] + l1 + v2[1]);
            var vd0 = (h[2] + s0) * ud0 - (h[0] + l0 + v2[0]);

            // 8
            return new(new(new List<FPInt> { ud0, ud1, FPInt.One(this) }, this), new(new List<FPInt> { vd0, vd1 }, this));
        }
        private PrimeDivisorFP CO2C(PrimeDivisorFP a, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // g == 2, deg u1 == 2
            //Console.WriteLine("Case 2c");

            // 0 
            var u1 = a.u[1];
            var u0 = a.u[0];
            var v1 = a.v[1];
            var v0 = a.v[0];

            // 1
            var vr1 = h[1] + 2 * v1 - h[2] * u1;
            var vr0 = h[0] + 2 * v0 - h[2] * u0;

            // 2
            var w0 = v1 * v1;
            var w1 = u1 * u1;
            var w2 = vr1 * vr1;
            var w3 = u1 * vr1;
            var r = u0 * w2 + vr0 * (vr0 - w3);

            // 3
            var invd1 = -vr1;
            var invd0 = vr0 - w3;

            // 4
            w3 = f[3] + w1;
            var w4 = 2 * u0;
            var td1 = 2 * (w1 - f[4] * u1) + w3 - w4 - h[2] * v1;
            var td0 = u1 * (2 * w4 - w3 + f[4] * u1 + h[2] * v1) + f[2] - w0 - 2 * f[4] * u0 - h[1] * v1 - h[2] * v0;

            // 5
            w0 = td0 * invd0;
            w1 = td1 * invd1;
            var sd1 = (invd0 + invd1) * (td0 + td1) - w0 - w1 * (FPInt.One(this) + u1);
            var sd0 =  w0 - u0 * w1;

            if (sd1.value == 0)
            {
                // 6
                w1 = r.Inverse();
                var s0 = sd0 * w1;
                w2 = u0 * s0 + v0 + h[0];

                // 7
                var ud0 = f[4] - s0 * s0 - s0 * h[2] - 2 * u1;

                // 8
                w1 = s0 * (u1 - ud0) - h[2] * ud0 + v1 + h[1];
                var vd0 = ud0 * w1 - w2;

                // 10
                return new(new(new List<FPInt> { ud0, FPInt.One(this) }, this), new(new List<FPInt> { vd0 }, this));
            }
            else
            {
                // 6
                w1 = (r * sd1).Inverse();
                w2 = r * w1;
                w3 = sd1 * sd1 * w1;
                w4 = r * w2;
                var w5 = w4 * w4;
                var sdd0 = sd0 * w2;

                // 7
                var ld2 = u1 + sdd0;
                var ld1 = u1 * sdd0 + u0;
                var ld0 = u0 * sdd0;

                // 8
                var ud0 = sdd0 * sdd0 + w4 * (h[2] * (sdd0 - u1) + 2 * v1 + h[1]) + w5 * (2 * u1 - f[4]);
                var ud1 = 2 * sdd0 + h[2] * w4 - w5;

                // 9
                w1 = ld2 - ud1;
                w2 = ud1 * w1 + ud0 - ld1;
                var vd1 = w2 * w3 - v1 - h[1] + h[2] * ud1;
                w2 = ud0 * w1 - ld0;
                var vd0 = w2 * w3 - v0 - h[0] + h[2] * ud0;

                // 10
                return new(new(new List<FPInt> { ud0, ud1, FPInt.One(this) }, this), new(new List<FPInt> { vd0, vd1 }, this));
            }
        }

        public ProjectiveDivisor ProjectiveAddition(ProjectiveDivisor a, ProjectiveDivisor b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // most common case

            // 0
            var d1 = a;
            var d2 = b;

            // 1
            var z = d1.z * d2.z;
            var ur21 = d1.z * d2.u1;
            var ur20 = d1.z * d2.u0;
            var vr21 = d1.z * d2.v1;
            var vr20 = d1.z * d2.v0;

            // 2
            var z1 = d1.u1 * d2.z - ur21;
            var z2 = ur20 - d1.u0 * d2.z;
            var z3 = d1.u1 * z1 + z2 * d1.z;
            var r = z2 * z3 + z1 * z1 * d1.u0;

            // 3
            var inv1 = z1;
            var inv0 = z3;

            // 4
            var w0 = d1.v0 * d2.z - vr20;
            var w1 = d1.v1 * d2.z - vr21;
            var w2 = inv0 * w0;
            var w3 = inv1 * w1;
            var s1 = (inv0 + d1.z * inv1) * (w0 + w1) - w2 - w3 * (d1.z + d1.u1);
            var s0 = w2 - d1.u0 * w3;

            if (s1.value == 0)
            {
                // s0 / (Z1^3*Z2^2) == sd0 at this point
                //var ud0 = (f[4] - d2.u1 - d1.u1) * r * r - s0 * s0 - s0 * h[2] * r;

                var r2 = r * r;
                var r3 = r * r2;
                var rz = r * z;

                var ud0 = (f[4] * z - d2.u1 * d1.z - d1.u1 * d2.z) * r2 - s0 * s0 * z - s0 * h[2] * rz; // r^2*z
                w1 = s0 * (d2.u1 * r2 * d1.z - ud0) + r3 * z * h[1] + d2.v1 * r3 * d1.z - r * h[2] * ud0; // r^3*z
                w2 = d2.u0 * d1.z * s0 + d2.v0 * d1.z * r + rz * h[0]; // r*z
                var vd0 = ud0 * w1 - w2 * r3 * rz; // r^5*z^2
                var zf = rz * rz * r3;

                return new(zf, ud0 * rz * r2, FPInt.Zero(this), vd0, zf, true);
            }
            else
            {
                // 5
                var R = z * r;
                s0 = s0 * z;
                var s3 = s1 * z;
                var Rr = R * s3;
                var S3 = s3 * s3;
                var S = s0 * s1;
                var Sr = s3 * s1;
                var Sa = s0 * s3;
                var Ra = Rr * Sr;

                // 6
                var l2 = Sr * ur21;
                var l0 = S * ur20;
                var l1 = (Sr + S) * (ur21 + ur20) - l2 - l0;
                l2 = l2 + Sa;

                // 7
                var ud0 = s0 * s0 + s1 * z1 * (s1 * (z1 + ur21) - 2 * s0) + z2 * Sr + R * (2 * s1 * vr21 + r * (z1 + 2 * ur21 - f[4] * z));
                var ud1 = 2 * Sa - Sr * z1 - R * R;

                // 8
                l2 = l2 - ud1;
                w0 = ud0 * l2 - S3 * l0;
                w1 = ud1 * l2 + S3 * (ud0 - l1);

                // 9
                var zd = Rr * S3;
                ud1 = Rr * ud1;
                ud0 = Rr * ud0;

                // 10
                var vd0 = w0 - Ra * vr20;
                var vd1 = w1 - Ra * vr21;

                // 11
                return new(ud1, ud0, vd1, vd0, zd, false);
            }
        }
        public ProjectiveDivisor ProjectiveDoubling(ProjectiveDivisor a, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // 0
            var u1 = a.u1;
            var u0 = a.u0;
            var v1 = a.v1;
            var v0 = a.v0;
            var z = a.z;

            // 1
            var z2 = z * z;
            var vr1 = 2 * v1;
            var vr0 = 2 * v0;
            var w0 = v1 * v1;
            var w1 = u1 * u1;
            var w2 = 4 * w0;
            var w3 = vr0 * z - u1 * vr1;
            var r = vr0 * w3 + w2 * u0;

            // 2
            var inv1 = -vr1;
            var inv0 = w3;

            // 3
            w3 = f[3] * z2 + w1;
            var w4 = 2 * u0;
            var t1 = 2 * w1 + w3 - z * (w4 + 2 * f[4] * u1);
            FPInt t0;
            if (f[4].value == 0)
                t0 = u1 * (2 * z * w4 - w3) + z * (f[2] * z2 - w0);
            else
                t0 = u1 * (z * (2 * w4 + f[4] * u1) - w3) + z * (z * (f[2] * z - 2 * f[4] * u0) - w0);

            // 4
            w0 = t0 * inv0;
            w1 = t1 * inv1;
            var s3 = (inv0 + inv1) * (t0 + t1) - w0 - (FPInt.One(this) + u1) * w1;
            var s1 = s3 * z;
            var s0 = w0 - z * u0 * w1;

            if (s1.value == 0)
            {
                // s0 / r == z^2 at this point
                var rz2 = r * z2;
                var rz3 = rz2 * z;
                var r2z3 = rz3 * r;
                var r3z5 = r2z3 * rz2;

                w2 = u0 * s0 + v0 * rz2 + h[0] * rz3;
                var ud0 = f[4] * r2z3 * z - s0 * s0 - s0 * h[2] * rz2 - 2 * u1 * r2z3;
                w1 = s0 * (u1 * r2z3 - ud0) - h[2] * ud0 * rz2 + v1 * r3z5 + h[1] * r3z5 * z;
                var vd0 = ud0 * w1 - w2 * r3z5 * rz2;

                var zf = r3z5 * rz3 * rz2;

                return new(zf, ud0 * r3z5 * z, FPInt.Zero(this), vd0, zf, true);
            }
            else
            {
                // 5
                var R = z2 * r;
                var Rr = R * s1;
                var S1 = s1 * s1;
                var S0 = s0 * s0;
                s1 = s1 * s3;
                s0 = s0 * s3;
                var S = s0 * z;
                var Ra = Rr * s1;

                // 6
                var l2 = u1 * s1;
                var l0 = u0 * s0;
                var l1 = (s1 + s0) * (u1 + u0) - l2 - l0;

                // 7
                var ud0 = S0 + R * (2 * s3 * v1 + z * r * (2 * u1 - f[4] * z));
                var ud1 = 2 * S - R * R;

                // 8
                l2 = l2 + S - ud1;
                w0 = ud0 * l2 - S1 * l0;
                w1 = ud1 * l2 + S1 * (ud0 - l1);

                // 9
                var zd = S1 * Rr;
                ud1 = Rr * ud1;
                ud0 = Rr * ud0;

                // 10
                var vd0 = w0 - Ra * v0;
                var vd1 = w1 - Ra * v1;

                // 11
                return new(ud1, ud0, vd1, vd0, zd, false);
            }
        }
        public NewDivisor NewAddition(NewDivisor a, NewDivisor b, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // f4 == 0
            // 0
            var d1 = a;
            var d2 = b;

            // 1
            var z13 = d1.Z1 * d1.Z2;
            var z23 = d2.Z1 * d2.Z2;
            var z14 = d1.z1 * z13;
            var z24 = d2.z1 * z23;
            var ur21 = d2.u1 * d1.z1;
            var ur20 = d2.u0 * d1.z1;
            var vr21 = d2.v1 * z14;
            var vr20 = d2.v0 * z14;

            // 2
            var y1 = d1.u1 * d2.z1 - ur21;
            var y2 = ur20 - d1.u0 * d2.z1;
            var y3 = d1.u1 * y1 + y2 * d1.z1;
            var r = y2 * y3 + y1 * y1 * d1.u0;
            var Zd2 = d1.Z1 * d2.Z1;
            var Zr2 = d1.Z2 * d2.Z2;
            var Zd1 = Zd2 * Zd2;
            Zr2 = Zr2 * Zd1 * r;
            Zd2 = Zd2 * Zr2;
            Zr2 = Zr2 * Zr2;
            var zd2 = Zd2 * Zd2;

            // 3
            var inv1 = y1;
            var inv0 = y3;

            // 4
            var w0 = d1.v0 * z24 - vr20;
            var w1 = d1.v1 * z24 - vr21;
            var w2 = inv0 * w0;
            var w3 = inv1 * w1;
            var s1 = (inv0 + d1.z1 * inv1) * (w0 + w1) - w2 - w3 * (d1.z1 + d1.u1);
            var s0 = w2 - d1.u0 * w3;

            if (s1.value == 0)
            {
                // s0 / (r*z13*z23) == sd0 at this point

                var s0r = r*z13*z23;
                var s0r2 = s0r * s0r;
                var s0r3 = s0r2 * s0r;
                var s0r4 = s0r2 * s0r2;
                //var s0rz = s0r * d2.z1 * z24;

                var ud0 = (f[4] - d2.u1 * d2.z1 - d1.u1 * d1.z1) * s0r2 - s0 * s0 - s0 * h[2] * s0r; // s0r2

                w1 = s0 * (d2.u1 * d2.z1 * s0r2 - ud0) + (h[1] + d2.v1 * z24) * s0r3 - h[2] * ud0 * s0r2; // s0r3
                w2 = d2.u0 *d2.z1 * s0 + (d2.v0 * z24 + h[0]) * s0r; // s0r
                var vd0 = ud0 * w1 - w2 * s0r4; // s0r5

                return new(s0r2, ud0, FPInt.Zero(this), vd0, s0r, s0r2, s0r2, s0r4, true);
            }
            else
            {
                // 5
                var S1 = s1 * s1;
                var S0 = s0 * Zd1;
                Zd1 = s1 * Zd1;
                var S = Zd1 * S0;
                S0 = S0 * S0;
                var R = r * Zd1;
                s0 = s0 * Zd1;
                s1 = s1 * Zd1;
                var zd1 = Zd1 * Zd1;

                // 6
                var l2 = s1 * ur21;
                var l0 = s0 * ur20;
                var l1 = (s0 + s1) * (ur20 + ur21) - l0 - l2;
                l2 = l2 + S;

                // 7
                var vd1 = R * vr21;
                var ud0 = S0 + y1 * (S1 * (y1 + ur21) - 2 * s0) + y2 * s1 + 2 * vd1 + (2 * ur21 + y1) * Zr2;
                var ud1 = 2 * S - y1 * s1 - zd2;

                // 8
                l2 = l2 - ud1;
                w0 = l2 * ud0;
                w1 = l2 * ud1;

                // 9
                vd1 = w1 - zd1 * (l1 + vd1 - ud0);
                var vd0 = w0 - zd1 * (l0 + R * vr20);

                // 10
                return new NewDivisor(ud1, ud0, vd1, vd0, Zd1, Zd2, zd1, zd2, false);
            }
        }
        public NewDivisor NewDoubling(NewDivisor a, PrimePolynomialFP h, PrimePolynomialFP f)
        {
            // f4 == 0

            // 0
            var u1 = a.u1;
            var u0 = a.u0;
            var v1 = a.v1;
            var v0 = a.v0;
            var Z1 = a.Z1;
            var Z2 = a.Z2;
            var z1 = a.z1;
            var z2 = a.z2;

            // 1
            var ur0 = u0 * z1;
            var w0 = v1 * v1;
            var w1 = u1 * u1;
            var w3 = v0 * z1 - u1 * v1;
            var r = w0 * u0 + v0 * w3;
            var Zr2 = z2 * r * z1;
            var Zd2 = 2 * Zr2 * Z1;
            Zr2 = Zr2 * Zr2;

            // 2
            var inv1 = -v1;
            var inv0 = w3;

            // 3
            var z3 = z1 * z1;
            w3 = f[3] * z3 + w1;
            var t1 = z2 * (2 * (w1 - ur0) + w3);
            z3 = z3 * z1;
            var t0 = z2 * (u1 * (4 * ur0 - w3) + z3 * f[2]) - w0;

            // 4
            w0 = t0 * inv0;
            w1 = t1 * inv1;
            var s1 = (inv0 + inv1) * (t0 + t1) - w0 - w1 * (FPInt.One(this) + u1);
            var s0 = w0 - w1 * ur0;

            if (s1.value == 0)
            {
                // s0 / r == Z1^2 here
                var rz1 = r * z1;
                var z21 = z1 * Z2;
                var z31 = z21 * Z1;
                var rz21 = r * z21;
                var rz31 = r * z31;
                var r2z42 = rz21 * rz21;
                var r2z62 = rz31 * rz31;
                var cube = r2z62 * rz31;
                var rZ1Z2 = r * Z1 * Z2;
                var zf = cube * r2z62;

                var w2 = u0 * s0 + v0 * rz1 + h[0] * rz1 * z31;
                var ud0 = f[4] * r2z62 - s0 * s0 - s0 * h[2] * rz31 - 2 * u1 * r2z42;
                w1 = s0 * (u1 * r2z42 - ud0) - h[2] * ud0 * rz31 + v1 * r2z62 * r + h[1] * cube;
                var vd0 = ud0 * w1 - w2 * rZ1Z2 * cube;

                return new(r2z62, ud0, FPInt.Zero(this), vd0, rz31, r2z62, r2z62, r2z62 * r2z62, true);
            }
            else
            {
                // 5
                var S0 = s0 * s0;
                var Zd1 = s1 * z1;
                var zd1 = Zd1 * Zd1;
                var S = s0 * Zd1;
                var R = r * Zd1;
                var zd2 = Zd2 * Zd2;
                s0 = s0 * s1;
                s1 = Zd1 * s1;

                // 6
                var l2 = s1 * u1;
                var l0 = s0 * u0;
                var l1 = (s0 + s1) * (u0 + u1) - l0 - l2;
                l2 = l2 + S;

                // 7
                var vd1 = R * v1;
                var ud0 = S0 + 4 * (vd1 + 2 * Zr2 * u1);
                var ud1 = 2 * S - zd2;

                // 8
                l2 = l2 - ud1;
                w0 = l2 * ud0;
                w1 = l2 * ud1;

                // 9
                vd1 = w1 - zd1 * (l1 + 2 * vd1 - ud0);
                var vd0 = w0 - zd1 * (l0 + 2 * R * v0);

                // 10
                return new NewDivisor(ud1, ud0, vd1, vd0, Zd1, Zd2, zd1, zd2, false);
            }
        }
    }
}
