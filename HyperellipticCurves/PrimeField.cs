﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace HyperellipticCurves
{
    public class PrimePolynomial
    {
        public List<int> poly { get; }
        public PrimeField field;
        public PrimePolynomial(List<int> poly, PrimeField field)
        {
            this.poly = new(poly);
            int degree = field.Degree(this.poly);
            this.poly.RemoveRange(degree + 1, this.poly.Count - degree - 1);

            this.field = field;
        }

        public static PrimePolynomial operator -(PrimePolynomial a)
        {
            return new PrimePolynomial(a.field.SubtractPoly(new List<int> { 0 }, a.poly), a.field);
        }

        public static PrimePolynomial operator +(PrimePolynomial a, PrimePolynomial b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new PrimePolynomial(a.field.AddPoly(a.poly, b.poly), a.field);
            else
                return null;
        }

        public static PrimePolynomial operator -(PrimePolynomial a, PrimePolynomial b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new PrimePolynomial(a.field.SubtractPoly(a.poly, b.poly), a.field);
            else
                return null;
        }

        public static PrimePolynomial operator *(PrimePolynomial a, PrimePolynomial b)
        {
            if (a.field.characteristic == b.field.characteristic)
                return new PrimePolynomial(a.field.MultiplyPoly(a.poly, b.poly), a.field);
            else
                return null;
        }

        public static PrimePolynomial operator %(PrimePolynomial a, PrimePolynomial b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                List<int> factor;
                return new PrimePolynomial(a.field.RemainderPoly(a.poly, b.poly, out factor), a.field);
            }
            else
                return null;
        }

        public static PrimePolynomial operator /(PrimePolynomial a, PrimePolynomial b)
        {
            if (a.field.characteristic == b.field.characteristic)
            {
                List<int> factor;
                a.field.RemainderPoly(a.poly, b.poly, out factor);
                return new PrimePolynomial(factor, a.field);
            }
            else
                return null;
        }
        public static PrimePolynomial Pow(PrimePolynomial a, int k)
        {
            return new PrimePolynomial(a.field.PowPoly(a.poly, k), a.field);
        }
    }

    public interface IField<T>
    {
        public T Add(T a, T b);
        public T Subtract(T a, T b);
        public T Multiply(T a, T b);
        public T Inverse(T a);
        public T Scalar(int a);
        public bool IsEqual(T a, T b);
        public T Clone(T a);
        public BigInteger Size();
        public T RandomNonZero(int seed = -1);
        public int Characteristic();
        public bool IsEqual(IField<T> b);
    }
    public class PrimeField : IField<int>
    {
        public int characteristic { get; }

        public PrimeField(int characteristic)
        {
            this.characteristic = characteristic;
        }

        public int Add(int a, int b)
        {
            return Methods.NumRemainder(a + b, characteristic);
        }
        public int Subtract(int a, int b)
        {
            return Methods.NumRemainder(a - b, characteristic);
        }
        public int Multiply(int a, int b)
        {
            return Methods.NumRemainder(a * b, characteristic);
        }
        public int Inverse(int a)
        {
            int inv, s;
            Methods.NumericalEuclid(a, characteristic, out inv, out s);
            inv = Methods.NumRemainder(inv, characteristic);
            return inv;
        }

        public List<int> AddPoly(List<int> a, List<int> b)
        {
            List<int> c = new List<int>();

            for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
                c.Add(Add(a[i], b[i]));

            if (a.Count > b.Count)
                for (int i = b.Count; i < a.Count; i++)
                    c.Add(a[i]);
            else
                for (int i = a.Count; i < b.Count; i++)
                    c.Add(b[i]);

            return c;
        }
        public List<int> SubtractPoly(List<int> a, List<int> b)
        {

            List<int> c = new List<int>();

            for (int i = 0; i < Math.Min(a.Count, b.Count); i++)
                c.Add(Subtract(a[i], b[i]));

            if (a.Count > b.Count)
                for (int i = b.Count; i < a.Count; i++)
                    c.Add(a[i]);
            else
                for (int i = a.Count; i < b.Count; i++)
                    c.Add(Subtract(0, b[i]));

            return c;
        }
        public List<int> MultiplyPoly(List<int> a, List<int> b)
        {
            List<int> c = new List<int>(a.Count + b.Count - 1);
            for (int i = 0; i < a.Count + b.Count - 1; i++)
                c.Add(0);
            for (int i = 0; i < a.Count; i++)
                for (int j = 0; j < b.Count; j++)
                    c[i + j] = Add(c[i + j], Multiply(a[i], b[j]));

            return c;
        }
        public List<int> PowPoly(List<int> a, int k)
        {
            var res = new List<int>(a);
            for (int i = 1; i < k; i++)
                res = MultiplyPoly(res, a);

            return res;
        }
        public List<int> RemainderPoly(List<int> a, List<int> b, out List<int> factor)
        {
            List<int> copy = new List<int>(a);
            factor = new List<int>(a);
            for (int i = 0; i < factor.Count; i++)
                factor[i] = 0;

            int cur = Degree(copy);
            int bf = Degree(b);

            int inv = Inverse(b[bf]);

            while (cur >= bf)
            {
                if (copy[cur] != 0)
                {
                    int diff = cur - bf;
                    for (int i = 0; i < bf; i++)
                        copy[diff + i] = Subtract(copy[diff + i], Multiply(Multiply(copy[cur], inv), b[i]));
                    factor[diff] = Multiply(copy[cur], inv);
                    copy[cur] = 0;
                }
                cur--;
            }

            return copy;
        }
        public PrimePolynomial EuclidPoly(PrimePolynomial a, PrimePolynomial b, out PrimePolynomial t, out PrimePolynomial s)
        {
            PrimePolynomial[] sl = { new(new List<int> { 0 }, this), new(new List<int> { 1 }, this) };
            PrimePolynomial[] tl = { new(new List<int> { 1 }, this), new(new List<int> { 0 }, this) };

            var ac = new PrimePolynomial(a.poly, this);
            var bc = new PrimePolynomial(b.poly, this);

            //Program.Print(MultiplyPoly(a.poly, tl[1].poly));
            //Program.Print(MultiplyPoly(b.poly, sl[1].poly));
            //Program.Print(bc.poly);
            //Console.WriteLine();

            while (LeadingCoeff(ac.poly) != 0 && LeadingCoeff(bc.poly) != 0)
            {
                List<int> q;
                List<int> c = ac.poly;
                ac = new(RemainderPoly(bc.poly, ac.poly, out q), this);

                //Program.Print(AddPoly(MultiplyPoly(c, q), ac.poly));
                //Program.Print(bc.poly);
                //Console.WriteLine();

                bc = new(c,this);

                sl = new PrimePolynomial[] { sl[1] - new PrimePolynomial(q, this) * sl[0], sl[0] };
                tl = new PrimePolynomial[] { tl[1] - new PrimePolynomial(q, this) * tl[0], tl[0] };

                //Program.Print(ac.poly);
                //Program.Print(q);
                //Program.Print(sl[0].poly);
                //Program.Print(tl[0].poly);
                //Console.WriteLine();

                //Program.Print(MultiplyPoly(a.poly, tl[0].poly));
                //Program.Print(MultiplyPoly(b.poly, sl[0].poly));
                //Program.Print(ac.poly);
                //Program.Print(bc.poly);
                //Console.WriteLine();
            }

            t = tl[1];
            s = sl[1];
            if (LeadingCoeff(bc.poly) != 0)
                return bc;
            else
                return new PrimePolynomial(new List<int> { 1 }, this);
        }
        public int Degree(List<int> polynomial)
        {
            int li = polynomial.FindLastIndex((int i) => i != 0);
            if (li < 0)
                return 0;
            else
                return li;
        }
        public virtual int LeadingCoeff(List<int> polynomial)
        {
            return polynomial[Degree(polynomial)];
        }

        public int Scalar(int a)
        {
            return Methods.NumRemainder(a, characteristic);
        }

        public bool IsEqual(int a, int b)
        {
            return a == b;
        }

        public int Clone(int a)
        {
            return a;
        }

        public BigInteger Size()
        {
            return characteristic;
        }

        public int RandomNonZero(int seed = -1)
        {
            var rand = new Random();
            return rand.Next() % characteristic;
        }

        public int Characteristic()
        {
            return characteristic;
        }

        public bool IsEqual(IField<int> b)
        {
            return characteristic == b.Characteristic();
        }

        public static bool operator ==(PrimeField a, PrimeField b)
        {
            return a.IsEqual(b);
        }
        public static bool operator !=(PrimeField a, PrimeField b)
        {
            return !(a == b);
        }
    }
}
