using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HyperellipticCurves
{
    class Program
    {
        // primes around 2^d, d is from 3 to 50
        static string[] primes = new string[] {
            "11",
            "17",
            "37",
            "71",
            "137",
            "263",
            "521",
            "1031",
            "2063",
            "4111",
            "8209",
            "16411",
            "32779",
            "65543",
            "131101",
            "262151",
            "524309",
            "1048583",
            "2097169",
            "4194319",
            "8388617",
            "16777259",
            "33554467",
            "67108879",
            "134217757",
            "268435463",
            "536870923",
            "1073741831",
            "2147483659",
            "4294967311",
            "8589934609",
            "17179869209",
            "34359738421",
            "68719476767",
            "137438953481",
            "274877906951",
            "549755813911",
            "1099511627791",
            "2199023255579",
            "4398046511119",
            "8796093022237",
            "17592186044423",
            "35184372088891",
            "70368744177679",
            "140737488355337",
            "281474976710677",
            "562949953421381",
            "1125899906842679",
            "2251799813685269",
            "4503599627370517"};

        public class Info
        {
            int multiplications;
            int inversions;
            double time;
            public Info(int multiplications, int inversions, double time)
            {
                this.multiplications = multiplications;
                this.inversions = inversions;
                this.time = time;
            }
            public override string ToString()
            {
                return multiplications + " " + inversions + " " + time;
            }
        }

        static void Main(string[] args)
        {
            //int t, s;
            //Methods.NumericalEuclid(123, 456, out t, out s);
            //Console.Read();

            //var rand = new Random();
            //var l = new List<int>(Enumerable.Repeat(0, 80));
            //l[0] = 1;
            //l[26] = 2;
            //l[79] = 1;
            //var field = new GaloisField(3, l);


            //var field = new GaloisField(3, new List<int> { 1, 2, 0, 0, 0, 1 });
            //var extension = new GaloisFieldExtension(field, new List<int> { 2, 1, 0, 0, 0, 0, 1 });

            //var field = new GaloisField(3, new List<int> { 2, 1, 0, 0, 0, 0, 1 });
            //var field = new GaloisField(3, new List<int> { 2, 1, 1 });
            //var field = new GaloisField(3, new List<int> { 1, 2, 0, 0, 0, 1 });
            //var field = new GaloisField(3, new List<int> { 2, 1});
            //var extension = new GaloisFieldExtension(field, new List<int> { 2, 1, 1 });

            //while (true)
            //{
            //    var re = extension.RandomNonZero();
            //    var root = GaloisFieldExtension.Root(re);
            //    if (root is not null)
            //        if (re != root[0] * root[0])
            //            Print(re);
            //}
            //var y = extension.Root(new(new List<GFElement> { field.Scalar(-1) }, extension));
            //var y = GaloisFieldExtension.SolveCubicEquation3L(1, extension);
            //Print(y);
            //Print(y * y * y - y);

            //var a = new GFElement(new List<int> { 2 }, field);
            //var b = new GFElement(new List<int> { 1 }, field);
            //EllipticCurve c = new EllipticCurve(a, b);

            //GFElement px;
            //GFElement py;

            //for (int i = 0; ; i++)
            //{
            //    px = c.field.RandomNonZero(i);
            //    var roots = c.field.Root(c.RightSide(px), 2);
            //    if (roots != null)
            //    {
            //        py = roots[0];
            //        break;
            //    }
            //    else
            //        continue;
            //}

            //px = new(new List<int> { 2, 0, 0 }, field);
            //py = new(new List<int> { 1, 0, 2 }, field);
            //var p = new ECP(px, py);

            //var pq = c.Mult(p, order);

            //var test1 = c.Mult(p, 5);
            //var p2 = c.Action(p, p);
            //var p4 = c.Action(p2, p2);
            //var test1 = c.Action(p, p4);

            //var test2 = p;
            //for (int i = 0; i < 4; i++)
            //{
            //    Console.WriteLine(i + 1);
            //    Print(test2.x);
            //    Print(test2.y);
            //    test2 = c.Action(p, test2);
            //}

            //Console.WriteLine(5);
            //Print(test2.x);
            //Print(test2.y);
            //Console.WriteLine();

            //Console.WriteLine(2);
            //Print(p2.x);
            //Print(p2.y);
            //Console.WriteLine(4);
            //Print(p4.x);
            //Print(p4.y);
            //Console.WriteLine(5);
            //Print(test1.x);
            //Print(test1.y);

            //var m = new byte[16];
            //curve.MapToGroup(m);

            //TestWeilPairing();

            // TODO: debug map to group
            // TOD0: field equality comparison
            SignatureScheme(79, false);
            //while (true)
            //    SignatureScheme(7, false);

            //Print(Methods.PrimeFactors(BigInteger.Parse("49269609804781974450852068861184694669")));
            //Print(Methods.PrimeFactors(BigInteger.Parse("4926960980478197445861184694669")));
            //Print(Methods.PrimeFactors(BigInteger.Parse("23")));


            //var before = DateTime.Now;
            //Print(FindIrreducible(3, 79));
            //Console.WriteLine((DateTime.Now - before).TotalSeconds);
        }

        public static void SignatureScheme(int l, bool positive)
        {
            var rand = new Random();

            // Step 1: keygen
            var a = new List<int> { 2 };
            var b = new List<int> { positive ? 1 : 2 };

            var cm = new EllipticCurveManager(a, b, l, positive);
            var test1 = cm.curve6.field.RandomNonZero();
            var test2 = cm.curve6.field.Inverse(test1);
            Console.ReadLine();

            // let's find P of order q
            // q is prime so just need qP = 0
            // TODO: change so we find primitive and raise it to the power of m/q
            var q = cm.q ?? throw new Exception("q wasn't initialized");
            var p = cm.RandomPointOfPrimeOrder(q);

            // "random" 0 < x < q
            BigInteger x = 0;
            while (x == 0)
            {
                var randBytes = new byte[q.GetByteCount() + 10];
                rand.NextBytes(randBytes);
                x = new BigInteger(randBytes);
                x = x < 0 ? -x : x;
                x %= q;
            }

            var r = cm.curve.Mult(p, x);

            // Step 2: signing
            var message = "Hello world!";
            var encodedMessage = Encoding.ASCII.GetBytes(message);
            var pm = cm.MapToGroup(encodedMessage, 1e-6f); // probability of a failure
            var sm = cm.curve.Mult(pm, x);

            // Step 3: verification
            // need a point with matching x-coordinate
            var roots2 = GaloisField.Root(cm.curve.RightSide(sm.x));
            if (roots2 == null)
                throw new Exception("No roots found: signature is invalid");

            var y = roots2[0];
            var s = new ECP(sm.x, y);
            //Print(s.y);
            //Print(sm.y);

            // convert everything to work in field6
            var s6 = new ECPExtension(new(new List<GFElement> { s.x }, cm.curve6.field), new(new List<GFElement> { s.y }, cm.curve6.field));
            var p6 = new ECPExtension(new(new List<GFElement> { p.x }, cm.curve6.field), new(new List<GFElement> { p.y }, cm.curve6.field));
            var r6 = new ECPExtension(new(new List<GFElement> { r.x }, cm.curve6.field), new(new List<GFElement> { r.y }, cm.curve6.field));
            var pm6 = new ECPExtension(new(new List<GFElement> { pm.x }, cm.curve6.field), new(new List<GFElement> { pm.y }, cm.curve6.field));

            //Console.WriteLine(s6 == pm6);
            //Console.WriteLine(r6 == p6);

            //Print(p.x);
            //Print(p.y);
            //Print(p6.x);
            //Print(p6.y);
            //Print(EllipticCurveManager.WeilPairing(p6, r6, q, cm.curve6));

            var u = EllipticCurveManager.WeilPairing(p6, cm.Auto6L(s6), q, cm.curve6);
            var v = EllipticCurveManager.WeilPairing(r6, cm.Auto6L(pm6), q, cm.curve6);

            if (u == v || u == cm.curve6.field.Inverse(v))
                Console.WriteLine("Congratulations!");
            else
                throw new Exception("Signature rejected");

        }
        
        public static void TestWeilPairing()
        {
            var field = new GaloisField(631, new List<int> { -2, 1 });
            var poly = new List<int>(Enumerable.Repeat(0, 7));
            poly[0] = 2;
            poly[1] = 1;
            poly[6] = 1;
            var field6 = new GaloisFieldExtension(field, poly);

            var a = new GFElement(new List<int> { 30 }, field);
            var b = new GFElement(new List<int> { 34 }, field);

            var a6 = new ExtensionElement(new List<GFElement> { a }, field6);
            var b6 = new ExtensionElement(new List<GFElement> { b }, field6);

            var curve = new EllipticCurveExtension(a6, b6);
            var p = new ECP(new(new List<int> { 36 }, field), new(new List<int> { 60 }, field));
            var q = new ECP(new(new List<int> { 121 }, field), new(new List<int> { 387 }, field));

            var p6 = new ECPExtension(new(new List<GFElement> { p.x }, field6), new(new List<GFElement> { p.y }, field6));
            var q6 = new ECPExtension(new(new List<GFElement> { q.x }, field6), new(new List<GFElement> { q.y }, field6));

            var res = EllipticCurveManager.WeilPairing(p6, q6, 5, curve);
            Print(res);
        }
        public static void fun1()
        {
            FP field = new FP(11);
            PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            PrimePolynomialFP f = new(new List<BigInteger> { 2, 1, 7, 3, 0, 1 }, field);

            List<Point> points = new()
            {
                field.p(1, 5),
                field.p(1, 6),
                field.p(2, 0),
                field.p(4, 5),
                field.p(4, 6),
                field.p(6, 4),
                field.p(6, 7),
                field.p(7, 4),
                field.p(7, 7),
                field.p(9, 4),
                field.p(9, 7),
                field.p(10, 2),
                field.p(10, 9),
            };

            var dt = new PrimeDivisorFP(new(new List<BigInteger> { 9, 3, 1 }, field), new(new List<BigInteger> { 3, 7 }, field));
            foreach (var p in points)
            {
                if (dt.u.Value(p.x).value == 0 && dt.v.Value(p.x).value == p.y.value)
                {
                    Console.WriteLine(p);
                }
            }
        }

        public static void fun2()
        {
            //1 5
            //1 6
            //2 0
            //4 5
            //4 6
            //6 4
            //6 7
            //7 4
            //7 7
            //9 4
            //9 7
            //10 2
            //10 9

            FP field = new FP(11);
            PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            PrimePolynomialFP f = new(new List<BigInteger> { 2, 1, 7, 3, 0, 1 }, field);
            var one = FPInt.One(field);
            var two = new FPInt(2, field);
            var three = new FPInt(3, field);

            // h == 0
            //for (int i = 0; i < field.characteristic; i++)
            //    for (int j = 0; j < field.characteristic; j++)
            //        if (new FPInt(j, field) * new FPInt(j, field) == f.Value(new FPInt(i, field)))
            //            Console.WriteLine(i + " " + j);

            //var d2 = new PrimeDivisorFP(new(new List<BigInteger> { 7, 1 }, field), new(new List<BigInteger> { 5 }, field));
            //var d1 = new PrimeDivisorFP(new(new List<BigInteger> { 10, 7, 1 }, field), new(new List<BigInteger> { 9, 1 }, field));
            var d2 = new PrimeDivisorFP(new(new List<BigInteger> { 10, 0, 1 }, field), new(new List<BigInteger> { 9, 7 }, field));
            //var d2 = new PrimeDivisorFP(new(new List<BigInteger> { 8, 7, 1 }, field), new(new List<BigInteger> { 2 }, field));
            //var d2 = new PrimeDivisorFP(new(new List<BigInteger> { 2, 8, 1 }, field), new(new List<BigInteger> { 10, 6 }, field));
            var dc = field.Cantor(d2, d2, h, f);
            Console.WriteLine("correct");
            dc.u.Print();
            dc.v.Print();

            var dco = field.CantorOptimized(d2, d2, h, f);
            Console.WriteLine("optimized");
            dco.u.Print();
            dco.v.Print();

            //var P1 = new ProjectiveDivisor(d1.u[1], d1.u[0], d1.v[1], d1.v[0], one);
            //var P2 = new ProjectiveDivisor(d2.u[1] * two, d2.u[0] * two, d2.v[1] * two, d2.v[0] * two, two);
            //var P2 = new ProjectiveDivisor(d2.u[1] * three, d2.u[0] * three, d2.v[1] * three, d2.v[0] * three, three);
            //var P2 = new ProjectiveDivisor(d2.u[1], d2.u[0], d2.v[1], d2.v[0], one, true);
            //var dcp = field.ProjectiveAddition(P1, P2, h, f);
            //var dcp = field.ProjectiveDoubling(P2, h, f);

            //Console.WriteLine();
            //var inv = dcp.z.Inverse();
            //Console.WriteLine((dcp.u0 * inv).value + " " + (dcp.u1 * inv).value);
            //Console.WriteLine((dcp.v0 * inv).value + " " + (dcp.v1 * inv).value);

            //var P1n = new NewDivisor(d1.u[1], d1.u[0], d1.v[1], d1.v[0], one, one, one, one);
            ////var P2 = new NewDivisor(d2.u[1], d2.u[0], d2.v[1], d2.v[0], one, one, one, one);
            ////var dcn = field.NewAddition(P1, P2, h, f);
            //var dcn = field.NewDoubling(P1n, h, f);

            //Console.WriteLine();
            //var inv1 = dcn.z1.Inverse();
            //var inv2 = inv1 * (dcn.Z1 * dcn.Z2).Inverse();
            //Console.WriteLine((dcn.u0 * inv1).value + " " + (dcn.u1 * inv1).value);
            //Console.WriteLine((dcn.v0 * inv2).value + " " + (dcn.v1 * inv2).value);
        }

        public static void fun3()
        {
            FP field = new FP(11);
            PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            PrimePolynomialFP f = new(new List<BigInteger> { 2, 1, 7, 3, 0, 1 }, field);

            List<Point> points = new()
            {
                field.p(1, 5),
                field.p(1, 6),
                field.p(2, 0),
                field.p(4, 5),
                field.p(4, 6),
                field.p(6, 4),
                field.p(6, 7),
                field.p(7, 4),
                field.p(7, 7),
                field.p(9, 4),
                field.p(9, 7),
                field.p(10, 2),
                field.p(10, 9),
            };

            var rand = new Random();

            int iters = 10000;
            for (int i = 0; i < iters; i++)
            {
                PrimeDivisorFP[] div = new PrimeDivisorFP[]{ null, null };
                Point[] usedPoints = new Point[] { null, null, null, null };

                for (int k = 0; k < 2; k++)
                {
                    int j1 = rand.Next(points.Count), j2 = j1;
                    while (points[j1].x == points[j2].x)
                        j2 = rand.Next(points.Count);

                    PrimeDivisorFP newDiv = null;

                    if (rand.Next(10) > 2)
                    {
                        // 2 points
                        var p1 = points[j1];
                        var p2 = points[j2];
                        usedPoints[k * 2] = p1;
                        usedPoints[k * 2 + 1] = p2;

                        var u = new PrimePolynomialFP(new List<BigInteger> { (p1.x * p2.x).value, (-p1.x - p2.x).value, 1 }, field);

                        var v0 = (p2.x * p1.y - p1.x * p2.y) * (p2.x - p1.x).Inverse();
                        var v1 = (p2.y - p1.y) * (p2.x - p1.x).Inverse();
                        var v = new PrimePolynomialFP(new List<FPInt> { v0, v1}, field);

                        newDiv = new PrimeDivisorFP(u, v);
                    }
                    else
                    {
                        // 1 point
                        var p = points[j1];
                        usedPoints[k * 2] = p;

                        var u = new PrimePolynomialFP(new List<BigInteger> { -p.x.value, 1 }, field);
                        var v = new PrimePolynomialFP(new List<BigInteger> { p.y.value }, field);

                        newDiv = new PrimeDivisorFP(u, v);
                    }

                    div[k] = newDiv;
                }

                var dc = field.Cantor(div[0], div[1], h, f);
                var dco = field.CantorOptimized(div[0], div[1], h, f);

                if (dc.u != dco.u || dc.v != dco.v)
                {
                    Console.WriteLine("error");

                    Console.WriteLine("results");
                    dc.u.Print();
                    dc.v.Print();
                    dco.u.Print();
                    dco.v.Print();
                    Console.WriteLine("divisors");
                    div[0].u.Print();
                    div[0].v.Print();
                    div[1].u.Print();
                    div[1].v.Print();
                    Console.WriteLine("points");
                    Console.WriteLine(usedPoints[0] + " | " + usedPoints[1]);
                    Console.WriteLine(usedPoints[2] + " | " + usedPoints[3]);
                    Console.Read();
                }
            }
        }
        public static void fun4()
        {
            FP field = new FP(11);
            PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            PrimePolynomialFP f = new(new List<BigInteger> { 2, 1, 7, 3, 0, 1 }, field);

            //var d1 = new PrimeDivisorFP(new(new List<BigInteger> { 8, 4, 1 }, field), new(new List<BigInteger> { 0, 5 }, field));
            //var d2 = new NewDivisor(new(4, field), new(8, field), new(10, field), new(0, field), new(1, field), new(2, field), new(1, field), new(4, field), false);

            var d1 = new PrimeDivisorFP(new(new List<BigInteger> { 3, 10, 1 }, field), new(new List<BigInteger> { 1, 1 }, field));
            var d2 = new NewDivisor(new(10, field), new(3, field), new(2, field), new(2, field), new(1, field), new(2, field), new(1, field), new(4, field), false);

            var dc = field.Cantor(d1, d1, h, f);
            Console.WriteLine("correct");
            dc.Print();

            var dco = field.NewDoubling(d2, h, f);
            Console.WriteLine("new");
            new PrimeDivisorFP(dco).Print();
        }
        public static void fun5()
        {
            //FP field = new FP(11);
            //PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            //PrimePolynomialFP f = new(new List<BigInteger> { 2, 1, 7, 3, 0, 1 }, field);

            //var d2 = new PrimeDivisorFP(new(new List<BigInteger> { 1, 3, 1 }, field), new(new List<BigInteger> { 2, 10 }, field));

            //var t = DateTime.Now;
            //field.method = FP.AdditionMethod.cantor;
            //var r = field.ScalarMult(d2, 10000, h, f);
            ////r.Print();
            //Console.WriteLine((DateTime.Now - t).TotalSeconds);

            var rand = new Random();

            FP field = new FP(BigInteger.Parse("129899216730745422747379980509"));
            PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            PrimePolynomialFP f = new(new List<BigInteger> {
                BigInteger.Parse("84554758087342364918400819975"),
                BigInteger.Parse("43842353021798749401773327333"),
                BigInteger.Parse("81016668528840831602117585943"),
                BigInteger.Parse("52690279977948565928475983553"),
                BigInteger.Parse("110854065858994061391078560211"),
                1 }, field);

            var p1 = field.RandomPoint(f, rand);
            var p2 = field.RandomPoint(f, rand);

            // 2 points
            var u = new PrimePolynomialFP(new List<BigInteger> { (p1.x * p2.x).value, (-p1.x - p2.x).value, 1 }, field);

            var v0 = (p2.x * p1.y - p1.x * p2.y) * (p2.x - p1.x).Inverse();
            var v1 = (p2.y - p1.y) * (p2.x - p1.x).Inverse();
            var v = new PrimePolynomialFP(new List<FPInt> { v0, v1 }, field);

            var d = new PrimeDivisorFP(u, v);
            d.Print();

            var t = DateTime.Now;
            field.method = FP.AdditionMethod.cantor;
            var r = field.ScalarMult(d, BigInteger.Parse("10000000000000000000000000000"), h, f);
            Console.WriteLine((DateTime.Now - t).TotalSeconds);
            r.Print();

            t = DateTime.Now;
            field.method = FP.AdditionMethod.cantorOptimized;
            r = field.ScalarMult(d, BigInteger.Parse("10000000000000000000000000000"), h, f);
            Console.WriteLine((DateTime.Now - t).TotalSeconds);
            r.Print();

            t = DateTime.Now;
            field.method = FP.AdditionMethod.cantorProjective;
            r = field.ScalarMult(d, BigInteger.Parse("10000000000000000000000000000"), h, f);
            Console.WriteLine((DateTime.Now - t).TotalSeconds);
            r.Print();

            t = DateTime.Now;
            field.method = FP.AdditionMethod.cantorNew;
            r = field.ScalarMult(d, BigInteger.Parse("10000000000000000000000000000"), h, f);
            Console.WriteLine((DateTime.Now - t).TotalSeconds);
            r.Print();
        }
        public static void fun6()
        {
            //FP field = new FP(11);
            //PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            //PrimePolynomialFP f = new(new List<BigInteger> { 2, 1, 7, 3, 0, 1 }, field);

            //FP field = new FP(BigInteger.Parse("129899216730745422747379980509"));
            //PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);
            //PrimePolynomialFP f = new(new List<BigInteger> {
            //    BigInteger.Parse("84554758087342364918400819975"),
            //    BigInteger.Parse("43842353021798749401773327333"),
            //    BigInteger.Parse("81016668528840831602117585943"),
            //    BigInteger.Parse("52690279977948565928475983553"),
            //    BigInteger.Parse("110854065858994061391078560211"),
            //    1 }, field);

            Dictionary<FP.AdditionMethod, List<Info>> addInfo = new();
            addInfo[FP.AdditionMethod.cantor] = new();
            addInfo[FP.AdditionMethod.cantorOptimized] = new();
            addInfo[FP.AdditionMethod.cantorProjective] = new();
            addInfo[FP.AdditionMethod.cantorNew] = new();

            Dictionary<FP.AdditionMethod, List<Info>> doubleInfo = new();
            doubleInfo[FP.AdditionMethod.cantor] = new();
            doubleInfo[FP.AdditionMethod.cantorOptimized] = new();
            doubleInfo[FP.AdditionMethod.cantorProjective] = new();
            doubleInfo[FP.AdditionMethod.cantorNew] = new();

            Dictionary<FP.AdditionMethod, List<Info>> scalarInfo = new();
            scalarInfo[FP.AdditionMethod.cantor] = new();
            scalarInfo[FP.AdditionMethod.cantorOptimized] = new();
            scalarInfo[FP.AdditionMethod.cantorProjective] = new();
            scalarInfo[FP.AdditionMethod.cantorNew] = new();

            var rand = new Random();

            int iter = 3;
            foreach (var prime in primes)
            {
                Console.WriteLine(iter);
                iter += 1;

                FP field = new FP(BigInteger.Parse(prime));
                PrimePolynomialFP h = new(new List<BigInteger> { 0 }, field);

                int divisorsIters = 30;
                int fIters = 30;

                var fs = new List<PrimePolynomialFP>();
                for (int i = 0; i < fIters; i++)
                    fs.Add(GenerateF(field, rand));

                var cases = new List<List<PrimeDivisorFP[]>>();
                var scalars = new List<List<BigInteger>>();

                for (int i = 0; i < fIters; i++)
                {
                    cases.Add(new());
                    scalars.Add(new());
                    var f = fs[i];

                    for (int j = 0; j < divisorsIters; j++)
                    {

                        PrimeDivisorFP[] div = new PrimeDivisorFP[] { null, null };

                        for (int k = 0; k < 2; k++)
                        {
                            var p1 = field.RandomPoint(f, rand);
                            var p2 = field.RandomPoint(f, rand);

                            // 2 points
                            var u = new PrimePolynomialFP(new List<BigInteger> { (p1.x * p2.x).value, (-p1.x - p2.x).value, 1 }, field);

                            var v0 = (p2.x * p1.y - p1.x * p2.y) * (p2.x - p1.x).Inverse();
                            var v1 = (p2.y - p1.y) * (p2.x - p1.x).Inverse();
                            var v = new PrimePolynomialFP(new List<FPInt> { v0, v1 }, field);

                            var newDiv = new PrimeDivisorFP(u, v);

                            div[k] = newDiv;
                        }

                        cases[i].Add(div);
                        scalars[i].Add(rand.Next());
                    }
                }

                //addInfo[FP.AdditionMethod.cantor].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantor, false));
                //addInfo[FP.AdditionMethod.cantorOptimized].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantorOptimized, false));
                //addInfo[FP.AdditionMethod.cantorProjective].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantorProjective, false));
                //addInfo[FP.AdditionMethod.cantorNew].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantorNew, false));

                //doubleInfo[FP.AdditionMethod.cantor].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantor, true));
                //doubleInfo[FP.AdditionMethod.cantorOptimized].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantorOptimized, true));
                //doubleInfo[FP.AdditionMethod.cantorProjective].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantorProjective, true));
                //doubleInfo[FP.AdditionMethod.cantorNew].Add(TestAlgorithm(field, h, fs, cases, FP.AdditionMethod.cantorNew, true));

                scalarInfo[FP.AdditionMethod.cantor].Add(TestScalar(field, h, fs, cases, scalars, FP.AdditionMethod.cantor));
                scalarInfo[FP.AdditionMethod.cantorOptimized].Add(TestScalar(field, h, fs, cases, scalars, FP.AdditionMethod.cantorOptimized));
                scalarInfo[FP.AdditionMethod.cantorProjective].Add(TestScalar(field, h, fs, cases, scalars, FP.AdditionMethod.cantorProjective));
                scalarInfo[FP.AdditionMethod.cantorNew].Add(TestScalar(field, h, fs, cases, scalars, FP.AdditionMethod.cantorNew));
            }

            //Console.WriteLine("cantor add:");
            //foreach (var info in addInfo[FP.AdditionMethod.cantor])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            //Console.WriteLine("cantorOptimized add:");
            //foreach (var info in addInfo[FP.AdditionMethod.cantorOptimized])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            //Console.WriteLine("cantorProjective add:");
            //foreach (var info in addInfo[FP.AdditionMethod.cantorProjective])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            //Console.WriteLine("cantorNew add:");
            //foreach (var info in addInfo[FP.AdditionMethod.cantorNew])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            /////////////////////////////////////////////////////////

            //Console.WriteLine("cantor double:");
            //foreach (var info in doubleInfo[FP.AdditionMethod.cantor])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            //Console.WriteLine("cantorOptimized double:");
            //foreach (var info in doubleInfo[FP.AdditionMethod.cantorOptimized])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            //Console.WriteLine("cantorProjective double:");
            //foreach (var info in doubleInfo[FP.AdditionMethod.cantorProjective])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            //Console.WriteLine("cantorNew double:");
            //foreach (var info in doubleInfo[FP.AdditionMethod.cantorNew])
            //    Console.WriteLine(info);
            //Console.WriteLine("");

            /////////////////////////////////////////////////////////

            Console.WriteLine("cantor scalar:");
            foreach (var info in scalarInfo[FP.AdditionMethod.cantor])
                Console.WriteLine(info);
            Console.WriteLine("");

            Console.WriteLine("cantorOptimized scalar:");
            foreach (var info in scalarInfo[FP.AdditionMethod.cantorOptimized])
                Console.WriteLine(info);
            Console.WriteLine("");

            Console.WriteLine("cantorProjective scalar:");
            foreach (var info in scalarInfo[FP.AdditionMethod.cantorProjective])
                Console.WriteLine(info);
            Console.WriteLine("");

            Console.WriteLine("cantorNew scalar:");
            foreach (var info in scalarInfo[FP.AdditionMethod.cantorNew])
                Console.WriteLine(info);
            Console.WriteLine("");
        }
        public static Info TestAlgorithm(FP field, PrimePolynomialFP h, List<PrimePolynomialFP> fs, List<List<PrimeDivisorFP[]>> cases, FP.AdditionMethod method, bool doubling)
        {
            var time = DateTime.Now;
            FPInt.multiplications = 0; FPInt.inversions = 0;

            for (int i = 0; i < fs.Count; i++)
            {
                var f = fs[i];

                foreach (var pair in cases[i])
                {
                    var a = pair[0];
                    var b = pair[1];
                    if (doubling)
                        b = pair[0];

                    switch (method)
                    {
                        case FP.AdditionMethod.cantor:
                            field.Cantor(a, b, h, f);
                            break;
                        case FP.AdditionMethod.cantorOptimized:
                            field.CantorOptimized(a, b, h, f);
                            break;
                        case FP.AdditionMethod.cantorProjective:
                            field.CantorProjective(new(a), new(b), h, f);
                            break;
                        case FP.AdditionMethod.cantorNew:
                            field.CantorNew(new(a), new(b), h, f);
                            break;
                    }
                }
            }

            //Console.WriteLine(FPInt.multiplications + " " + FPInt.inversions + " " + (DateTime.Now - time).TotalSeconds);
            return new(FPInt.multiplications, FPInt.inversions, (DateTime.Now - time).TotalSeconds);
        }
        public static Info TestScalar(FP field, PrimePolynomialFP h, List<PrimePolynomialFP> fs, List<List<PrimeDivisorFP[]>> cases, List<List<BigInteger>> scalars, FP.AdditionMethod method)
        {
            var time = DateTime.Now;
            FPInt.multiplications = 0; FPInt.inversions = 0;
            field.method = method;

            for (int i = 0; i < fs.Count; i++)
            {
                var f = fs[i];

                for (int j = 0; j < cases[i].Count; j++)
                {
                    var a = cases[i][j][0];
                    var scalar = scalars[i][j];
                    field.ScalarMult(a, scalar, h, f);
                }
            }

            //Console.WriteLine(FPInt.multiplications + " " + FPInt.inversions + " " + (DateTime.Now - time).TotalSeconds);
            return new(FPInt.multiplications, FPInt.inversions, (DateTime.Now - time).TotalSeconds);
        }
        public static BigInteger FindPrime(BigInteger lowerBound)
        {
            lowerBound += 12 - lowerBound % 6;
            int[] ks = new int[] { -1, 1};

            while (true)
            {
                foreach (var k1 in ks)
                {
                    var test = lowerBound + k1;
                    BigInteger cp = 6;
                    bool prime = true;

                    //Console.WriteLine("test " + test);

                    while ((cp - 1) * (cp - 1) <= test && prime != false) {
                        //Console.WriteLine("cp " + cp);
                        foreach (var k2 in ks)
                        {
                            if (test % (cp + k2) == 0)
                                prime = false;
                        }
                        cp += 6;
                    }

                    if (prime == true)
                        return test;
                }

                lowerBound += 6;
            }
        }
        public static PrimePolynomialFP GenerateF(FP field, Random rand)
        {
            var coeffs = new List<FPInt> { null, null, null, null, null, FPInt.One(field) };

            while (true)
            {
                for (int i = 0; i < 5; i++)
                    coeffs[i] = new(rand.Next(), field);

                var f = new PrimePolynomialFP(coeffs, field);
                var fd = f.Derivative();

                PrimePolynomialFP t, s;
                if (field.EuclidPoly(f, fd, out t, out s).Degree() == 0)
                    return f;
            }
        }

        public static void Print(List<int> polynomial)
        {
            for (int i = 0; i < polynomial.Count; i++)
                Console.Write(polynomial[i] + " ");
            Console.WriteLine();
        }
        public static void Print(List<BigInteger> polynomial)
        {
            for (int i = 0; i < polynomial.Count; i++)
                Console.Write(polynomial[i] + " ");
            Console.WriteLine();
        }
        
        public static void Print(GFElement element)
        {
            Print(element.p);
        }
        public static void Print(ExtensionElement element)
        {
            foreach (GFElement el in element.p)
                Print(el);
            Console.WriteLine();
        }
    }
}
