using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace HyperellipticCurves
{
    class Program
    {
        static void Main(string[] args)
        {
            //var encodedMessage = File.ReadAllBytes("Short Signatures from the Weil Pairing.pdf");
            //SignatureSchemeUpdated(79, false, encodedMessage);

            // 79-, 97+, 149+, 163+, 163-, 167+

            float failureP = 1e-6f;
            var message = Encoding.ASCII.GetBytes("Hello world!");

            //SignatureScheme(79, false, message, failureP);
            SignatureSchemeUpdated(79, false, message, failureP);
        }

        // l - dimension of the field F_3^l
        // if positive then curve is y^2 = x^3 + 2x + 1, else y^2 = x^3 + 2x + 2
        public static void SignatureScheme(int l, bool positive, byte[] message, float failureP)
        {
            (EllipticCurveManager cm, ECP<int> p, ECP<int> sm, ECP<int> r) = Sign_v1(l, positive, message, failureP);
            Verify_v1(cm, p, sm, r, message, failureP);
        }
        public static void SignatureSchemeUpdated(int l, bool positive, byte[] message, float failureP)
        {
            (EllipticCurveManager cm, ECP<int> p, ECP<int> sm, ECP<int> r) = Sign_v2(l, positive, message, failureP);
            Verify_v2(cm, p, sm, r, message, failureP);
        }

        public static (EllipticCurveManager, ECP<int>, ECP<int>, ECP<int>) Sign_v1(int l, bool positive, byte[] message, float failurePl)
        {
            var rand = new Random();
            var curTime = DateTime.Now;

            // Step 1: keygen
            var a = new List<int> { 2 };
            var b = new List<int> { positive ? 1 : 2 };
            var cm = new EllipticCurveManager(a, b, l, positive, true);

            // let's find P of order q
            var q = cm.q;
            var p = cm.pq;

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

            Console.WriteLine("Signature scheme v1: generated keys");
            Console.WriteLine($"Last step took: {(DateTime.Now - curTime).TotalSeconds}");
            curTime = DateTime.Now;

            // Step 2: signing
            var pm = cm.MapToGroup(message, 1e-6f); // probability of a failure
            var sm = cm.curve.Mult(pm, x);

            Console.WriteLine("Signature scheme v1: signed the message");
            Console.WriteLine($"Last step took: {(DateTime.Now - curTime).TotalSeconds}");
            curTime = DateTime.Now;
            return (cm, p, sm, r);
        }
        public static void Verify_v1(EllipticCurveManager cm, ECP<int> p, ECP<int> sm, ECP<int> r, byte[] message, float failureP)
        {
            var curTime = DateTime.Now;
            // Step 3: verification
            // need a point with matching x-coordinate
            var y = GaloisField<int>.Root(cm.curve.RightSide(sm.x));
            if (y is null)
                throw new Exception("No roots found: signature is invalid");
            var s = new ECP<int>(sm.x, y);

            var pm = cm.MapToGroup(message, failureP);

            // convert everything to work in F_3^6l
            var s6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { s.x }, cm.curve6.field), new(new List<GFElement<int>> { s.y }, cm.curve6.field));
            var p6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { p.x }, cm.curve6.field), new(new List<GFElement<int>> { p.y }, cm.curve6.field));
            var r6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { r.x }, cm.curve6.field), new(new List<GFElement<int>> { r.y }, cm.curve6.field));
            var pm6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { pm.x }, cm.curve6.field), new(new List<GFElement<int>> { pm.y }, cm.curve6.field));

            var u = EllipticCurveManager.WeilPairing(p6, cm.Auto6L(s6), cm.q, cm.curve6, cm.s);
            var v = EllipticCurveManager.WeilPairing(r6, cm.Auto6L(pm6), cm.q, cm.curve6, cm.s);

            if (u == v || u == cm.curve6.field.Inverse(v))
            {
                Console.WriteLine("Signature scheme v1: verified the message");
                Console.WriteLine($"Last step took: {(DateTime.Now - curTime).TotalSeconds}");
                Console.WriteLine();
            }
            else
                throw new Exception("Signature rejected");
        }
        public static (EllipticCurveManager, ECP<int>, ECP<int>, ECP<int>) Sign_v2(int l, bool positive, byte[] message, float failurePl)
        {
            var rand = new Random();
            var curTime = DateTime.Now;

            // Step 1: keygen
            var a = new List<int> { 2 };
            var b = new List<int> { positive ? 1 : 2 };

            var cm = new EllipticCurveManager(a, b, l, positive, true);

            // let's find P of order q
            // q is prime so just need qP = 0
            var q = cm.q;
            var p = cm.pq;

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

            Console.WriteLine("Signature scheme v2: generated keys");
            Console.WriteLine($"Last step took: {(DateTime.Now - curTime).TotalSeconds}");
            curTime = DateTime.Now;

            // Step 2: signing
            var pm = cm.MapToGroupUpdated(message, 1e-6f); // probability of a failure
            var sm = cm.curve.Mult(pm, x);

            Console.WriteLine("Signature scheme v2: signed the message");
            Console.WriteLine($"Last step took: {(DateTime.Now - curTime).TotalSeconds}");
            curTime = DateTime.Now;
            return (cm, p, sm, r);
        }
        public static void Verify_v2(EllipticCurveManager cm, ECP<int> p, ECP<int> sm, ECP<int> r, byte[] message, float failureP)
        {
            var curTime = DateTime.Now;
            // Step 3: verification
            // need a point with matching x-coordinate
            var leftSide = sm.y * sm.y - cm.curve.b;
            if (cm.solver.Trace(leftSide) != 0)
                throw new Exception("No roots found: signature is invalid");
            var root = cm.solver.Solve(leftSide);
            var rootPoly = root.p;
            rootPoly[0] = sm.x.p[0];
            var xt = new GFElement<int>(rootPoly, cm.curve.field); // TODO check order q

            var s = new ECP<int>(xt, sm.y);
            var pm = cm.MapToGroupUpdated(message, failureP);

            // convert everything to work in F_3^6l
            var s6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { s.x }, cm.curve6.field), new(new List<GFElement<int>> { s.y }, cm.curve6.field));
            var p6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { p.x }, cm.curve6.field), new(new List<GFElement<int>> { p.y }, cm.curve6.field));
            var r6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { r.x }, cm.curve6.field), new(new List<GFElement<int>> { r.y }, cm.curve6.field));
            var pm6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { pm.x }, cm.curve6.field), new(new List<GFElement<int>> { pm.y }, cm.curve6.field));

            var u = EllipticCurveManager.WeilPairing(p6, cm.Auto6L(s6), cm.q, cm.curve6, cm.s);
            var v = EllipticCurveManager.WeilPairing(r6, cm.Auto6L(pm6), cm.q, cm.curve6, cm.s);

            if (u == v)
            {
                Console.WriteLine("Signature scheme v2: verified the message");
                Console.WriteLine($"Last step took: {(DateTime.Now - curTime).TotalSeconds}");
                Console.WriteLine();
            }
            else
                throw new Exception("Signature rejected");
        }
        public static void TestWeilPairing()
        {
            var field = new GaloisField<int>(new PrimeField(631), new List<int> { -2, 1 });
            var poly = new List<int>(Enumerable.Repeat(0, 7));
            poly[0] = 2;
            poly[1] = 1;
            poly[6] = 1;
            var field6 = new GaloisField<GFElement<int>>(field, poly);

            var a = new GFElement<int>(new List<int> { 30 }, field);
            var b = new GFElement<int>(new List<int> { 34 }, field);

            var a6 = new GFElement<GFElement<int>>(new List<GFElement<int>> { a }, field6);
            var b6 = new GFElement<GFElement<int>>(new List<GFElement<int>> { b }, field6);

            var curve = new EllipticCurve<GFElement<int>>(a6, b6);
            var p = new ECP<GFElement<int>>(
                new(new List<GFElement<int>> { field.Scalar(36) }, field6),
                new(new List<GFElement<int>> { field.Scalar(60) }, field6));
            var q = new ECP<GFElement<int>>(
                new(new List<GFElement<int>> { field.Scalar(121) }, field6),
                new(new List<GFElement<int>> { field.Scalar(387) }, field6));

            //var s = new(new(new List<GFElement> { curve.field.baseField.Scalar(0) }, curve.field), new(new List<GFElement> { curve.field.baseField.Scalar(36) }, curve.field));
            //var res = EllipticCurveManager.WeilPairing(p, q, 5, curve);
            //Print(res);
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
        
        public static void Print(GFElement<int> element)
        {
            Print(element.p);
        }
        public static void Print(GFElement<GFElement<int>> element)
        {
            foreach (var el in element.p)
                Print(el);
            Console.WriteLine();
        }
    }
}
