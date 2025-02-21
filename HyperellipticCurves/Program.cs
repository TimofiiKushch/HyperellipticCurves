using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace HyperellipticCurves
{
    class Program
    {
        static void Main(string[] args)
        {
            //TestWeilPairing();

            // TODO: debug map to group
            SignatureScheme(79, false);

            //while (true)
            //    SignatureScheme(7, false);

        }

        public static void SignatureScheme(int l, bool positive)
        {
            var rand = new Random();

            // Step 1: keygen
            var a = new List<int> { 2 };
            var b = new List<int> { positive ? 1 : 2 };

            var cm = new EllipticCurveManager(a, b, l, positive, true);

            // let's find P of order q
            // q is prime so just need qP = 0
            // TODO: change so we find primitive and raise it to the power of m/q
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

            x = 2; // test
            var r = cm.curve.Mult(p, x);

            // Step 2: signing
            var message = "Hello world!";
            var encodedMessage = Encoding.ASCII.GetBytes(message);
            var pm = cm.MapToGroup(encodedMessage, 1e-6f); // probability of a failure
            var sm = cm.curve.Mult(pm, x);

            // Step 3: verification
            // need a point with matching x-coordinate
            var y = GaloisField<int>.Root(cm.curve.RightSide(sm.x));
            if (y is null)
                throw new Exception("No roots found: signature is invalid");
            var s = new ECP<int>(sm.x, y);

            // convert everything to work in field6
            var s6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { s.x }, cm.curve6.field), new(new List<GFElement<int>> { s.y }, cm.curve6.field));
            var p6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { p.x }, cm.curve6.field), new(new List<GFElement<int>> { p.y }, cm.curve6.field));
            var r6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { r.x }, cm.curve6.field), new(new List<GFElement<int>> { r.y }, cm.curve6.field));
            var pm6 = new ECP<GFElement<int>>(new(new List<GFElement<int>> { pm.x }, cm.curve6.field), new(new List<GFElement<int>> { pm.y }, cm.curve6.field));

            var u = EllipticCurveManager.WeilPairing(p6, cm.Auto6L(s6), q, cm.curve6, cm.s);
            var v = EllipticCurveManager.WeilPairing(r6, cm.Auto6L(pm6), q, cm.curve6, cm.s);

            if (u == v || u == cm.curve6.field.Inverse(v))
                Console.WriteLine("Congratulations!");
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

            var res = EllipticCurveManager.WeilPairing(p, q, 5, curve);
            Print(res);
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
