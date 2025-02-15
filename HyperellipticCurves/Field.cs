using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HyperellipticCurves
{
    //public static List<int> FindPrimitivePolynomial(int p, int n)
    //{
    //    List<int> primes_pn; 
    //    List<int> degrees_pn;
    //    Factorize(((int)Math.Pow(p, n) - 1) / (p - 1), out primes_pn, out degrees_pn);

    //    List<int> primes_p1;
    //    List<int> degrees_p1;
    //    Factorize((p - 1), out primes_p1, out degrees_p1);

    //    var rand = new Random();

    //    List<int> l = new List<int>(n+1);
    //    for (int i = 0; i < n + 1; i++)
    //        l.Add(0);
    //    l[n] = 1;

    //    while (true)
    //    {
    //        bool bad = false;
    //        // 1
    //        for (int i = 0; i < n; i++)
    //            l[i] = rand.Next(0, p);

    //        // 2
    //        if (l[0] == 0)
    //            break;

    //        foreach (int prime in primes_p1)
    //        {
    //            int a = l[0];
    //            for (int i = 1; i < (p - 1) / prime; i++)
    //                a = (a * l[0]) % p;
    //            if (a == 1)
    //            {
    //                bad = true;
    //                break;
    //            }
    //        }

    //        if (bad)
    //            continue;

    //        // 3

    //    }

    //}

    //public static void Factorize(int number, out List<int> primes, out List<int> degrees)
    //{
    //    primes = new List<int>();
    //    degrees = new List<int>();

    //    while (number > 1)
    //    {
    //        bool found = false;
    //        for (int i = 2; i < 1000; i++)
    //            if (number % i == 0)
    //            {
    //                if (primes.Contains(i))
    //                    degrees[primes.IndexOf(i)] += 1;
    //                else
    //                {
    //                    primes.Add(i);
    //                    degrees.Add(1);
    //                }

    //                number /= i;
    //                found = true;
    //                break;
    //            }

    //        if (found)
    //            continue;

    //        // pollard
    //        int a = 2;
    //        int b = 2;
    //        while (true)
    //        {
    //            a = ((int)Math.Pow(a,2) + 1) % number;
    //            b = ((int)Math.Pow(b,2) + 1) % number;
    //            b = ((int)Math.Pow(b,2) + 1) % number;

    //            int d = PrimeField.NumericalEuclid((a - b) % number, number);
    //            if (d > 1 && d < number)
    //            {
    //                if (primes.Contains(d))
    //                    degrees[primes.IndexOf(d)] += 1;
    //                else
    //                {
    //                    primes.Add(d);
    //                    degrees.Add(d);
    //                }

    //                number /= d;
    //                break;
    //            }
    //        }
    //    }

    //}
}
