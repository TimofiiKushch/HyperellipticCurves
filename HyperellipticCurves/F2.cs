using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace HyperellipticCurves
{
    public class F2 : PrimeField
    {
        public F2() : base(2)
        {

        }

        override public int Add(int a, int b)
        {
            return a ^ b;
        }

        override public int Subtract(int a, int b)
        {
            return a ^ b;
        }

        override public int Multiply(int a, int b)
        {
            return a & b;
        }
        override public int Inverse(int a)
        {
            return 1;
        }
        override public int LeadingCoeff(List<int> polynomial)
        {
            return 1;
        }

    }
}
