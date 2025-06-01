# include <cmath>
# include <cstdlib>
# include <cstring>
# include <fstream>
# include <iomanip>
# include <iostream>

using namespace std;

# include "tetrahedron.hpp"

namespace tetrahedron {

//****************************************************************************80

int i4_gcd ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    i4_gcd() finds the greatest common divisor of two I4's.
//
//  Discussion:
//
//    Note that only the absolute values of I and J are
//    considered, so that the result is always nonnegative.
//
//    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
//
//    If I and J have no common factor, I4_GCD is returned as 1.
//
//    Otherwise, using the Euclidean algorithm, I4_GCD is the
//    greatest common divisor of I and J.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, two numbers whose GCD is desired.
//
//    Output, int I4_GCD, the greatest common divisor of I and J.
//
{
  int p;
  int q;
  int r;
//
//  Return immediately if either I or J is zero.
//
  if ( i == 0 )
  {
    q = i4_max ( 1, abs ( j ) );
    return q;
  }
  else if ( j == 0 )
  {
    q = i4_max ( 1, abs ( i ) );
    return q;
  }
//
//  Set IP to the larger of I and J, IQ to the smaller.
//  This way, we can alter IP and IQ as we go.
//
  p = i4_max ( abs ( i ), abs ( j ) );
  q = i4_min ( abs ( i ), abs ( j ) );
//
//  Carry out the Euclidean algorithm.
//
  for ( ; ; )
  {
    r = p % q;

    if ( r == 0 )
    {
      break;
    }
    p = q;
    q = r;
  }

  return q;
}
//****************************************************************************80

int i4_lcm ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    i4_lcm() computes the least common multiple of two I4's.
//
//  Discussion:
//
//    The least common multiple may be defined as
//
//      LCM(I,J) = ABS( I * J ) / GCF(I,J)
//
//    where GCF(I,J) is the greatest common factor of I and J.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the integers whose LCM is desired.
//
//    Output, int I4_LCM, the least common multiple of I and J.
//    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
//
{
  int value;

  value = abs ( i * ( j / i4_gcd ( i, j ) ) );

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    i4_max() returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    i4_min() returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    i4vec_copy() copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

int i4vec_lcm ( int n, int v[] )

//****************************************************************************80
//
//  Purpose:
//
//    i4vec_lcm() returns the least common multiple of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The value LCM returned has the property that it is the smallest integer
//    which is evenly divisible by every element of V.
//
//    The entries in V may be negative.
//
//    If any entry of V is 0, then LCM is 0.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    02 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of V.
//
//    Input, int V[N], the vector.
//
//    Output, int I4VEC_LCM, the least common multiple of V.
//
{
  int i;
  int lcm;

  lcm = 1;

  for ( i = 0; i < n; i++ )
  {
    if ( v[i] == 0 )
    {
      lcm = 0;
      break;
    }
    lcm = i4_lcm ( lcm, v[i] );
  }
  return lcm;
}
//****************************************************************************80

double polygon_area_3d ( int n, double v[], double normal[] )

//****************************************************************************80
//
//  Purpose:
//
//    polygon_area_3d() computes the area of a polygon in 3D.
//
//  Discussion:
//
//    The computation is not valid unless the vertices really do lie
//    in a plane, so that the polygon that is defined is "flat".
//
//    The polygon does not have to be "regular", that is, neither its
//    sides nor its angles need to be equal.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allen Van Gelder,
//    Efficient Computation of Polygon Area and Polyhedron Volume,
//    Graphics Gems V, edited by Alan Paeth,
//    AP Professional, 1995, T385.G6975.
//
//  Parameters:
//
//    Input, int N, the number of vertices.
//
//    Input, double V[3*N], the coordinates of the vertices.
//    The vertices should be listed in neighboring order.
//
//    Output, double NORMAL[3], the unit normal vector to the polygon.
//
//    Output, double POLYGON_AREA_3D, the area of the polygon.
//
{
# define DIM_NUM 3

  double area;
  int i;
  int ip1;

  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;

  for ( i = 0; i < n; i++ )
  {
    if ( i < n - 1 )
    {
      ip1 = i + 1;
    }
    else
    {
      ip1 = 0;
    }
//
//  Compute the cross product and add it to NORMAL.
//
    normal[0] = normal[0] + v[1+i*3] * v[2+ip1*3] - v[2+i*3] * v[1+ip1*3];
    normal[1] = normal[1] + v[2+i*3] * v[0+ip1*3] - v[0+i*3] * v[2+ip1*3];
    normal[2] = normal[2] + v[0+i*3] * v[1+ip1*3] - v[1+i*3] * v[0+ip1*3];
  }

  area = r8vec_norm ( DIM_NUM, normal );

  if ( area != 0.0 )
  {
    normal[0] = normal[0] / area;
    normal[1] = normal[1] / area;
    normal[2] = normal[2] / area;
  }
  else
  {
    normal[0] = 1.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
  }

  area = 0.5 * area;

  return area;
# undef DIM_NUM
}
//****************************************************************************80

double r8_acos ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    r8_acos() computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double R8_ACOS, an angle whose cosine is C.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( c <= -1.0 )
  {
    value = r8_pi;
  }
  else if ( 1.0 <= c )
  {
    value = 0.0;
  }
  else
  {
    value = acos ( c );
  }
  return value;
}
//****************************************************************************80

void r8_swap ( double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    r8_swap() switches two R8's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = *x;
  *x = *y;
  *y = z;

  return;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    r8_uniform_01() returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  const int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_copy ( int m, int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8mat_copy() copies one R8MAT to another.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
}
//****************************************************************************80

double r8mat_det_4d ( double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8mat_det_4d() computes the determinant of a 4 by 4 R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
//****************************************************************************80

int r8mat_solve ( int n, int rhs_num, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8mat_solve() uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R8MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a pivot row.
//
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( fabs ( apivot ) < fabs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    r8mat_transpose_print() prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, an optional title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    r8mat_transpose_print_some() prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double r8vec_angle_3d ( double u[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_angle_3d() computes the angle between two vectors in 3D.
//
//  Modified:
//
//    07 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double U[3], V[3], the vectors.
//
//    Output, double ANGLE, the angle between the two vectors.
//
{
  double angle;
  double angle_cos;
  double u_norm;
  double uv_dot;
  double v_norm;

  uv_dot = r8vec_dot ( 3, u, v );

  u_norm = sqrt ( r8vec_dot ( 3, u, u ) );

  v_norm = sqrt ( r8vec_dot ( 3, v, v ) );

  angle_cos = uv_dot / u_norm / v_norm;

  angle = r8_acos ( angle_cos );

  return angle;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_copy() copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double *r8vec_cross_3d ( double v1[3], double v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_cross_3d() computes the cross product of two R8VEC's in 3D.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_3D[3], the cross product vector.
//
{
  double *v3;

  v3 = new double[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
//****************************************************************************80

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_dot() computes the dot product of a pair of R8VEC's in ND.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }

  return value;
}
//****************************************************************************80

double r8vec_length ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_length() returns the Euclidean length of an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the vector.
//
//    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = value + pow ( x[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_max() returns the value of the maximum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = - HUGE_VAL;

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_norm ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_norm() returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM, the L2 norm of A.
//
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_print() prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_transpose_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_transpose_print() prints an R8VEC "transposed".
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Example:
//
//    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
//    TITLE = 'My vector:  '
//
//    My vector:   1.0    2.1    3.2    4.3    5.4
//                 6.5    7.6    8.7    9.8   10.9
//                11.0
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    11 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;

  cout << title << "\n";

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    cout << "  ";
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    r8vec_zero() zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

void shape_print ( int point_num, int face_num, int face_order_max,
  double point_coord[], int face_order[], int face_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    shape_print() prints information about a polyhedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_ORDER_MAX, the number of vertices per face.
//
//    Input, double POINT_COORD[DIM_NUM*POINT_NUM], the point coordinates.
//
//    Input, int FACE_ORDER[FACE_NUM], the number of vertices per face.
//
//    Input, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
//    contains the index of the I-th point in the J-th face.  The
//    points are listed in the counter clockwise direction defined
//    by the outward normal at the face.
//
{
# define DIM_NUM 3

  int i;
  int j;

  cout << "\n";
  cout << "shape_print()\n";
  cout << "  Information about a polytope.\n";
  cout << "\n";
  cout << "  The number of vertices is " << point_num << "\n";
  cout << "\n";
  cout << "  Vertices:\n";
  cout << "\n";
  cout << "     Index          X               Y               Z\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(8) << j + 1 << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setprecision(8) << setw(16) << point_coord[i+j*DIM_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The number of faces is " << face_num << "\n";
  cout << "  The maximum order of any face is " << face_order_max << "\n";
  cout << "\n";
  cout << "     Index     Order         Indices of Nodes in Face\n";
  for ( j = 1; j <= face_order_max; j++ )
  {
    cout << setw(8) << j;
  }
  cout << "\n";
  cout << "                      ";
  cout << "\n";

  for ( j = 0; j < face_num; j++ )
  {
    cout << "  " << setw(8) << j + 1
         << "  " << setw(8) << face_order[j]
         << "  ";
    for ( i = 0; i < face_order[j]; i++ )
    {
      cout << setw(8) << face_point[i+j*face_order_max];
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

double *tetrahedron_barycentric ( double tetra[3*4], double p[3] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_barycentric() returns the barycentric coordinates of a point.
//
//  Discussion:
//
//    The barycentric coordinates of a point P with respect to
//    a tetrahedron are a set of four values C(1:4), each associated
//    with a vertex of the tetrahedron.  The values must sum to 1.
//    If all the values are between 0 and 1, the point is contained
//    within the tetrahedron.
//
//    The barycentric coordinate of point X related to vertex A can be
//    interpreted as the ratio of the volume of the tetrahedron with
//    vertex A replaced by vertex X to the volume of the original
//    tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Input, double P[3], the point to be checked.
//
//    Output, double C[4], the barycentric coordinates of the point with
//    respect to the tetrahedron.
//
{
# define N 3
# define RHS_NUM 1

  double a[N*(N+RHS_NUM)];
  double *c;
  int info;
//
//  Set up the linear system
//
//    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
//    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
//    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
//
//  which is satisfied by the barycentric coordinates.
//

  a[0+0*N] = tetra[0+1*3] - tetra[0+0*3];
  a[1+0*N] = tetra[1+1*3] - tetra[1+0*3];
  a[2+0*N] = tetra[2+1*3] - tetra[2+0*3];

  a[0+1*N] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*N] = tetra[1+2*3] - tetra[1+0*3];
  a[2+1*N] = tetra[2+2*3] - tetra[2+0*3];

  a[0+2*N] = tetra[0+3*3] - tetra[0+0*3];
  a[1+2*N] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*N] = tetra[2+3*3] - tetra[2+0*3];

  a[0+3*N] = p[0]         - tetra[0+0*3];
  a[1+3*N] = p[1]         - tetra[1+0*3];
  a[2+3*N] = p[2]         - tetra[2+0*3];
//
//  Solve the linear system.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "tetrahedron_barycentric(): Fatal error!\n";
    cerr << "  The linear system is singular.\n";
    cerr << "  The input data does not form a proper tetrahedron.\n";
    exit ( 1 );
  }

  c = new double[4];

  c[1] = a[0+3*N];
  c[2] = a[1+3*N];
  c[3] = a[2+3*N];

  c[0] = 1.0 - c[1] - c[2] - c[3];

  return c;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

double *tetrahedron_centroid ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_centroid() computes the centroid of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_CENTROID[3], the coordinates of the centroid.
//
{
  double *centroid;
  int i;

  centroid = new double[3];

  centroid[0] = ( tetra[0+0*3] + tetra[0+1*3]
                + tetra[0+2*3] + tetra[0+3*3] );
  centroid[1] = ( tetra[1+0*3] + tetra[1+1*3]
                + tetra[1+2*3] + tetra[1+3*3] );
  centroid[2] = ( tetra[2+0*3] + tetra[2+1*3]
                + tetra[2+2*3] + tetra[2+3*3] );

  for ( i = 0; i < 3; i++ )
  {
    centroid[i] = centroid[i] / 4.0;
  }

  return centroid;
}
//****************************************************************************80

void tetrahedron_circumsphere ( double tetra[3*4], double &r, double pc[3] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_circumsphere() computes the circumsphere of a tetrahedron.
//
//  Discussion:
//
//    The circumsphere, or circumscribed sphere, of a tetrahedron is the
//    sphere that passes through the four vertices.  The circumsphere is not
//    necessarily the smallest sphere that contains the tetrahedron.
//
//    Surprisingly, the diameter of the sphere can be found by solving
//    a 3 by 3 linear system.  This is because the vectors P2 - P1,
//    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
//    right triangle with the diameter through P1.  Hence, the dot product of
//    P2 - P1 with that diameter is equal to the square of the length
//    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
//    the diameter vector originating at P1, and hence the radius and
//    center.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    10 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double &R, PC[3], the coordinates of the center of the
//    circumscribed sphere, and its radius.  If the linear system is
//    singular, then R = -1, PC[] = 0.
//
{
  double a[3*4];
  int info;
//
//  Set up the linear system.
//
  a[0+0*3] = tetra[0+1*3] - tetra[0+0*3];
  a[0+1*3] = tetra[1+1*3] - tetra[1+0*3];
  a[0+2*3] = tetra[2+1*3] - tetra[2+0*3];
  a[0+3*3] = pow ( tetra[0+1*3] - tetra[0+0*3], 2 )
           + pow ( tetra[1+1*3] - tetra[1+0*3], 2 )
           + pow ( tetra[2+1*3] - tetra[2+0*3], 2 );

  a[1+0*3] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*3] = tetra[1+2*3] - tetra[1+0*3];
  a[1+2*3] = tetra[2+2*3] - tetra[2+0*3];
  a[1+3*3] = pow ( tetra[0+2*3] - tetra[0+0*3], 2 )
           + pow ( tetra[1+2*3] - tetra[1+0*3], 2 )
           + pow ( tetra[2+2*3] - tetra[2+0*3], 2 );

  a[2+0*3] = tetra[0+3*3] - tetra[0+0*3];
  a[2+1*3] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*3] = tetra[2+3*3] - tetra[2+0*3];
  a[2+3*3] = pow ( tetra[0+3*3] - tetra[0+0*3], 2 )
           + pow ( tetra[1+3*3] - tetra[1+0*3], 2 )
           + pow ( tetra[2+3*3] - tetra[2+0*3], 2 );
//
//  Solve the linear system.
//
  info = r8mat_solve ( 3, 1, a );
//
//  If the system was singular, return a consolation prize.
//
  if ( info != 0 )
  {
    r = -1.0;
    r8vec_zero ( 3, pc );
    return;
  }
//
//  Compute the radius and center.
//
  r = 0.5 * sqrt
    ( a[0+3*3] * a[0+3*3]
    + a[1+3*3] * a[1+3*3]
    + a[2+3*3] * a[2+3*3] );

  pc[0] = tetra[0+0*3] + 0.5 * a[0+3*3];
  pc[1] = tetra[1+0*3] + 0.5 * a[1+3*3];
  pc[2] = tetra[2+0*3] + 0.5 * a[2+3*3];

  return;
}
//****************************************************************************80

bool tetrahedron_contains_point ( double tetra[3*4], double p[3] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_contains_point(): a tetrahedron contains a point.
//
//  Discussion:
//
//    Thanks to Saiful Akbar for pointing out that the array of barycentric
//    coordinated was not being deleted!  29 January 2006
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    29 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Input, double P[3], the point to be checked.
//
//    Output, bool tetrahedron_contains_point, is TRUE if the point is inside
//    the tetrahedron or on its boundary, and FALSE otherwise.
//
{
  double *c;
  bool value;

  c = tetrahedron_barycentric ( tetra, p );
//
//  If the point is in the tetrahedron, its barycentric coordinates
//  must be nonnegative.
//
  value = 0.0 <= c[0] &&
          0.0 <= c[1] &&
          0.0 <= c[2] &&
          0.0 <= c[3];

  delete [] c;

  return value;
}
//****************************************************************************80

double *tetrahedron_dihedral_angles ( double tetra[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_dihedral_angles() computes dihedral angles of a tetrahedron.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron,
//    which can be labeled as A, B, C and D.
//
//    Output, double TETRAHEDRON_DIHEDRAL_ANGLES[6], the dihedral angles
//    along the axes AB, AC, AD, BC, BD and CD, respectively.
//
{
  double ab[3];
  double *abc_normal;
  double *abd_normal;
  double ac[3];
  double *acd_normal;
  double ad[3];
  double *angle;
  double bc[3];
  double *bcd_normal;
  double bd[3];
  double cd[3];
  int i;
  const double r8_pi = 3.141592653589793;

  tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd );

  abc_normal = r8vec_cross_3d ( ac, ab );
  abd_normal = r8vec_cross_3d ( ab, ad );
  acd_normal = r8vec_cross_3d ( ad, ac );
  bcd_normal = r8vec_cross_3d ( bc, bd );

  angle = new double[6];

  angle[0] = r8vec_angle_3d ( abc_normal, abd_normal );
  angle[1] = r8vec_angle_3d ( abc_normal, acd_normal );
  angle[2] = r8vec_angle_3d ( abd_normal, acd_normal );
  angle[3] = r8vec_angle_3d ( abc_normal, bcd_normal );
  angle[4] = r8vec_angle_3d ( abd_normal, bcd_normal );
  angle[5] = r8vec_angle_3d ( acd_normal, bcd_normal );

  for ( i = 0; i < 6; i++ )
  {
    angle[i] = r8_pi - angle[i];
  }

  delete [] abc_normal;
  delete [] abd_normal;
  delete [] acd_normal;
  delete [] bcd_normal;

  return angle;
}
//****************************************************************************80

double *tetrahedron_edge_length ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_edge_length() returns edge lengths of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    10 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double EDGE_LENGTH[6], the length of the edges.
//
{
  double *edge_length;
  int i;
  int j1;
  int j2;
  int k;
  double v[3];

  edge_length = new double[6];

  k = 0;
  for ( j1 = 0; j1 < 3; j1++ )
  {
    for ( j2 = j1 + 1; j2 < 4; j2++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        v[i] = tetra[i+j2*3] - tetra[i+j1*3];
      }
      edge_length[k] = r8vec_length ( 3, v );
      k = k + 1;
    }
  }

  return edge_length;
}
//****************************************************************************80

void tetrahedron_edges ( double tetra[3*4], double ab[], double ac[],
  double ad[], double bc[], double bd[], double cd[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_edges() returns the edges of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    11 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double AB[3], AC[3], AD[3], BC[3], BD[3], CD[3], the edges.
//
{
  int i;
//
//  Compute the vectors that represent the sides.
//
  for ( i = 0; i < 3; i++ )
  {
    ab[i] = tetra[i+1*3] - tetra[i+0*3];
    ac[i] = tetra[i+2*3] - tetra[i+0*3];
    ad[i] = tetra[i+3*3] - tetra[i+0*3];
    bc[i] = tetra[i+2*3] - tetra[i+1*3];
    bd[i] = tetra[i+3*3] - tetra[i+1*3];
    cd[i] = tetra[i+3*3] - tetra[i+2*3];
  }

  return;
}
//****************************************************************************80

void tetrahedron_face_angles ( double tetra[], double angles[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_face_angles() returns the 12 face angles of a tetrahedron.
//
//  Discussion:
//
//    The tetrahedron has 4 triangular faces.  This routine computes the
//    3 planar angles associated with each face.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    03 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4] the tetrahedron vertices.
//
//    Output, double ANGLES[3*4], the face angles.
//
{
  double *tri;

  tri = new double[3*3];
//
//  Face 123
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+2*3];
  tri[1+2*3] = tetra[1+2*3];
  tri[2+2*3] = tetra[2+2*3];

  triangle_angles_3d ( tri, angles );
//
//  Face 124
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  triangle_angles_3d ( tri, angles+3 );
//
//  Face 134
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  triangle_angles_3d ( tri, angles+6 );
//
//  Face 234
//
  tri[0+0*3] = tetra[0+1*3];
  tri[1+0*3] = tetra[1+1*3];
  tri[2+0*3] = tetra[2+1*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  triangle_angles_3d ( tri, angles+9 );

  delete [] tri;

  return;
}
//****************************************************************************80

void tetrahedron_face_areas ( double tetra[], double areas[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_face_areas() returns the 4 face areas of a tetrahedron.
//
//  Discussion:
//
//    The tetrahedron has 4 triangular faces.  This routine computes the
//    areas associated with each face.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4] the tetrahedron vertices.
//
//    Output, double AREAS[4], the face areas.
//
{
  double *tri;

  tri = new double[3*3];
//
//  Face 123
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+2*3];
  tri[1+2*3] = tetra[1+2*3];
  tri[2+2*3] = tetra[2+2*3];

  areas[0] = triangle_area_3d ( tri );
//
//  Face 124
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  areas[1] = triangle_area_3d ( tri );
//
//  Face 134
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  areas[2] = triangle_area_3d ( tri );
//
//  Face 234
//
  tri[0+0*3] = tetra[0+1*3];
  tri[1+0*3] = tetra[1+1*3];
  tri[2+0*3] = tetra[2+1*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  areas[3] = triangle_area_3d ( tri );

  delete [] tri;

  return;
}
//****************************************************************************80

void tetrahedron_insphere ( double tetra[3*4], double &r, double pc[3] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_insphere() finds the insphere of a tetrahedron.
//
//  Discussion:
//
//    The insphere of a tetrahedron is the inscribed sphere, which touches
//    each face of the tetrahedron at a single point.
//
//    The points of contact are the centroids of the triangular faces
//    of the tetrahedron.  Therefore, the point of contact for a face
//    can be computed as the average of the vertices of that face.
//
//    The sphere can then be determined as the unique sphere through
//    the four given centroids.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Schneider, David Eberly,
//    Geometric Tools for Computer Graphics,
//    Elsevier, 2002,
//    ISBN: 1558605940,
//    LC: T385.G6974.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double &R, PC[3], the radius and the center
//    of the sphere.
//
{
  double b[4*4];
  double gamma;
  int i;
  int j;
  double l123;
  double l124;
  double l134;
  double l234;
  double *n123;
  double *n124;
  double *n134;
  double *n234;
  double v21[3];
  double v31[3];
  double v41[3];
  double v32[3];
  double v42[3];
  double v43[3];

  tetrahedron_edges ( tetra, v21, v31, v41, v32, v42, v43 );

  n123 = r8vec_cross_3d ( v21, v31 );
  n124 = r8vec_cross_3d ( v41, v21 );
  n134 = r8vec_cross_3d ( v31, v41 );
  n234 = r8vec_cross_3d ( v42, v32 );

  l123 = r8vec_length ( 3, n123 );
  l124 = r8vec_length ( 3, n124 );
  l134 = r8vec_length ( 3, n134 );
  l234 = r8vec_length ( 3, n234 );

  delete [] n123;
  delete [] n124;
  delete [] n134;
  delete [] n234;

  for ( i = 0; i < 3; i++ )
  {
    pc[i] = ( l234 * tetra[i+0*3]
            + l134 * tetra[i+1*3]
            + l124 * tetra[i+2*3]
            + l123 * tetra[i+3*3] )
            / ( l234 + l134 + l124 + l123 );
  }

  for ( j = 0; j < 4; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      b[i+j*4] = tetra[i+j*3];
    }
    b[3+j*4] = 1.0;
  }

  gamma = fabs ( r8mat_det_4d ( b ) );

  r = gamma / ( l234 + l134 + l124 + l123 );

  return;
}
//****************************************************************************80

void tetrahedron_lattice_layer_point_next ( int c[], int v[], bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_lattice_layer_point_next(): next tetrahedron lattice layer point.
//
//  Discussion:
//
//    The tetrahedron lattice layer L is bounded by the lines
//
//      0 <= X,
//      0 <= Y,
//      0 <= Z,
//      L - 1 < X / C[0] + Y / C[1] + Z/C[2] <= L.
//
//    In particular, layer L = 0 always contains the single point (0,0).
//
//    This function returns, one at a time, the points that lie within
//    a given tetrahedron lattice layer.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int C[4], coefficients defining the
//    lattice layer in the first 3 entries, and the laver index in C[3].
//    The coefficients should be positive, and C[3] must be nonnegative.
//
//    Input/output, int V[3].  On first call for a given layer,
//    the input value of V is not important.  On a repeated call for the same
//    layer, the input value of V should be the output value from the previous
//    call.  On output, V contains the next lattice layer point.
//
//    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
//    that this is the first call for a given layer.  Thereafter, the input
//    value should be the output value from the previous call.  On output,
//    MORE is TRUE if the returned value V is a new point.
//    If the output value is FALSE, then no more points were found,
//    and V was reset to 0, and the lattice layer has been exhausted.
//
{
  int c1n;
  int lhs;
  int n = 3;
  int rhs1;
  int rhs2;
//
//  Treat layer C[3] = 0 specially.
//
  if ( c[3] == 0 )
  {
    if ( !(*more) )
    {
      v[0] = 0;
      v[1] = 0;
      v[2] = 0;
      *more = true;
    }
    else
    {
      *more = false;
    }
    return;
  }
//
//  Compute the first point.
//
  if ( !(*more) )
  {
    v[0] = ( c[n] - 1 ) * c[0] + 1;
    v[1] = 0;
    v[2] = 0;
    *more = true;
  }
  else
  {
    c1n = i4vec_lcm ( n, c );

    rhs1 = c1n * ( c[n] - 1 );
    rhs2 = c1n *   c[n];
//
//  Can we simply increase X?
//
    v[0] = v[0] + 1;

    lhs = ( c1n / c[0] ) * v[0]
        + ( c1n / c[1] ) * v[1]
        + ( c1n / c[2] ) * v[2];

    if ( lhs <= rhs2 )
    {
    }
//
//  No.  Increase Y, and set X so we just exceed RHS1, if possible.
//
    else
    {
      v[1] = v[1] + 1;

      v[0] = ( c[0] * ( rhs1 - ( c1n / c[1] ) * v[1]
                             - ( c1n / c[2] ) * v[2] ) ) / c1n;
      v[0] = i4_max ( v[0], 0 );

      lhs = ( c1n / c[0] ) * v[0]
          + ( c1n / c[1] ) * v[1]
          + ( c1n / c[2] ) * v[2];

      if ( lhs <= rhs1 )
      {
        v[0] = v[0] + 1;
        lhs = lhs + c1n / c[0];
      }
//
//  We have increased Y by 1.  Have we stayed below the upper bound?
//
      if ( lhs <= rhs2 )
      {
      }
//
//  No.  Increase Z, and set X so we just exceed RHS1, if possible.
//
      else
      {
        v[2] = v[2] + 1;
        v[1] = 0;
        v[0] = ( c[0] * ( rhs1 - ( c1n / c[1] ) * v[1]
                               - ( c1n / c[2] ) * v[2] ) ) / c1n;
        v[0] = i4_max ( v[0], 0 );

        lhs = ( c1n / c[0] ) * v[0]
            + ( c1n / c[1] ) * v[1]
            + ( c1n / c[2] ) * v[2];

        if ( lhs <= rhs1 )
        {
          v[0] = v[0] + 1;
          lhs = lhs + c1n / c[0];
        }

        if ( lhs <= rhs2 )
        {
        }
        else
        {
          *more = false;
          v[0] = 0;
          v[1] = 0;
          v[2] = 0;
        }
      }
    }
  }
  return;
}
//****************************************************************************80

void tetrahedron_lattice_point_next ( int c[], int v[], bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_lattice_point_next() returns the next tetrahedron lattice point.
//
//  Discussion:
//
//    The lattice tetrahedron is defined by the vertices:
//
//      (0,0,0), (C[3]/C[0],0,0), (0,C[3]/C[1],0) and (0,0,C[3]/C[2])
//
//    The lattice tetrahedron is bounded by the lines
//
//      0 <= X,
//      0 <= Y
//      0 <= Z,
//      X / C[0] + Y / C[1] + Z / C[2] <= C[3]
//
//    Lattice points are listed one at a time, starting at the origin,
//    with X increasing first.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int C[4], coefficients defining the
//    lattice tetrahedron.  These should be positive.
//
//    Input/output, int V[3].  On first call, the input
//    value is not important.  On a repeated call, the input value should
//    be the output value from the previous call.  On output, V contains
//    the next lattice point.
//
//    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
//    that this is the first call for a given tetrahedron.  Thereafter, the input
//    value should be the output value from the previous call.  On output,
//    MORE is TRUE if not only is the returned value V a lattice point,
//    but the routine can be called again for another lattice point.
//    If the output value is FALSE, then no more lattice points were found,
//    and V was reset to 0, and the routine should not be called further
//    for this tetrahedron.
//
{
  int c1n;
  int lhs;
  int n = 3;
  int rhs;

  if ( !(*more) )
  {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
    *more = true;
  }
  else
  {
    c1n = i4vec_lcm ( n, c );

    rhs = c1n * c[n];

    lhs =        c[1] * c[2] * v[0]
        + c[0] *        c[2] * v[1]
        + c[0] * c[1]        * v[2];

    if ( lhs + c1n / c[0] <= rhs )
    {
      v[0] = v[0] + 1;
    }
    else
    {
      lhs = lhs - c1n * v[0] / c[0];
      v[0] = 0;
      if ( lhs + c1n / c[1] <= rhs )
      {
        v[1] = v[1] + 1;
      }
      else
      {
        lhs = lhs - c1n * v[1] / c[1];
        v[1] = 0;
        if ( lhs + c1n / c[2] <= rhs )
        {
          v[2] = v[2] + 1;
        }
        else
        {
          v[2] = 0;
          *more = false;
        }
      }
    }
  }
  return;
}
//****************************************************************************80

double tetrahedron_quality1 ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_quality1(): "quality" of a tetrahedron.
//
//  Discussion:
//
//    The quality of a tetrahedron is 3.0 times the ratio of the radius of
//    the inscribed sphere divided by that of the circumscribed sphere.
//
//    An equilateral tetrahredron achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double TETRAHEDRON_QUALITY1, the quality of the tetrahedron.
//
{
  double pc[3];
  double quality;
  double r_in;
  double r_out;

  tetrahedron_circumsphere ( tetra, r_out, pc );

  tetrahedron_insphere ( tetra, r_in, pc );

  quality = 3.0 * r_in / r_out;

  return quality;
}
//****************************************************************************80

double tetrahedron_quality2 ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_quality2(): "quality" of a tetrahedron.
//
//  Discussion:
//
//    The quality measure #2 of a tetrahedron is:
//
//      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
//
//    where
//
//      RIN = radius of the inscribed sphere;
//      LMAX = length of longest side of the tetrahedron.
//
//    An equilateral tetrahredron achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Desheng Wang,
//    The Optimal Centroidal Voronoi Tesselations and the Gersho's
//    Conjecture in the Three-Dimensional Space,
//    Computers and Mathematics with Applications,
//    Volume 49, 2005, pages 1355-1373.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double TETRAHEDRON_QUALITY2, the quality of the tetrahedron.
//
{
  double *edge_length;
  double l_max;
  double pc[3];
  double quality2;
  double r_in;

  edge_length = tetrahedron_edge_length ( tetra );

  l_max = r8vec_max ( 6, edge_length );

  tetrahedron_insphere ( tetra, r_in, pc );

  quality2 = 2.0 * sqrt ( 6.0 ) * r_in / l_max;

  delete [] edge_length;

  return quality2;
}
//****************************************************************************80

double tetrahedron_quality3 ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_quality3() computes the mean ratio of a tetrahedron.
//
//  Discussion:
//
//    This routine computes QUALITY3, the eigenvalue or mean ratio of
//    a tetrahedron.
//
//      QUALITY3 = 12 * ( 3 * volume )^(2/3) / (sum of square of edge lengths).
//
//    This value may be used as a shape quality measure for the tetrahedron.
//
//    For an equilateral tetrahedron, the value of this quality measure
//    will be 1.  For any other tetrahedron, the value will be between
//    0 and 1.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    This version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_QUALITY3, the mean ratio of the tetrahedron.
//
{
  double ab[3];
  double ac[3];
  double ad[3];
  double bc[3];
  double bd[3];
  double cd[3];
  double denom;
  double lab;
  double lac;
  double lad;
  double lbc;
  double lbd;
  double lcd;
  double quality3;
  double volume;
//
//  Compute the vectors representing the sides of the tetrahedron.
//
  tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd );
//
//  Compute the squares of the lengths of the sides.
//
  lab = pow ( ab[0], 2 ) + pow ( ab[1], 2 ) + pow ( ab[2], 2 );
  lac = pow ( ac[0], 2 ) + pow ( ac[1], 2 ) + pow ( ac[2], 2 );
  lad = pow ( ad[0], 2 ) + pow ( ad[1], 2 ) + pow ( ad[2], 2 );
  lbc = pow ( bc[0], 2 ) + pow ( bc[1], 2 ) + pow ( bc[2], 2 );
  lbd = pow ( bd[0], 2 ) + pow ( bd[1], 2 ) + pow ( bd[2], 2 );
  lcd = pow ( cd[0], 2 ) + pow ( cd[1], 2 ) + pow ( cd[2], 2 );
//
//  Compute the volume.
//
  volume = fabs (
      ab[0] * ( ac[1] * ad[2] - ac[2] * ad[1] )
    + ab[1] * ( ac[2] * ad[0] - ac[0] * ad[2] )
    + ab[2] * ( ac[0] * ad[1] - ac[1] * ad[0] ) ) / 6.0;

  denom = lab + lac + lad + lbc + lbd + lcd;

  if ( denom == 0.0 )
  {
    quality3 = 0.0;
  }
  else
  {
    quality3 = 12.0 * pow ( 3.0 * volume, 2.0 / 3.0 ) / denom;
  }

  return quality3;
}
//****************************************************************************80

double tetrahedron_quality4 ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_quality4() computes the minimum solid angle of a tetrahedron.
//
//  Discussion:
//
//    This routine computes a quality measure for a tetrahedron, based
//    on the sine of half the minimum of the four solid angles.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    This version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double QUALITY4, the value of the quality measure.
//
{
  double ab[3];
  double ac[3];
  double ad[3];
  double bc[3];
  double bd[3];
  double cd[3];
  double denom;
  double l1;
  double l2;
  double l3;
  double lab;
  double lac;
  double lad;
  double lbc;
  double lbd;
  double lcd;
  double quality4;
  double volume;
//
//  Compute the edges.
//
  tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd );
//
//  Compute the lengths of the sides.
//
  lab = r8vec_length ( 3, ab );
  lac = r8vec_length ( 3, ac );
  lad = r8vec_length ( 3, ad );
  lbc = r8vec_length ( 3, bc );
  lbd = r8vec_length ( 3, bd );
  lcd = r8vec_length ( 3, cd );
//
//  Compute the volume.
//
  volume = fabs (
      ab[0] * ( ac[1] * ad[2] - ac[2] * ad[1] )
    + ab[1] * ( ac[2] * ad[0] - ac[0] * ad[2] )
    + ab[2] * ( ac[0] * ad[1] - ac[1] * ad[0] ) ) / 6.0;

  quality4 = 1.0;

  l1 = lab + lac;
  l2 = lab + lad;
  l3 = lac + lad;

  denom = ( l1 + lbc ) * ( l1 - lbc )
        * ( l2 + lbd ) * ( l2 - lbd )
        * ( l3 + lcd ) * ( l3 - lcd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = fmin ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lab + lbc;
  l2 = lab + lbd;
  l3 = lbc + lbd;

  denom = ( l1 + lac ) * ( l1 - lac )
        * ( l2 + lad ) * ( l2 - lad )
        * ( l3 + lcd ) * ( l3 - lcd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = fmin ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lac + lbc;
  l2 = lac + lcd;
  l3 = lbc + lcd;

  denom = ( l1 + lab ) * ( l1 - lab )
        * ( l2 + lad ) * ( l2 - lad )
        * ( l3 + lbd ) * ( l3 - lbd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = fmin ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lad + lbd;
  l2 = lad + lcd;
  l3 = lbd + lcd;

  denom = ( l1 + lab ) * ( l1 - lab )
        * ( l2 + lac ) * ( l2 - lac )
        * ( l3 + lbc ) * ( l3 - lbc );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = fmin ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  quality4 = quality4 * 1.5 * sqrt ( 6.0 );

  return quality4;
}
//****************************************************************************80

void tetrahedron_rhombic_shape ( int point_num, int face_num,
  int face_order_max, double point_coord[], int face_order[],
  int face_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_rhombic_shape() describes a rhombic tetrahedron.
//
//  Discussion:
//
//    Call tetrahedron_rhombic_size() first, to get dimension information.
//
//    The tetrahedron is described using 10 nodes.  If we label the vertices
//    P0, P1, P2 and P3, then the extra nodes lie halfway between vertices,
//    and have the labels P01, P02, P03, P12, P13 and P23.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    17 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anwei Liu, Barry Joe,
//    Quality Local Refinement of Tetrahedral Meshes Based
//    on 8-Subtetrahedron Subdivision,
//    Mathematics of Computation,
//    Volume 65, Number 215, July 1996, pages 1183-1200.
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_ORDER_MAX, the maximum number of vertices per face.
//
//    Output, double POINT_COORD[3*POINT_NUM], the vertices.
//
//    Output, int FACE_ORDER[FACE_NUM], the number of vertices
//    for each face.
//
//    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
//    contains the index of the I-th point in the J-th face.  The
//    points are listed in the counter clockwise direction defined
//    by the outward normal at the face.
//
{
  double a;
  double b;
  double c;
  double d;
  int face;
  double z = 0.0;

  a =        1.0   / sqrt ( 3.0 );
  b = sqrt ( 2.0 ) / sqrt ( 3.0 );
  c = sqrt ( 3.0 ) /        6.0;
  d =        1.0   / sqrt ( 6.0 );
//
//  Set the point coordinates.
//
  point_coord[0+0*3] = -b;
  point_coord[1+0*3] =  z;
  point_coord[2+0*3] =  z;

  point_coord[0+1*3] =  z;
  point_coord[1+1*3] = -a;
  point_coord[2+1*3] =  z;

  point_coord[0+2*3] =  z;
  point_coord[1+2*3] =  a;
  point_coord[2+2*3] =  z;

  point_coord[0+3*3] =  z;
  point_coord[1+3*3] =  z;
  point_coord[2+3*3] =  b;

  point_coord[0+4*3] = -d;
  point_coord[1+4*3] = -c;
  point_coord[2+4*3] =  z;

  point_coord[0+5*3] = -d;
  point_coord[1+5*3] =  c;
  point_coord[2+5*3] =  z;

  point_coord[0+6*3] = -d;
  point_coord[1+6*3] =  z;
  point_coord[2+6*3] =  d;

  point_coord[0+7*3] =  z;
  point_coord[1+7*3] =  z;
  point_coord[2+7*3] =  z;

  point_coord[0+8*3] =  z;
  point_coord[1+8*3] = -c;
  point_coord[2+8*3] =  d;

  point_coord[0+9*3] =  z;
  point_coord[1+9*3] =  c;
  point_coord[2+9*3] =  d;
//
//  Set the face orders.
//
  for ( face = 0; face < face_num; face++ )
  {
    face_order[face] = 6;
  }
//
//  Set faces.
//
  face_point[0+0*face_order_max] =  1;
  face_point[1+0*face_order_max] =  5;
  face_point[2+0*face_order_max] =  2;
  face_point[3+0*face_order_max] =  9;
  face_point[4+0*face_order_max] =  4;
  face_point[5+0*face_order_max] =  7;

  face_point[0+1*face_order_max] =  2;
  face_point[1+1*face_order_max] =  8;
  face_point[2+1*face_order_max] =  3;
  face_point[3+1*face_order_max] = 10;
  face_point[4+1*face_order_max] =  4;
  face_point[5+1*face_order_max] =  9;

  face_point[0+2*face_order_max] =  3;
  face_point[1+2*face_order_max] =  6;
  face_point[2+2*face_order_max] =  1;
  face_point[3+2*face_order_max] =  7;
  face_point[4+2*face_order_max] =  4;
  face_point[5+2*face_order_max] = 10;

  face_point[0+3*face_order_max] =  1;
  face_point[1+3*face_order_max] =  6;
  face_point[2+3*face_order_max] =  3;
  face_point[3+3*face_order_max] =  8;
  face_point[4+3*face_order_max] =  2;
  face_point[5+3*face_order_max] =  5;

  return;
}
//****************************************************************************80

void tetrahedron_rhombic_size ( int *point_num, int *edge_num,
  int *face_num, int *face_order_max )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_rhombic_size() gives "sizes" for a rhombic tetrahedron.
//
//  Discussion:
//
//    Call this routine first, in order to learn the required dimensions
//    of arrays to be set up by tetrahedron_rhombic_shape().
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *POINT_NUM, the number of vertices.
//
//    Output, int *EDGE_NUM, the number of edges.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *FACE_ORDER_MAX, the maximum order of any face.
//
{
  *point_num = 10;
  *edge_num = 6;
  *face_num = 4;
  *face_order_max = 6;

  return;
}
//****************************************************************************80

void tetrahedron_sample ( double tetra[3*4], int n, int &seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_sample() returns random points in a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double P[3*N], random points in the tetrahedron.
//
{
# define DIM_NUM 3

  double alpha;
  double beta;
  double gamma;
  int i;
  int j;
  int k;
  double *p12;
  double *p13;
  double r;
  double *t;

  p12 = new double[DIM_NUM];
  p13 = new double[DIM_NUM];
  t = new double[DIM_NUM*3];

  for ( k = 0; k < n; k++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the tetrahedron's volume.
//
//  Imagine a plane, parallel to face 1, so that the volume between
//  vertex 1 and the plane is R percent of the full tetrahedron volume.
//
//  The plane will intersect sides 12, 13, and 14 at a fraction
//  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
//
    alpha = pow ( r, 1.0 / 3.0 );
//
//  Determine the coordinates of the points on sides 12, 13 and 14 intersected
//  by the plane, which form a triangle TR.
//
    for ( i = 0; i < DIM_NUM; i++ )
    {
      for ( j = 0; j < 3; j++ )
      {
        t[i+j*3] = ( 1.0 - alpha ) * tetra[i+0*3]
                 +         alpha   * tetra[i+(j+1)*3];
      }
    }
//
//  Now choose, uniformly at random, a point in this triangle.
//
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    beta = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
    for ( i = 0; i < DIM_NUM; i++ )
    {
      p12[i] = ( 1.0 - beta ) * t[i+0*3]
             +         beta   * t[i+1*3];

      p13[i] = ( 1.0 - beta ) * t[i+0*3]
             +         beta   * t[i+2*3];
    }
//
//  Now choose, uniformly at random, a point on the line L.
//
    gamma = r8_uniform_01 ( seed );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      p[i+k*3] = gamma * p12[i] + ( 1.0 - gamma ) * p13[i];
    }
  }

  delete [] p12;
  delete [] p13;
  delete [] t;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_shape ( int point_num, int face_num, int face_order_max,
  double point_coord[], int face_order[], int face_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_shape() describes a tetrahedron.
//
//  Discussion:
//
//    The vertices lie on the unit sphere.
//
//    The dual of the tetrahedron is the tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    07 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_ORDER_MAX, the maximum number of vertices per face.
//
//    Output, double POINT_COORD[3*POINT_NUM], the point coordinates.
//
//    Output, int FACE_ORDER[FACE_NUM], the number of vertices
//    for each face.
//
//    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
//    contains the index of the I-th point in the J-th face.  The
//    points are listed in the counter clockwise direction defined
//    by the outward normal at the face.
//
{
# define DIM_NUM 3

  static int face_order_save[4] = {
    3, 3, 3, 3 };
  static int face_point_save[3*4] = {
       1, 3, 2,
       1, 2, 4,
       1, 4, 3,
       2, 3, 4 };
  static double point_coord_save[3*4] = {
        0.942809,    0.000000,   -0.333333,
       -0.471405,    0.816497,   -0.333333,
       -0.471405,   -0.816497,   -0.333333,
        0.000000,    0.000000,    1.000000 };

  i4vec_copy ( face_num, face_order_save, face_order );
  i4vec_copy ( face_order_max*face_num, face_point_save, face_point );
  r8vec_copy ( DIM_NUM*point_num, point_coord_save, point_coord );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_size ( int *point_num, int *edge_num, int *face_num,
  int *face_order_max )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_size() gives "sizes" for a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *POINT_NUM, the number of points.
//
//    Output, int *EDGE_NUM, the number of edges.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *FACE_ORDER_MAX, the maximum order of any face.
//
{
  *point_num = 4;
  *edge_num = 12;
  *face_num = 4;
  *face_order_max = 3;

  return;
}
//****************************************************************************80

double *tetrahedron_solid_angles ( double tetra[] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_solid_angles() computes solid angles of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    07 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_SOLID_ANGLES[4], the solid angles.
//
{
  double *angle;
  double *dihedral_angles;
  const double r8_pi = 3.141592653589793;

  dihedral_angles = tetrahedron_dihedral_angles ( tetra );

  angle = new double[4];

  angle[0] = dihedral_angles[0]
           + dihedral_angles[1]
           + dihedral_angles[2] - r8_pi;

  angle[1] = dihedral_angles[0]
           + dihedral_angles[3]
           + dihedral_angles[4] - r8_pi;

  angle[2] = dihedral_angles[1]
           + dihedral_angles[3]
           + dihedral_angles[5] - r8_pi;

  angle[3] = dihedral_angles[2]
           + dihedral_angles[4]
           + dihedral_angles[5] - r8_pi;

  delete [] dihedral_angles;

  return angle;
}
//****************************************************************************80

int tetrahedron_unit_lattice_point_num ( int s )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_unit_lattice_point_num(): count lattice points.
//
//  Discussion:
//
//    The tetrahedron is assumed to be the unit tetrahedron:
//
//    ( (0,0,0), (1,0,0), (0,1,0), (0,0,1) )
//
//    or a copy of this tetrahedron scaled by an integer S:
//
//    ( (0,0,0), (S,0,0), (0,S,0), (0,0,S) ).
//
//    The routine returns the number of integer lattice points that appear
//    inside the tetrahedron, or on its faces, edges or vertices.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    03 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Matthias Beck, Sinai Robins,
//    Computing the Continuous Discretely,
//    Springer, 2006,
//    ISBN13: 978-0387291390,
//    LC: QA640.7.B43.
//
//  Parameters:
//
//    Input, int S, the scale factor.
//
//    Output, int TETRAHEDRON_UNIT_LATTICE_POINT_NUM, the number of lattice points.
//
{
  int n;

  n = ( ( s + 3 ) * ( s + 2 ) * ( s + 1 ) ) / 6;

  return n;
}
//****************************************************************************80

double tetrahedron_volume ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    tetrahedron_volume() computes the volume of a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
//
{
  double a[4*4];
  int i;
  int j;
  double volume;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      a[i+j*4] = tetra[i+j*3];
    }
  }

  i = 3;
  for ( j = 0; j < 4; j++ )
  {
    a[i+j*4] = 1.0;
  }

  volume = fabs ( r8mat_det_4d ( a ) ) / 6.0;

  return volume;
}
//****************************************************************************80

void triangle_angles_3d ( double t[3*3], double angle[3] )

//****************************************************************************80
//
//  Purpose:
//
//    triangle_angles_3d() computes the angles of a triangle in 3D.
//
//  Discussion:
//
//    The law of cosines is used:
//
//      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
//
//    where GAMMA is the angle opposite side C.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[3*3], the triangle vertices.
//
//    Output, double ANGLE[3], the angles opposite
//    sides P1-P2, P2-P3 and P3-P1, in radians.
//
{
  double a;
  double b;
  double c;
  const double r8_pi = 3.141592653589793;

  a = sqrt ( pow ( t[0+1*3] - t[0+0*3], 2 )
           + pow ( t[1+1*3] - t[1+0*3], 2 )
           + pow ( t[2+1*3] - t[2+0*3], 2 ) );

  b = sqrt ( pow ( t[0+2*3] - t[0+1*3], 2 )
           + pow ( t[1+2*3] - t[1+1*3], 2 )
           + pow ( t[2+2*3] - t[2+1*3], 2 ) );

  c = sqrt ( pow ( t[0+0*3] - t[0+2*3], 2 )
           + pow ( t[1+0*3] - t[1+2*3], 2 )
           + pow ( t[2+0*3] - t[2+2*3], 2 ) );
//
//  Take care of a ridiculous special case.
//
  if ( a == 0.0 && b == 0.0 && c == 0.0 )
  {
    angle[0] = 2.0 * r8_pi / 3.0;
    angle[1] = 2.0 * r8_pi / 3.0;
    angle[2] = 2.0 * r8_pi / 3.0;
    return;
  }

  if ( c == 0.0 || a == 0.0 )
  {
    angle[0] = r8_pi;
  }
  else
  {
    angle[0] = r8_acos ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) );
  }

  if ( a == 0.0 || b == 0.0 )
  {
    angle[1] = r8_pi;
  }
  else
  {
    angle[1] = r8_acos ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) );
  }

  if ( b == 0.0 || c == 0.0 )
  {
    angle[2] = r8_pi;
  }
  else
  {
    angle[2] = r8_acos ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
  }

  return;
}
//****************************************************************************80

double triangle_area_3d ( double t[3*3] )

//****************************************************************************80
//
//  Purpose:
//
//    triangle_area_3d() computes the area of a triangle in 3D.
//
//  Discussion:
//
//    This routine uses the fact that the norm of the cross product vector
//    is the area of the parallelogram they form.  The triangle they
//    form has half that area.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[3*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_3D, the area of the triangle.
//
{
  double area;
  double *cross;
  int i;
//
//  Compute the cross product vector.
//
  cross = new double[3];

  cross[0] = ( t[1+1*3] - t[1+0*3] )
           * ( t[2+2*3] - t[2+0*3] )
           - ( t[2+1*3] - t[2+0*3] )
           * ( t[1+2*3] - t[1+0*3] );

  cross[1] = ( t[2+1*3] - t[2+0*3] )
           * ( t[0+2*3] - t[0+0*3] )
           - ( t[0+1*3] - t[0+0*3] )
           * ( t[2+2*3] - t[2+0*3] );

  cross[2] = ( t[0+1*3] - t[0+0*3] )
           * ( t[1+2*3] - t[1+0*3] )
           - ( t[1+1*3] - t[1+0*3] )
           * ( t[0+2*3] - t[0+0*3] );

  area = 0.0;
  for ( i = 0; i < 3; i++ )
  {
    area = area + pow ( cross[i], 2 );
  }

  area = 0.5 * sqrt ( area );

  delete [] cross;

  return area;
}

}
