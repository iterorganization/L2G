#include "bicubic.hpp"
#include <cmath>
#include <cstdlib>
#include<rkf45.hpp>

RKF45::RKF45(){
    m_abserr_save = -1.0;
    m_flag_save = -1000;
    m_h = -1.0;
    m_init = -1000;
    m_kflag = -1000;
    m_kop = -1;
    m_nfe = -1;
    m_relerr_save = -1.0;
    m_remin = 1.0E-12;
    m_r8epsilon = 2.220446049250313E-016;
    return;
}

double RKF45::r8_max(double x, double y)
//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if (y < x)
  {
    return x;
  }
  else
  {
    return y;
  }
}

double RKF45::r8_min (double x, double y)
//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if (y < x)
  {
    return y;
  }
  else
  {
    return x;
  }
}

double RKF45::r8_abs(double x)
//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  if (0.0 <= x)
  {
    return x;
  }
  else
  {
    return (-x);
  }
}

double RKF45::r8_epsilon()

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{

  return m_r8epsilon;
}

double RKF45::r8_sign (double x)
//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  if ( x < 0.0 )
  {
    return ( -1.0 );
  }
  else
  {
    return ( +1.0 );
  }
}

void RKF45::r8_flt(double y[2], double yp[2]){
    // Derivative function for the RKF45. It evaluates the ratio between the
    // (R, Z) magnetic components and the toroidal magnetic component in the
    // cylindrical coordinate system.
    double factor;

    BI_DATA *context = new BI_DATA;
    m_interp_psi->populateContext(context);
    context->r = y[0] - m_r_move;
    context->z = y[1] - m_z_move;
    m_interp_psi->getFirstDerivativeValues(context);

    factor = y[0] / m_vacuum_fpol;
    yp[0] = - context->valdy * factor;
    yp[1] =   context->valdx * factor;
    delete context;
}

void RKF45::r8_fehl(double y[2], double h, double yp[2],
                    double f1[2], double f2[2], double f3[2], double f4[2],
                    double f5[2], double s[2])
//****************************************************************************80
//
//  Purpose:
//
//    R8_FEHL takes one Fehlberg fourth-fifth order step.
//
//  Discussion:
//
//    This version of the routine uses DOUBLE real arithemtic.
//
//    This routine integrates a system of 3 first order ordinary differential
//    equations of the form
//      dY(i)/dT = F(T,Y(1:3))
//    where the initial values Y and the initial derivatives
//    YP are specified at the starting point T.
//
//    The routine advances the solution over the fixed step H and returns
//    the fifth order (sixth order accurate locally) solution
//    approximation at T+H in array S.
//
//    The formulas have been grouped to control loss of significance.
//    The routine should be called with an H not smaller than 13 units of
//    roundoff in T so that the various independent arguments can be
//    distinguished.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Erwin Fehlberg,
//    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
//    NASA Technical Report R-315, 1969.
//
//    Lawrence Shampine, Herman Watts, S Davenport,
//    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
//    SIAM Review,
//    Volume 18, pages 376-411, 1976.
//
//  Parameters:
//
//    Input, external F, a user-supplied subroutine to evaluate the
//    derivatives Y'(T), of the form:
//
//      void f ( double t, double y[], double yp[] )
//
//    Input, double Y[2], the current value of the dependent variable.
//
//    Input, double T, the current value of the independent variable.
//
//    Input, double H, the step size to take.
//
//    Input, double YP[2], the current value of the derivative of the
//    dependent variable.
//
//    Output, double F1[2], F2[2], F3[2], F4[2], F5[2], derivative
//    values needed for the computation.
//
//    Output, double S[2], the estimate of the solution at T+H.
//
{
  double ch;
  int i;

  ch = h / 4.0;

  for ( i = 0; i < 2; i++ )
  {
    f5[i] = y[i] + ch * yp[i];
  }

  r8_flt ( f5, f1);

  ch = 3.0 * h / 32.0;

  for ( i = 0; i < 2; i++ )
  {
    f5[i] = y[i] + ch * ( yp[i] + 3.0 * f1[i] );
  }

  r8_flt ( f5, f2);

  ch = h / 2197.0;

  for ( i = 0; i < 2; i++ )
  {
    f5[i] = y[i] + ch *
    ( 1932.0 * yp[i]
    + ( 7296.0 * f2[i] - 7200.0 * f1[i] )
    );
  }

  r8_flt ( f5, f3);

  ch = h / 4104.0;

  for ( i = 0; i < 2; i++ )
  {
    f5[i] = y[i] + ch *
    (
      ( 8341.0 * yp[i] - 845.0 * f3[i] )
    + ( 29440.0 * f2[i] - 32832.0 * f1[i] )
    );
  }

  r8_flt ( f5, f4);

  ch = h / 20520.0;

  for ( i = 0; i < 2; i++ )
  {
    f1[i] = y[i] + ch *
    (
      ( -6080.0 * yp[i]
      + ( 9295.0 * f3[i] - 5643.0 * f4[i] )
      )
    + ( 41040.0 * f1[i] - 28352.0 * f2[i] )
    );
  }

  r8_flt ( f1, f5 );
//
//  Ready to compute the approximate solution at T+H.
//
  ch = h / 7618050.0;

  for ( i = 0; i < 2; i++ )
  {
    s[i] = y[i] + ch *
    (
      ( 902880.0 * yp[i]
      + ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) )
    + ( 3953664.0 * f2[i] + 277020.0 * f5[i] )
    );
  }

  return;
}

int RKF45::r8_rkf45(double y[2], double yp[2], double *t,
                    double tout, double *relerr, double abserr, int flag)

//****************************************************************************80
//
//  Purpose:
//
//    R8_RKF45 carries out the Runge-Kutta-Fehlberg method.
//
//  Discussion:
//
//    This version of the routine uses DOUBLE real arithmetic.
//
//    This routine is primarily designed to solve non-stiff and mildly stiff
//    differential equations when derivative evaluations are inexpensive.
//    It should generally not be used when the user is demanding
//    high accuracy.
//
//    This routine integrates a system of NEQN first-order ordinary differential
//    equations of the form:
//
//      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
//
//    where the Y(1:NEQN) are given at T.
//
//    Typically the subroutine is used to integrate from T to TOUT but it
//    can be used as a one-step integrator to advance the solution a
//    single step in the direction of TOUT.  On return, the parameters in
//    the call list are set for continuing the integration.  The user has
//    only to call again (and perhaps define a new value for TOUT).
//
//    Before the first call, the user must
//
//    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
//      and declare F in an EXTERNAL statement;
//
//    * initialize the parameters:
//      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
//      In particular, T should initially be the starting point for integration,
//      Y should be the value of the initial conditions, and FLAG should
//      normally be +1.
//
//    Normally, the user only sets the value of FLAG before the first call, and
//    thereafter, the program manages the value.  On the first call, FLAG should
//    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
//    have been reset by the program to the value of 2 (or -2 in single
//    step mode), and the user can continue to call the routine with that
//    value of FLAG.
//
//    (When the input magnitude of FLAG is 1, this indicates to the program
//    that it is necessary to do some initialization work.  An input magnitude
//    of 2 lets the program know that that initialization can be skipped,
//    and that useful information was computed earlier.)
//
//    The routine returns with all the information needed to continue
//    the integration.  If the integration reached TOUT, the user need only
//    define a new TOUT and call again.  In the one-step integrator
//    mode, returning with FLAG = -2, the user must keep in mind that
//    each step taken is in the direction of the current TOUT.  Upon
//    reaching TOUT, indicated by the output value of FLAG switching to 2,
//    the user must define a new TOUT and reset FLAG to -2 to continue
//    in the one-step integrator mode.
//
//    In some cases, an error or difficulty occurs during a call.  In that case,
//    the output value of FLAG is used to indicate that there is a problem
//    that the user must address.  These values include:
//
//    * 3, integration was not completed because the input value of RELERR, the
//      relative error tolerance, was too small.  RELERR has been increased
//      appropriately for continuing.  If the user accepts the output value of
//      RELERR, then simply reset FLAG to 2 and continue.
//
//    * 4, integration was not completed because more than MAXNFE derivative
//      evaluations were needed.  This is approximately (MAXNFE/6) steps.
//      The user may continue by simply calling again.  The function counter
//      will be reset to 0, and another MAXNFE function evaluations are allowed.
//
//    * 5, integration was not completed because the solution vanished,
//      making a pure relative error test impossible.  The user must use
//      a non-zero ABSERR to continue.  Using the one-step integration mode
//      for one step is a good way to proceed.
//
//    * 6, integration was not completed because the requested accuracy
//      could not be achieved, even using the smallest allowable stepsize.
//      The user must increase the error tolerances ABSERR or RELERR before
//      continuing.  It is also necessary to reset FLAG to 2 (or -2 when
//      the one-step integration mode is being used).  The occurrence of
//      FLAG = 6 indicates a trouble spot.  The solution is changing
//      rapidly, or a singularity may be present.  It often is inadvisable
//      to continue.
//
//    * 7, it is likely that this routine is inefficient for solving
//      this problem.  Too much output is restricting the natural stepsize
//      choice.  The user should use the one-step integration mode with
//      the stepsize determined by the code.  If the user insists upon
//      continuing the integration, reset FLAG to 2 before calling
//      again.  Otherwise, execution will be terminated.
//
//    * 8, invalid input parameters, indicates one of the following:
//      NEQN <= 0;
//      T = TOUT and |FLAG| /= 1;
//      RELERR < 0 or ABSERR < 0;
//      FLAG == 0  or FLAG < -2 or 8 < FLAG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2012
//
//  Author:
//
//    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Erwin Fehlberg,
//    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
//    NASA Technical Report R-315, 1969.
//
//    Lawrence Shampine, Herman Watts, S Davenport,
//    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
//    SIAM Review,
//    Volume 18, pages 376-411, 1976.
//
//  Parameters:
//
//    Input, external F, a user-supplied subroutine to evaluate the
//    derivatives Y'(T), of the form:
//
//      void f ( double t, double y[], double yp[] )
//
//    Input, int NEQN, the number of equations to be integrated.
//
//    Input/output, double Y[NEQN], the current solution vector at T.
//
//    Input/output, double YP[NEQN], the derivative of the current solution
//    vector at T.  The user should not set or alter this information!
//
//    Input/output, double *T, the current value of the independent variable.
//
//    Input, double TOUT, the output point at which solution is desired.
//    TOUT = T is allowed on the first call only, in which case the routine
//    returns with FLAG = 2 if continuation is possible.
//
//    Input, double *RELERR, ABSERR, the relative and absolute error tolerances
//    for the local error test.  At each step the code requires:
//      abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
//    for each component of the local error and the solution vector Y.
//    RELERR cannot be "too small".  If the routine believes RELERR has been
//    set too small, it will reset RELERR to an acceptable value and return
//    immediately for user action.
//
//    Input, int FLAG, indicator for status of integration. On the first call,
//    set FLAG to +1 for normal use, or to -1 for single step mode.  On
//    subsequent continuation steps, FLAG should be +2, or -2 for single
//    step mode.
//
//    Output, int RKF45_D, indicator for status of integration.  A value of 2
//    or -2 indicates normal progress, while any other value indicates a
//    problem that should be addressed.
//
{
# define MAXNFE 3000


  double ae;
  double dt;
  double ee;
  double eeoet;
  double eps;
  double esttol;
  double et;
  double *f1;
  double *f2;
  double *f3;
  double *f4;
  double *f5;
  int flag_return;


  bool hfaild;
  double hmin;
  int i;

  int k;


  int mflag;

  bool output;
  double relerr_min;


  double s;
  double scale;
  double tol;
  double toln;
  double ypk;

  flag_return = flag;
//
//  Check the input parameters.
//
  eps = r8_epsilon ( );

  // if ( neqn < 1 )
  // {
  //   flag_return = 8;
    // std::cerr << "\n";
    // std::cerr << "R8_RKF45 - Fatal error!\n";
    // std::cerr << "  Invalid input value of NEQN.\n";
  //   return flag_return;
  // }

  if ( (*relerr) < 0.0 )
  {
    flag_return = 8;
    // std::cerr << "\n";
    // std::cerr << "R8_RKF45 - Fatal error!\n";
    // std::cerr << "  Invalid input value of RELERR.\n";
    return flag_return;
  }

  if ( abserr < 0.0 )
  {
    flag_return = 8;
    // std::cerr << "\n";
    // std::cerr << "R8_RKF45 - Fatal error!\n";
    // std::cerr << "  Invalid input value of ABSERR.\n";
    return flag_return;
  }

  if ( flag_return == 0 || 8 < flag_return  || flag_return < -2 )
  {
    flag_return = 8;
    // std::cerr << "\n";
    // std::cerr << "R8_RKF45 - Fatal error!\n";
    // std::cerr << "  Invalid input.\n";
    return flag_return;
  }

  mflag = abs ( flag_return );
//
//  Is this a continuation call?
//
  if ( mflag != 1 )
  {
    if ( *t == tout && m_kflag != 3 )
    {
      flag_return = 8;
      return flag_return;
    }
//
//  FLAG = -2 or +2:
//
    if ( mflag == 2 )
    {
      if ( m_kflag == 3 )
      {
        flag_return = m_flag_save;
        mflag = abs ( flag_return );
      }
      else if ( m_init == 0 )
      {
        flag_return = m_flag_save;
      }
      else if ( m_kflag == 4 )
      {
        m_nfe = 0;
      }
      else if ( m_kflag == 5 && abserr == 0.0 )
      {
        // std::cerr << "\n";
        // std::cerr << "R8_RKF45 - Fatal error!\n";
        // std::cerr << "  KFLAG = 5 and ABSERR = 0.0\n";
        exit ( 1 );
      }
      else if ( m_kflag == 6 && (*relerr) <= m_relerr_save && abserr <= m_abserr_save )
      {
        // std::cerr << "\n";
        // std::cerr << "R8_RKF45 - Fatal error!\n";
        // std::cerr << "  KFLAG = 6 and\n";
        // std::cerr << "  RELERR <= RELERR_SAVE and\n";
        // std::cerr << "  " << (*relerr) << " <= " << m_relerr_save << "\n";
        // std::cerr << "  ABSERR <= ABSERR_SAVE\n";
        // std::cerr << "  " << abserr << " <= " << m_abserr_save << "\n";
        exit ( 2 );
      }
    }
//
//  FLAG = 3, 4, 5, 6, 7 or 8.
//
    else
    {
      if ( flag_return == 3 )
      {
        flag_return = m_flag_save;
        if ( m_kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
      else if ( flag_return == 4 )
      {
        m_nfe = 0;
        flag_return = m_flag_save;
        if ( m_kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
      else if ( flag_return == 5 && 0.0 < abserr )
      {
        flag_return = m_flag_save;
        if ( m_kflag == 3 )
        {
          mflag = abs ( flag_return );
        }
      }
//
//  Integration cannot be continued because the user did not respond to
//  the instructions pertaining to FLAG = 5, 6, 7 or 8.
//
      else
      {
        // std::cerr << "\n";
        // std::cerr << "R8_RKF45 - Fatal error!\n";
        // std::cerr << "  Integration cannot be continued.\n";
        // std::cerr << "  The user did not respond to the output\n";
        // std::cerr << "  value FLAG = 5, 6, 7, or 8.\n";
        exit ( 3 );
        //flag_return = 2;
        //flag_return = 1;
        //flag_return = -2;
        //mflag = abs ( flag_return );
      }
    }
  }
//
//  Save the input value of FLAG.
//  Set the continuation flag KFLAG for subsequent input checking.
//
  m_flag_save = flag_return;
  m_kflag = 0;
//
//  Save RELERR and ABSERR for checking input on subsequent calls.
//
  m_relerr_save = (*relerr);
  m_abserr_save = abserr;
//
//  Restrict the relative error tolerance to be at least
//
//    2*EPS+REMIN
//
//  to avoid limiting precision difficulties arising from impossible
//  accuracy requests.
//
  relerr_min = 2.0 * r8_epsilon ( ) + m_remin;
//
//  Is the relative error tolerance too small?
//
  if ( (*relerr) < relerr_min )
  {
    (*relerr) = relerr_min;
    m_kflag = 3;
    flag_return = 3;
    return flag_return;
  }

  dt = tout - *t;
//
//  Initialization:
//
//  Set the initialization completion indicator, INIT;
//  set the indicator for too many output points, KOP;
//  evaluate the initial derivatives
//  set the counter for function evaluations, NFE;
//  estimate the starting stepsize.
//
  f1 = new double[2];
  f2 = new double[2];
  f3 = new double[2];
  f4 = new double[2];
  f5 = new double[2];

  if ( mflag == 1 )
  {
    m_init = 0;
    m_kop = 0;
    r8_flt ( y, yp);
    m_nfe = 1;

    if ( *t == tout )
    {
      flag_return = 2;
      return flag_return;
    }

  }

  if ( m_init == 0 )
  {
    m_init = 1;
    m_h = r8_abs ( dt );
    toln = 0.0;

    for ( k = 0; k < 2; k++ )
    {
      tol = (*relerr) * r8_abs ( y[k] ) + abserr;
      if ( 0.0 < tol )
      {
        toln = tol;
        ypk = r8_abs ( yp[k] );
        if ( tol < ypk * std::pow ( m_h, 5 ) )
        {
          m_h = std::pow ( ( tol / ypk ), 0.2 );
        }
      }
    }

    if ( toln <= 0.0 )
    {
      m_h = 0.0;
    }

    m_h = r8_max ( m_h, 26.0 * eps * r8_max ( r8_abs ( *t ), r8_abs ( dt ) ) );

    if ( flag_return < 0 )
    {
      m_flag_save = -2;
    }
    else
    {
      m_flag_save = 2;
    }
  }
//
//  Set stepsize for integration in the direction from T to TOUT.
//
  m_h = r8_sign ( dt ) * r8_abs ( m_h );
//
//  Test to see if too may output points are being requested.
//
  if ( 2.0 * r8_abs ( dt ) <= r8_abs ( m_h ) )
  {
    m_kop = m_kop + 1;
  }
//
//  Unnecessary frequency of output.
//
  if ( m_kop == 100 )
  {
    m_kop = 0;
    delete [] f1;
    delete [] f2;
    delete [] f3;
    delete [] f4;
    delete [] f5;
    flag_return = 7;
    return flag_return;
  }
//
//  If we are too close to the output point, then simply extrapolate and return.
//
  if ( r8_abs ( dt ) <= 26.0 * eps * r8_abs ( *t ) )
  {
    *t = tout;
    for ( i = 0; i < 2; i++ )
    {
      y[i] = y[i] + dt * yp[i];
    }
    r8_flt ( y, yp );
    m_nfe = m_nfe + 1;

    delete [] f1;
    delete [] f2;
    delete [] f3;
    delete [] f4;
    delete [] f5;
    flag_return = 2;
    return flag_return;
  }
//
//  Initialize the output point indicator.
//
  output = false;
//
//  To avoid premature underflow in the error tolerance function,
//  scale the error tolerances.
//
  scale = 2.0 / (*relerr);
  ae = scale * abserr;
//
//  Step by step integration.
//
  for ( ; ; )
  {
    hfaild = false;
//
//  Set the smallest allowable stepsize.
//
    hmin = 26.0 * eps * r8_abs ( *t );
//
//  Adjust the stepsize if necessary to hit the output point.
//
//  Look ahead two steps to avoid drastic changes in the stepsize and
//  thus lessen the impact of output points on the code.
//
    dt = tout - *t;

    if ( 2.0 * r8_abs ( m_h ) <= r8_abs ( dt ) )
    {
    }
    else
//
//  Will the next successful step complete the integration to the output point?
//
    {
      if ( r8_abs ( dt ) <= r8_abs ( m_h ) )
      {
        output = true;
        m_h = dt;
      }
      else
      {
        m_h = 0.5 * dt;
      }

    }
//
//  Here begins the core integrator for taking a single step.
//
//  The tolerances have been scaled to avoid premature underflow in
//  computing the error tolerance function ET.
//  To avoid problems with zero crossings, relative error is measured
//  using the average of the magnitudes of the solution at the
//  beginning and end of a step.
//  The error estimate formula has been grouped to control loss of
//  significance.
//
//  To distinguish the various arguments, H is not permitted
//  to become smaller than 26 units of roundoff in T.
//  Practical limits on the change in the stepsize are enforced to
//  smooth the stepsize selection process and to avoid excessive
//  chattering on problems having discontinuities.
//  To prevent unnecessary failures, the code uses 9/10 the stepsize
//  it estimates will succeed.
//
//  After a step failure, the stepsize is not allowed to increase for
//  the next attempted step.  This makes the code more efficient on
//  problems having discontinuities and more effective in general
//  since local extrapolation is being used and extra caution seems
//  warranted.
//
//  Test the number of derivative function evaluations.
//  If okay, try to advance the integration from T to T+H.
//
    for ( ; ; )
    {
//
//  Have we done too much work?
//
      if ( MAXNFE < m_nfe )
      {
        m_kflag = 4;
        delete [] f1;
        delete [] f2;
        delete [] f3;
        delete [] f4;
        delete [] f5;
        flag_return = 4;
        return flag_return;
      }
//
//  Advance an approximate solution over one step of length H.
//
      r8_fehl ( y, m_h, yp, f1, f2, f3, f4, f5, f1 );
      m_nfe = m_nfe + 5;
//
//  Compute and test allowable tolerances versus local error estimates
//  and remove scaling of tolerances.  The relative error is
//  measured with respect to the average of the magnitudes of the
//  solution at the beginning and end of the step.
//
      eeoet = 0.0;

      for ( k = 0; k < 2; k++ )
      {
        et = r8_abs ( y[k] ) + r8_abs ( f1[k] ) + ae;

        if ( et <= 0.0 )
        {
          delete [] f1;
          delete [] f2;
          delete [] f3;
          delete [] f4;
          delete [] f5;
          flag_return = 5;
          return flag_return;
        }

        ee = r8_abs
        ( ( -2090.0 * yp[k]
          + ( 21970.0 * f3[k] - 15048.0 * f4[k] )
          )
        + ( 22528.0 * f2[k] - 27360.0 * f5[k] )
        );

        eeoet = r8_max ( eeoet, ee / et );

      }

      esttol = r8_abs ( m_h ) * eeoet * scale / 752400.0;

      if ( esttol <= 1.0 )
      {
        break;
      }
//
//  Unsuccessful step.  Reduce the stepsize, try again.
//  The decrease is limited to a factor of 1/10.
//
      hfaild = true;
      output = false;

      if ( esttol < 59049.0 )
      {
        s = 0.9 / std::pow ( esttol, 0.2 );
      }
      else
      {
        s = 0.1;
      }

      m_h = s * m_h;

      if ( r8_abs ( m_h ) < hmin )
      {
        m_kflag = 6;
        delete [] f1;
        delete [] f2;
        delete [] f3;
        delete [] f4;
        delete [] f5;
        flag_return = 6;
        return flag_return;
      }

    }
//
//  We exited the loop because we took a successful step.
//  Store the solution for T+H, and evaluate the derivative there.
//
    *t = *t + m_h;
    for ( i = 0; i < 2; i++ )
    {
      y[i] = f1[i];
    }
    r8_flt ( y, yp );
    m_nfe = m_nfe + 1;
//
//  Choose the next stepsize.  The increase is limited to a factor of 5.
//  If the step failed, the next stepsize is not allowed to increase.
//
    if ( 0.0001889568 < esttol )
    {
      s = 0.9 / std::pow ( esttol, 0.2 );
    }
    else
    {
      s = 5.0;
    }

    if ( hfaild )
    {
      s = r8_min ( s, 1.0 );
    }

    m_h = r8_sign ( m_h ) * r8_max ( s * r8_abs ( m_h ), hmin );
//
//  End of core integrator
//
//  Should we take another step?
//
    if ( output )
    {
      *t = tout;
      delete [] f1;
      delete [] f2;
      delete [] f3;
      delete [] f4;
      delete [] f5;
      flag_return = 2;
      return flag_return;
    }

    if ( flag_return <= 0 )
    {
      delete [] f1;
      delete [] f2;
      delete [] f3;
      delete [] f4;
      delete [] f5;
      flag_return = -2;
      return flag_return;
    }

  }
# undef MAXNFE
}