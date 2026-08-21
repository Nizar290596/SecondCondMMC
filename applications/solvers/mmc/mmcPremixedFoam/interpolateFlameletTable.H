/*---------------------------------------------------------------------------*\
|                          Interpolation Routine                              |
|-----------------------------------------------------------------------------|
| Interpolation between two values for the the flamelet tables                |
\*---------------------------------------------------------------------------*/

scalar interpolate2Table(scalar ci, scalar c0, scalar c1, scalar Y0, scalar Y1)
{
    scalar intValue;
    intValue = Y0 + (Y1 -Y0)*(ci - c0)/(c1 - c0);
    return intValue;
}
