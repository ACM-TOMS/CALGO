// Application includes
#include    <cubUtil.h>
#include    <elrint3d.h>

static const double scale = 80; // A conservative factor to control error estimation

// Static functions
static bool is_infinite_boundary(double a, double b)
// Checks if either boundary value is infinite
{
    return ((a <= -INFINITY) || (a >= INFINITY) ||
            (b <= -INFINITY) || (b >= INFINITY));
}

static bool twoInRange(double a, double b, double scope)
// Checks if a and b are within the specified range - used to assess whether or not
// the algorithm can produce better results by continuing sequence advancement        .
{
    if (fabs(a) <= fabs(b) * scope && fabs(b) <= fabs(a) * scope)
    {
        return true;
    }
    return false;
}

static bool threeInRange(double a, double b, double c, double scope = 15.0)
// Checks if a, b and c are within the specified range
{
    if (twoInRange(a, b, scope) && twoInRange(b, c, scope) && twoInRange(a, c, scope))
    {
        return true;
    }
    return false;
}

Elrint3d::BoundaryDef::BoundaryDef(double x1, double x2, double y1, double y2, double z1, double z2, INFINITE_MAP infMap) :
cx1(x1), cx2(x2),
cy1(y1), cy2(y2),
cz1(z1), cz2(z2),
fy1(NULL), fy2(NULL),
fz1(NULL), fz2(NULL),
byCons(true), bzCons(true),
bxInfty(is_infinite_boundary(x1, x2)),
byInfty(is_infinite_boundary(y1, y2)),
bzInfty(is_infinite_boundary(z1, z2)),
inftyMap(infMap)
{
}

Elrint3d::BoundaryDef::BoundaryDef(double x1, double x2, BOUND_FUNC_1 y1, BOUND_FUNC_1 y2, double z1, double z2, INFINITE_MAP infMap) :
cx1(x1), cx2(x2),
cy1(DBL_MAX), cy2(DBL_MAX),
cz1(z1), cz2(z2),
fy1(y1), fy2(y2),
fz1(NULL), fz2(NULL),
byCons(false), bzCons(true),
bxInfty(is_infinite_boundary(x1, x2)),
byInfty(false),
bzInfty(is_infinite_boundary(z1, z2)),
inftyMap(infMap)
{
}

Elrint3d::BoundaryDef::BoundaryDef(double x1, double x2, double y1, double y2, BOUND_FUNC_2 z1, BOUND_FUNC_2 z2, INFINITE_MAP infMap) :
cx1(x1), cx2(x2),
cy1(y1), cy2(y2),
cz1(DBL_MAX), cz2(DBL_MAX),
fy1(NULL), fy2(NULL),
fz1(z1), fz2(z2),
byCons(true), bzCons(false),
bxInfty(is_infinite_boundary(x1, x2)),
byInfty(is_infinite_boundary(y1, y2)),
bzInfty(false),
inftyMap(infMap)
{
}

Elrint3d::BoundaryDef::BoundaryDef(double x1, double x2, BOUND_FUNC_1 y1, BOUND_FUNC_1 y2, BOUND_FUNC_2 z1, BOUND_FUNC_2 z2, INFINITE_MAP infMap) :
cx1(x1), cx2(x2),
cy1(DBL_MAX), cy2(DBL_MAX),
cz1(DBL_MAX), cz2(DBL_MAX),
fy1(y1), fy2(y2),
fz1(z1), fz2(z2),
byCons(false), bzCons(false),
bxInfty(is_infinite_boundary(x1, x2)),
byInfty(false),
bzInfty(false),
inftyMap(infMap)
{
}

void Elrint3d::reset()
// Resets the routine so it calls 'evaluate' multiple
// times without creating a new instance each time
{
    EmbeddedCubRule<double>::reset();
    estError = DBL_MAX;
    eflag = EF_LOW;
    bStop = false;
    bFirstNoImprovement = false;
    evalSum = 0.0;
}

double Elrint3d::BoundaryDef::trans(double* px, double* py, double* pz) const
// Transforms the coordinates from the unit cube to the original domain
// and returns the weight associated with the point (px, py, pz)
{
    double coeff_x, coeff_y, coeff_z;
    if (!bxInfty)
    {
        *px = (*px)*(cx2 - cx1) + cx1; // x coordinate
        coeff_x = (cx2 - cx1);
    }
    else
    {
        coeff_x = inftyMap(cx1, cx2, *px);
    }

    if (!byInfty)
    {
        double ty1, ty2;
        if (byCons)
            // Constant boundary for y
        {
            ty1 = cy1;
            ty2 = cy2;
        }
        else
            // y boundary is variable
        {
            ty1 = fy1(*px);
            ty2 = fy2(*px);
        }

        *py = (*py) * (ty2 - ty1) + ty1; // y coordinate
        coeff_y = (ty2 - ty1);
    }
    else
    {
        coeff_y = inftyMap(cy1, cy2, *py);
    }

    if (!bzInfty)
    {
        double tz1, tz2;
        if (bzCons)
            // Constant boundary for z
        {
            tz1 = cz1;
            tz2 = cz2;
        }
        else
            // z boundary is variable
        {
            tz1 = fz1(*px, *py);
            tz2 = fz2(*px, *py);
        }

        *pz = (*pz) *(tz2 - tz1) + tz1; // z coordinate
        coeff_z = (tz2 - tz1);
    }
    else
    {
        coeff_z = inftyMap(cz1, cz2, *pz);
    }

    return coeff_x * coeff_y * coeff_z; // Associated weight
}

double Elrint3d::evaluate()
{
    if (bound.bInfinite())
    {
        double res1 = EmbeddedCubRule<double>::evaluate();
        int ef1 = errFlag();
        double ee1 = estErr();
        nEvals = EmbeddedCubRule<double>::evals();

        if (ef1 != NORMAL)
            // Abnormal termination, indicating it is likely that the current
            // infinite mapping is not suitable.  Switch to alternative mapping.
        {
            reset();
            bound.setInfMapping(logarithmicInfMap<double>);
            double res2 = EmbeddedCubRule<double>::evaluate();
            int ef2 = errFlag();
            double ee2 = estErr();
            nEvals += EmbeddedCubRule<double>::evals();

            if (ef2 == NORMAL || ee2 < ee1)
                // Second choice of infinite mapping is successful
                // or returns a better result than the first mapping
            {
                estError = ee2;
                eflag = ef2;
                return res2;
            }
            else
                // Otherwise, use the result based on the original mapping
            {
                estError = ee1;
                eflag = ef1;
                return res1;
            }
        }
        else
            // A satisfactory result is achieved with the first mapping
        {
            estError = ee1;
            eflag = ef1;
            return res1;
        }
    }
    else
    {
        return EmbeddedCubRule<double>::evaluate();
    }
}

void Elrint3d::postAdvancement(int iLat)
// Computes the error estimate for the latest cubature and tests for termination (either
// because the required accuracy has been achieved or because it cannot be achieved).
// (Note that several heuristics are used in this function; it is not thread safe)
{
    static double preEstErr = DBL_MAX;
    preEstErr = estError; // Use of static variable => not thread safe

    int curLat = iCurLat(); // Current lattice index

    if (curLat >= 3)
        // An error estimate can only be computed after one augmentation cycle
    {
        double e[3];
        double diff_e; // Error estimate based on difference between cubatures

        e[0] = rel_err(interCubs[curLat], interCubs[curLat - 3]);
        e[1] = rel_err(interCubs[curLat], interCubs[curLat - 2]);
        e[2] = rel_err(interCubs[curLat], interCubs[curLat - 1]);
        estError = diff_e = errPerCycle(e);

        if (curLat >= 6)
            // The extrapolated value can be computed after two augmentation cycles have been completed
        {
            int start_i = curLat % 6;
            double ev = epsilon_extrapolation(interCubs[start_i], interCubs[start_i + 3], interCubs[start_i + 6]);
            estError = rel_err(interCubs[curLat], ev);
        }

        if ((curLat == 6) && e[2] > 1.0e-3)
            // After two augmentation cycles, if the difference between successive cubatures from
            // the augmented sequence is not small enough, the algorithm will be terminated
        {
            bStop = true;
            estError = max3(e[1], e[2], e[0]) * pow(scale, 2);
            eflag = NOT_SUITABLE;
            return;
        }

        if (curLat >= 3 && max3(e[1], e[2], e[0]) > 1.0e-2)
            // If, after the first augmentation, there is little improvement
            // in the cubature, use a very conservative error estimate
        {
            estError = diff_e * pow(scale, 5.0);

        }

        if (curLat >= 7 && estError > 5.0e-6)
            // Terminate the algorithm if the accuracy obtained beyond
            // two augmentation cycles is deemed unsatisfactory
        {
            eflag = NOT_IMPROVEMENT;
            bStop = true;
            return;
        }

        // Test for the presence of rounding error
        if (e[2] * 30 > e[0] && curLat >= 6 && e[2] > 1.0e-14)
        {
            if (bFirstNoImprovement)
            {
                eflag = NOT_IMPROVEMENT;
                bStop = true;
                estError = diff_e*scale;
                return;
            }
            else
            {
                estError = diff_e*scale;
            }
            bFirstNoImprovement = true;
        }
        else
        {
            if (bFirstNoImprovement)
            {
                estError = diff_e*scale;
            }
            bFirstNoImprovement = false;
        }

        // From the second augmentation cycle onwards, check the behaviour
        // of the error estimates for the four most recent cubatures
        if (e[2] > 5.0e-13 && curLat >= 6 &&
                threeInRange(fabs(interCubs[curLat] - interCubs[curLat - 1]),
                fabs(interCubs[curLat - 1] - interCubs[curLat - 2]),
                fabs(interCubs[curLat - 2] - interCubs[curLat - 3])))
        {
            estError = diff_e*scale;
        }

        // Choose betwen the conservative error estimate and the error
        // estimate based on the extrapolated cubature value 
        if (curLat == 6 && estError < 5.0e-15 && preEstErr > 1.0e-12)
        {
            estError = max3(e[0], e[1], e[2]);
        }

        // See if the requested accuracy has been achieved
        if (estError <= max(relReqTol, absReqTol * interCubs[curLat]))
        {
            bStop = true;
            eflag = NORMAL;
            return;
        }

        // See if the maximum number of points has been reached
        if (curLat >= 9 || maxEvals <= sequence()->count(curLat))
        {
            bStop = true;
            eflag = REACH_MAX_NUMBER;
            return;
        }
    }
}

bool Elrint3d::bTerminated() const
// Checks if the algorithm should be terminated
{
    return bStop || iCurLat() >= 9;
}

double Elrint3d::funcEval(const double* pa) const
// Returns the value of the integrand function at the point 'pa'
{
    return rIntegrand.fun(pa);
}

double Elrint3d::trans(double * const pa) const
// Transforms coordinates in the unit cube to coordinates in the original integration
// domain and returns the weight associated with the point.
// (This function overrides the virtual function in 'EmbeddedCubRule'.)
{
    return bound.trans(&pa[0], &pa[1], &pa[2]) * pa[3];
}

double Elrint3d::sequenceAdvance(int iLat)
// Overrides virtual function in EmbeddedCubRule 
{
    double sum = EmbeddedCubRule<double>::sequenceAdvance(iLat);
    evalSum += sum;
    interCubs[iLat] = evalSum / sequence()->count(iLat); // Intermediate cubatures
    return sum;
}

double Elrint3d::errPerCycle(const double e[]) const
{
    if (min3(e[0], e[1], e[2]) > 1.0e-3)
        // It seems the algorithm can't do much, so return a very conservative error estimate
    {
        return max3(e[0], e[1], e[2]) * scale;
    }

    if (e[2] != 0.0)
    {
        if (e[1] / e[2] < 5.0)
            // Little improvement after the first augmentation
        {
            return e[0];
        }
    }

    switch (iCurLat() % 3)
    {
        case 0:
            return e[2] * 4.5; // Large improvement at end of augmentation cycle; use least conservative error estimate
        case 1:
            return e[1] * 3.0; // Always use a conservative error estimate for the first intermediate lattice after an augmentation cycle
        case 2:
            return e[1] * 1.5; // Use a less conservative estimate for the second intermediate lattice
        default:
            // Cannot happen
            throw runtime_error("Abnormal condition");
    }
}

double Elrint3d::sum(const double* pa, int n) const
// Uses Kahan summation
{
    return kahan_sum(pa, n);
}

