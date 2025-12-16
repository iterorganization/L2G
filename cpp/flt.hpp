#ifndef FLT_H
#define FLT_H

#if defined(_WIN32)
    #ifdef flt_STATIC
        #define FLT_API
    #elif flt_EXPORTS
        #define FLT_API __declspec(dllexport)
    #else
        #define FLT_API __declspec(dllimport)
    #endif
#else
    #ifdef flt_EXPORTS
        #define FLT_API __attribute__ ((visibility ("default")))
    #else
        #define FLT_API
    #endif
#endif

// Main FLT code here. Instead of writing what is here in short text, here are
// the things the code does not do. It does not check and/or raises exceptions
// if the input data does not match. It is out of scope of the code. The main
// point of this code is that all the data is expected to be fully and
// correctly provided. Failing to do so will result in UB or SF. The focus of
// the code is to have as less code as possible in order to focus on precision
// and performance of the code. That does not mean it has high(est) performance
// possible, but just that it was developed at the time as best as it could be.

// Update: Exceptions are provided but not in hot code.

// Bi-cubic interpolation method
#include <bicubic.hpp>

#include <cmath>
#include <tlas.hpp>
#include <vector>

/// Class that solves the FL equation to obtain FL segments and checking if
/// there is any intersection with the shadowing geometry. The FL equation is
/// solved in a Cylindrical coordinate system and by using an axisymmetric
/// plasma equilibrium input data.
class FLT_API FLT
{
private:

    BICUBIC_INTERP *m_interp_psi;

    /// Array of R points describing the radial coordinates in unit of meter.\n
    /// Together the R and Z points create the (R, Z) plane in units (m, m) if
    /// one make pairs for every value in one array with the values of the other
    /// array. Or a meshgrid.
    std::vector<double> m_r_points;

    /// Left border of major radius R value that defines the computation
    /// domain. Along with m_r_max, m_z_min and m_z_max this is used to catch
    /// FLs that might escape the tokamak area. The values are set when the
    /// R and Z points are delivered via setRARR and setZARR methods. In
    /// meters.
    double m_r_min=3.0;
    /// Right border of major radius R value that defines the computation
    /// domain. In meters.
    double m_r_max=9.0;

    /// Bottom border of the Z value that defines the computation domain. In
    /// meters.
    double m_z_min=-6.0;
    /// Top border of the Z value that defines the computation domain. In
    /// meters.
    double m_z_max=6.0;

    /// Array of Z points, describing the vertical coordinates in unit of
    /// meter.\n
    /// Together the R and Z points create the (R, Z) plane in units (m, m) if
    /// one make pairs for every value in one array with the values of the other
    /// array. Or a meshgrid.
    std::vector<double> m_z_points;

    /// 1D array of poloidal flux values in unit of Webb/rad. The size of this
    /// array must be, m_NDIMR * m_NDIMZ. In other words the rows of the PSI
    /// matrix are concatenated into one array.
    std::vector<double> m_psi_values;
    std::vector<std::vector<double>> m_reshaped_psi;

    /// 1D array of poloidal flux values in units of Webb/rad. This defines the
    /// definition space of the FPOL function, which gives us the toroidal
    /// component of the magnetic field times the major radius. NOTE: Not sure
    /// if this is the correct way of obtaining the toroidal component of the
    /// magnetic field outside of the plasma.\n
    /// Note PSI == FLUX in terms of units. Different names are used to
    /// distinguish the role, as one acts as function values and other as
    /// function points.
    std::vector<double> m_flux_points; //FPOL domain

    /// 1D array of poloidal current function values in unit of m T for the
    /// toroidal component interpolator
    std::vector<double> m_fpol_values;

    /// The poloidal current function F = B_t R in the vacuum.
    double m_vacuum_fpol;

    /// Number of threads to use when running OpenMP
    int m_number_of_threads=1;

    /// Initial values in (R, Z, Phi) in units (m, m, rad).\n
    /// A vector is used regardless of parallel or sequential.
    std::vector<double> m_initial_y;

    // Solver parameters
    // The following parameters are passed to the RKF45 solver.

    /// Absolute error
    double m_abserr=1e-5;
    /// Relative error
    double m_relerr=1e-5;
    /// Time resolution. This basically tells the solver, we want points at
    /// these steps. In radians.
    double m_t_step=1e-4;

    /// Self intersection avoidance length. Length at which FLT does not check
    /// for intersections.
    double m_self_intersection_avoidance_length=5e-3;
    /// m_max_fieldline_length states the maximum length to follow a FL. In meters.
    double m_max_fieldline_length=1.0;

    /// Starting toroidal direction of a fieldline. By default to determine the
    /// direction of a fieldline, we first take the sign of the FPol or the
    /// toroidal component of the magnetic field. Then, depending on the normal
    /// of our point, by checking the sign of the dot product between the
    /// magnetic field vector and the target normal, we flip the direction so
    /// that we are following the FL on the right side of our element.
    std::vector<int> m_directions;

    /// As name implied, starting points of the field lines
    std::vector<double> m_origin_points;

    /// The radial displacement moves in unit of meter. This shift corresponds
    /// to the shift of the plasma. Instead of applying this shift before
    /// interpolation, it is added negatively to the location points when
    /// solving a fieldline.
    double m_r_move = 0.0;
    double m_z_move = 0.0;

    /// Pointer to the TLAS object that loads tinybvh::BVH objects, used for
    /// checking for intersection tests with any loaded geometry inside the
    /// tlas.

    TLAS* m_tlas;

    /// Boolean variable that tells us if magnetic data was loaded
    bool m_equilibrium_loaded=false;

    /// This is the RKF4(5) method for getting the next approximate solution to
    /// a 2D system of PDEs (non-time-dependent). The factors are calculated
    /// based on FORMULA 2 Table III in Fehlberg NASA Technical Report 287.
    /// @param[in] y Is the initial [X, Y] (or [R, Z] in Cylindrical
    ///                 coordinate system) point.
    /// @param[in] h is the current parametric time step or in the actual sense
    ///              the current toroidal angle step
    /// @param[in, out] yp is the derivative value in [R, Z] direction at y point
    /// @param[in, out] k2 is the slope factor.
    /// @param[in, out] k3 is the slope factor.
    /// @param[in, out] k4 is the slope factor.
    /// @param[in, out] k5 is the slope factor.
    /// @param[in, out] k6 is the slope factor.
    /// @param[in, out] s is the solution.
    /// @param[in, out] context is the struct necessary for running
    ///                 interpolation functions
    void fehl_step(double y[2], double h, double yp[2], double k2[2],
                   double k3[2], double k4[2], double k5[2], double k6[2],
                   double s[2], BI_DATA *context) noexcept;

    /// Function that calculates the partial derivatives of the axis-symmetric
    /// non-time dependent PDE's for fieldlines in a Cylindrical coordinate
    /// system.
    /// @param[in] y is the (R, Z) position where we want to calculate the
    ///            partial derivative values.
    /// @param[in, out] yp holds the values of the partial derivatives at
    ///                 points (R, Z)
    /// @param[in, out] context is the struct necessary for running
    ///                 interpolation functions
    void flt_pde(double y[2], double yp[2], BI_DATA *context) noexcept;

public:
    FLT();
    ~FLT();

    /// Output vector. Resized when the origin points are specified. This one
    /// contains geometry IDs which are >= 0 of stored geometries in tlas
    /// object. -1 means that the particular fieldline reach its maximum
    /// fieldline length. -2 means that the fieldline left the defined
    /// computational area of the equilibrium grid (outside the min/max of [R,
    /// Z])
    std::vector<int> m_out_geom_hit_ids;
    /// Output vector. Resized when the origin points are specified. This one
    /// contains primite IDs which are >= 0 of triangles of a stored geometry
    /// in the tlas object. -1 means that the particular fieldline reached
    /// its maximum fieldline length. -2 means that the fieldline left the
    /// defined computational area of the equilibrium grid (outside the
    /// min/max of [R, Z])
    std::vector<int> m_out_prim_hit_ids;
    /// Output vector. Resized when the origin points are specified. This one
    /// stores the lengths of the fieldlines.
    std::vector<double> m_out_fieldline_lengths;

    /// Sets the poloidal magnetic flux array and the r, z arrays describing
    /// the [R, Z] points. Vectors r and z must be sorted from smallest to
    /// largest.
    /// @param[in] r 1D vector that holds radial R points.
    /// @param[in] z 1D vector that holds vertical Z points.
    /// @param[in] psi flattened 2D to 1D vector that holds the flux values.
    ///                Size of psi is [dim r x dim z] in row major
    ///                configuration.
    void setPoloidalMagneticFlux(std::vector<double> r, std::vector<double> z,
                                 std::vector<double> psi);

    /// Sets the poloidal current value (Bt R) in vacuum. Units (m T)
    /// @param[in] vacuum_fpol is the F=Bt*R value in vacuum.
    void setVacuumFPOL(double vacuum_fpol){m_vacuum_fpol = vacuum_fpol;};

    /// Sets the plasma shift in (m, m)
    /// @param[in] r_move is the radial shift.
    /// @param[in] z_move is the vertical shift.
    void setShift(double r_move, double z_move);

    /// Sets the desired absolute error for the RKF45 solver
    /// @param[in] abserr is the RKF45 parameter for absolute accuracy.
    void setAbsError(double abserr){m_abserr = abserr;};

    /// Sets the desired relative error for the RKF45 solver
    /// @param[in] relerr is the RKF45 parameter for relative accuracy.
    void setRelError(double relerr){m_relerr = relerr;};

    /// Sets the desired toroidal angle step (the parametric time of the
    /// partial differential equations). Since the solver is an adaptive rkf45
    /// desired mean that at each time the solver tries to find the next
    /// solution it will use this step as it's starting move.
    /// @param[in] t_step is the toroidal angle step at which we wish to gather
    ///                   FLT points, used for FLT or basically tracking.
    void setDesiredStep(double t_step){
        m_t_step = t_step;
    }

    /// Sets the maximum fieldline length in meters to follow a FL.
    /// @param[in] max_fieldline_length is the maximum connection length to 
    ///                                 which we follow a FL, if the 
    ///                                 termination is determined by its 
    ///                                 length. Otherwise ignored if following 
    ///                                 till aspecified toroidal angle.
    void setMaximumFieldlineLength(double max_fieldline_length){
        m_max_fieldline_length = max_fieldline_length;
    };

    /// Calculates the barycenter of a set of triangle points. Orientation is
    /// determined from the order of the ID list of a triangle. :
    ///                      p3 - - p2
    ///                       \     /
    ///                        \   /
    ///                          p1
    /// @param[in] p1 is the first point of the triangle
    /// @param[in] p2 is the second point of the triangle
    /// @param[in] p3 is the third point of the triangle
    void getBaryCenter(double p1[3], double p2[3], double p3[3]);

    /// Sets the self-intersection avoidance length in order to avoid getting
    /// false positive self intersection results. In meters.
    ///
    /// Note that when this is set to a positive value, the RKF45 solver steps
    /// is shortened so that we do not jump out of this length in the first
    /// solver step.
    ///
    /// @param[in] value is the value at which we would like to avoid
    ///                  intersection test when following a FL. This value
    ///                  should be small, in the range of millimeters. Default
    ///                  unit is meters.
    void setSelfIntersectionAvoidanceLength(double value) {m_self_intersection_avoidance_length = value;};

    /// Sets the number of threads to use when running in parallel with OpenMP.
    /// The functions doesn't check the maximum value but insures the value
    /// is not below or equal 0.
    /// @param[in] n is the number of desired CPU to run in parallel.
    void setNumberOfThreads(int n){
        if (n <= 0){
            n = 1;
        }
        m_number_of_threads=n;
    }

    // USER BEWARE
    //     WHEN SETTING THE INPUT ARRAYS, IT'S ON YOU TO MAKE SURE THAT THE
    //     INPUT ARRAYS SIZE MATCH IN THE CONTEXT OF HOW MANY POINTS ARE.
    // USER BEWARE

    /// Set the origin points or starting points to run FLT form
    /// @param[in] origin_points is a vector of (3xN) size containing origin points.
    void setPoints(std::vector<double> origin_points);

    /// Set the toroidal direction (positive (-1) or negative (1) orientation)
    /// for the input data.
    /// @param[in] directions is a vector of size N containing the direction
    void setStartingFLDirection(std::vector<int> directions){m_directions=directions;};

    // ONLY WHEN THE FUNCTION SETPOINTS AND SETSTARTINGFLDIRECTION WHERE CALLED
    // MAY YOU CALL RUNFLT

    /// The main function that performs the FLT on the input data set via
    /// functions setPoints and setStartingFLDirection

    /// Start the FLT. When finished (successfully) the results will be stored
    /// in the m_out_* attributes.
    void runFLT();

    /// Get points of a fieldline starting from a (r, z, Phi) and in a selected
    /// toroidal direction. Points are stored in the provided vector.
    /// @param[in] r Starting R
    /// @param[in] z Starting Z
    /// @param[in] phi Starting toroidal angle
    /// @param[in] direction Starting toroidal direction. (Positive means along
    ///            the toroidal direction of the magnetic field)
    /// @param[in, out] storage vector for storing points (continuously, as in
    ///                         (r1, z1, phi1, r2, z2, phi2, ...))
    /// @param[in] with_flt Enable intersection with the shadowing geometry.
    void getFL(const double r, const double z, const double phi,
                const int direction, std::vector<double>& storage,
                const bool with_flt);

    /// For a given point in the (R, Z) space in units of (m, m)
    /// return the poloidal and toroidal component of the magnetic field.
    /// @param[in] r is the radial position of a point. In meters.
    /// @param[in] z is the vertical position of a point. In meters.
    /// @param[in, out] out holds the value of the poloidal and toroidal
    ///                     component of the magnetic field. The parameter out
    ///                     must be a vector of size 2.
    void getBCyln(double r, double z, std::vector<double> &out);

    /// For a given point in the (R, Z, Phi) space in units of (m, m, rad) return
    /// the magnetic field vector in Cartesian coordinate system. The values
    /// are written in the out variable.
    /// @param[in] r is the radial position of a point. In meters.
    /// @param[in] z is the vertical position of a point. In meters.
    /// @param[in] phi is the toroidal position of the point. In radians.
    /// @param[in, out] out holds the value of the magnetic field vector in the
    ///                     Cartesian coordinate system. Out must be a vector
    ///                     of at least size 3.
    void getBCart(double r, double z, double phi, std::vector<double> &out);

    /// For a given point in the (R, Z) space in units of (m, m) return the
    /// poloidal current function in units of m T (Bt R).
    /// @returns fpol which is the value of the poloidal current function
    /// @param[in] r is the radial position. In meters.
    /// @param[in] z is the vertical position. In meters.
    double getFPol(double r, double z);

    /// Get's the Poloidal current function in vacuum
    /// @returns  the value of the poloidal current function in the vacuum. 
    ///           It is constant for whole vacuum.
    double getVacuumFPOL(){return m_vacuum_fpol;};

    /// For a given point in the (R, Z) space in units of (m, m) return the
    /// poloidal flux value in Webb/rad
    /// @returns flux which is the poloidal magnetic flux value
    /// @param[in] r is the radial position
    /// @param[in] z is the vertical position
    double getPoloidalFlux(double r, double z);

    /// Sets the pointer to the TLAS object responsible for the tracing.
    /// @param[in] tlas_obj is the TLAS accelerated structure used for
    ///                      performing the intersection tests.
    void setTLAS(TLAS* tlas_obj);
};

#endif //FLT_H