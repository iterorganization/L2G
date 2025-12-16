#ifndef BICUBIC_H
#define BICUBIC_H

#if defined(_WIN32)
    #ifdef flt_STATIC
        #define BICUBIC_API
    #elif flt_EXPORTS
        #define BICUBIC_API __declspec(dllexport)
    #else
        #define BICUBIC_API __declspec(dllimport)
    #endif
#else
    #ifdef flt_EXPORTS
        #define BICUBIC_API __attribute__ ((visibility ("default")))
    #else
        #define BICUBIC_API
    #endif
#endif

#include <vector>

// BI_DATA struct is used for providing context when obtaining interpolation
// values. Also more OpenMP friendly than plain shared class attributes.

struct BICUBIC_API BI_DATA
{
    double r=0;
    double z=0;
    double val=0;
    double valdx=0;
    double valdy=0;
    double dx=0.0;
    double dy=0.0;
    double inv_dx=1.0;
    double inv_dy=1.0;
    double minx=0.0;
    double miny=0.0;
    double maxx=0.0;
    double maxy=0.0;
    int nx=0;
    int ny=0;
};

/// Class that implements the 2D bicubic interpolation algorithm. The aim is to
/// have as light-as-possible code that is accurate and most importantly very
/// fast. The test directory contains the testing programs to see if the code
/// is accurate enough, fast enough and if possible scalable with more cores
/// via OpenMP. OpenMP scalability is ensured with the help if BI_DATA struct.
/// It acts as a context for each call request, reducing the number of OpenMP
/// false-sharing that might happen if the variables in the struct were class
/// attributes.
class BICUBIC_API BICUBIC_INTERP
{

private:
    /// Rectilinear weights. Obtained from the input X- and Y-axis arrays and
    /// are needed when obtaining the partial derivatives.
    double m_dx=1.0, m_dy=1.0;

    /// Domain borders.
    double m_minx, m_maxx, m_miny, m_maxy;

    /// Input 1D arrays of the X- and Y-axis. These two arrays form the cell
    /// grid or a rectilinear mesh grid.
    std::vector<double> m_x, m_y;

    /// Size of the X- and Y-axis arrays.
    int m_nx, m_ny;

    /// 1D vector that holds the constants. The Size is (Nx*Ny*16);
    std::vector<double> m_a;

public:

    BICUBIC_INTERP(){};
    ~BICUBIC_INTERP(){};

    /// 2D vectors for storing the discrete values of the function we wish to
    /// interpolate, the first partial derivatives and the mixed partial
    /// derivative.
    std::vector<std::vector<double>> m_f;
    std::vector<std::vector<double>> m_fdx;
    std::vector<std::vector<double>> m_fdy;
    std::vector<std::vector<double>> m_fdxdy;

    /// Function that takes the X-axis and Y-axis 1D arrays and the 2D array
    /// which we wish to interpolate. For the bicubic interpolation the first
    /// partial derivatives are required, as well as the second partial mixed
    /// derivative. This function obtains the derivative values by using higher
    /// finite difference table (combination of central, forward and backward,
    /// depending where are we.).
    ///
    /// WARNING! If you do not set the arrays or the input data, then the other
    /// functions will cause a segfault.
    ///
    /// @param[in] x is the X-axis 1D vector array of size Nx - Columns
    /// @param[in] y is the Y-axis 1D vector array of size Ny - Rows
    /// @param[in] f is the X-axis 2D vector array (vector of vectors) of size
    ///              (Ny, Nx)
    void setArrays(std::vector<double> x, std::vector<double> y,
                   std::vector<std::vector<double>> f);

    /// Function that sets the constant values to a BI_DATA context. Used once
    /// on a BI_DATA construct after the setArrays was used and before calling
    /// gevValues.
    /// @param[in] context is the structure containing the information on the
    ///                    interpolated data, as well as holding the input and
    ///                    output variables for interpolation queries.
    void populateContext(BI_DATA *context) noexcept {
        context->dx = m_dx;
        context->dy = m_dy;
        context->inv_dx = 1.0 / m_dx;
        context->inv_dy = 1.0 / m_dy;
        context->minx = m_minx;
        context->miny = m_miny;
        context->nx = m_nx;
        context->ny = m_ny;
        context->maxx = m_maxx;
        context->maxy = m_maxy;
    }

    /// Obtain the function value and its first partial derivatives at point
    /// (x, y). **DOES NOT EXPTRAPOLATE**. If the Query point points outside
    /// the numerical domain defined by the function domain, it is
    /// automatically clipped to the border. The BI_DATA construct holds the
    /// arguments and return values altogether. You need to set the BI_DATA->r,
    /// and BI_DATA->z and obtain the values in variables BI_DATA->val,
    /// BI_DATA->valdx and BI_DATA->valdy.
    /// @param[in] context is the structure containing the information on the
    ///                    interpolated data, as well as holding the input and
    ///                    output variables for interpolation queries.
    void getValues(BI_DATA *context) noexcept;

    /// Obtain first derivative values at (x, y)
    /// @param[in] context is the structure containing the information on the
    ///                    interpolated data, as well as holding the input and
    ///                    output variables for interpolation queries.
    void getFirstDerivativeValues(BI_DATA *context) noexcept;

    /// Obtain secondary derivative values at (x, y)
    /// @param[in] context is the structure containing the information on the
    ///                    interpolated data, as well as holding the input and
    ///                    output variables for interpolation queries.
    void getSecondDerivativeValues(BI_DATA *context) noexcept;
};
#endif /*BICUBIC_H*/
