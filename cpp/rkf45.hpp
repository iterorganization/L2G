#ifndef RKF45_H
#define RKF45_H

#if defined(_WIN32)
    #ifdef flt_STATIC
        #define RKF45_API
    #elif flt_EXPORTS
        #define RKF45_API __declspec(dllexport)
    #else
        #define RKF45_API __declspec(dllimport)
    #endif
#else
    #ifdef flt_EXPORTS
        #define RKF45_API __attribute__ ((visibility ("default")))
    #else
        #define RKF45_API
    #endif
#endif

#include <bicubic.hpp>

class RKF45_API RKF45
{
private:
    double m_abserr_save;
    int m_flag_save;
    double m_h;
    int m_init;
    int m_kflag;
    int m_kop;
    int m_nfe;
    double m_relerr_save;
    double m_remin;
    double m_r8epsilon;

    // Externally set variables
    BICUBIC_INTERP *m_interp_psi;
    bool m_interpolator_set=false;
    int m_omp_thread=0;
    double m_r_move = 0.0;
    double m_z_move = 0.0;
    double m_vacuum_fpol;
    void r8_flt(double y[2], double yp[2]);
    double r8_abs(double x);
    double r8_epsilon();
    double r8_max(double x, double y);
    double r8_min(double x, double y);
    void r8_fehl(double y[2], double h, double yp[2], double f1[2],
                 double f2[2], double f3[2], double f4[2], double f5[2],
                 double s[2]);
    double r8_sign(double x);

public:
    RKF45();
    ~RKF45(){if (m_interpolator_set){delete m_interp_psi;}};
    void set_r_move(double r_move){m_r_move = r_move;};
    void set_z_move(double z_move){m_z_move = z_move;};
    void set_vacuum_fpol(double vacuum_fpol){m_vacuum_fpol = vacuum_fpol;};
    void set_interpolator(BICUBIC_INTERP *interp){m_interp_psi = interp;m_interpolator_set=true;};
    BICUBIC_INTERP* get_interpolator(){return m_interp_psi;};
    void set_omp_thread(int omp_thread){m_omp_thread=omp_thread;};
    int r8_rkf45(double y[2], double yp[2], double *t, double tout,
                 double *relerr, double abserr, int flag);

};

#endif /*RKF45_H*/