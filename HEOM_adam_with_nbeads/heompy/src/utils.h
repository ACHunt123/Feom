// File: utils.h
#ifndef UTILS_H
#define UTILS_H
#include <string>
#include "general.h"
#include "timer.h"
//*****************************************************************************
// Signal handling
//*****************************************************************************
void signal2exitfile(int signum);
//*****************************************************************************
// check_file
//*****************************************************************************
bool check_file(const String &filename);
bool check_exit();
//*****************************************************************************
// check_stability
//*****************************************************************************
bool check_stability(const Real &threshold, const Real &position);
//*****************************************************************************
// Write-out functions
//*****************************************************************************
template <typename T>
void write_out(const String & filename, const T &item1);
template <typename T>
void write_out_default(const String & filename, const T &item1);
template <typename T1, typename T2>
void write_out(const String & filename, const T1 &item1, const T2 &item2);
template <typename T>
void write_out(const String & filename, const T &item1, const Vec &vec);
template <typename T>
void write_out(const String & filename, const std::initializer_list<T> &numlist);
void write_out_counted(const String & filename,
                                const std::initializer_list<Real> &numlist);
void write_out_counted(const String & filename, const Vec &numlist);
template <typename T1, typename T2>
void write_out(const String & filename, const T1 &item1, const T2 &item2,
        const Vec &vec);
template <typename T>
void write_out(const String & filename, const T &item1, const Vec &vec);
void write_out(const String & filename, const Vec &vec);
void write_out_columns(const String & filename,
                                    const std::initializer_list<Mat> &veclist);
//*****************************************************************************
// Text extraction utilities
//*****************************************************************************
String extract_string (const String & original);
CppVec<String> extract_string_to_vector (const String & original);
void extract_line (const String & line, int &output);
void extract_line (const String & line, long &output);
void extract_line (const String & line, double &output);
void extract_line (const String & line, IntVec &output);
void extract_line (const String & line, Vec &output);
void extract_line (const String & line, String &output);
//*****************************************************************************
// test_precision
//*****************************************************************************
void test_precision();
//*****************************************************************************
// Progress_reporter
//*****************************************************************************
class Progress_reporter
{
private:
    Timer timer;
    Int percent;
    Real total;
public:
    Progress_reporter (const Real &inp_total);
    void print_progress(const Int &iteration);
};
//*****************************************************************************
// Recorder
//*****************************************************************************
class Recorder
{
public:
    bool used{false};
    Int i_count{};
    bool first_trajectory{true};
    bool first_record{true};
    bool first_recordA{true};
    Int allocation_length{1000};
    bool autocorrelation{true};
    CplxVec cplx_init_val;
    CplxMat cplxtraj;
    Vec outvec;
    bool collect_data{false};
    Int n_steps{}, n_horizontal_samples{}, n_trajs{0}, index{-1},
                                    n_discarded{0}, n_data_samples{0}, n_aggr{};
    Int n_quantities{};
    bool horizontal{false};
    bool is_data_complex{false};
    Int n_data{};
    Mat data, data_tmp, data_stdev;
    Real partfun{0}, partfun_tmp{0};
    Real discarded_t{0}, discarded_t_stdev{0};
    Recorder() {}
    void recordA(const std::initializer_list<Cplx> &numlist);
    template <typename T>
    void recordA(const T &numlist);
    void record(const std::initializer_list<Cplx> &numlist);
    template <typename T>
    void record(const T &numlist);
    void discard_trajectory();
    //void average_samples();
    bool check_convergence(const Real & error_threshold);
    void add_constant(const Real &constant);
    void write_data(const String &name, const String &kind,
                                const Real &timestep, const Int &quantity=-1);
    void write_log(const String & logfile,const Real &timestep);
};
//*****************************************************************************
// TCF_recorder
//*****************************************************************************
class TCF_recorder: public Recorder
{
public:
    TCF_recorder() {}
    void end_trajectory(const Cplx &phase=Cplx(1,0), const Real &partfun_sample=1);
};
//*****************************************************************************
// Average_recorder
//*****************************************************************************
class Average_recorder: public Recorder
{
public:
    Average_recorder() {}
    void end_trajectory(const Cplx &phase=Cplx(1,0), const Real &partfun_sample=1);
};

//template <typename T>
//void aggregate_recorders(T & tot, T & src, const int &rank, const int &n_procs);

#endif
