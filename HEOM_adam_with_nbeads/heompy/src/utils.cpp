// File: utils.cpp
#include "general.h"
#include <iomanip>
#include <iostream>
// For extract_string
#include <string>
#include <vector>
// For extraction
#include <sstream>
#include <algorithm>
#include <iterator>
// For isnan:
#include <cmath>
// For file manipulation and C-style output
#include <cstdio>
// For exit
#include <cstdlib>

#include <initializer_list>
#include "utils.h"

int precision{17};
int width{26};

bool check_file(const String &filename)
{
    std::ifstream f(filename.c_str());
    return f.good();
}
bool check_exit()
{
    String filename{"exit"};
    if (check_file(filename))
    {
        return true;
    }
    else
    {
        return false;
    }
}
//*****************************************************************************
// Signal handling
//*****************************************************************************
void signal2exitfile(int signum)
{
    std::cout << "\nInterrupt signal (" << signum << ") received.\n";
    std::ofstream file("exit",std::ios::app);
    if (!file.is_open())
    {
        std::cerr << "Uh oh, exitfile  could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << "exit";
        file.close();
    }
}
//*****************************************************************************
// check_stability
//*****************************************************************************
bool check_stability(const Real &threshold, const Real &position)
{
    if (!std::isnormal(position))
    {
        std::cerr << "Position is NaN. Maybe the time step is too big.\n";
        return false;
    }
    else if (std::abs(position)>threshold)
    {
        std::cerr << "\nPosition is too big\n";
        std::cerr << "An instability has occurred\n";
        std::cerr << "Last position is:\n"<< position << "\n";
        return false;
    }
    return true;
}

//*****************************************************************************
// Write-out functions
//*****************************************************************************
template <typename T>
void write_out_default(const String & filename, const T &item1)
{
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << item1;
        file << "\n";
        file.close();
    }
}
template void write_out_default<>(const String & filename, const Real &item1);
template void write_out_default<>(const String & filename, const Int &item1);
template void write_out_default<>(const String & filename, const String & item1);
//*****************************************************************************
template <typename T>
void write_out(const String & filename, const T & item1)
{
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        file << item1;
        file << "\n";
        file.close();
    }
}
template void write_out<>(const String & filename, const Real &item1);
template void write_out<>(const String & filename, const Int &item1);
template void write_out<>(const String & filename, const String & item1);
//*****************************************************************************
template <typename T1, typename T2>
void write_out(const String & filename, const T1 & item1, const T2 &item2)
{
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        file << item1;
        file << "     ";
        file << item2;
        file << "\n";
        file.close();
    }
}
template void write_out<>(const String & filename, const Real &item1,
                                                            const Real &item2);
//*****************************************************************************
void write_out_counted(const String & filename,
                                const std::initializer_list<Real> &numlist)
{
    static Int count{0};
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        file << static_cast<Real>(count);
        for(auto elem: numlist)
        {
            file << "     ";
            file << elem;
        }
        file << "\n";
        file.close();
        count++;
    }
}
//*****************************************************************************
void write_out_counted(const String & filename, const Vec &numlist)
{
    static Int count{0};
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        file << static_cast<Real>(count);
        for(auto elem: numlist)
        {
            file << "     ";
            file << elem;
        }
        file << "\n";
        file.close();
        count++;
    }
}
//*****************************************************************************
template <typename T>
void write_out (const String & filename, const std::initializer_list<T> &numlist)
{
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        for(auto elem: numlist)
        {
            file << elem;
            file << "     ";
        }
        file << "\n";
        file.close();
    }
}
template void write_out (const String & filename, const std::initializer_list<Real> &numlist);
//*****************************************************************************
template <typename T1, typename T2>
void write_out(const String & filename, const T1 &item1, const T2 &item2,
                                                                const Vec &vec)
{
    Int length{static_cast<Int>(vec.n_rows)};
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        file << item1;
        file << "     ";
        file << item2;
        for (int i=0;i<length;i++)
        {
            file << vec(i);
        }
        file << "\n";
        file.close();
    }
}
template void write_out<>(const String & filename, const Real &item1,
                                            const Real &item2, const Vec &vec);
//*****************************************************************************
template <typename T>
void write_out(const String & filename, const T &item1, const Vec &vec)
{
    Int length{static_cast<Int>(vec.n_rows)};
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        file << item1;
        for (int i=0;i<length;i++)
        {
            file << "     ";
            file << vec(i);
        }
        file << "\n";
        file.close();
    }
}
template void write_out<>(const String & filename, const Real &item1,
                                                            const Vec &vec);
//*****************************************************************************

void write_out (const String & filename, const Vec &vec)
{
    Int length{static_cast<Int>(vec.n_rows)};
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        file << std::setprecision(precision);
        file << std::setw(width);
        file << std::scientific;
        for (int i=0;i<length;i++)
        {
            file << vec(i);
            file << "     ";
        }
        file << "\n";
        file.close();
    }
}

//******************************************************************************

void write_out_columns (const String & filename,
                                   const std::initializer_list<Mat> &veclist)
{
    // Opening a file
    std::ofstream file(filename,std::ios::app);
    if (!file.is_open())
    {
        std::cerr << "Uh oh, "<< filename <<" could not be opened for writing!"
                                                                << std::endl;
    }
    else
    {
        // Checking that columns are of the same length
        Int length{0};
        bool first{true}, equal{true};
        for(auto elem: veclist)
        {
            if (first)
            {
                length = elem.n_rows;
                first = false;
            }
            else
            {
                if (Int(elem.n_rows) != length)
                {
                    equal = false;
                }
            }
        }
        if (equal)
        {
            file << std::setprecision(precision);
            file << std::setw(width);
            file << std::scientific;
            for (Int i=0; i<length; i++)
            {
                for(auto elem: veclist)
                {
                    for (Int j=0; j<elem.n_cols; j++)
                    {
                        file << elem(i,j);
                        file << "     ";
                    }
                }
                file << "\n";
            }
        }
        else
        {
            std::cerr << "Uh oh, vectors that were supposed to be written into "
                        << filename <<" were of different lengths" << std::endl;
        }
        file.close();
    }
}
//*****************************************************************************
// Text extraction utilities
//*****************************************************************************
String extract_string (const String & original)
{
    std::istringstream iss(original);
    std::vector<String> results(std::istream_iterator<String>{iss},
                                         std::istream_iterator<String>());
    return results[1];
}
CppVec<String> extract_string_to_vector (const String & original)
{
    std::istringstream iss(original);
    CppVec<String> results(std::istream_iterator<String>{iss},
                                         std::istream_iterator<String>());
    results.erase(results.begin());
    return results;
}

void extract_line (const String & line, int &output)
{
    String current_value{};
    current_value = extract_string(line);
    output = std::stoi(current_value);
}
void extract_line (const String & line, long &output)
{
    String current_value{};
    current_value = extract_string(line);
    output =  std::stoi(current_value);
}

void extract_line (const String & line, double &output)
{
    String current_value{};
    current_value = extract_string(line);
    output =  std::stod(current_value);
}

void extract_line (const String & line, IntVec &output)
{
    CppVec<String> temp_vec;
    temp_vec = extract_string_to_vector(line);
    output = IntVec(temp_vec.size());
    for (Int i=0; i<temp_vec.size(); i++)
    {
        output(i) = std::stoi(temp_vec[i]);
    }
}
void extract_line (const String & line, Vec &output)
{
    CppVec<String> temp_vec;
    temp_vec = extract_string_to_vector(line);
    output = Vec(temp_vec.size());
    for (Int i=0; i<temp_vec.size(); i++)
    {
        output(i) = std::stod(temp_vec[i]);
    }
}
void extract_line (const String & line, String &output)
{
    String current_value{};
    current_value = extract_string(line);
    output = current_value;
}
//*****************************************************************************
// test_precision
//*****************************************************************************
void test_precision ()
{
    std::cout << "Running test_precision function" << std::endl;
    float PI0;
    double PI1;
    long double PI2;
    PI0 = 3.14159265358979323846264338327950288419716939937510L;
    PI1 = 3.14159265358979323846264338327950288419716939937510L;
    PI2 = 3.14159265358979323846264338327950288419716939937510L;
    std::cout.precision(50);
    std::cout << "PI0 = " << PI0 << std::endl;
    std::cout << "PI1 = " << PI1 << std::endl;
    std::cout << "PI2 = " << PI2 << std::endl;
    std::cout << "PIx = 3.14159265358979323846264338327950288419716939937510" << std::endl;
    std::cout << "sizeof(float) = " << sizeof(float) << std::endl;
    std::cout << "sizeof(double) = " << sizeof(double) << std::endl;
    std::cout << "sizeof(long double) = " << sizeof(long double) << std::endl;
}

//*****************************************************************************
// Progress_reporter
//*****************************************************************************
Progress_reporter::Progress_reporter(const Real &inp_total)
{
    timer.start();
    total = inp_total;
    percent = static_cast<Int>(total/100);
    if (percent==0) percent=1;
}

void Progress_reporter::print_progress(const Int &iteration)
{
    if (iteration%percent == 0 || iteration < 5)
    {
        Real total_time, time_remaining;
        Real time = timer.end();
        total_time = time*total/iteration;
        time_remaining = total_time-time;
        if (time_remaining<60.0)
        {
            std::printf("\rprogress: %2ld %%, expected time remaining,     %6.4f s",
                static_cast<long>(iteration*100/total), time_remaining);
        }
        else if (time_remaining<3600.0)
        {
            time_remaining /= 60.0;
            std::printf("\rprogress: %2ld %%, expected time remaining,   %6.4f min",
                static_cast<long>(iteration*100/total), time_remaining);
        }
        else if (time_remaining<(3600*24))
        {
            time_remaining /= 3600.0;
            std::printf("\rprogress: %2ld %%, expected time remaining, %6.4f hours",
                static_cast<long>(iteration*100/total), time_remaining);
        }
        else
        {
            time_remaining /= 3600.0*24;
            std::printf("\rprogress: %2ld %%, expected time remaining,  %6.4f days",
                static_cast<long>(iteration*100/total), time_remaining);
        }

        fflush(stdout);
    }
}
//*****************************************************************************
// Recorder
//*****************************************************************************
// If the TCF is not an autocorrelation function, the first of the two quantities
// (assuming that this is a two-point correlation function) is recorded with
// recordA.
void Recorder::recordA(const std::initializer_list<Cplx> &numlist)
{
    if (!used) used = true;
    if (first_recordA)
    {
        autocorrelation = false;
        first_recordA = false;
        n_quantities = numlist.size();
        cplx_init_val = CplxVec(n_quantities, a::fill::zeros);
    }
    if (Int(numlist.size()) != n_quantities)
    {
        std::cerr <<
           "\nRecorder: The number of recorded quantities has changed.\n";
        std::exit(1);
    }
    i_count=0;
    for(auto elem: numlist)
    {
        cplx_init_val(i_count) = elem;
        i_count++;
    }
}
//*****************************************************************************
template <typename T>
void Recorder::recordA(const T &numlist)
{
    if (!used) used = true;
    if (first_recordA)
    {
        autocorrelation = false;
        first_recordA = false;
        n_quantities = numlist.n_rows;
        cplx_init_val = CplxVec(n_quantities, a::fill::zeros);
    }
    if (static_cast<Int>(numlist.n_rows) != n_quantities)
    {
        std::cerr <<
           "\nRecorder: The number of recorded quantities has changed.\n";
        std::exit(1);
    }
    for (int i=0; i<n_quantities; i++)
    {
        cplx_init_val(i) = numlist(i);
    }
}
template void Recorder::recordA<>(const Vec &numlist);
template void Recorder::recordA<>(const CplxVec &numlist);
//*****************************************************************************
void Recorder::record(const std::initializer_list<Cplx> &numlist)
{
    if (!used) used = true;
    if (first_record)
    {
        first_record=false;
        n_quantities = numlist.size();
        cplxtraj = CplxMat(allocation_length, n_quantities, a::fill::zeros);
    }
    if (Int(numlist.size()) != n_quantities)
    {
        std::cerr <<
           "\nRecorder: The number of recorded quantities has changed.\n";
        std::exit(1);
    }
    index += 1;
    i_count=0;
    for(auto elem: numlist)
    {
        cplxtraj(index,i_count) = elem;
        i_count++;
    }
    // If vector is full, increase the size
    if (first_trajectory && (index == static_cast<Int>(cplxtraj.n_rows-1)))
    {
        cplxtraj.resize(cplxtraj.n_rows+allocation_length, n_quantities);
    }
}
//*****************************************************************************
template <typename T>
void Recorder::record(const T &numlist)
{
    if (!used) used = true;
    if (first_record)
    {
        first_record=false;
        n_quantities = numlist.n_rows;
        cplxtraj = CplxMat(allocation_length, n_quantities, a::fill::zeros);
    }
    if (static_cast<Int>(numlist.n_rows) != n_quantities)
    {
        std::cerr <<
           "\nRecorder: The number of recorded quantities has changed.\n";
        std::exit(1);
    }
    index += 1;
    for (int i=0; i<n_quantities; i++)
    {
        cplxtraj(index,i) = numlist(i);
    }
    // If vector is full, increase the size
    if (first_trajectory && (index == static_cast<Int>(cplxtraj.n_rows-1)))
    {
        cplxtraj.resize(cplxtraj.n_rows+allocation_length, n_quantities);
    }
}
template void Recorder::record<>(const Vec &numlist);
template void Recorder::record<>(const CplxVec &numlist);
//*****************************************************************************
void Recorder::discard_trajectory()
{
    if (used)
    {
        discarded_t += index;
        discarded_t_stdev += index*index;
        index = -1;
        n_discarded += 1;
    }
}
//*****************************************************************************
bool Recorder::check_convergence(const Real & error_threshold)
{
    Real error_max, range;
    bool is_converged{false};
    if (data.size()!=0)
    {
        error_max = a::max(data_stdev.col(0));
        range = a::max(data.col(0))-a::min(data.col(0));
        std::cout << "\nRange " << range << "\n";
        std::cout << "\nMaximal error " << error_max << "\n";
        std::cout << "Threshold value " << range*error_threshold << "\n\n";
        if (error_max<range*error_threshold) is_converged = true;
        // If there is NaN in the data, the calculation is not converged
        for (int i=0; i<n_data; i++)
        {
            if (std::isnan(data_stdev(i,0)))
            {
                is_converged = false;
                break;
            }
        }
    }
    return is_converged;
}
//*****************************************************************************
void Recorder::add_constant(const Real &constant)
{
    data = data + constant;
}
//*****************************************************************************
void Recorder::write_data(const String &name, const String &kind,
                                    const Real &timestep, const Int &quantity)
{
    if (n_trajs<1||data.size()==0)
    {
        std::cout << "No TCF written since no trajectories were recorded.\n";
        std::cout << "(at least one on each CPU is required)\n";
    }
    else
    {
        String file_real_au_dat{get_file_real_au_dat(name,kind)};
        String file_real_au_err{get_file_real_au_err(name,kind)};
        //String file_real_fs_dat{get_file_real_fs_dat(name,kind)};
        //String file_real_fs_err{get_file_real_fs_err(name,kind)};
        String file_data{};
        Vec indeces(n_data+1);
        for (int i_step=0;i_step<=n_data;i_step++)
        {
            indeces(i_step) = i_step;
        }
        if (quantity==-1) // i.e. all quantities
        {
            file_data = file_real_au_dat;
            std::remove(file_data.data());
            write_out_columns(file_data,{indeces*timestep,
                                        data
                                        });
            file_data = file_real_au_err;
            std::remove(file_data.data());
            write_out_columns(file_data,{indeces*timestep,
                                        data_stdev
                                        });

            //file_data = file_real_fs_dat;
            //std::remove(file_data.data());
            //write_out_columns(file_data,{indeces*timestep*au2fs,
            //                            data
            //                            });
            //file_data = file_real_fs_err;
            //std::remove(file_data.data());
            //write_out_columns(file_data,{indeces*timestep*au2fs,
            //                            data_stdev
            //                            });

        }
        else
        {
            file_data = file_real_au_dat;
            std::remove(file_data.data());
            write_out_columns(file_data,{indeces*timestep,
                                        data.col(quantity)
                                        });
            file_data = file_real_au_err;
            std::remove(file_data.data());
            write_out_columns(file_data,{indeces*timestep,
                                        data_stdev.col(quantity)
                                        });

            //file_data = file_real_fs_dat;
            //std::remove(file_data.data());
            //write_out_columns(file_data,{indeces*timestep*au2fs,
            //                            data.col(quantity)
            //                            });
            //file_data = file_real_fs_err;
            //std::remove(file_data.data());
            //write_out_columns(file_data,{indeces*timestep*au2fs,
            //                            data_stdev.col(quantity)
            //                            });

        }
    }
}
//*****************************************************************************
void Recorder::write_log(const String & logfile,const Real &timestep)
{
    std::ofstream log(logfile,std::ios::app);
    if (!log.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< logfile <<" could not be opened for writing!"
                                                                << std::endl;
        std::exit(1);
    }
    if (n_trajs==0)
    {
        log << "Recorder log:\n"
            << n_trajs << " trajectories averaged\n"
            << n_data_samples
            << " trajectory samples (including horizontal sampling)\n"
            << n_discarded << " trajectories discarded\n"
            << "\n";

    }
    else
    {
        log << "Recorder log:\n"
            << n_trajs << " trajectories averaged\n"
            << n_data_samples
            << " trajectory samples (including horizontal sampling)\n"
            << n_discarded << " trajectories discarded\n"
            << 100*n_discarded/(n_discarded+n_trajs) << " % trajectories discarded\n"
            << discarded_t*timestep << " average length of a discarded trajectory\n"
            << discarded_t_stdev*timestep << "standard deviation of above\n"
            << n_aggr << " independent threads used to get errors in TCFs\n"
            << partfun << " average phase (for Matsubara dynamics)\n"
            << "\n";
    }
    log.close();
}
//*****************************************************************************
// TCF_recorder
//*****************************************************************************
void TCF_recorder::end_trajectory(const Cplx &phase, const Real &partfun_sample)
{
    if (used)
    {
        partfun += partfun_sample;
        // Use first trajectory obtain the trajectory and TCF length
        if (horizontal && autocorrelation)
        {
            if (first_trajectory)
            {
                is_data_complex = true;
                n_steps = index;
                if (n_steps%2==0)
                {
                    n_data = n_steps/2;
                    n_horizontal_samples = n_data+1;
                }
                else
                {
                    n_data = (n_steps-1)/2;
                    n_horizontal_samples = n_data+2;
                }
                data = Mat(n_data+1, n_quantities, a::fill::zeros);
            }
            if (index == n_steps)
            {
                for (int step=0; step<=n_data; step++)
                {
                    data.row(step) += a::real(phase
                     *a::sum(cplxtraj.rows(0,n_horizontal_samples-1)
                             %cplxtraj.rows(step,step+n_horizontal_samples-1)));
                }
            }
            else
            {
                std::cerr << "Incomplete trajectory => Trajectory discarded\n";
            }
        }
        else // not horizontal
        {
            if (first_trajectory)
            {
                is_data_complex = true;
                n_steps = index;
                n_data = n_steps;
                n_horizontal_samples = 1;
                data = Mat(n_data+1, n_quantities, a::fill::zeros);
            }
            if (index == n_steps)
            {
                if (autocorrelation)
                {
                    for (int step=0; step<=n_data; step++)
                    {
                        data.row(step) +=
                            a::real(phase*cplxtraj.row(0)%cplxtraj.row(step));
                    }
                }
                else
                {
                    for (int step=0; step<=n_data; step++)
                    {
                        data.row(step) +=
                            a::real(phase*cplx_init_val%cplxtraj.row(step));
                    }
                }
            }
            else
            {
                std::cerr << "Incomplete trajectory => Trajectory discarded\n";
            }

        }
        if (first_trajectory)
        {
            first_trajectory = false;
            // Set size of the trajectory vector to the actual size
            // Done after extraction of data such that elements can be lost
            cplxtraj.set_size(n_steps+1, n_quantities);
        }
        n_trajs += 1;
        index = -1;
        n_data_samples = n_trajs*n_horizontal_samples;
    }
}
//*****************************************************************************
// Average_recorder
//*****************************************************************************
void Average_recorder::end_trajectory(const Cplx &phase, const Real &partfun_sample)
{
    if (used)
    {
        partfun += partfun_sample;
        // Use first trajectory obtain the trajectory and TCF length
        if (first_trajectory)
        {
            first_trajectory = false;
            is_data_complex = true;
            n_steps = index;
            n_data = n_steps;
            data = Mat(n_data+1, n_quantities, a::fill::zeros);
            cplxtraj.resize(n_steps+1, n_quantities);
            n_horizontal_samples = 1;
        }
        if (index == n_steps)
        {
            data += a::real(phase*cplxtraj);
        }
        else
        {
            std::cerr << "Incomplete trajectory => Trajectory discarded\n";
        }
        n_trajs += 1;
        index = -1;
        n_data_samples = n_trajs*n_horizontal_samples;
    }
}
//*****************************************************************************
// Helper functions
//*****************************************************************************
//template <typename T>
//void aggregate_recorders(T & tot, T & src, const int &rank, const int &n_procs)
//{
//    if (rank==0)
//    {
//        tot.is_data_complex = src.is_data_complex;
//        tot.n_quantities = src.n_quantities;
//        tot.n_data = src.n_data;
//        tot.n_trajs = 0;
//        tot.n_discarded = 0;
//        tot.discarded_t = 0;
//        tot.discarded_t_stdev = 0;
//        tot.n_data_samples = 0;
//        tot.partfun = 0;
//        tot.n_aggr = n_procs;
//        tot.data.copy_size(src.data);
//        tot.data.fill(0);
//        tot.data_stdev.copy_size(src.data);
//        tot.data_stdev.fill(0);
//    }
//    MPI_Reduce(&src.n_data_samples, &tot.n_data_samples, 1,
//                                    MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&src.n_trajs, &tot.n_trajs, 1,
//                                    MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&src.n_discarded, &tot.n_discarded, 1,
//                                    MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&src.discarded_t, &tot.discarded_t, 1,
//                                    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&src.discarded_t_stdev, &tot.discarded_t_stdev, 1,
//                                    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    // Note that intermediate data is first calculated in each process
//    // before gathering it in the main process
//
//    // Deciding on whether data should be collected
//    if (src.n_trajs==0)
//    {
//        src.collect_data=false;
//    }
//    else
//    {
//        src.collect_data=true;
//    }
//    MPI_Allreduce(MPI_IN_PLACE, &src.collect_data, 1,
//            MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
//
//    if (src.collect_data)
//    {
//        // Calculating and summing TCF data
//        src.data_tmp = src.data
//                        /(src.n_data_samples*src.partfun/src.n_trajs);
//        MPI_Reduce(src.data_tmp.memptr(), tot.data.memptr(),
//                src.data_tmp.n_elem, MPI_DOUBLE,
//                MPI_SUM, 0, MPI_COMM_WORLD);
//        // Calculating and summing squares of TFC data
//        src.data_stdev = a::square(src.data_tmp);
//        MPI_Reduce(src.data_stdev.memptr(), tot.data_stdev.memptr(),
//                src.data_stdev.n_elem, MPI_DOUBLE,
//                MPI_SUM, 0, MPI_COMM_WORLD);
//        // Calculating and summing prtition function
//        src.partfun_tmp = src.partfun/src.n_trajs;
//        MPI_Reduce(&src.partfun_tmp, &tot.partfun,
//                1, MPI_DOUBLE,
//                MPI_SUM, 0, MPI_COMM_WORLD);
//    }
//    // Calculating errors
//    if (rank==0)
//    {
//        if (src.collect_data)
//        {
//            tot.data /= tot.n_aggr;
//            tot.data_stdev /= tot.n_aggr;
//            tot.partfun /= tot.n_aggr;
//            if (tot.n_aggr>1)
//            {
//                tot.data_stdev = std::sqrt(1.0/(tot.n_aggr-1))
//                        *a::sqrt(tot.data_stdev-a::square(tot.data));
//            }
//            else tot.data_stdev.fill(0);
//        }
//        tot.discarded_t /= tot.n_discarded;
//        tot.discarded_t_stdev /= tot.n_discarded;
//        if (tot.n_discarded>1)
//        {
//            tot.discarded_t_stdev = std::sqrt(1.0/(tot.n_discarded-1))
//                *std::sqrt(tot.discarded_t_stdev-tot.discarded_t*tot.discarded_t);
//        }
//        else tot.discarded_t_stdev = 0;
//    }
//}
//template void aggregate_recorders<>(TCF_recorder & tot, TCF_recorder & src, const int &rank, const int &n_procs);
//template void aggregate_recorders<>(Average_recorder & tot, Average_recorder & src, const int &rank, const int &n_procs);
