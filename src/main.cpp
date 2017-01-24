#include <memory>
#include <iostream>
#include <exception>
#include "appcross.h"

#ifdef _OPENMP
extern "C" {
    int omp_get_num_procs();
    int omp_get_max_threads();
    void omp_set_num_threads(int);
}
#endif

int main(int argc, char *argv[])
{
    std::cerr << "Cross v1.0 (Built on " __DATE__ " at " __TIME__ ")\n"
                 "  https://github.com/njau-sri/cross\n";

#ifdef _OPENMP
    int num_procs = omp_get_num_procs();
    if (num_procs < omp_get_max_threads() * 2)
        omp_set_num_threads(num_procs < 4 ? 1 : num_procs / 2);
#endif

    try {
        return std::make_shared<AppCross>()->run(argc, argv);
    }
    catch (const std::exception &e) {
        std::cerr << "FATAL: exception caught: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "FATAL: unknown exception caught\n";
        return 1;
    }

    return 0;
}
