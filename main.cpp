#include <iostream>
#include <exception>

int cross(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    try {
        return cross(argc, argv);
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
